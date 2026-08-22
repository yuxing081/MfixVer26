#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
project_root="$(cd -- "$script_dir/../../.." && pwd -P)"
case_dir="$project_root/validation/cases/restart_pair_pfr"
results_root="$project_root/validation/results"
run_root="$results_root/restart_pair_pfr"
solver="${1:-$project_root/code/build/mfixsolver}"

if [[ ! -x "$solver" ]]; then
  echo "solver is not executable: $solver" >&2
  exit 1
fi
solver="$(cd -- "$(dirname -- "$solver")" && pwd -P)/$(basename -- "$solver")"

if [[ ! -f "$case_dir/mfix.dat" || ! -f "$case_dir/particle_input.dat" ]]; then
  echo "validation input is incomplete: $case_dir" >&2
  exit 1
fi

mkdir -p "$results_root"
if [[ "$run_root" != "$results_root/restart_pair_pfr" ]]; then
  echo "refusing unexpected results path: $run_root" >&2
  exit 1
fi
rm -rf -- "$run_root"
mkdir -p "$run_root/continuous" "$run_root/stage1" \
  "$run_root/restart1" "$run_root/restart2"

prepare_case() {
  local destination="$1"
  local run_name="$2"
  local run_type="$3"
  local stop_time="$4"
  cp -- "$case_dir/mfix.dat" "$destination/mfix.dat"
  cp -- "$case_dir/particle_input.dat" "$destination/particle_input.dat"
  sed -i \
    -e "s/^run_name = .*/run_name = '$run_name'/" \
    -e "s/^run_type = .*/run_type = '$run_type'/" \
    -e "s/^tstop = .*/tstop = $stop_time/" \
    "$destination/mfix.dat"
}

run_case() {
  local destination="$1"
  (cd "$destination" && "$solver" > solver.log 2>&1)
}

prepare_case "$run_root/continuous" "JKR_R26_PFR_CONT" "new" "8.0e-8"
run_case "$run_root/continuous"

prepare_case "$run_root/stage1" "JKR_R26_PFR" "new" "5.0e-9"
run_case "$run_root/stage1"

for branch in restart1 restart2; do
  prepare_case "$run_root/$branch" "JKR_R26_PFR" "RESTART_1" "8.0e-8"
  cp -- "$run_root/stage1/JKR_R26_PFR.RES" \
    "$run_root/$branch/JKR_R26_PFR.RES"
  cp -- "$run_root/stage1"/JKR_R26_PFR.SP* "$run_root/$branch/"
  cp -- "$run_root/stage1/JKR_R26_PFR.TXT" \
    "$run_root/$branch/JKR_R26_PFR.TXT"
  if [[ -f "$run_root/stage1/JKR_R26_PFR_DES.RES" ]]; then
    cp -- "$run_root/stage1/JKR_R26_PFR_DES.RES" \
      "$run_root/$branch/JKR_R26_PFR_DES.RES"
  fi
done

run_case "$run_root/restart1"

sed -i "s/^run_type = .*/run_type = 'RESTART_2'/" \
  "$run_root/restart2/mfix.dat"
rm -f -- "$run_root/restart2"/JKR_R26_PFR.SP* \
  "$run_root/restart2/JKR_R26_PFR.TXT"
run_case "$run_root/restart2"

for branch in continuous restart1 restart2; do
  output="$run_root/$branch/particle_output.csv"
  if [[ ! -s "$output" ]]; then
    echo "missing output: $output" >&2
    exit 1
  fi
done

continuous_sha="$(sha256sum "$run_root/continuous/particle_output.csv" | awk '{print $1}')"
restart1_sha="$(sha256sum "$run_root/restart1/particle_output.csv" | awk '{print $1}')"
restart2_sha="$(sha256sum "$run_root/restart2/particle_output.csv" | awk '{print $1}')"
printf 'continuous %s\nrestart1 %s\nrestart2 %s\n' \
  "$continuous_sha" "$restart1_sha" "$restart2_sha" \
  | tee "$run_root/sha256.txt"

if [[ "$continuous_sha" != "$restart1_sha" ]]; then
  echo "RESTART_1 particle output differs from continuous run" >&2
  exit 1
fi
if [[ "$continuous_sha" == "$restart2_sha" ]]; then
  echo "RESTART_2 matches continuous"
else
  echo "RESTART_2 differs from continuous (sensitivity branch)"
fi

echo "restart pair PFR regression passed: $run_root"
