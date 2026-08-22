#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
project_root="$(cd -- "$script_dir/../../.." && pwd -P)"
case_dir="$project_root/validation/cases/restart_wall_oblique"
results_root="$project_root/validation/results"
run_root="$results_root/restart_wall_oblique"

solver_input="${1:-$project_root/code/build/mfixsolver}"
if [[ "$solver_input" != /* ]]; then
  solver_input="$project_root/${solver_input#./}"
fi
if [[ ! -x "$solver_input" ]]; then
  echo "solver is not executable: $solver_input" >&2
  exit 1
fi
solver="$(cd -- "$(dirname -- "$solver_input")" && pwd -P)/$(basename -- "$solver_input")"

if [[ ! -f "$case_dir/mfix.dat" || ! -f "$case_dir/particle_input.dat" ]]; then
  echo "validation input is incomplete: $case_dir" >&2
  exit 1
fi

mkdir -p "$results_root"
expected_root="$project_root/validation/results/restart_wall_oblique"
if [[ "$run_root" != "$expected_root" ]]; then
  echo "refusing unexpected results path: $run_root" >&2
  exit 1
fi
rm -rf -- "$run_root"
mkdir -p "$run_root/continuous" "$run_root/stage1" "$run_root/restart1"

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

prepare_case "$run_root/continuous" "JKR_R26_WALL_CONT" "new" "6.0e-7"
run_case "$run_root/continuous"

prepare_case "$run_root/stage1" "JKR_R26_WALL_SPLIT" "new" "2.0e-7"
run_case "$run_root/stage1"

prepare_case "$run_root/restart1" "JKR_R26_WALL_SPLIT" "RESTART_1" "6.0e-7"
cp -- "$run_root/stage1/JKR_R26_WALL_SPLIT.RES" \
  "$run_root/restart1/JKR_R26_WALL_SPLIT.RES"
cp -- "$run_root/stage1"/JKR_R26_WALL_SPLIT.SP* "$run_root/restart1/"
cp -- "$run_root/stage1/JKR_R26_WALL_SPLIT.TXT" \
  "$run_root/restart1/JKR_R26_WALL_SPLIT.TXT"
if [[ -f "$run_root/stage1/JKR_R26_WALL_SPLIT_DES.RES" ]]; then
  cp -- "$run_root/stage1/JKR_R26_WALL_SPLIT_DES.RES" \
    "$run_root/restart1/JKR_R26_WALL_SPLIT_DES.RES"
fi
run_case "$run_root/restart1"

for branch in continuous restart1; do
  output="$run_root/$branch/particle_output.csv"
  if [[ ! -s "$output" ]]; then
    echo "missing output: $output" >&2
    exit 1
  fi
done

continuous_sha="$(sha256sum "$run_root/continuous/particle_output.csv" | awk '{print $1}')"
restart1_sha="$(sha256sum "$run_root/restart1/particle_output.csv" | awk '{print $1}')"
printf 'continuous %s\nrestart1 %s\n' "$continuous_sha" "$restart1_sha" \
  | tee "$run_root/sha256.txt"

if [[ "$continuous_sha" != "$restart1_sha" ]]; then
  echo "RESTART_1 particle output differs from continuous run" >&2
  exit 1
fi

echo "restart wall oblique regression passed: $run_root"
