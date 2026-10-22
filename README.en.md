# JKR Contact Model for Spherical Particles in MFIX 26.1.2

[简体中文](README.md) | English

This project is based on the official MFIX 26.1.2 source and implements JKR adhesive contact, sliding, rolling, and torsion for spherical particles. It also preserves contact history across neighbor-list rebuilds, parallel migration, and restart runs.

The current local source baseline is `9b5d639`, corresponding to the GitHub source release [`e72b226`](https://github.com/yuxing081/MfixVer26/commit/e72b22660868c2c4ac9d9ed16013f87659fadf35). The `code/` directory contains same-name override files relative to the official version; it is not a complete MFIX installation package.

## Main changes relative to official MFIX

| Change | Current implementation | Main code entry points |
|---|---|---|
| Spherical-particle JKR contact and resistance | Normal adhesion, tangential sliding, elastic rolling, rolling damping, and torsion for particle–particle and particle–wall contacts; supports phase, wall, and per-particle surface energies, plus reduced-stiffness parameter mapping for actual and simulated Young's moduli | [Particle pairs](code/calc_force_dem.f), [walls](code/calc_collision_wall_mod.f), [shared kernels](code/usr_mod.f), [parameter checks](code/check_solids_dem.f) |
| Contact-history direction | Transports tangential history `PFT` and rolling history `PFR` as the contact normal rotates and during common spin, before applying increments and limits; `PFTW` stores the relative twist angle; removes a false stop condition triggered by round-off residuals in near-zero histories | [History transport](code/usr_mod.f), [history variables](code/discretelement_mod.f) |
| MPI migration and neighbor rebuilds | Retains MFIX global particle IDs and contact slots, and completes history transfer when either endpoint migrates, including corner forwarding and restoration; an ownership migration triggers a neighbor rebuild so history is available before the next force evaluation | [Migration communication](code/mpi_funs_des_mod.f), [packing](code/mpi_pack_des_mod.f), [unpacking](code/mpi_unpack_des_mod.f), [neighbor management](code/neighbour.f) |
| Wall and neighbor geometry | Handles duplicate contacts at coplanar STL seams and safe wall-history transfer in the JKR path, and corrects the half-cell neighbor-pruning threshold | [Wall contact](code/calc_collision_wall_mod.f), [grid search](code/desgrid_mod.f) |
| Restart and stop-time saving | Saves JKR contact state, sliding/rolling/torsional history, the previous normal, and history time; the current JKR restart format is 2.2; a pure DEM run that receives a stop request writes a matching pair of restart files | [Restart read](code/read_res0_des.f), [restart write](code/write_res0_des.f), [stop handling](code/mfix.f) |
| Force-chain output and GUI input | Selects six contact fields by output region: contact state, contact radius, effective surface energy, tangential displacement, rolling displacement, and elastic twist angle; normalization uses a copy so multi-region output does not overwrite the original force magnitude; project keywords are registered | [VTP output](code/vtp_mod.f), [output parameters](code/output_control.inc), [GUI keywords](code/usr_init_namelist.f) |

Enable JKR with `USE_JKR=.TRUE.` while using the ordinary spherical Hertz contact path (`DES_COLL_MODEL='HERTZIAN'`). Non-spherical SuperDEM/GSP and LSD are unsupported, and JKR cannot be combined with the native `USE_COHESION`. JKR-specific physics and history extensions are controlled by the switch; the multi-region output and pure-DEM stop-saving fixes also apply when JKR is disabled.

Existing validation records cover contact response, history direction, wall geometry, selected MPI decompositions and migration cases, same-version restart, and short shear runs with one thousand particles. These results do not establish complete validation for arbitrary process counts, restarts after changing the partition, long-duration agglomerate breakup, or macroscopic rheology. The solver was not rerun for this release.

## Build and run

In a WSL/Linux environment with MFIX 26.1.2 installed, build from this repository's `code/` directory:

```sh
conda activate mfix-26.1.2
cd /path/to/MfixVer26/code
build_mfixsolver -j        # serial
build_mfixsolver --dmp -j  # MPI/DMP
```

Prepare `mfix.dat` and the required particle inputs in a case directory, then invoke the solver built above:

```sh
cd /path/to/case
/path/to/MfixVer26/code/mfixsolver
# DMP example: the product of the partition counts must match the process count
mpiexec --bind-to none -np 8 /path/to/MfixVer26/code/mfixsolver_dmp
```

The build uses the same-name override mechanism with the flat `code/` source layout. Do not modify the official installed source or use a leftover solver from an old build directory. For a restart that preserves contact history, use `RESTART_1` and keep both `RUN_NAME.RES` and `RUN_NAME_DES.RES` from the same run. The current JKR 2.2 read path is incompatible with older 2.0/2.1 history files.

## Documentation entry points

The table below describes the documentation in the complete local workspace. The GitHub version includes the published source and the Chinese and English READMEs; the other documents retain the historical versions already on the remote, and local documents that have not been synchronized are listed as paths. `validation/` is an independent repository; cloning this repository does not include its validation records or cases.

| Topic | Document |
|---|---|
| Current capabilities, validated scope, and known limitations | `docs/STATUS.md` |
| Formulas, parameters, reduced stiffness, rolling, torsion, and history direction | `docs/CONTACT_MODEL.md` |
| Contact identity, MPI migration, wall history, and restart format | `docs/CONTACT_HISTORY.md` |
| Build, time step, stop/restart, force-chain fields, and diagnostics | `docs/RUNNING_AND_OUTPUT.md` |
| Original evidence for upstream limitations and the handling at the time | [MFIX upstream limitations](docs/MFIX_UPSTREAM_LIMITATIONS.md) |
| Literature and external implementation references (the local index is newer than the remote copy) | [References](docs/references/README.md) |
| Migration process, conclusions from earlier versions, and original-text recovery methods | `docs/archive/CHANGE_RECORD.md` |

The upstream-limitations document retains its historical record; this page and the local `docs/STATUS.md` provide the current state. Detailed validation methods and tolerances are in the local `validation/docs/VALIDATION.md`; stage progress is in `validation/docs/TEST_PROGRESS.md`.

## Directory responsibilities

The table describes the complete local workspace. Except for `code/` and documents already present on the remote, the listed content is not necessarily published with this repository.

| Directory | Purpose |
|---|---|
| `code/` | Production Fortran files modified relative to the official version, kept flat to use the native same-name override mechanism |
| `model/` | Read-only snapshot of the official source; excluded from production repository version control |
| `docs/` | Current technical notes, reference index, and historical records |
| `tests/force_chain/` | Regression tests for the production output interface |
| `validation/` | Independent Git repository containing cases, reference formulas, validation entry points, and results |
| `UsrVer26/` | Independent Git repository containing user-customized versions and synchronization records |
| `FlowVer26/` | Flow cases and post-processing workspace |

Maintain each responsibility in its corresponding directory. Production code does not contain case-specific history seeds or energy statistics. The existing double-precision diagnostics are enabled explicitly with `MFIX_JKR_DIAGNOSTICS` and are off by default.

## Baseline and runtime conventions

- Official environment: `/home/yuxing/miniforge3/envs/mfix-26.1.2`, running in WSL Ubuntu-22.04.
- Official source: `/home/yuxing/miniforge3/envs/mfix-26.1.2/share/mfix/src`; do not build or modify files there.
- Legacy JKR source: MFIX 23.4.1, commit `ea55027498a4b939a2982cbb77d137201b034ffc`; only model and history semantics that were checked are migrated.
- The current JKR restart format is 2.2; do not try to make old history compatible by rewriting the file header.

Complete build and validation instructions are maintained in the local `docs/RUNNING_AND_OUTPUT.md`. Before starting new work, check the Git status of the production, validation, and user-version repositories separately, and preserve existing changes.
