# WENO Standalone

A standalone OpenMP/OpenACC offload harness for the WENO reconstruction
kernels from MFC (`src/simulation/m_weno.fpp`). It is not a physics-faithful
mini-app; it allocates representative WENO work arrays with synthetic
non-zero data and exercises every device kernel family used by MFC's WENO
path so they can be benchmarked and profiled in isolation.

The harness runs the following kernel families:

1. Reshape kernels from `s_initialize_weno()` for x/y/z
2. Order-1 copy kernels from `s_weno()` for x/y/z
3. WENO3 reconstruction
4. WENO5 reconstruction
5. WENO7 reconstruction
6. Monotonicity-preserving kernel from `s_preserve_monotonicity()`

## Loading modules

`wenoss.sh` doubles as the module loader. Source it (don't execute) with
the `load` subcommand to load the module set for a compiler combination:

```bash
source wenoss.sh load cce_mp     # CCE with OpenMP offload
source wenoss.sh load amd_mp     # AMD Flang with OpenMP offload
```

This runs `module purge`, loads the right modules, exports
`WENO_COMPILER_COMBO`, and auto-activates the rocprof-compute venv if it
exists. `load` is the only subcommand accepted when sourcing — everything
else (build, run, profile, …) requires executing the script.

## Building

Builds go through CMake (`CMakeLists.txt`). The convenience wrapper
`wenoss.sh` calls CMake for you and handles run/profile/sweep on top of
it. Either path works.

Remember to `source wenoss.sh load <combo>` first so the right compiler
modules are loaded.

### Via `wenoss.sh` (recommended)

```bash
./wenoss.sh build cce_mp       # -> build/cce_mp/weno_standalone_cce_mp
./wenoss.sh build amd_mp       # -> build/amd_mp/weno_standalone_amd_mp
```

`COMPILER` defaults to `$WENO_COMPILER_COMBO` (set by sourcing `wenoss.sh`),
then falls back to `amd_mp`. Run `./wenoss.sh help` for the full list of
subcommands.

### Directly with CMake

```bash
cmake -S . -B build/amd_mp -DCOMPILER=amd_mp
cmake --build build/amd_mp
```

Per-combo flag baselines (encoded in `CMakeLists.txt`):

| Combo     | Offload flags                                                           | Extra flags                                     |
|-----------|-------------------------------------------------------------------------|-------------------------------------------------|
| `cce_mp`  | `-fopenmp`                                                              | `-h keepfiles`                                  |
| `cce_acc` | *not supported — source uses OpenMP target, not OpenACC*                | —                                               |
| `amd_mp`  | `-fopenmp -fopenmp-targets=amdgcn-amd-amdhsa --offload-arch=gfx90a`     | `-flto-partitions=8` in `OPT_FLAGS`             |

Common overrides (CMake cache variables — pass with `-D` at configure
time, or as env-vars to `wenoss.sh`):

```bash
cmake -S . -B build/amd_mp -DCOMPILER=amd_mp -DOPT_FLAGS="-O2"
cmake -S . -B build/amd_mp -DCOMPILER=amd_mp -DEXTRA_FLAGS="-g"
cmake -S . -B build/amd_mp -DCOMPILER=amd_mp -DAMD_GPU_ARCH=gfx942
./wenoss.sh build-save-temps amd_mp          # -save-temps / -fsave-loopmark variant
./wenoss.sh clean                            # remove build/ and results/
```

## Running

```bash
./wenoss.sh run amd_mp         # alias: benchmark
```

The `run` subcommand sets the standard offload environment and invokes
the binary. To run manually:

```bash
ROCR_VISIBLE_DEVICE=0 OMP_TARGET_OFFLOAD=MANDATORY ./build/amd_mp/weno_standalone_amd_mp
```

## Inputs

The harness takes inputs from two places:

### 1. `weno_input.nml` (Fortran namelist)

Lives at `inputs/weno_input.nml` and is passed to the binary automatically
by `wenoss.sh` (override with `INPUT_NML=...`). Manual invocations can
also supply any path as the first positional argument.

```fortran
&grid_params
  m_cells   = 199,       ! grid cells in x
  n_cells   = 199,       ! grid cells in y
  p_cells   = 199,       ! grid cells in z
  v_size    = 8,         ! number of conservative variables
  buff_size = 6,         ! halo/buffer width
  weno_eps  = 1e-16,     ! WENO regularization epsilon
  wenojs    = T,         ! use WENO-JS weights (vs. TENO)
  teno_CT   = 1e-5       ! TENO cutoff threshold (unused when wenojs=T)
/
```

### 2. Environment variables

| Variable                       | Default | Meaning                                                                 |
|--------------------------------|---------|-------------------------------------------------------------------------|
| `WENO_WARMUP_ITERS`            | `1`     | Warmup iterations before timing                                         |
| `WENO_BENCH_ITERS`             | `5`     | Timed benchmark iterations                                              |
| `WENO_KERNEL_MODE`             | `all`   | Which kernel family to run (`all`, `weno3`, `weno5`, `weno7`, ...)      |
| `WENO5_BENCH_SCOPE`            | `full`  | WENO5-only: `full` (with reshape/copy) or `kernel` (reconstruction only)|
| `WENO5_SPLIT_KERNELS`          | `0`     | WENO5: split combined kernel into discrete kernels                      |
| `WENO5_SPECIALIZED_COMBINED`   | `0`     | WENO5: use specialized combined kernel                                  |
| `OMP_TARGET_OFFLOAD`           | `MANDATORY` | Fail hard if offload unavailable (set by `wenoss.sh run`)           |
| `ROCR_VISIBLE_DEVICE[S]`       | `0`     | Which AMD GPU to use                                                    |

`wenoss.sh` forwards all of the `WENO_*` / `WENO5_*` variables to the
binary, so they can be overridden on the command line:

```bash
WENO_KERNEL_MODE=weno5 WENO_BENCH_ITERS=20 ./wenoss.sh run amd_mp
```

## Profiling

```bash
./wenoss.sh profile-weno5 amd_mp      # rocprofv3 trace for WENO5
./wenoss.sh list-pmc                  # list available PMC counters
```

Output lands in `results/` by default (override with `RESULTS_DIR=...`).

```bash
./wenoss.sh setup-venv        # Python venv for rocprof-compute (login node — needs internet)
./wenoss.sh profile-compute   # rocprof-compute profile
./wenoss.sh roofline          # --roof-only profile
./wenoss.sh sweep             # drive toolchain/optimize_amdflang.sh
```

## Correctness testing

The harness prints `Checksum X/Y/Z/T :` at the end of every run — four
reductions (`sum(abs(arr))`) over the final kernel output arrays.
`./wenoss.sh test` parses those four values, diffs them against a
committed golden file, and fails if the relative error exceeds
`WENO_TEST_RTOL` (default `1e-10`).

```bash
./wenoss.sh test-update amd_mp    # write tests/golden/all.txt
./wenoss.sh test amd_mp           # re-run, diff against the golden
./wenoss.sh test cce_mp           # same golden — cce_mp must agree within rtol
WENO_TEST_RTOL=1e-8 ./wenoss.sh test cce_mp
WENO_TEST_KERNEL_MODE=weno5 ./wenoss.sh test amd_mp   # separate golden per mode
```

Golden files are keyed only by kernel mode (`tests/golden/<mode>.txt`) and
shared across compilers — cce_mp and amd_mp run the same Fortran so they
should agree within `WENO_TEST_RTOL`. Cross-compiler floating-point
reassociation (OpenMP reduction ordering, `-O3` tree reductions, different
math libs) can exceed the default `1e-10`; if the test fails on one
compiler but passes on the one you used to write the golden, loosen the
tolerance before assuming a correctness regression.

The run is pinned to `warmup=0`, `bench=1`, `WENO5_SPLIT_KERNELS=0`,
`WENO5_SPECIALIZED_COMBINED=0` for determinism; `test-update` refuses to
commit an all-zero or non-finite checksum (would indicate offload didn't
execute).

## Layout

```
CMakeLists.txt         CMake build configuration
wenoss.sh              module loader (when sourced) + build / run / profile driver (when executed)
inputs/                harness source + default grid/WENO parameters
  weno_standalone.f90
  weno_input.nml
toolchain/             profiling + sweep helpers, venv requirements, archived notes
  profile_weno5_rocprofv3.sh
  optimize_amdflang.sh
  requirements.txt     (rocprof-compute Python deps, used by setup-venv)
  old/                 (archived HANDOFF.md / PROGRESS.md)
tests/golden/          committed reference checksums (see Correctness testing)
```
