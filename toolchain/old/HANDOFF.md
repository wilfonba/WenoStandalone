# WENO Standalone Handoff

## Goal

Current goal is no longer broad kernel coverage. The harness already covers the intended kernel families.

The active work is to optimize `s_weno5_kernel` in [`weno_standalone.f90`](./weno_standalone.f90) using:

- `rocprofv3` for GPU metrics
- `amdflang -save-temps` for generated device IR/assembly inspection

Primary optimization concerns:

- register pressure / possible spilling
- inefficient memory access patterns

## Current State

- The standalone harness still includes:
  - reshape/init loops for x/y/z
  - order-1 copy loops for x/y/z
  - WENO3 kernels
  - WENO5 kernels
  - WENO7 kernels
  - monotonicity-preserving kernels
- Benchmark controls were added:
  - `WENO_WARMUP_ITERS`
  - `WENO_BENCH_ITERS`
- Kernel selection was added:
  - `WENO_KERNEL_MODE=all`
  - `WENO_KERNEL_MODE=weno5`
- There is now a focused execution path for just the WENO5 kernels.
- Build/profiling automation now exists in:
  - [`Makefile`](./Makefile)
  - [`scripts/optimize_amdflang.sh`](./scripts/optimize_amdflang.sh)
  - [`scripts/profile_weno5_rocprofv3.sh`](./scripts/profile_weno5_rocprofv3.sh)
  - `make list-pmc`
- Runtime thread-limit plumbing is now explicit in the harness tooling:
  - `OMP_TEAMS_THREAD_LIMIT` is passed through by `make run`
  - `OMP_TEAMS_THREAD_LIMIT` is passed through by `make profile-weno5`
  - `OMP_TEAMS_THREAD_LIMIT` is passed through by `scripts/optimize_amdflang.sh`
  - unset `OMP_TEAMS_THREAD_LIMIT` is no longer exported as an empty string
- There is now an experimental split toggle for WENO5:
  - `WENO5_SPLIT_KERNELS=0` uses the original combined WENO5 kernel
  - `WENO5_SPLIT_KERNELS=1` dispatches separate left/right WENO5 kernels
  - the split path keeps the existing target-data lifetime, so it does not add CPU/GPU transfer traffic between the two WENO5 kernels
- There is now an experimental specialized-combined toggle for WENO5:
  - `WENO5_SPECIALIZED_COMBINED=0` uses the generic combined WENO5 kernel path
  - `WENO5_SPECIALIZED_COMBINED=1` dispatches x/y/z-specialized combined WENO5 kernels that reference the harness module coefficient/input arrays directly
  - this is currently a diagnostic experiment for descriptor/live-state pressure, not a final MFC-portable recommendation
- A first WENO5 source optimization pass is now in place in [`weno_standalone.f90`](./weno_standalone.f90):
  - fixed-size temporaries in `s_weno5_kernel` were scalarized
  - left/right reconstruction reuses scalar state instead of boxed small arrays
  - the OpenMP loop/clause structure was preserved
- One more small WENO5 cleanup is now in place in [`weno_standalone.f90`](./weno_standalone.f90):
  - the combined WENO5 kernel no longer recomputes the left-side `poly0/poly1/poly2` block twice
  - this preserved the current kernel structure and OpenMP clauses
- Progress and TODO tracking now lives in:
  - [`PROGRESS.md`](./PROGRESS.md)

## Important User Constraint

- Do not change the clause structure of the original kernels.
- Do not apply harness-specific optimizations that would change behavior in the MFC application.
- Always run all GPU-enabled executions outside the Codex sandbox.
- Always run all profiling, including every `rocprofv3` invocation and every mandatory-offload validation, outside the Codex sandbox.
- Treat outside-sandbox execution for GPU/profiling work as a hard rule, not a fallback option.
- Current compiler of record:
  - `/sw/crusher/ums/compilers/afar/therock-23.1.0-gfx90a-7.12.0-bb5005b6/bin/amdflang`
- Current profiler of record:
  - `/sw/crusher/ums/compilers/afar/therock-23.1.0-gfx90a-7.12.0-bb5005b6/bin/rocprofv3`
- Current runtime convention:
  - direct execution, no `srun`
  - export `ROCR_VISIBLE_DEVICE=0`
  - export `ROCR_VISIBLE_DEVICES=0`

## Validation Last Run

The latest successful execution in this environment was host fallback for WENO5-only mode:

```bash
OMP_TARGET_OFFLOAD=DISABLED \
WENO_KERNEL_MODE=weno5 \
WENO_WARMUP_ITERS=0 \
WENO_BENCH_ITERS=1 \
./build/weno_standalone_amdflang
```

Observed output:

```text
Checksum X : 1.3253000244E+06
Checksum Y : 1.3244180297E+06
Checksum Z : 1.3231500395E+06
Checksum T : 3.9728680936E+06
Kernel time per iter     : 0.155236 s
All kernel families executed with finite non-zero output.
```

`amdflang` build with offload flags also succeeds, and `make build-save-temps` emits device temp files.

Latest successful execution after the WENO5 scalarization pass:

```bash
OMP_TARGET_OFFLOAD=DISABLED \
WENO_KERNEL_MODE=weno5 \
WENO_WARMUP_ITERS=0 \
WENO_BENCH_ITERS=1 \
./build/weno_standalone_amdflang
```

Observed output:

```text
Checksum X : 1.3253000244E+06
Checksum Y : 1.3244180297E+06
Checksum Z : 1.3231500395E+06
Checksum T : 3.9728680936E+06
Kernel time per iter     : 0.088077 s
All kernel families executed with finite non-zero output.
```

Latest successful execution after removing the duplicate left-side `poly*` block:

```bash
OMP_TARGET_OFFLOAD=DISABLED \
WENO_KERNEL_MODE=weno5 \
WENO_WARMUP_ITERS=0 \
WENO_BENCH_ITERS=1 \
./build/weno_standalone_amdflang
```

Observed output:

```text
Checksum X : 1.3253000244E+06
Checksum Y : 1.3244180297E+06
Checksum Z : 1.3231500395E+06
Checksum T : 3.9728680936E+06
Kernel time per iter     : 0.063990 s
All kernel families executed with finite non-zero output.
```

Latest successful execution on the current explicit-shape plus staged-poly source state:

```bash
OMP_TARGET_OFFLOAD=DISABLED \
WENO_KERNEL_MODE=weno5 \
WENO_WARMUP_ITERS=0 \
WENO_BENCH_ITERS=1 \
./build/weno_standalone_amdflang
```

Observed output:

```text
Checksum X : 1.3253000244E+06
Checksum Y : 1.3244180297E+06
Checksum Z : 1.3231500395E+06
Checksum T : 3.9728680936E+06
Kernel time per iter     : 0.059108 s
All kernel families executed with finite non-zero output.
```

## Key Findings

- WENO5 device kernel symbol in generated assembly:
  - baseline: `__omp_offloading_8116438_e30011e2__QMm_weno_standalonePs_weno5_kernel_l449`
  - current: `__omp_offloading_8116438_e30011e2__QMm_weno_standalonePs_weno5_kernel_l455`
- Baseline WENO5 assembly metadata in `weno_standalone_amdflang_save_temps.amdgcn.gfx90a.img.lto.s`:
  - `num_vgpr = 256`
  - `num_agpr = 23`
  - `numbered_sgpr = 96`
- Baseline MLIR showed boxed private temporaries for WENO5:
  - `dvd`
  - `poly`
  - `beta`
  - `alpha`
  - `omega`
  - `delta`

These are the current strongest clues for register-pressure work.

New source-side change already made:

- `s_weno5_kernel` no longer uses the boxed fixed-size private arrays that were present in source for:
  - `dvd`
  - `poly`
  - `beta`
  - `alpha`
  - `omega`
  - `delta`
- Regenerated device MLIR now shows scalar WENO5 privates instead of boxed array privates.
- Refreshed assembly now shows a materially improved WENO5 footprint:
  - `num_vgpr = 167`
  - `num_agpr = 0`
  - `numbered_sgpr = 96`
  - `sgpr_spill_count = 34`
- First real `rocprofv3` run also succeeded outside the Codex sandbox and produced:
  - `results/rocprofv3_weno5/weno5_profile_results.db`
  - WENO5 kernel calls: `33`
  - WENO5 total duration: `2105300 ns`
  - WENO5 average duration: `6.380e+04 ns`
  - WENO5 inclusive percent: `55.423138`
- The profiling wrapper now supports explicit machine-readable outputs:
  - `ROCPROFV3_OUTPUT_FORMAT` default: `csv json rocpd`
  - generated files include CSV, JSON, and SQLite DB outputs under `results/rocprofv3_weno5`
- Available hardware counters can now be dumped with:
  - `make list-pmc`
  - output file: `results/rocprofv3_avail_pmc.txt`
- Raw PMC arguments can be passed through:
  - `ROCPROFV3_PMC_ARGS='--pmc COUNTER1 COUNTER2' make profile-weno5`
  - or multi-pass groups via repeated raw `--pmc ...` in `ROCPROFV3_PMC_ARGS`
- First PMC baseline already collected:
  - output directory: `results/rocprofv3_weno5_pmc_baseline`
  - counter file: `results/rocprofv3_weno5_pmc_baseline/weno5_profile_counter_collection.csv`
  - counter set:
    - `OccupancyPercent`
    - `MeanOccupancyPerCU`
    - `Wavefronts`
    - `MemUnitBusy`
    - `MemUnitStalled`
    - `SQ_INSTS_VMEM_RD`
    - `SQ_INSTS_VMEM_WR`
    - `TCC_HIT`
    - `TCC_MISS`
  - initial reading:
    - occupancy is only about `23%`
    - memory units are busy about `55%` to `57%`
    - inferred `TCC` hit rate is roughly `67%`
    - this is consistent with the remaining low-occupancy / register-pressure hypothesis
- Second PMC baseline already collected:
  - output directory: `results/rocprofv3_weno5_pmc_memory_detail`
  - counter file: `results/rocprofv3_weno5_pmc_memory_detail/weno5_profile_counter_collection.csv`
  - counter set:
    - `TA_BUSY_avr`
    - `TA_DATA_STALLED_BY_TC_CYCLES`
    - `TCC_BUSY_avr`
    - `TCC_TAG_STALL`
    - `TCP_PENDING_STALL_CYCLES`
    - `TCP_TCC_READ_REQ_LATENCY`
    - `TCP_TCP_LATENCY`
    - `WriteUnitStalled`
  - initial reading:
    - `TCP_PENDING_STALL_CYCLES` is consistently around `1.7M` to `2.0M`
    - `TCP_TCC_READ_REQ_LATENCY` is consistently around `81M` to `86M`
    - `TCP_TCP_LATENCY` is consistently around `192M` to `200M`
    - `TCC_TAG_STALL` is consistently around `24k` to `25k`
    - this suggests the remaining bottleneck is not only occupancy; there is also substantial TCP/TCC-side memory latency and pending-stall behavior
- Follow-on stencil-value cache is now in source:
  - `s_weno5_kernel` caches `v_rs_ws(j-2:j+2)` neighbors in scalars before forming `dvd`
  - refreshed assembly remains at the improved state:
    - `num_vgpr = 167`
    - `sgpr_spill_count = 32`
  - short validation profile under `results/rocprofv3_weno5_post_stencil_cache` reports:
    - WENO5 average duration: `6.363e+04 ns`
  - this change is effectively neutral-to-slightly-positive and is currently kept
- One broader experiment was discarded:
  - caching `dL` / `dR` into scalars pushed WENO5 to:
    - `num_vgpr = 215`
    - `sgpr_spill_count = 41`
  - that tradeoff was rejected
- Two additional portable source experiments were also rejected:
  - adding `contiguous` to the `s_weno5_kernel` dummy arrays pushed WENO5 to:
    - `num_vgpr = 165`
    - `sgpr_spill_count = 31`
  - dropping the cached outer stencil values `vjp2` / `vjm2` and forming `dvd_p1` / `dvd_m2` directly pushed WENO5 to:
    - `num_vgpr = 167`
    - `sgpr_spill_count = 34`
  - both tradeoffs were rejected
- One more structural experiment was also rejected:
  - splitting the left/right reconstruction work into inner `block` scopes preserved the host checksum baseline
  - regenerated device MLIR introduced explicit block-local allocas
  - refreshed WENO5 assembly stayed unchanged at:
    - `num_vgpr = 167`
    - `sgpr_spill_count = 34`
  - that change was reverted
- One more portable interface experiment was also rejected:
  - removing the `target` attribute from the generic WENO5 dummy arrays compiled and preserved the host checksum baseline
  - regenerated device MLIR for `_QMm_weno_standalonePs_weno5_kernel` dropped the dummy `fir.target` markings
  - authoritative `rocprofv3` timing under `results/rocprofv3_weno5_no_target_dummy` reports:
    - WENO5 total duration: `2115056 ns`
    - WENO5 average duration: `6.409e+04 ns`
  - relative to the current generic combined baseline:
    - `results/rocprofv3_weno5_current_baseline_fixed`
    - WENO5 total duration regressed slightly from `2097308 ns`
  - current interpretation:
    - removing `target` alone does not materially relieve the combined-kernel descriptor/live-state bottleneck
    - that source change was reverted
- One additional source cleanup was kept:
  - the combined WENO5 kernel had a duplicate left-side `poly0` / `poly1` / `poly2` evaluation block
  - removing it preserved the host checksum baseline and improved host-fallback time
  - refreshed save-temps still reports:
    - `num_vgpr = 167`
    - `num_agpr = 0`
    - `numbered_sgpr = 96`
  - current interpretation:
    - this is a valid cleanup, but it does not change the combined kernel's top-line register allocation under the current compiler
- One new combined-vs-split codegen finding is now established:
  - combined WENO5 save-temps currently reports:
    - `num_vgpr = 167`
    - `numbered_sgpr = 96`
    - `sgpr_spill_count = 34`
  - split left/right WENO5 save-temps currently report:
    - `num_vgpr = 126`
    - `numbered_sgpr = 92`
    - `sgpr_spill_count = 0`
  - current interpretation:
    - part of the remaining combined-kernel SGPR pressure is plausibly coming from carrying both left and right array descriptors/live state in the same kernel
    - this is now a stronger hypothesis than another obvious arithmetic redundancy in source
- One new authoritative runtime recheck was completed outside the sandbox on the current source state:
  - clean default runtime output:
    - `results/rocprofv3_weno5_current_baseline_fixed/weno5_profile_kernel_stats.csv`
    - WENO5 total duration: `2097308 ns`
    - WENO5 average duration: `6.355e+04 ns`
  - `OMP_TEAMS_THREAD_LIMIT=192` output:
    - `results/rocprofv3_weno5_threadlimit192_current/weno5_profile_kernel_stats.csv`
    - WENO5 total duration: `2079062 ns`
    - WENO5 average duration: `6.300e+04 ns`
  - current interpretation:
    - `OMP_TEAMS_THREAD_LIMIT=192` is still the best low-risk runtime knob tested so far for the current combined kernel
    - on the current source state, its measured WENO5 gain is about `0.9%`
- One narrow runtime sweep was also completed outside the sandbox on the current source state:
  - `OMP_TEAMS_THREAD_LIMIT=160`:
    - trace confirms `Workgroup_Size_X = 160`
    - output: `results/rocprofv3_weno5_threadlimit160_narrow/weno5_profile_kernel_stats.csv`
    - WENO5 average duration: `7.097e+04 ns`
  - `OMP_TEAMS_THREAD_LIMIT=192`:
    - trace confirms `Workgroup_Size_X = 192`
    - output: `results/rocprofv3_weno5_threadlimit192_narrow/weno5_profile_kernel_stats.csv`
    - WENO5 average duration: `6.255e+04 ns`
  - `OMP_TEAMS_THREAD_LIMIT=224`:
    - trace confirms `Workgroup_Size_X = 224`
    - output: `results/rocprofv3_weno5_threadlimit224_narrow/weno5_profile_kernel_stats.csv`
    - WENO5 average duration: `6.501e+04 ns`
  - current interpretation:
    - `192` remains the best team-size setting tested so far
    - `224` is slightly worse than `192`
    - `160` is clearly worse and should not be preferred
- One PMC comparison was also completed outside the sandbox for current default vs `OMP_TEAMS_THREAD_LIMIT=192`:
  - baseline PMC output:
    - `results/rocprofv3_weno5_current_baseline_pmc/weno5_profile_counter_collection.csv`
    - `OccupancyPercent`: about `23.24`
    - `MeanOccupancyPerCU`: about `7.44`
    - `Wavefronts`: `1512`
    - `MemUnitBusy`: about `55.25`
    - `MemUnitStalled`: about `0.513`
    - `SQ_INSTS_VMEM_RD`: `290304`
    - `SQ_INSTS_VMEM_WR`: `18144`
    - inferred `TCC` hit rate: about `66.4%`
  - `192`-thread PMC output:
    - `results/rocprofv3_weno5_threadlimit192_current_pmc/weno5_profile_counter_collection.csv`
    - `OccupancyPercent`: about `25.05`
    - `MeanOccupancyPerCU`: about `8.02`
    - `Wavefronts`: `1320`
    - `MemUnitBusy`: about `55.22`
    - `MemUnitStalled`: about `0.554`
    - `SQ_INSTS_VMEM_RD`: `290304`
    - `SQ_INSTS_VMEM_WR`: `18144`
    - inferred `TCC` hit rate: about `66.7%`
  - current interpretation:
    - the `192`-thread benefit is mainly an occupancy/work-distribution improvement
    - it does not reduce measured VMEM traffic
    - it does not materially change cache behavior
    - it does not solve the deeper memory-latency bottleneck
- One clause-based thread-limit experiment is now in source for the generic combined kernel:
  - source change:
    - added `thread_limit(192)` directly to the generic combined `s_weno5_kernel` OpenMP directive
  - host-fallback validation preserved the checksum baseline:
    - `Checksum T : 3.9728680936E+06`
  - short trace validation under `results/rocprofv3_weno5_clause_threadlimit192` confirms:
    - WENO5 `Workgroup_Size_X = 192`
    - WENO5 `Grid_Size_X = 84480`
  - important codegen note:
    - refreshed save-temps still reports the same top-line generic combined kernel metadata:
      - `num_vgpr = 167`
      - `numbered_sgpr = 96`
      - `sgpr_spill_count = 34`
    - the runtime trace is the authoritative confirmation that the clause changed the actual dispatched workgroup size
  - authoritative timing outside the sandbox:
    - clause-based output:
      - `results/rocprofv3_weno5_clause_threadlimit192_full/weno5_profile_kernel_stats.csv`
      - WENO5 total duration: `2051215 ns`
      - WENO5 average duration: `6.216e+04 ns`
    - relative to current generic baseline:
      - `results/rocprofv3_weno5_current_baseline_fixed/weno5_profile_kernel_stats.csv`
      - WENO5 total duration improves from `2097308 ns` to `2051215 ns`
      - measured gain is about `2.2%`
    - relative to prior runtime-only `192` tuning:
      - `results/rocprofv3_weno5_threadlimit192_current/weno5_profile_kernel_stats.csv`
      - WENO5 total duration improves from `2079062 ns` to `2051215 ns`
      - measured gain is about `1.3%`
  - clause-based PMC output:
    - `results/rocprofv3_weno5_clause_threadlimit192_pmc/weno5_profile_counter_collection.csv`
    - aggregated WENO5 reading:
      - `OccupancyPercent`: about `25.18`
      - `MeanOccupancyPerCU`: about `8.06`
      - `Wavefronts`: `1320`
      - `MemUnitBusy`: about `55.84`
      - `MemUnitStalled`: about `0.569`
      - `SQ_INSTS_VMEM_RD`: `290304`
      - `SQ_INSTS_VMEM_WR`: `18144`
      - inferred `TCC` hit rate: about `66.6%`
  - current interpretation:
    - this behaves similarly to the earlier runtime-only `192` tuning in terms of memory traffic and cache behavior
    - the gain still appears to come mainly from occupancy/work-distribution
    - within the generic combined-kernel variants tested so far, this is the best measured WENO5 timing
- One explicit-shape dummy experiment is now in source for the generic combined kernel:
  - source change:
    - converted the generic combined `s_weno5_kernel` dummy arrays from assumed-shape to explicit-shape bounds matching the harness allocations
    - retained the accepted `thread_limit(192)` clause in the same kernel
  - host-fallback validation preserved the checksum baseline:
    - `Checksum T : 3.9728680936E+06`
  - regenerated device MLIR for `_QMm_weno_standalonePs_weno5_kernel` now shows plain pointer dummy arguments without the earlier descriptor-style `fir.target` markings
  - authoritative runtime trace under `results/rocprofv3_weno5_explicit_shape_combined` confirms:
    - WENO5 `Workgroup_Size_X = 192`
    - WENO5 `Grid_Size_X = 84480`
    - trace-reported launch registers:
      - `VGPR_Count = 36`
      - `Accum_VGPR_Count = 132`
      - `SGPR_Count = 112`
  - refreshed save-temps for the generic combined kernel currently reports:
    - `num_vgpr = 187`
    - `numbered_sgpr = 96`
    - `sgpr_spill_count = 45`
    - `.max_flat_workgroup_size = 192`
  - important interpretation:
    - the top-line static register metadata looks worse than the prior clause-based generic baseline
    - but that metadata did not predict runtime correctly for this experiment
  - authoritative timing outside the sandbox:
    - explicit-shape output:
      - `results/rocprofv3_weno5_explicit_shape_combined/weno5_profile_kernel_stats.csv`
      - WENO5 total duration: `1502409 ns`
      - WENO5 average duration: `4.553e+04 ns`
    - relative to the clause-based generic combined baseline:
      - `results/rocprofv3_weno5_clause_threadlimit192_full/weno5_profile_kernel_stats.csv`
      - WENO5 total duration improves from `2051215 ns` to `1502409 ns`
      - measured gain is about `26.8%`
    - relative to the earlier generic default baseline:
      - `results/rocprofv3_weno5_current_baseline_fixed/weno5_profile_kernel_stats.csv`
      - WENO5 total duration improves from `2097308 ns` to `1502409 ns`
      - measured gain is about `28.4%`
  - explicit-shape PMC output:
    - `results/rocprofv3_weno5_explicit_shape_combined_pmc/weno5_profile_counter_collection.csv`
    - aggregated WENO5 reading:
      - `OccupancyPercent`: about `21.16`
      - `MeanOccupancyPerCU`: about `6.77`
      - `Wavefronts`: `1320`
      - `MemUnitBusy`: about `32.45`
      - `MemUnitStalled`: about `0.712`
      - `SQ_INSTS_VMEM_RD`: `65016`
      - `SQ_INSTS_VMEM_WR`: `18144`
      - inferred `TCC` hit rate: about `55.0%`
  - current interpretation:
    - this portable interface change behaves much more like the specialized-combined diagnostic path than like the earlier runtime-only `192` tuning
    - occupancy is lower than the clause-based generic baseline, but VMEM reads collapse from `290304` to `65016`
    - the strongest current hypothesis is that removing the generic assumed-shape descriptor machinery cuts effective memory traffic enough to dominate the worse top-line register metadata
    - this is now the best measured portable WENO5 source result in the tree so far
- One follow-on latency-staging change is now also in source for the explicit-shape generic combined kernel:
  - source change:
    - stages the `polyL(j, :, :)` / `polyR(j, :, :)` coefficient row values into reused just-in-time scalars before each left/right polynomial reconstruction
    - does not introduce broad persistent left/right weight caches like the rejected `dL` / `dR` staging
  - host-fallback validation preserved the checksum baseline:
    - `Checksum T : 3.9728680936E+06`
  - refreshed save-temps for the generic combined kernel reports:
    - `num_vgpr = 164`
    - `numbered_sgpr = 96`
    - `sgpr_spill_count = 5`
    - `.max_flat_workgroup_size = 192`
  - authoritative timing outside the sandbox:
    - staged-poly output:
      - `results/rocprofv3_weno5_explicit_shape_poly_stage/weno5_profile_kernel_stats.csv`
      - WENO5 total duration: `1482574 ns`
      - WENO5 average duration: `4.493e+04 ns`
    - relative to the plain explicit-shape generic baseline:
      - `results/rocprofv3_weno5_explicit_shape_combined/weno5_profile_kernel_stats.csv`
      - WENO5 total duration improves from `1502409 ns` to `1482574 ns`
      - measured gain is about `1.3%`
  - staged-poly PMC output:
    - `results/rocprofv3_weno5_explicit_shape_poly_stage_pmc/weno5_profile_counter_collection.csv`
    - aggregated WENO5 reading:
      - `OccupancyPercent`: about `20.83`
      - `MeanOccupancyPerCU`: about `6.67`
      - `Wavefronts`: `1320`
      - `MemUnitBusy`: about `32.09`
      - `MemUnitStalled`: about `0.696`
      - `SQ_INSTS_VMEM_RD`: `65016`
      - `SQ_INSTS_VMEM_WR`: `18144`
      - inferred `TCC` hit rate: about `55.0%`
  - current interpretation:
    - the extra gain does not come from another reduction in VMEM traffic
    - the likely benefit is lower instruction/addressing overhead plus improved live-state behavior around the polynomial coefficient loads
    - this is now the current best measured portable WENO5 source result in the tree
- One source experiment aimed specifically at combined-kernel descriptor/live-state pressure is now in the tree behind `WENO5_SPECIALIZED_COMBINED=1`:
  - it preserves the WENO5 checksum baseline in host fallback
  - save-temps for the specialized x/y/z combined kernels currently reports:
    - x kernel `__omp_offloading_8116438_e30011e2__QMm_weno_standalonePs_weno5_kernel_x_l624`
    - y kernel `__omp_offloading_8116438_e30011e2__QMm_weno_standalonePs_weno5_kernel_y_l782`
    - z kernel `__omp_offloading_8116438_e30011e2__QMm_weno_standalonePs_weno5_kernel_z_l940`
    - each specialized kernel currently shows:
      - `vgpr_count = 161`
      - `sgpr_count = 100`
      - `sgpr_spill_count = 4`
  - relative to the generic combined kernel:
    - `vgpr_count` improves from `167` to `161`
    - `sgpr_spill_count` improves from `34` to `4`
  - authoritative timing outside the sandbox:
    - generic combined baseline:
      - `results/rocprofv3_weno5_current_baseline_fixed/weno5_profile_kernel_stats.csv`
      - total duration: `2097308 ns`
    - specialized combined:
      - `results/rocprofv3_weno5_specialized_combined/weno5_profile_kernel_stats.csv`
      - total specialized WENO5 duration: `1525452 ns`
      - effective improvement vs generic combined is about `27%`
    - specialized combined with `OMP_TEAMS_THREAD_LIMIT=192`:
      - `results/rocprofv3_weno5_specialized_combined_threadlimit192/weno5_profile_kernel_stats.csv`
      - total specialized WENO5 duration: `1524487 ns`
      - interpretation:
        - the specialized source change dominates
        - `192` adds essentially no meaningful gain on top of it
  - specialized PMC output:
    - `results/rocprofv3_weno5_specialized_combined_pmc/weno5_profile_counter_collection.csv`
    - aggregated reading:
      - `OccupancyPercent`: about `20.79`
      - `MeanOccupancyPerCU`: about `6.66`
      - `Wavefronts`: `1512`
      - `MemUnitBusy`: about `34.27`
      - `MemUnitStalled`: about `0.600`
      - `SQ_INSTS_VMEM_RD`: `74088`
      - `SQ_INSTS_VMEM_WR`: `18144`
      - inferred `TCC` hit rate: about `52.5%`
  - current interpretation:
    - this experiment is much faster despite lower occupancy
    - the strongest measured shift is a large drop in VMEM reads plus much lower SGPR spilling
    - this is strong evidence that combined-kernel interface/live-descriptor pressure is a real bottleneck
    - because the implementation is harness-specialized, the result should currently be treated as a diagnostic clue for an MFC-valid refactor, not as a final accepted optimization
- One additional specialized-path experiment was rejected:
  - additionally changing the specialized x/y/z kernels to explicit-shape `vL` / `vR` dummy arguments preserved the host checksum baseline
  - refreshed save-temps for the specialized kernels remained effectively unchanged at:
    - `num_vgpr = 161`
    - `numbered_sgpr = 96`
    - `sgpr_spill_count = 4`
  - authoritative timing outside the sandbox under `results/rocprofv3_weno5_specialized_combined_explicit_shape` reports:
    - x kernel total duration: `672163 ns`
    - y kernel total duration: `664166 ns`
    - z kernel total duration: `512803 ns`
    - total specialized WENO5 duration: `1849132 ns`
  - relative to the earlier specialized-combined baseline:
    - `results/rocprofv3_weno5_specialized_combined/weno5_profile_kernel_stats.csv`
    - total specialized WENO5 duration regressed from `1525452 ns`
  - current interpretation:
    - explicit-shape `vL` / `vR` dummies do not help on top of the already-specialized x/y/z path
    - that source change was reverted
- Split-kernel experiment outcome:
  - the experimental split path lowers the per-kernel compiler footprint to:
    - left kernel: `num_vgpr = 126`, `sgpr_spill_count = 0`
    - right kernel: `num_vgpr = 126`, `sgpr_spill_count = 0`
  - but the current combined kernel still wins end-to-end:
    - plain profile baseline: `results/rocprofv3_weno5_split_eval_combined`
      - combined WENO5 total duration: `2085939 ns`
      - combined WENO5 average duration: `6.321e+04 ns`
    - split profile: `results/rocprofv3_weno5_split_eval_split`
      - left+right total duration: `2937464 ns`
      - average per left/right pair: about `8.901e+04 ns`
  - split PMC result:
    - combined PMC file: `results/rocprofv3_weno5_split_eval_combined_pmc/weno5_profile_counter_collection.csv`
    - split PMC file: `results/rocprofv3_weno5_split_eval_split_pmc/weno5_profile_counter_collection.csv`
    - combined WENO5:
      - occupancy about `23.21%`
      - `SQ_INSTS_VMEM_RD = 290304`
      - `SQ_INSTS_VMEM_WR = 18144`
      - inferred `TCC` hit rate about `67%`
    - split WENO5, summed over left+right:
      - occupancy averages only about `21.49%`
      - `SQ_INSTS_VMEM_RD = 417312`
      - `SQ_INSTS_VMEM_WR = 18144`
      - inferred `TCC` hit rate drops to about `50%`
  - conclusion:
    - the split path removes spilling but loses on launch overhead and duplicated memory traffic
    - current best-known direction remains the combined WENO5 kernel
- Runtime threads-per-team sweep outcome for the combined kernel:
  - this was explored with `OMP_TEAMS_THREAD_LIMIT`, so the kernel clause structure was not changed
  - the runtime knob does change actual launch shape:
    - baseline profile uses `Workgroup_Size_X = 256`
    - `OMP_TEAMS_THREAD_LIMIT=192` produces `Workgroup_Size_X = 192`
    - `OMP_TEAMS_THREAD_LIMIT=128` produces `Workgroup_Size_X = 128`
    - `OMP_TEAMS_THREAD_LIMIT=64` produces `Workgroup_Size_X = 64`
  - important constraint:
    - this does not change the compile-time WENO5 spill metadata
    - it is a runtime occupancy/scheduling lever, not a source-side spill reduction
  - timing sweep artifacts:
    - baseline: `results/rocprofv3_weno5_split_eval_combined/weno5_profile_kernel_stats.csv`
      - WENO5 total duration: `2085939 ns`
      - WENO5 average duration: `6.321e+04 ns`
    - `192` threads/team: `results/rocprofv3_weno5_threadlimit192/weno5_profile_kernel_stats.csv`
      - WENO5 total duration: `2040338 ns`
      - WENO5 average duration: `6.183e+04 ns`
    - `128` threads/team: `results/rocprofv3_weno5_threadlimit128/weno5_profile_kernel_stats.csv`
      - WENO5 total duration: `2092174 ns`
      - WENO5 average duration: `6.340e+04 ns`
    - `64` threads/team: `results/rocprofv3_weno5_threadlimit64/weno5_profile_kernel_stats.csv`
      - WENO5 total duration: `2457625 ns`
      - WENO5 average duration: `7.447e+04 ns`
  - PMC comparison for the best smaller team size:
    - baseline PMC file: `results/rocprofv3_weno5_split_eval_combined_pmc/weno5_profile_counter_collection.csv`
      - occupancy about `23.21%`
      - `MeanOccupancyPerCU` about `7.43`
      - `Wavefronts = 1512`
      - `SQ_INSTS_VMEM_RD = 290304`
      - `SQ_INSTS_VMEM_WR = 18144`
    - `192`-thread PMC file: `results/rocprofv3_weno5_threadlimit192_pmc/weno5_profile_counter_collection.csv`
      - occupancy about `24.98%`
      - `MeanOccupancyPerCU` about `7.99`
      - `Wavefronts = 1320`
      - `SQ_INSTS_VMEM_RD = 290304`
      - `SQ_INSTS_VMEM_WR = 18144`
  - conclusion:
    - `OMP_TEAMS_THREAD_LIMIT=192` is the only smaller team size tested so far that helps
    - the gain is modest, about `2%`
    - the main effect appears to be slightly better occupancy; memory-traffic counts stay essentially unchanged
    - `128` is neutral-to-slightly-worse, and `64` is clearly worse
- Historical note:
  - `results/rocprofv3_weno5_post_beta_simplify` and `results/rocprofv3_weno5_pmc_post_beta` came from an earlier harness-specific beta-term simplification pass
  - under the current user constraint, that shortcut is not a valid direction for future work and is not the basis for new source changes
- Updated occupancy/memory baseline for the current best version:
  - output directory: `results/rocprofv3_weno5_pmc_post_beta`
  - counter file: `results/rocprofv3_weno5_pmc_post_beta/weno5_profile_counter_collection.csv`
  - initial reading:
    - occupancy remains about `22%` to `23%`
    - `SQ_INSTS_VMEM_RD` dropped from `290304` to `263088`
    - `SQ_INSTS_VMEM_WR` stayed at `18144`
    - `TCC_MISS` stayed about `130k`
    - `TCC_HIT` dropped to about `177k` to `180k`
  - interpretation:
    - the latest accepted change reduced measured VMEM read traffic without changing occupancy much

## Current Blocker

GPU profiling is blocked in the Codex sandboxed shell, not by the source code or this node.

Confirmed diagnostics:

```text
TARGET AMDGPU RTL --> Failed to initialize AMDGPU's HSA library
Registered plugin AMDGPU with 0 visible device(s)
Skipping plugin AMDGPU with no visible devices
```

And at the HSA layer:

```text
Unable to open /dev/kfd read-write: No such file or directory
Failed to get user name to check for video group membership
```

Implication:

- the node and user environment do have working ROCm access outside the sandbox
- but the Codex sandbox hides `/dev/kfd` and `/dev/dri`
- so GPU commands must always be run with outside-sandbox escalation from this session

## Recommended Next Steps

1. Use outside-sandbox escalation for every GPU-enabled run and every profiling run from Codex, without exception.
2. Inspect:
   - `weno_standalone_amdflang_save_temps.amdgcn.gfx90a.img.lto.s`
   - `weno_standalone-openmp-amdgcn-amd-amdhsa-gfx90a-llvmir.mlir`
   - `results/rocprofv3_weno5/weno5_profile_results.db`
   - `results/rocprofv3_weno5/weno5_profile_kernel_stats.csv`
   - `results/rocprofv3_weno5/weno5_profile_results.json`
   - `results/rocprofv3_avail_pmc.txt`
   - `results/rocprofv3_weno5_pmc_baseline/weno5_profile_counter_collection.csv`
   - `results/rocprofv3_weno5_pmc_memory_detail/weno5_profile_counter_collection.csv`
3. Keep the current best portable generic-combined source state:
   - explicit-shape dummy arrays in `s_weno5_kernel`
   - staged `polyL` / `polyR` row scalars
   - explicit `thread_limit(192)` clause
4. Use `results/rocprofv3_weno5_explicit_shape_poly_stage` as the baseline for any new comparison, not the older assumed-shape generic kernel results.
5. Bias the next source change toward further portable reductions in descriptor/live-state or addressing overhead, using the specialized-combined path only as a diagnostic upper bound.
6. Re-profile after the next code change and keep only changes that improve target-GPU behavior.

For the latest detailed state and task tracking, see [`PROGRESS.md`](./PROGRESS.md).
