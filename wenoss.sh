#!/usr/bin/env bash
#
# wenoss.sh — weno_standalone build + run + profile driver.
#
# Two usages:
#
#   1. Source to load compiler modules (replaces the old modules.sh):
#        source wenoss.sh load <cce_mp|amd_mp>
#      This runs `module purge`, loads the module set for that combo, and
#      exports WENO_COMPILER_COMBO so later invocations pick it up.
#
#   2. Execute to drive builds / runs / profiling via CMake:
#        ./wenoss.sh <subcommand> [COMPILER]
#
# Subcommands (execute mode):
#   build               configure + compile with CMake
#   build-save-temps    build a second binary with compiler-temps flags
#   run | benchmark     build, then launch the binary with WENO_* env vars
#   profile-weno5       rocprofv3 kernel / stats profile
#   profile-compute     rocprof-compute (Omniperf) profile
#   roofline            rocprof-compute --roof-only profile
#   list-pmc            dump available PMC counters via rocprofv3-avail
#   sweep               run toolchain/optimize_amdflang.sh
#   setup-venv          create the rocprof-compute Python venv (login node!)
#   test | check        run the binary and diff checksums against the golden file
#   test-update         regenerate tests/golden/<mode>.txt (shared across compilers)
#   clean               remove build/ and results/
#   help                show this message
#
# COMPILER defaults to $WENO_COMPILER_COMBO (set by sourcing this script)
# and falls back to amd_mp. Every knob below is env-overridable.

# ---------- module loading (used when sourced) ----------
# Defined up front so it's available to the sourced branch below. Avoids
# `set -u` / `exit`: sourcing must not nuke the caller's shell.
wenoss_load() {
    local combo="${1:-}"
    if [[ -z "$combo" ]]; then
        echo "Usage: source wenoss.sh load <cce_mp|amd_mp>" >&2
        return 1
    fi

    # Activate the rocprof-compute venv first (if it exists). Doing this
    # after `module purge` / `module load` drops the venv's PATH entries,
    # so activate up front and let the subsequent module ops layer on top.
    # Run `./wenoss.sh setup-venv` on a login node once to create the venv.
    local _wenoss_dir
    _wenoss_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    if [[ -f "${_wenoss_dir}/.venv-rocprof-compute/bin/activate" ]]; then
        # shellcheck disable=SC1091
        . "${_wenoss_dir}/.venv-rocprof-compute/bin/activate"
        echo "Activated venv: ${_wenoss_dir}/.venv-rocprof-compute"
    fi

    module purge
    # shellcheck disable=SC1091
    source /opt/cray/pe/cpe/25.09/restore_lmod_system_defaults.sh

    # NOTE: do NOT wrap `module load` in `( … )` — parentheses create a
    # subshell and the module state is lost when it exits. Use `{ … }` or
    # plain statements so the loads persist in the caller's shell.
    case "$combo" in
        cce_mp)
            module load PrgEnv-cray cpe/25.09 rocm/6.4.2 craype-accel-amd-gfx90a cray-python cce/20.0.0
            ;;
        cce_acc)
            module load cpe/25.09 rocm/6.4.2 craype-accel-amd-gfx90a cray-python PrgEnv-cray
            ;;
        amd_mp)
            module load amd/7.1.1 rocm/7.1.1 cray-python
            ;;
        *)
            echo "Unknown compiler combo: $combo" >&2
            echo "Supported: cce_mp, amd_mp (cce_acc module set also provided)" >&2
            return 1
            ;;
    esac

    echo "Loaded modules for $combo"
    export WENO_COMPILER_COMBO="$combo"
}

# When sourced, dispatch on the first arg. Only `load <combo>` is valid;
# anything else is rejected so execute-only subcommands can't silently run
# against the caller's shell. ${BASH_SOURCE[0]} != $0 signals sourcing.
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
    _wenoss_sub="${1:-}"
    case "${_wenoss_sub}" in
        load)
            shift
            wenoss_load "$@"
            _wenoss_rc=$?
            unset _wenoss_sub
            return ${_wenoss_rc}
            ;;
        "")
            echo "Usage: source wenoss.sh load <cce_mp|amd_mp>" >&2
            unset _wenoss_sub
            return 1
            ;;
        *)
            echo "source wenoss.sh: unknown subcommand '${_wenoss_sub}' (only 'load' is valid when sourcing)" >&2
            echo "For build/run/profile, execute: ./wenoss.sh ${_wenoss_sub} ..." >&2
            unset _wenoss_sub
            return 1
            ;;
    esac
fi

# ---------- executed (driver) mode ----------
set -euo pipefail

# ---------- project layout ----------
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_ROOT="${BUILD_ROOT:-${PROJECT_DIR}/build}"
RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/results}"

# ---------- compiler selection ----------
COMPILER="${COMPILER:-${2:-${WENO_COMPILER_COMBO:-amd_mp}}}"
case "${COMPILER}" in
    cce_mp|amd_mp) ;;
    cce_acc)
        echo "COMPILER=cce_acc is not supported: weno_standalone.f90 uses OpenMP target, not OpenACC." >&2
        exit 1 ;;
    *)
        echo "Unknown COMPILER \"${COMPILER}\". Supported: cce_mp, amd_mp" >&2
        exit 1 ;;
esac

BUILD_DIR="${BUILD_ROOT}/${COMPILER}"
BIN="${BUILD_DIR}/weno_standalone_${COMPILER}"
SAVE_TEMPS_BIN="${BUILD_DIR}/weno_standalone_${COMPILER}_save_temps"

# Namelist input file. Passed as the binary's first positional arg so it
# works from any cwd (the Fortran default is relative and brittle).
INPUT_NML="${INPUT_NML:-${PROJECT_DIR}/inputs/weno_input.nml}"

# ---------- runtime / profiling knobs (mirror old Makefile defaults) ----------
RUNNER="${RUNNER:-}"
ROCR_VISIBLE_DEVICE="${ROCR_VISIBLE_DEVICE:-0}"
OMP_TARGET_OFFLOAD_MODE="${OMP_TARGET_OFFLOAD_MODE:-MANDATORY}"
OMP_TEAMS_THREAD_LIMIT="${OMP_TEAMS_THREAD_LIMIT:-}"
WENO_KERNEL_MODE="${WENO_KERNEL_MODE:-weno5}"
WENO5_BENCH_SCOPE="${WENO5_BENCH_SCOPE:-full}"
WENO5_SPLIT_KERNELS="${WENO5_SPLIT_KERNELS:-0}"
WENO5_SPECIALIZED_COMBINED="${WENO5_SPECIALIZED_COMBINED:-0}"
WENO_WARMUP_ITERS="${WENO_WARMUP_ITERS:-5}"
WENO_BENCH_ITERS="${WENO_BENCH_ITERS:-10}"

# correctness test knobs — see cmd_test / cmd_test_update
TEST_GOLDEN_DIR="${TEST_GOLDEN_DIR:-${PROJECT_DIR}/tests/golden}"
WENO_TEST_RTOL="${WENO_TEST_RTOL:-1e-10}"
WENO_TEST_KERNEL_MODE="${WENO_TEST_KERNEL_MODE:-weno5}"
WENO_TEST_WARMUP_ITERS="${WENO_TEST_WARMUP_ITERS:-0}"
WENO_TEST_BENCH_ITERS="${WENO_TEST_BENCH_ITERS:-1}"

# rocprofv3
ROCPROFV3="${ROCPROFV3:-/sw/crusher/ums/compilers/afar/therock-23.1.0-gfx90a-7.12.0-bb5005b6/bin/rocprofv3}"
ROCPROFV3_AVAIL="${ROCPROFV3_AVAIL:-rocprofv3-avail}"
ROCPROFV3_ARGS="${ROCPROFV3_ARGS:--r --stats}"
ROCPROFV3_OUTPUT_FORMAT="${ROCPROFV3_OUTPUT_FORMAT:-csv}"
ROCPROFV3_OUTPUT_DIR="${ROCPROFV3_OUTPUT_DIR:-${RESULTS_DIR}/rocprofv3_weno5}"
ROCPROFV3_PMC_ARGS="${ROCPROFV3_PMC_ARGS:-}"

# rocprof-compute (Omniperf)
ROCPROF_COMPUTE="${ROCPROF_COMPUTE:-rocprof-compute}"
ROCPROF_COMPUTE_RUN_NAME="${ROCPROF_COMPUTE_RUN_NAME:-weno_standalone}"
ROCPROF_COMPUTE_OUTPUT_DIR="${ROCPROF_COMPUTE_OUTPUT_DIR:-${RESULTS_DIR}/rocprof_compute_weno5}"
ROCPROF_COMPUTE_PROFILE_ARGS="${ROCPROF_COMPUTE_PROFILE_ARGS:-}"

# Python venv for rocprof-compute
VENV_DIR="${VENV_DIR:-${PROJECT_DIR}/.venv-rocprof-compute}"
ROCPROF_COMPUTE_REQUIREMENTS="${ROCPROF_COMPUTE_REQUIREMENTS:-${PROJECT_DIR}/toolchain/requirements.txt}"

# ---------- helpers ----------
log() { printf '>> %s\n' "$*"; }
die() { printf '!! %s\n' "$*" >&2; exit 1; }

# Run env-vars that the binary + profiler shims read.
run_env=(
    "ROCR_VISIBLE_DEVICE=${ROCR_VISIBLE_DEVICE}"
    "ROCR_VISIBLE_DEVICES=${ROCR_VISIBLE_DEVICE}"
    "OMP_TARGET_OFFLOAD=${OMP_TARGET_OFFLOAD_MODE}"
    "WENO_KERNEL_MODE=${WENO_KERNEL_MODE}"
    "WENO5_BENCH_SCOPE=${WENO5_BENCH_SCOPE}"
    "WENO5_SPLIT_KERNELS=${WENO5_SPLIT_KERNELS}"
    "WENO5_SPECIALIZED_COMBINED=${WENO5_SPECIALIZED_COMBINED}"
    "WENO_WARMUP_ITERS=${WENO_WARMUP_ITERS}"
    "WENO_BENCH_ITERS=${WENO_BENCH_ITERS}"
)
[[ -n "${OMP_TEAMS_THREAD_LIMIT}" ]] && run_env+=("OMP_TEAMS_THREAD_LIMIT=${OMP_TEAMS_THREAD_LIMIT}")

# Shell-tokenize $RUNNER in case the caller passed e.g. "srun -n 1".
# shellcheck disable=SC2206
runner_argv=( ${RUNNER} )

cmake_configure() {
    mkdir -p "${BUILD_DIR}"
    cmake -S "${PROJECT_DIR}" -B "${BUILD_DIR}" -DCOMPILER="${COMPILER}"
}

cmd_build() {
    cmake_configure
    cmake --build "${BUILD_DIR}" --target "weno_standalone_${COMPILER}"
}

cmd_build_save_temps() {
    cmake_configure
    cmake --build "${BUILD_DIR}" --target build-save-temps
}

cmd_run() {
    cmd_build
    log "running ${BIN}"
    env "${run_env[@]}" "${runner_argv[@]}" "${BIN}" "${INPUT_NML}"
}

cmd_profile_weno5() {
    cmd_build
    mkdir -p "${RESULTS_DIR}"
    log "rocprofv3 profile -> ${ROCPROFV3_OUTPUT_DIR}"
    env "${run_env[@]}" \
        ROCPROFV3="${ROCPROFV3}" \
        ROCPROFV3_ARGS="${ROCPROFV3_ARGS}" \
        ROCPROFV3_OUTPUT_FORMAT="${ROCPROFV3_OUTPUT_FORMAT}" \
        ROCPROFV3_OUTPUT_DIR="${ROCPROFV3_OUTPUT_DIR}" \
        ROCPROFV3_PMC_ARGS="${ROCPROFV3_PMC_ARGS}" \
        OMP_TARGET_OFFLOAD_MODE="${OMP_TARGET_OFFLOAD_MODE}" \
        RUNNER="${RUNNER}" \
        WENO_INPUT_NML="${INPUT_NML}" \
        "${PROJECT_DIR}/toolchain/profile_weno5_rocprofv3.sh" "${BIN}"
}

check_venv() {
    [[ -f "${VENV_DIR}/bin/activate" ]] \
        || die "venv not found at ${VENV_DIR}. Run '$0 setup-venv' on a login node first."
}

activate_venv() {
    check_venv
    # shellcheck disable=SC1091
    . "${VENV_DIR}/bin/activate"
}

cmd_setup_venv() {
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        log "WARNING: looks like a compute node (SLURM_JOB_ID=${SLURM_JOB_ID})."
        log "         pip will likely fail without network. Run on a login node."
    fi
    [[ -f "${ROCPROF_COMPUTE_REQUIREMENTS}" ]] \
        || die "requirements file not found: ${ROCPROF_COMPUTE_REQUIREMENTS}"
    if [[ ! -f "${VENV_DIR}/bin/activate" ]]; then
        log "creating venv at ${VENV_DIR}"
        python -m venv "${VENV_DIR}"
    fi
    # shellcheck disable=SC1091
    . "${VENV_DIR}/bin/activate"
    pip install --upgrade pip
    pip install -r "${ROCPROF_COMPUTE_REQUIREMENTS}"
}

cmd_profile_compute() {
    cmd_build
    mkdir -p "${RESULTS_DIR}" "${ROCPROF_COMPUTE_OUTPUT_DIR}"
    activate_venv
    log "rocprof-compute profile -> ${ROCPROF_COMPUTE_OUTPUT_DIR}"
    # shellcheck disable=SC2086
    env "${run_env[@]}" \
        "${ROCPROF_COMPUTE}" profile \
            -n "${ROCPROF_COMPUTE_RUN_NAME}" \
            --path "${ROCPROF_COMPUTE_OUTPUT_DIR}" \
            ${ROCPROF_COMPUTE_PROFILE_ARGS} \
            -- "${runner_argv[@]}" "${BIN}" "${INPUT_NML}"
}

cmd_roofline() {
    cmd_build
    mkdir -p "${RESULTS_DIR}" "${ROCPROF_COMPUTE_OUTPUT_DIR}"
    activate_venv
    log "rocprof-compute --roof-only -> ${ROCPROF_COMPUTE_OUTPUT_DIR}"
    # shellcheck disable=SC2086
    env "${run_env[@]}" \
        "${ROCPROF_COMPUTE}" profile \
            -n "${ROCPROF_COMPUTE_RUN_NAME}" \
            --path "${ROCPROF_COMPUTE_OUTPUT_DIR}" \
            --roof-only \
            ${ROCPROF_COMPUTE_PROFILE_ARGS} \
            -- "${runner_argv[@]}" "${BIN}" "${INPUT_NML}"
}

cmd_list_pmc() {
    mkdir -p "${RESULTS_DIR}"
    bash -lc "module load rocm/7.2.0 && \"${ROCPROFV3_AVAIL}\" list --pmc" \
        | tee "${RESULTS_DIR}/rocprofv3_avail_pmc.txt"
}

cmd_sweep() {
    mkdir -p "${RESULTS_DIR}"
    env "${run_env[@]}" \
        OMP_TARGET_OFFLOAD_MODE="${OMP_TARGET_OFFLOAD_MODE}" \
        RUNNER="${RUNNER}" \
        "${PROJECT_DIR}/toolchain/optimize_amdflang.sh"
}

cmd_clean() {
    rm -rf "${BUILD_ROOT}" "${RESULTS_DIR}"
}

# ---------- correctness test helpers ----------
# The harness already prints `Checksum X/Y/Z/T : <value>` at the end of
# each run (sum(abs(arr)) over each output dimension + total). The test
# tool captures those four numbers, compares against a committed golden
# file, and fails if the relative error exceeds WENO_TEST_RTOL.
#
# Settings are pinned so the run is deterministic:
#   WENO_KERNEL_MODE=<WENO_TEST_KERNEL_MODE>, WENO5_BENCH_SCOPE=full,
#   WENO5_SPLIT_KERNELS=0, WENO5_SPECIALIZED_COMBINED=0,
#   warmup=${WENO_TEST_WARMUP_ITERS}, bench=${WENO_TEST_BENCH_ITERS}.

_test_golden_file() {
    # Shared across compilers: cce_mp and amd_mp run the same Fortran and
    # should agree on the checksum within WENO_TEST_RTOL. Keyed by kernel
    # mode only.
    echo "${TEST_GOLDEN_DIR}/${WENO_TEST_KERNEL_MODE}.txt"
}

# Run the binary with pinned env vars; binary's full stdout goes to stdout.
# Build output / log lines go to stderr so the caller can redirect stdout
# cleanly into a file for parsing.
_test_run_pinned() {
    cmd_build >&2
    log "running test: kernel_mode=${WENO_TEST_KERNEL_MODE} warmup=${WENO_TEST_WARMUP_ITERS} bench=${WENO_TEST_BENCH_ITERS}" >&2
    env \
        ROCR_VISIBLE_DEVICE="${ROCR_VISIBLE_DEVICE}" \
        ROCR_VISIBLE_DEVICES="${ROCR_VISIBLE_DEVICE}" \
        OMP_TARGET_OFFLOAD="${OMP_TARGET_OFFLOAD_MODE}" \
        WENO_KERNEL_MODE="${WENO_TEST_KERNEL_MODE}" \
        WENO5_BENCH_SCOPE=full \
        WENO5_SPLIT_KERNELS=0 \
        WENO5_SPECIALIZED_COMBINED=0 \
        WENO_WARMUP_ITERS="${WENO_TEST_WARMUP_ITERS}" \
        WENO_BENCH_ITERS="${WENO_TEST_BENCH_ITERS}" \
        "${runner_argv[@]}" "${BIN}" "${INPUT_NML}"
}

# From harness stdout, extract "<label> <value>" pairs for X/Y/Z/T.
_test_extract_checksums() {
    awk '/^  Checksum [XYZT] :/ { print $2, $4 }' "$1"
}

cmd_test() {
    local golden; golden="$(_test_golden_file)"
    [[ -f "${golden}" ]] \
        || die "golden file not found: ${golden}. Generate with: ./wenoss.sh test-update ${COMPILER}"

    local out; out="$(mktemp)"
    local actual; actual="$(mktemp)"
    # shellcheck disable=SC2064
    trap "rm -f '${out}' '${actual}'" RETURN

    _test_run_pinned > "${out}"
    _test_extract_checksums "${out}" > "${actual}"

    log "comparing against ${golden} (rtol=${WENO_TEST_RTOL})"
    python3 - "${actual}" "${golden}" "${WENO_TEST_RTOL}" <<'PYEOF'
import sys
actual_path, expected_path, rtol = sys.argv[1], sys.argv[2], float(sys.argv[3])

def read(path):
    with open(path) as f:
        return {k: float(v) for k, v in (line.split() for line in f if line.strip())}

actual, expected = read(actual_path), read(expected_path)
fail = False
for key in ('X', 'Y', 'Z', 'T'):
    a = actual.get(key, float('nan'))
    e = expected.get(key, float('nan'))
    denom = max(abs(e), 1e-300)
    rel = abs(a - e) / denom
    status = 'OK  ' if rel <= rtol else 'FAIL'
    print(f"  {key}  expected={e: .10e}  actual={a: .10e}  rel_err={rel:.3e}  {status}")
    if rel > rtol or a != a:  # NaN => fail
        fail = True
sys.exit(1 if fail else 0)
PYEOF
    log "test passed"
}

cmd_test_update() {
    mkdir -p "${TEST_GOLDEN_DIR}"
    local golden; golden="$(_test_golden_file)"

    local out; out="$(mktemp)"
    # shellcheck disable=SC2064
    trap "rm -f '${out}'" RETURN

    _test_run_pinned > "${out}"

    # Refuse to overwrite the golden with obviously-bad output (all zeros
    # = offload didn't execute; non-finite = NaN/Inf blew up).
    local chk_t
    chk_t="$(awk '/^  Checksum T :/ { print $4 }' "${out}")"
    [[ -n "${chk_t}" ]] || die "no 'Checksum T' line in binary output — did the run succeed?"
    python3 - "${chk_t}" <<'PYEOF'
import math, sys
v = float(sys.argv[1])
if not math.isfinite(v) or v == 0.0:
    sys.exit(f"refusing to write golden: Checksum T = {v!r}")
PYEOF

    _test_extract_checksums "${out}" > "${golden}"
    log "wrote golden: ${golden}"
    cat "${golden}"
}

cmd_help() { sed -n '3,31p' "${BASH_SOURCE[0]}"; }

# ---------- dispatch ----------
subcmd="${1:-help}"
case "${subcmd}" in
    build)              cmd_build ;;
    build-save-temps)   cmd_build_save_temps ;;
    run|benchmark)      cmd_run ;;
    profile-weno5)      cmd_profile_weno5 ;;
    profile-compute)    cmd_profile_compute ;;
    roofline)           cmd_roofline ;;
    list-pmc)           cmd_list_pmc ;;
    sweep)              cmd_sweep ;;
    setup-venv)         cmd_setup_venv ;;
    test|check)         cmd_test ;;
    test-update)        cmd_test_update ;;
    clean)              cmd_clean ;;
    help|-h|--help)     cmd_help ;;
    *) die "unknown subcommand: ${subcmd}. Try '$0 help'." ;;
esac
