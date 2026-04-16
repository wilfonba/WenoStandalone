#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 <application>" >&2
    exit 1
fi

APP="$1"
WENO_INPUT_NML="${WENO_INPUT_NML:-}"
ROCPROFV3="${ROCPROFV3:-/sw/crusher/ums/compilers/afar/therock-23.1.0-gfx90a-7.12.0-bb5005b6/bin/rocprofv3}"
ROCPROFV3_ARGS="${ROCPROFV3_ARGS:---kernel-trace --stats --summary}"
ROCPROFV3_OUTPUT_FORMAT="${ROCPROFV3_OUTPUT_FORMAT:-csv json rocpd}"
ROCPROFV3_OUTPUT_DIR="${ROCPROFV3_OUTPUT_DIR:-results/rocprofv3_weno5}"
ROCPROFV3_PMC_ARGS="${ROCPROFV3_PMC_ARGS:-}"
RUNNER="${RUNNER:-}"
OMP_TARGET_OFFLOAD_MODE="${OMP_TARGET_OFFLOAD_MODE:-MANDATORY}"
OMP_TEAMS_THREAD_LIMIT="${OMP_TEAMS_THREAD_LIMIT:-}"
ROCR_VISIBLE_DEVICE="${ROCR_VISIBLE_DEVICE:-0}"
WENO5_BENCH_SCOPE="${WENO5_BENCH_SCOPE:-full}"
WENO5_SPLIT_KERNELS="${WENO5_SPLIT_KERNELS:-0}"
WENO5_SPECIALIZED_COMBINED="${WENO5_SPECIALIZED_COMBINED:-0}"
WENO_WARMUP_ITERS="${WENO_WARMUP_ITERS:-1}"
WENO_BENCH_ITERS="${WENO_BENCH_ITERS:-5}"

mkdir -p "${ROCPROFV3_OUTPUT_DIR}"

if [[ ! -x "${ROCPROFV3}" ]]; then
    echo "rocprofv3 not found: ${ROCPROFV3}" >&2
    exit 1
fi

runner_cmd=()
if [[ -n "${RUNNER}" ]]; then
    # shellcheck disable=SC2206
    runner_cmd=(${RUNNER})
fi

# shellcheck disable=SC2206
rocprof_cmd=(${ROCPROFV3_ARGS})
output_format_cmd=()
if [[ -n "${ROCPROFV3_OUTPUT_FORMAT}" ]]; then
    # shellcheck disable=SC2206
    output_format_cmd=(-f ${ROCPROFV3_OUTPUT_FORMAT})
fi

pmc_cmd=()
if [[ -n "${ROCPROFV3_PMC_ARGS}" ]]; then
    # shellcheck disable=SC2206
    pmc_cmd=(${ROCPROFV3_PMC_ARGS})
fi

runtime_env=(
    ROCR_VISIBLE_DEVICE="${ROCR_VISIBLE_DEVICE}"
    ROCR_VISIBLE_DEVICES="${ROCR_VISIBLE_DEVICE}"
    OMP_TARGET_OFFLOAD="${OMP_TARGET_OFFLOAD_MODE}"
    WENO_KERNEL_MODE=weno5
)
if [[ -n "${OMP_TEAMS_THREAD_LIMIT}" ]]; then
    runtime_env+=(OMP_TEAMS_THREAD_LIMIT="${OMP_TEAMS_THREAD_LIMIT}")
fi

app_args=()
if [[ -n "${WENO_INPUT_NML}" ]]; then
    app_args+=("${WENO_INPUT_NML}")
fi

(set -x;
env "${runtime_env[@]}" \
    WENO5_BENCH_SCOPE="${WENO5_BENCH_SCOPE}" WENO5_SPLIT_KERNELS="${WENO5_SPLIT_KERNELS}" WENO5_SPECIALIZED_COMBINED="${WENO5_SPECIALIZED_COMBINED}" WENO_WARMUP_ITERS="${WENO_WARMUP_ITERS}" WENO_BENCH_ITERS="${WENO_BENCH_ITERS}" \
    "${runner_cmd[@]}" "${ROCPROFV3}" -d "${ROCPROFV3_OUTPUT_DIR}" -o weno5_profile "${output_format_cmd[@]}" "${rocprof_cmd[@]}" "${pmc_cmd[@]}" -- "${APP}" "${app_args[@]}" \
)
