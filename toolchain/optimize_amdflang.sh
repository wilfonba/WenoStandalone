#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC_FILE="${SRC_FILE:-${ROOT_DIR}/inputs/weno_standalone.f90}"
BUILD_DIR="${BUILD_DIR:-${ROOT_DIR}/build}"
RESULTS_DIR="${RESULTS_DIR:-${ROOT_DIR}/results}"
AMD_FLANG="${AMD_FLANG:-/sw/crusher/ums/compilers/afar/therock-23.1.0-gfx90a-7.12.0-bb5005b6/bin/amdflang}"
AMD_GPU_ARCH="${AMD_GPU_ARCH:-gfx90a}"
COMMON_FLAGS=(${COMMON_FLAGS:-})
OMP_FLAGS=(${OMP_FLAGS:--fopenmp -fopenmp-targets=amdgcn-amd-amdhsa --offload-arch=${AMD_GPU_ARCH}})
WARMUP_ITERS="${WENO_WARMUP_ITERS:-1}"
BENCH_ITERS="${WENO_BENCH_ITERS:-5}"
RUNNER="${RUNNER:-}"
OMP_TARGET_OFFLOAD_MODE="${OMP_TARGET_OFFLOAD_MODE:-MANDATORY}"
OMP_TEAMS_THREAD_LIMIT="${OMP_TEAMS_THREAD_LIMIT:-}"
WENO_KERNEL_MODE="${WENO_KERNEL_MODE:-all}"
WENO5_BENCH_SCOPE="${WENO5_BENCH_SCOPE:-full}"
WENO5_SPLIT_KERNELS="${WENO5_SPLIT_KERNELS:-0}"
WENO5_SPECIALIZED_COMBINED="${WENO5_SPECIALIZED_COMBINED:-0}"
ROCR_VISIBLE_DEVICE="${ROCR_VISIBLE_DEVICE:-0}"
RESULT_CSV="${RESULTS_DIR}/amdflang_optimization_results.csv"

mkdir -p "${BUILD_DIR}" "${RESULTS_DIR}"

if [[ ! -x "${AMD_FLANG}" ]]; then
    echo "amdflang not found: ${AMD_FLANG}" >&2
    exit 1
fi

runner_cmd=()
if [[ -n "${RUNNER}" ]]; then
    # shellcheck disable=SC2206
    runner_cmd=(${RUNNER})
fi

runtime_env=(
    ROCR_VISIBLE_DEVICE="${ROCR_VISIBLE_DEVICE}"
    ROCR_VISIBLE_DEVICES="${ROCR_VISIBLE_DEVICE}"
    OMP_TARGET_OFFLOAD="${OMP_TARGET_OFFLOAD_MODE}"
    WENO_KERNEL_MODE="${WENO_KERNEL_MODE}"
    WENO5_BENCH_SCOPE="${WENO5_BENCH_SCOPE}"
    WENO5_SPLIT_KERNELS="${WENO5_SPLIT_KERNELS}"
    WENO5_SPECIALIZED_COMBINED="${WENO5_SPECIALIZED_COMBINED}"
    WENO_WARMUP_ITERS="${WARMUP_ITERS}"
    WENO_BENCH_ITERS="${BENCH_ITERS}"
)
if [[ -n "${OMP_TEAMS_THREAD_LIMIT}" ]]; then
    runtime_env+=(OMP_TEAMS_THREAD_LIMIT="${OMP_TEAMS_THREAD_LIMIT}")
fi

variants=(
    "baseline|-O3"
    "ofast|-Ofast"
    "o3_amdopt|-O3 -famd-opt"
    "ofast_amdopt|-Ofast -famd-opt"
    "o3_amdopt_loop|-O3 -famd-opt -floop-interchange -fexperimental-loop-fusion"
    "ofast_amdopt_fastmath|-Ofast -famd-opt -ffast-math"
)

extract_metric() {
    local key="$1"
    awk -F':' -v name="${key}" '$1 ~ name {gsub(/^[ \t]+/, "", $2); sub(/[ \t]+s$/, "", $2); print $2}' | tail -n 1
}

extract_checksum() {
    awk -F':' '$1 ~ /Checksum T/ {gsub(/^[ \t]+/, "", $2); print $2}' | tail -n 1
}

printf 'variant,binary,compile_status,run_status,checksum_total,kernel_time_total_s,kernel_time_per_iter_s,flags\n' > "${RESULT_CSV}"

best_variant=""
best_time=""

for entry in "${variants[@]}"; do
    IFS='|' read -r variant_name variant_flags <<< "${entry}"
    binary_path="${BUILD_DIR}/weno_${variant_name}"
    compile_log="${RESULTS_DIR}/${variant_name}_compile.log"
    run_log="${RESULTS_DIR}/${variant_name}_run.log"

    # shellcheck disable=SC2206
    variant_flag_array=(${variant_flags})

    echo "==> Building ${variant_name}"
    if "${AMD_FLANG}" "${COMMON_FLAGS[@]}" "${OMP_FLAGS[@]}" "${variant_flag_array[@]}" "${SRC_FILE}" -o "${binary_path}" \
        >"${compile_log}" 2>&1; then
        compile_status="ok"
    else
        compile_status="fail"
        printf '%s,%s,%s,%s,%s,%s,%s,"%s"\n' \
            "${variant_name}" "${binary_path}" "${compile_status}" "not-run" "" "" "" "${variant_flags}" >> "${RESULT_CSV}"
        continue
    fi

    echo "==> Running ${variant_name}"
    if env "${runtime_env[@]}" "${runner_cmd[@]}" "${binary_path}" >"${run_log}" 2>&1; then
        run_status="ok"
    else
        run_status="fail"
    fi

    checksum_total="$(extract_checksum < "${run_log}")"
    kernel_total="$(extract_metric "Kernel time total" < "${run_log}")"
    kernel_iter="$(extract_metric "Kernel time per iter" < "${run_log}")"

    printf '%s,%s,%s,%s,%s,%s,%s,"%s"\n' \
        "${variant_name}" "${binary_path}" "${compile_status}" "${run_status}" "${checksum_total}" "${kernel_total}" "${kernel_iter}" "${variant_flags}" \
        >> "${RESULT_CSV}"

    if [[ "${run_status}" == "ok" && -n "${kernel_iter}" ]]; then
        if [[ -z "${best_time}" ]] || awk "BEGIN {exit !(${kernel_iter} < ${best_time})}"; then
            best_time="${kernel_iter}"
            best_variant="${variant_name}"
        fi
    fi
done

echo
echo "Results written to ${RESULT_CSV}"
if [[ -n "${best_variant}" ]]; then
    echo "Best variant: ${best_variant} (${best_time} s/iter)"
else
    echo "No successful benchmark runs were recorded."
fi
