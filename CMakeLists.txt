cmake_minimum_required(VERSION 3.20)

# Compiler selection: cce_mp, cce_acc, or amd_mp
#   cce_mp  - CCE (ftn) with OpenMP offload
#   cce_acc - CCE (ftn) with OpenACC offload (not supported — source uses OpenMP target)
#   amd_mp  - AMD Flang with OpenMP offload
# Load matching modules first via: source modules.sh <combo>
#
# Typical use:
#   source modules.sh amd_mp
#   cmake -S . -B build -DCOMPILER=amd_mp
#   cmake --build build
#
# `modules.sh` exports WENO_COMPILER_COMBO, which is used as the default if
# -DCOMPILER is not passed.
if(DEFINED ENV{WENO_COMPILER_COMBO} AND NOT DEFINED COMPILER)
    set(_default_compiler "$ENV{WENO_COMPILER_COMBO}")
else()
    set(_default_compiler "amd_mp")
endif()
set(COMPILER "${_default_compiler}" CACHE STRING
    "Compiler combo: cce_mp, cce_acc, amd_mp")
set_property(CACHE COMPILER PROPERTY STRINGS cce_mp cce_acc amd_mp)

set(AMD_FLANG
    "/sw/crusher/ums/compilers/afar/therock-23.1.0-gfx90a-7.12.0-bb5005b6/bin/amdflang"
    CACHE FILEPATH "Path to AMD Flang driver")
set(AMD_GPU_ARCH "gfx90a" CACHE STRING "AMD GPU architecture")
set(CCE_FLANG "ftn" CACHE FILEPATH "CCE Fortran driver")

# Select compiler and flags BEFORE project() so CMAKE_Fortran_COMPILER sticks.
if(COMPILER STREQUAL "cce_mp")
    set(CMAKE_Fortran_COMPILER "${CCE_FLANG}" CACHE FILEPATH "" FORCE)
    set(_offload_flags  "-fopenmp")
    set(_extra_flags    "-h;keepfiles")
    set(_opt_extra_flags "")
    set(_bin_name       "weno_standalone_cce_mp")
    set(_save_temps_flags "-fsave-loopmark")
    set(_save_temps_bin "weno_standalone_cce_mp_save_temps")
elseif(COMPILER STREQUAL "cce_acc")
    message(FATAL_ERROR
        "COMPILER=cce_acc is not supported: weno_standalone.f90 uses raw "
        "OpenMP target directives, not OpenACC. Use cce_mp or amd_mp.")
elseif(COMPILER STREQUAL "amd_mp")
    set(CMAKE_Fortran_COMPILER "${AMD_FLANG}" CACHE FILEPATH "" FORCE)
    set(_offload_flags  "-fopenmp;-fopenmp-targets=amdgcn-amd-amdhsa;--offload-arch=${AMD_GPU_ARCH}")
    set(_extra_flags    "")
    set(_opt_extra_flags "-flto-partitions=8")
    set(_bin_name       "weno_standalone_amd_mp")
    set(_save_temps_flags "-save-temps")
    set(_save_temps_bin "weno_standalone_amd_mp_save_temps")
else()
    message(FATAL_ERROR
        "Unknown COMPILER \"${COMPILER}\". Supported values: cce_mp, cce_acc, amd_mp")
endif()

project(weno_standalone LANGUAGES Fortran)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(OPT_FLAGS "-O3" CACHE STRING "Base optimization flags")
set(EXTRA_FLAGS "" CACHE STRING "Extra user flags appended to compile line")
separate_arguments(_opt_flags_list      NATIVE_COMMAND "${OPT_FLAGS}")
separate_arguments(_user_extra_list     NATIVE_COMMAND "${EXTRA_FLAGS}")
separate_arguments(_opt_extra_list      NATIVE_COMMAND "${_opt_extra_flags}")
separate_arguments(_save_temps_list     NATIVE_COMMAND "${_save_temps_flags}")

# Matches the Makefile: COMMON_FLAGS OFFLOAD_FLAGS OPT_FLAGS EXTRA_FLAGS
set(_common_compile_flags
    ${_offload_flags}
    ${_opt_flags_list}
    ${_opt_extra_list}
    ${_extra_flags}
    ${_user_extra_list})

set(_src ${CMAKE_CURRENT_SOURCE_DIR}/inputs/weno_standalone.f90)

# Primary build target
add_executable(${_bin_name} ${_src})
target_compile_options(${_bin_name} PRIVATE ${_common_compile_flags})
target_link_options(${_bin_name}    PRIVATE ${_offload_flags})

# `build-save-temps` equivalent — separate binary with extra flags.
# Opt-in via: cmake --build build --target build-save-temps
add_executable(${_save_temps_bin} EXCLUDE_FROM_ALL ${_src})
target_compile_options(${_save_temps_bin} PRIVATE
    ${_common_compile_flags} ${_save_temps_list})
target_link_options(${_save_temps_bin} PRIVATE
    ${_offload_flags} ${_save_temps_list})

add_custom_target(build           DEPENDS ${_bin_name})
add_custom_target(build-save-temps DEPENDS ${_save_temps_bin})

# `make run` equivalent — launch the freshly built binary with the same env
# vars the Makefile sets. Overridable at configure time with -D<var>=<value>.
set(RUNNER                        "" CACHE STRING "Optional launcher (e.g. srun) prefixing the binary")
set(ROCR_VISIBLE_DEVICE           "0"          CACHE STRING "")
set(OMP_TARGET_OFFLOAD_MODE       "MANDATORY"  CACHE STRING "")
set(OMP_TEAMS_THREAD_LIMIT        ""           CACHE STRING "")
set(WENO_KERNEL_MODE              "weno5"      CACHE STRING "")
set(WENO5_BENCH_SCOPE             "full"       CACHE STRING "")
set(WENO5_SPLIT_KERNELS           "0"          CACHE STRING "")
set(WENO5_SPECIALIZED_COMBINED    "0"          CACHE STRING "")
set(WENO_WARMUP_ITERS             "5"          CACHE STRING "")
set(WENO_BENCH_ITERS              "10"         CACHE STRING "")

set(_run_env
    ROCR_VISIBLE_DEVICE=${ROCR_VISIBLE_DEVICE}
    ROCR_VISIBLE_DEVICES=${ROCR_VISIBLE_DEVICE}
    OMP_TARGET_OFFLOAD=${OMP_TARGET_OFFLOAD_MODE}
    WENO_KERNEL_MODE=${WENO_KERNEL_MODE}
    WENO5_BENCH_SCOPE=${WENO5_BENCH_SCOPE}
    WENO5_SPLIT_KERNELS=${WENO5_SPLIT_KERNELS}
    WENO5_SPECIALIZED_COMBINED=${WENO5_SPECIALIZED_COMBINED}
    WENO_WARMUP_ITERS=${WENO_WARMUP_ITERS}
    WENO_BENCH_ITERS=${WENO_BENCH_ITERS})
if(OMP_TEAMS_THREAD_LIMIT)
    list(APPEND _run_env OMP_TEAMS_THREAD_LIMIT=${OMP_TEAMS_THREAD_LIMIT})
endif()

separate_arguments(_runner_list NATIVE_COMMAND "${RUNNER}")

add_custom_target(run
    COMMAND ${CMAKE_COMMAND} -E env ${_run_env}
            ${_runner_list} $<TARGET_FILE:${_bin_name}>
    DEPENDS ${_bin_name}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    USES_TERMINAL
    VERBATIM)

add_custom_target(benchmark DEPENDS run)

message(STATUS "weno_standalone: COMPILER=${COMPILER}")
message(STATUS "weno_standalone: Fortran compiler=${CMAKE_Fortran_COMPILER}")
message(STATUS "weno_standalone: primary target=${_bin_name}")
