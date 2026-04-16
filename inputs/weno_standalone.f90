!> Standalone OpenMP offload harness for the kernels in src/simulation/m_weno.fpp
!!
!! This file is not a physics-faithful mini-app. Its job is narrower:
!! allocate representative WENO work arrays, populate them with smooth
!! non-trivial data, and execute every OpenMP target kernel family used by
!! m_weno.fpp:
!!   1. reshape kernels from s_initialize_weno() for x/y/z
!!   2. order-1 copy kernels from s_weno() for x/y/z
!!   3. WENO3 reconstruction kernel
!!   4. WENO5 reconstruction kernel
!!   5. WENO7 reconstruction kernel
!!   6. monotonicity-preserving kernel from s_preserve_monotonicity()
!!
!! The coefficient values are synthetic but finite and non-zero. They are
!! chosen to drive all stencil accesses and arithmetic paths needed to
!! exercise the device loops.
!!
!! Frontier build notes
!! --------------------
!! GCC toolchain requested for validation:
!!   /sw/crusher/ums/compilers/gcc/15.2.1-20260306/bin/gfortran \
!!     -fopenmp -fopenmp-targets=amdgcn-amdhsa \
!!     -foffload-options=amdgcn-amdhsa=-march=gfx90a \
!!     -ffree-line-length-none -O3 \
!!     -o weno_standalone_gcc weno_standalone.f90
!!
!! AMD LLVM / amdflang toolchain:
!!   /sw/crusher/ums/compilers/afar/therock-23.1.0-gfx90a-7.12.0-bb5005b6/bin/amdflang \
!!     -fopenmp -fopenmp-targets=amdgcn-amd-amdhsa --offload-arch=gfx90a \
!!     -O3 \
!!     -o weno_standalone_amdflang weno_standalone.f90
!!
!! CPU-only fallback with the GCC toolchain above:
!!   /sw/crusher/ums/compilers/gcc/15.2.1-20260306/bin/gfortran \
!!     -fopenmp -foffload=disable -ffree-line-length-none -O3 \
!!     -o weno_standalone_cpu weno_standalone.f90
!!
!! Example Frontier run:
!!   ROCR_VISIBLE_DEVICE=0 ROCR_VISIBLE_DEVICES=0 ./weno_standalone_gcc
!! or
!!   ROCR_VISIBLE_DEVICE=0 ROCR_VISIBLE_DEVICES=0 ./weno_standalone_amdflang
!!
!! Optional OpenMP environment settings for debugging:
!!   export OMP_TARGET_OFFLOAD=MANDATORY
!!   export LIBOMPTARGET_INFO=4
!!   export OMP_DISPLAY_ENV=VERBOSE

module m_weno_standalone

    use omp_lib
    use, intrinsic :: ieee_arithmetic
    implicit none

    integer, parameter :: wp = kind(1.0d0)

    ! MFC Non-parameter variables (always)
    integer :: m_cells = 199
    integer :: n_cells = 199
    integer :: p_cells = 199
    integer :: v_size  = 8
    integer :: buff_size = 6

    ! MFC Non-parameter variables (when using --case-optimization)
    real(wp) :: weno_eps = 1.0e-16_wp
    logical :: wenojs = .true.
    real(wp) :: teno_CT = 1.0e-5_wp

    ! MFC Parameter variables (when using --case-optimization)
    integer, parameter :: weno_polyn_max   = 3
    integer, parameter :: weno_stencils_max = 4
    integer, parameter :: beta_terms_max   = 6
    logical, parameter :: mapped_weno = .false.
    logical, parameter :: wenoz = .false.
    logical, parameter :: teno = .false.
    integer, parameter :: wenoz_q = 2

    ! Standalone specific options
    logical :: weno5_split_kernels = .false.
    logical :: weno5_specialized_combined = .false.

    ! Bounds variables
    integer :: is1_beg, is1_end
    integer :: is2_beg, is2_end
    integer :: is3_beg, is3_end
    integer :: core1_beg, core1_end

    real(wp), allocatable, target :: v_vf(:,:,:,:)
    real(wp), allocatable, target :: v_rs_ws_x(:,:,:,:), v_rs_ws_y(:,:,:,:), v_rs_ws_z(:,:,:,:)

    real(wp), allocatable, target :: poly_coef_cbL_x(:,:,:), poly_coef_cbR_x(:,:,:)
    real(wp), allocatable, target :: poly_coef_cbL_y(:,:,:), poly_coef_cbR_y(:,:,:)
    real(wp), allocatable, target :: poly_coef_cbL_z(:,:,:), poly_coef_cbR_z(:,:,:)

    real(wp), allocatable, target :: d_cbL_x(:,:), d_cbR_x(:,:)
    real(wp), allocatable, target :: d_cbL_y(:,:), d_cbR_y(:,:)
    real(wp), allocatable, target :: d_cbL_z(:,:), d_cbR_z(:,:)

    real(wp), allocatable, target :: beta_coef_x(:,:,:)
    real(wp), allocatable, target :: beta_coef_y(:,:,:)
    real(wp), allocatable, target :: beta_coef_z(:,:,:)

    !$omp declare target(v_vf, v_rs_ws_x, v_rs_ws_y, v_rs_ws_z, &
    !$omp&              poly_coef_cbL_x, poly_coef_cbR_x, poly_coef_cbL_y, poly_coef_cbR_y, &
    !$omp&              poly_coef_cbL_z, poly_coef_cbR_z, d_cbL_x, d_cbR_x, d_cbL_y, d_cbR_y, &
    !$omp&              d_cbL_z, d_cbR_z, beta_coef_x, beta_coef_y, beta_coef_z, &
    !$omp&              m_cells, n_cells, p_cells, v_size, &
    !$omp&              is1_beg, is1_end, is2_beg, is2_end, is3_beg, is3_end, core1_beg, core1_end)

contains

    subroutine s_read_input(input_file)

        character(len=*), intent(in) :: input_file
        integer :: iu, ios
        logical :: file_exists

        namelist /grid_params/ m_cells, n_cells, p_cells, v_size, &
            buff_size, weno_eps, wenojs, teno_CT

        inquire(file=trim(input_file), exist=file_exists)
        if (.not. file_exists) then
            write(*,'(a,a,a)') '  Input file "', trim(input_file), '" not found; using defaults.'
            return
        end if

        open(newunit=iu, file=trim(input_file), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,'(a,a,a)') '  WARNING: could not open "', trim(input_file), '"; using defaults.'
            return
        end if

        read(iu, nml=grid_params, iostat=ios)
        if (ios /= 0) then
            write(*,'(a)') '  WARNING: error reading namelist grid_params; using defaults.'
        end if

        close(iu)

        write(*,'(a,i0)') '  m_cells (from input)     : ', m_cells
        write(*,'(a,i0)') '  n_cells (from input)     : ', n_cells
        write(*,'(a,i0)') '  p_cells (from input)     : ', p_cells
        write(*,'(a,i0)') '  v_size  (from input)     : ', v_size

    end subroutine s_read_input

    subroutine s_initialize_data()

        integer :: ibeg, iend

        is1_beg = -buff_size
        is1_end = m_cells + buff_size
        is2_beg = -buff_size
        is2_end = n_cells + buff_size
        is3_beg = -buff_size
        is3_end = p_cells + buff_size

        core1_beg = is1_beg + weno_polyn_max
        core1_end = is1_end - weno_polyn_max

        ibeg = core1_beg
        iend = core1_end

        allocate(v_vf(is1_beg-weno_polyn_max:is1_end+weno_polyn_max, &
                      is2_beg-weno_polyn_max:is2_end+weno_polyn_max, &
                      is3_beg-weno_polyn_max:is3_end+weno_polyn_max, 1:v_size))

        allocate(v_rs_ws_x(is1_beg-weno_polyn_max:is1_end+weno_polyn_max, is2_beg:is2_end, is3_beg:is3_end, 1:v_size))
        allocate(v_rs_ws_y(is1_beg-weno_polyn_max:is1_end+weno_polyn_max, is2_beg:is2_end, is3_beg:is3_end, 1:v_size))
        allocate(v_rs_ws_z(is1_beg-weno_polyn_max:is1_end+weno_polyn_max, is2_beg:is2_end, is3_beg:is3_end, 1:v_size))

        allocate(poly_coef_cbL_x(ibeg:iend, 0:weno_stencils_max, 0:weno_polyn_max-1))
        allocate(poly_coef_cbR_x(ibeg:iend, 0:weno_stencils_max, 0:weno_polyn_max-1))
        allocate(poly_coef_cbL_y(ibeg:iend, 0:weno_stencils_max, 0:weno_polyn_max-1))
        allocate(poly_coef_cbR_y(ibeg:iend, 0:weno_stencils_max, 0:weno_polyn_max-1))
        allocate(poly_coef_cbL_z(ibeg:iend, 0:weno_stencils_max, 0:weno_polyn_max-1))
        allocate(poly_coef_cbR_z(ibeg:iend, 0:weno_stencils_max, 0:weno_polyn_max-1))

        allocate(d_cbL_x(0:weno_stencils_max, ibeg:iend))
        allocate(d_cbR_x(0:weno_stencils_max, ibeg:iend))
        allocate(d_cbL_y(0:weno_stencils_max, ibeg:iend))
        allocate(d_cbR_y(0:weno_stencils_max, ibeg:iend))
        allocate(d_cbL_z(0:weno_stencils_max, ibeg:iend))
        allocate(d_cbR_z(0:weno_stencils_max, ibeg:iend))

        allocate(beta_coef_x(ibeg:iend, 0:weno_stencils_max, 0:beta_terms_max-1))
        allocate(beta_coef_y(ibeg:iend, 0:weno_stencils_max, 0:beta_terms_max-1))
        allocate(beta_coef_z(ibeg:iend, 0:weno_stencils_max, 0:beta_terms_max-1))

        call s_fill_input_field()
        call s_fill_coefficients(poly_coef_cbL_x, poly_coef_cbR_x, d_cbL_x, d_cbR_x, beta_coef_x, 0.03_wp)
        call s_fill_coefficients(poly_coef_cbL_y, poly_coef_cbR_y, d_cbL_y, d_cbR_y, beta_coef_y, 0.05_wp)
        call s_fill_coefficients(poly_coef_cbL_z, poly_coef_cbR_z, d_cbL_z, d_cbR_z, beta_coef_z, 0.07_wp)

        !$omp target enter data map(alloc:v_vf, v_rs_ws_x, v_rs_ws_y, v_rs_ws_z)
        !$omp target enter data map(alloc:poly_coef_cbL_x, poly_coef_cbR_x, poly_coef_cbL_y, poly_coef_cbR_y)
        !$omp target enter data map(alloc:poly_coef_cbL_z, poly_coef_cbR_z, d_cbL_x, d_cbR_x, d_cbL_y, d_cbR_y)
        !$omp target enter data map(alloc:d_cbL_z, d_cbR_z, beta_coef_x, beta_coef_y, beta_coef_z)

        !$omp target update to(v_vf, poly_coef_cbL_x, poly_coef_cbR_x, poly_coef_cbL_y, poly_coef_cbR_y)
        !$omp target update to(poly_coef_cbL_z, poly_coef_cbR_z, d_cbL_x, d_cbR_x, d_cbL_y, d_cbR_y)
        !$omp target update to(d_cbL_z, d_cbR_z, beta_coef_x, beta_coef_y, beta_coef_z)
        !$omp target update to(m_cells, n_cells, p_cells, v_size)
        !$omp target update to(is1_beg, is1_end, is2_beg, is2_end, is3_beg, is3_end, core1_beg, core1_end)

    end subroutine s_initialize_data

    subroutine s_finalize_data()

        !$omp target exit data map(release:v_vf, v_rs_ws_x, v_rs_ws_y, v_rs_ws_z)
        !$omp target exit data map(release:poly_coef_cbL_x, poly_coef_cbR_x, poly_coef_cbL_y, poly_coef_cbR_y)
        !$omp target exit data map(release:poly_coef_cbL_z, poly_coef_cbR_z, d_cbL_x, d_cbR_x, d_cbL_y, d_cbR_y)
        !$omp target exit data map(release:d_cbL_z, d_cbR_z, beta_coef_x, beta_coef_y, beta_coef_z)

        deallocate(v_vf, v_rs_ws_x, v_rs_ws_y, v_rs_ws_z)
        deallocate(poly_coef_cbL_x, poly_coef_cbR_x, poly_coef_cbL_y, poly_coef_cbR_y)
        deallocate(poly_coef_cbL_z, poly_coef_cbR_z, d_cbL_x, d_cbR_x, d_cbL_y, d_cbR_y)
        deallocate(d_cbL_z, d_cbR_z, beta_coef_x, beta_coef_y, beta_coef_z)

    end subroutine s_finalize_data

    subroutine s_fill_input_field()

        integer :: i, j, k, l
        real(wp) :: x, y, z

        do l = 1, v_size
            do k = lbound(v_vf, 3), ubound(v_vf, 3)
                z = real(k - is3_beg + 1, wp)/real(max(1, is3_end - is3_beg + 3), wp)
                do j = lbound(v_vf, 2), ubound(v_vf, 2)
                    y = real(j - is2_beg + 1, wp)/real(max(1, is2_end - is2_beg + 3), wp)
                    do i = lbound(v_vf, 1), ubound(v_vf, 1)
                        x = real(i - is1_beg + 1, wp)/real(max(1, is1_end - is1_beg + 3), wp)
                        v_vf(i, j, k, l) = sin(2.1_wp*x*real(l, wp)) + cos(3.7_wp*y) + 0.25_wp*sin(4.3_wp*z) &
                                           + 0.02_wp*real((i + 2*j - k + l), wp)
                    end do
                end do
            end do
        end do

        v_rs_ws_x = 0.0_wp
        v_rs_ws_y = 0.0_wp
        v_rs_ws_z = 0.0_wp

    end subroutine s_fill_input_field

    subroutine s_fill_coefficients(polyL, polyR, dL, dR, beta_coef, offset)

        real(wp), intent(out) :: polyL(core1_beg:, 0:, 0:)
        real(wp), intent(out) :: polyR(core1_beg:, 0:, 0:)
        real(wp), intent(out) :: dL(0:, core1_beg:)
        real(wp), intent(out) :: dR(0:, core1_beg:)
        real(wp), intent(out) :: beta_coef(core1_beg:, 0:, 0:)
        real(wp), intent(in) :: offset

        integer :: j, r, s
        real(wp) :: baseL(0:weno_stencils_max), baseR(0:weno_stencils_max)

        baseL = [0.30_wp, 0.27_wp, 0.20_wp, 0.14_wp, 0.09_wp]
        baseR = [0.09_wp, 0.14_wp, 0.20_wp, 0.27_wp, 0.30_wp]

        do j = core1_beg, core1_end
            do r = 0, weno_stencils_max
                dL(r, j) = baseL(r) + offset*0.10_wp
                dR(r, j) = baseR(r) + offset*0.10_wp

                do s = 0, weno_polyn_max-1
                    polyL(j, r, s) = offset + 0.015_wp*real((r + 1)*(s + 1), wp)
                    polyR(j, r, s) = -offset + 0.012_wp*real((2*r + s + 1), wp)
                end do

                beta_coef(j, r, :) = 0.0_wp
                beta_coef(j, r, 0) = 0.8_wp + offset + 0.05_wp*real(r + 1, wp)
                beta_coef(j, r, 2) = 0.9_wp + offset + 0.04_wp*real(r + 1, wp)
                beta_coef(j, r, 3) = 1.0_wp + offset + 0.03_wp*real(r + 1, wp)
                beta_coef(j, r, 5) = 1.1_wp + offset + 0.02_wp*real(r + 1, wp)
            end do
        end do

    end subroutine s_fill_coefficients

    subroutine s_initialize_weno_x()
        integer :: j, k, l, q

        ! Source: src/simulation/m_weno.fpp:1211 ($:GPU_PARALLEL_LOOP)
        !$omp target teams distribute parallel do collapse(4)
        do j = 1, v_size
            do q = is3_beg, is3_end
                do l = is2_beg, is2_end
                    do k = is1_beg - weno_polyn_max, is1_end + weno_polyn_max
                        v_rs_ws_x(k, l, q, j) = v_vf(k, l, q, j)
                    end do
                end do
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine s_initialize_weno_x

    subroutine s_initialize_weno_y()
        integer :: j, k, l, q

        ! Source: src/simulation/m_weno.fpp:1228 ($:GPU_PARALLEL_LOOP)
        !$omp target teams distribute parallel do collapse(4)
        do j = 1, v_size
            do q = is3_beg, is3_end
                do l = is2_beg, is2_end
                    do k = is1_beg - weno_polyn_max, is1_end + weno_polyn_max
                        v_rs_ws_y(k, l, q, j) = v_vf(l, k, q, j)
                    end do
                end do
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine s_initialize_weno_y

    subroutine s_initialize_weno_z()
        integer :: j, k, l, q

        ! Source: src/simulation/m_weno.fpp:1245 ($:GPU_PARALLEL_LOOP)
        !$omp target teams distribute parallel do collapse(4)
        do j = 1, v_size
            do q = is3_beg, is3_end
                do l = is2_beg, is2_end
                    do k = is1_beg - weno_polyn_max, is1_end + weno_polyn_max
                        v_rs_ws_z(k, l, q, j) = v_vf(q, l, k, j)
                    end do
                end do
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine s_initialize_weno_z

    subroutine s_order1_copy_x(vL_x, vR_x)
        real(wp), intent(inout), target :: vL_x(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR_x(is1_beg:, is2_beg:, is3_beg:, 1:)
        integer :: i, j, k, l

        ! Source: src/simulation/m_weno.fpp:672 ($:GPU_PARALLEL_LOOP)
        !$omp target teams distribute parallel do collapse(4)
        do i = 1, v_size
            do l = is3_beg, is3_end
                do k = is2_beg, is2_end
                    do j = is1_beg, is1_end
                        vL_x(j, k, l, i) = v_vf(j, k, l, i)
                        vR_x(j, k, l, i) = v_vf(j, k, l, i)
                    end do
                end do
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine s_order1_copy_x

    subroutine s_order1_copy_y(vL_y, vR_y)
        real(wp), intent(inout), target :: vL_y(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR_y(is1_beg:, is2_beg:, is3_beg:, 1:)
        integer :: i, j, k, l

        ! Source: src/simulation/m_weno.fpp:685 ($:GPU_PARALLEL_LOOP)
        !$omp target teams distribute parallel do collapse(4)
        do i = 1, v_size
            do l = is3_beg, is3_end
                do k = is2_beg, is2_end
                    do j = is1_beg, is1_end
                        vL_y(j, k, l, i) = v_vf(k, j, l, i)
                        vR_y(j, k, l, i) = v_vf(k, j, l, i)
                    end do
                end do
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine s_order1_copy_y

    subroutine s_order1_copy_z(vL_z, vR_z)
        real(wp), intent(inout), target :: vL_z(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR_z(is1_beg:, is2_beg:, is3_beg:, 1:)
        integer :: i, j, k, l

        ! Source: src/simulation/m_weno.fpp:698 ($:GPU_PARALLEL_LOOP)
        !$omp target teams distribute parallel do collapse(4)
        do i = 1, v_size
            do l = is3_beg, is3_end
                do k = is2_beg, is2_end
                    do j = is1_beg, is1_end
                        vL_z(j, k, l, i) = v_vf(l, k, j, i)
                        vR_z(j, k, l, i) = v_vf(l, k, j, i)
                    end do
                end do
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine s_order1_copy_z

    subroutine s_weno3_kernel(v_rs_ws, polyL, polyR, dL, dR, beta_coef, vL, vR)

        real(wp), intent(in), target :: v_rs_ws(is1_beg-weno_polyn_max:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(in), target :: polyL(core1_beg:, 0:, 0:)
        real(wp), intent(in), target :: polyR(core1_beg:, 0:, 0:)
        real(wp), intent(in), target :: dL(0:, core1_beg:)
        real(wp), intent(in), target :: dR(0:, core1_beg:)
        real(wp), intent(in), target :: beta_coef(core1_beg:, 0:, 0:)
        real(wp), intent(inout), target :: vL(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR(is1_beg:, is2_beg:, is3_beg:, 1:)

        real(wp) :: dvd(-1:0), poly(0:1), alpha(0:1), omega(0:1), beta(0:1)
        real(wp) :: tau
        integer, parameter :: weno_num_stencils = 1
        integer :: i, j, k, l

        ! Source: src/simulation/m_weno.fpp:715 ($:GPU_PARALLEL_LOOP)
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& private(beta, dvd, poly, omega, alpha, tau)
        do l = is3_beg, is3_end
            do k = is2_beg, is2_end
                do j = is1_beg, is1_end
                    do i = 1, v_size
                        ! reconstruct from left side

                        alpha(:) = 0._wp
                        omega(:) = 0._wp
                        beta(:) = weno_eps

                        dvd(0) = v_rs_ws(j + 1, k, l, i) - v_rs_ws(j, k, l, i)
                        dvd(-1) = v_rs_ws(j, k, l, i) - v_rs_ws(j - 1, k, l, i)

                        poly(0) = v_rs_ws(j, k, l, i) + polyL(j, 0, 0)*dvd(0)
                        poly(1) = v_rs_ws(j, k, l, i) + polyL(j, 1, 0)*dvd(-1)

                        beta(0) = beta_coef(j, 0, 0)*dvd(0)*dvd(0) + weno_eps
                        beta(1) = beta_coef(j, 1, 0)*dvd(-1)*dvd(-1) + weno_eps

                        if (wenojs) then
                            alpha(0:weno_num_stencils) = dL(0:weno_num_stencils, j)/(beta(0:weno_num_stencils)**2._wp)

                        elseif (mapped_weno) then
                            alpha(0:weno_num_stencils) = dL(0:weno_num_stencils, j)/(beta(0:weno_num_stencils)**2._wp)
                            omega = alpha/sum(alpha)
                            alpha(0:weno_num_stencils) = (dL(0:weno_num_stencils, j)*(1._wp + dL(0:weno_num_stencils, j) - 3._wp*omega(0:weno_num_stencils)) + omega(0:weno_num_stencils)**2._wp) &
                                                         *(omega(0:weno_num_stencils)/(dL(0:weno_num_stencils, j)**2._wp + omega(0:weno_num_stencils)*(1._wp - 2._wp*dL(0:weno_num_stencils, j))))

                        elseif (wenoz) then
                            tau = abs(beta(1) - beta(0))
                            alpha(0:weno_num_stencils) = dL(0:weno_num_stencils, j)*(1._wp + tau/beta(0:weno_num_stencils))

                        end if

                        omega = alpha/sum(alpha)
                        vL(j, k, l, i) = omega(0)*poly(0) + omega(1)*poly(1)

                        ! reconstruct from right side

                        poly(0) = v_rs_ws(j, k, l, i) + polyR(j, 0, 0)*dvd(0)
                        poly(1) = v_rs_ws(j, k, l, i) + polyR(j, 1, 0)*dvd(-1)

                        if (wenojs) then
                            alpha(0:weno_num_stencils) = dR(0:weno_num_stencils, j)/(beta(0:weno_num_stencils)**2._wp)

                        elseif (mapped_weno) then
                            alpha(0:weno_num_stencils) = dR(0:weno_num_stencils, j)/(beta(0:weno_num_stencils)**2._wp)
                            omega = alpha/sum(alpha)
                            alpha(0:weno_num_stencils) = (dR(0:weno_num_stencils, j)*(1._wp + dR(0:weno_num_stencils, j) - 3._wp*omega(0:weno_num_stencils)) + omega(0:weno_num_stencils)**2._wp) &
                                                         *(omega(0:weno_num_stencils)/(dR(0:weno_num_stencils, j)**2._wp + omega(0:weno_num_stencils)*(1._wp - 2._wp*dR(0:weno_num_stencils, j))))

                        elseif (wenoz) then
                            alpha(0:weno_num_stencils) = dR(0:weno_num_stencils, j)*(1._wp + tau/beta(0:weno_num_stencils))

                        end if

                        omega = alpha/sum(alpha)
                        vR(j, k, l, i) = omega(0)*poly(0) + omega(1)*poly(1)
                    end do
                end do
            end do
        end do
        !$omp end target teams distribute parallel do

    end subroutine s_weno3_kernel

    subroutine s_weno5_kernel(v_rs_ws, polyL, polyR, dL, dR, beta_coef, vL, vR)

        real(wp), intent(in) :: v_rs_ws(is1_beg-weno_polyn_max:is1_end+weno_polyn_max, &
                                        is2_beg:is2_end, is3_beg:is3_end, 1:v_size)
        real(wp), intent(in) :: polyL(core1_beg:core1_end, 0:weno_stencils_max, 0:weno_polyn_max-1)
        real(wp), intent(in) :: polyR(core1_beg:core1_end, 0:weno_stencils_max, 0:weno_polyn_max-1)
        real(wp), intent(in) :: dL(0:weno_stencils_max, core1_beg:core1_end)
        real(wp), intent(in) :: dR(0:weno_stencils_max, core1_beg:core1_end)
        real(wp), intent(in) :: beta_coef(core1_beg:core1_end, 0:weno_stencils_max, 0:weno_polyn_max-1)
        real(wp), intent(inout) :: vL(is1_beg:is1_end, is2_beg:is2_end, is3_beg:is3_end, 1:v_size)
        real(wp), intent(inout) :: vR(is1_beg:is1_end, is2_beg:is2_end, is3_beg:is3_end, 1:v_size)

        real(wp) :: dvd_m2, dvd_m1, dvd_0, dvd_p1
        real(wp) :: poly0, poly1, poly2
        real(wp) :: poly_c00, poly_c01, poly_c10, poly_c11, poly_c20, poly_c21
        real(wp) :: beta0, beta1, beta2
        real(wp) :: alpha0, alpha1, alpha2
        real(wp) :: omega0, omega1, omega2
        real(wp) :: delta0, delta1, delta2
        real(wp) :: tau, sum_alpha, vj, vjm2, vjm1, vjp1, vjp2
        integer, parameter :: weno_num_stencils = 2
        integer :: i, j, k, l

        ! Source: src/simulation/m_weno.fpp:800 ($:GPU_PARALLEL_LOOP)
        ! Inner sequential GPU_LOOP sites in source: 804, 859, 868, 875, 913, 919
        !$omp target teams distribute parallel do collapse(3) thread_limit(192) &
        !$omp& private(dvd_m2, dvd_m1, dvd_0, dvd_p1, poly0, poly1, poly2, poly_c00, poly_c01, poly_c10, &
        !$omp&         poly_c11, poly_c20, poly_c21, beta0, beta1, beta2, &
        !$omp&         alpha0, alpha1, alpha2, omega0, omega1, omega2, tau, delta0, delta1, delta2, &
        !$omp&         sum_alpha, vj, vjm2, vjm1, vjp1, vjp2)
        do l = is3_beg, is3_end
            do k = is2_beg, is2_end
                do j = core1_beg, core1_end
                    do i = 1, v_size
                        ! reconstruct from left side

                        alpha0 = 0._wp
                        alpha1 = 0._wp
                        alpha2 = 0._wp
                        omega0 = 0._wp
                        omega1 = 0._wp
                        omega2 = 0._wp
                        delta0 = 0._wp
                        delta1 = 0._wp
                        delta2 = 0._wp

                        vj = v_rs_ws(j, k, l, i)
                        vjp2 = v_rs_ws(j + 2, k, l, i)
                        vjp1 = v_rs_ws(j + 1, k, l, i)
                        vjm1 = v_rs_ws(j - 1, k, l, i)
                        vjm2 = v_rs_ws(j - 2, k, l, i)

                        dvd_p1 = vjp2 - vjp1
                        dvd_0  = vjp1 - vj
                        dvd_m1 = vj - vjm1
                        dvd_m2 = vjm1 - vjm2

                        poly_c00 = polyL(j, 0, 0)
                        poly_c01 = polyL(j, 0, 1)
                        poly_c10 = polyL(j, 1, 0)
                        poly_c11 = polyL(j, 1, 1)
                        poly_c20 = polyL(j, 2, 0)
                        poly_c21 = polyL(j, 2, 1)

                        poly0 = vj + poly_c00*dvd_p1 + poly_c01*dvd_0
                        poly1 = vj + poly_c10*dvd_0 + poly_c11*dvd_m1
                        poly2 = vj + poly_c20*dvd_m1 + poly_c21*dvd_m2

                        beta0 = beta_coef(j, 0, 0)*dvd_p1*dvd_p1 + beta_coef(j, 0, 1)*dvd_p1*dvd_0 &
                              + beta_coef(j, 0, 2)*dvd_0*dvd_0 + weno_eps
                        beta1 = beta_coef(j, 1, 0)*dvd_0*dvd_0 + beta_coef(j, 1, 1)*dvd_0*dvd_m1 &
                              + beta_coef(j, 1, 2)*dvd_m1*dvd_m1 + weno_eps
                        beta2 = beta_coef(j, 2, 0)*dvd_m1*dvd_m1 + beta_coef(j, 2, 1)*dvd_m1*dvd_m2 &
                              + beta_coef(j, 2, 2)*dvd_m2*dvd_m2 + weno_eps

                        if (wenojs) then
                            alpha0 = dL(0, j)/(beta0**2._wp)
                            alpha1 = dL(1, j)/(beta1**2._wp)
                            alpha2 = dL(2, j)/(beta2**2._wp)

                        elseif (mapped_weno) then
                            alpha0 = dL(0, j)/(beta0**2._wp)
                            alpha1 = dL(1, j)/(beta1**2._wp)
                            alpha2 = dL(2, j)/(beta2**2._wp)
                            sum_alpha = alpha0 + alpha1 + alpha2
                            omega0 = alpha0/sum_alpha
                            omega1 = alpha1/sum_alpha
                            omega2 = alpha2/sum_alpha
                            alpha0 = (dL(0, j)*(1._wp + dL(0, j) - 3._wp*omega0) + omega0**2._wp) &
                                   *(omega0/(dL(0, j)**2._wp + omega0*(1._wp - 2._wp*dL(0, j))))
                            alpha1 = (dL(1, j)*(1._wp + dL(1, j) - 3._wp*omega1) + omega1**2._wp) &
                                   *(omega1/(dL(1, j)**2._wp + omega1*(1._wp - 2._wp*dL(1, j))))
                            alpha2 = (dL(2, j)*(1._wp + dL(2, j) - 3._wp*omega2) + omega2**2._wp) &
                                   *(omega2/(dL(2, j)**2._wp + omega2*(1._wp - 2._wp*dL(2, j))))

                        elseif (wenoz) then
                            tau = abs(beta2 - beta0)
                            alpha0 = dL(0, j)*(1._wp + (tau/beta0))
                            alpha1 = dL(1, j)*(1._wp + (tau/beta1))
                            alpha2 = dL(2, j)*(1._wp + (tau/beta2))

                        elseif (teno) then
                            tau = abs(beta2 - beta0)
                            alpha0 = (1._wp + tau/beta0)**6._wp
                            alpha1 = (1._wp + tau/beta1)**6._wp
                            alpha2 = (1._wp + tau/beta2)**6._wp
                            sum_alpha = alpha0 + alpha1 + alpha2
                            omega0 = alpha0/sum_alpha
                            omega1 = alpha1/sum_alpha
                            omega2 = alpha2/sum_alpha

                            if (omega0 < teno_CT) then
                                delta0 = 0._wp
                            else
                                delta0 = 1._wp
                            end if
                            if (omega1 < teno_CT) then
                                delta1 = 0._wp
                            else
                                delta1 = 1._wp
                            end if
                            if (omega2 < teno_CT) then
                                delta2 = 0._wp
                            else
                                delta2 = 1._wp
                            end if
                            alpha0 = delta0*dL(0, j)
                            alpha1 = delta1*dL(1, j)
                            alpha2 = delta2*dL(2, j)
                        end if

                        sum_alpha = alpha0 + alpha1 + alpha2
                        omega0 = alpha0/sum_alpha
                        omega1 = alpha1/sum_alpha
                        omega2 = alpha2/sum_alpha
                        vL(j, k, l, i) = omega0*poly0 + omega1*poly1 + omega2*poly2

                        ! reconstruct from right side

                        poly_c00 = polyR(j, 0, 0)
                        poly_c01 = polyR(j, 0, 1)
                        poly_c10 = polyR(j, 1, 0)
                        poly_c11 = polyR(j, 1, 1)
                        poly_c20 = polyR(j, 2, 0)
                        poly_c21 = polyR(j, 2, 1)

                        poly0 = vj + poly_c00*dvd_p1 + poly_c01*dvd_0
                        poly1 = vj + poly_c10*dvd_0 + poly_c11*dvd_m1
                        poly2 = vj + poly_c20*dvd_m1 + poly_c21*dvd_m2

                        if (wenojs) then
                            alpha0 = dR(0, j)/(beta0**2._wp)
                            alpha1 = dR(1, j)/(beta1**2._wp)
                            alpha2 = dR(2, j)/(beta2**2._wp)

                        elseif (mapped_weno) then
                            alpha0 = dR(0, j)/(beta0**2._wp)
                            alpha1 = dR(1, j)/(beta1**2._wp)
                            alpha2 = dR(2, j)/(beta2**2._wp)
                            sum_alpha = alpha0 + alpha1 + alpha2
                            omega0 = alpha0/sum_alpha
                            omega1 = alpha1/sum_alpha
                            omega2 = alpha2/sum_alpha
                            alpha0 = (dR(0, j)*(1._wp + dR(0, j) - 3._wp*omega0) + omega0**2._wp) &
                                   *(omega0/(dR(0, j)**2._wp + omega0*(1._wp - 2._wp*dR(0, j))))
                            alpha1 = (dR(1, j)*(1._wp + dR(1, j) - 3._wp*omega1) + omega1**2._wp) &
                                   *(omega1/(dR(1, j)**2._wp + omega1*(1._wp - 2._wp*dR(1, j))))
                            alpha2 = (dR(2, j)*(1._wp + dR(2, j) - 3._wp*omega2) + omega2**2._wp) &
                                   *(omega2/(dR(2, j)**2._wp + omega2*(1._wp - 2._wp*dR(2, j))))

                        elseif (wenoz) then
                            alpha0 = dR(0, j)*(1._wp + (tau/beta0))
                            alpha1 = dR(1, j)*(1._wp + (tau/beta1))
                            alpha2 = dR(2, j)*(1._wp + (tau/beta2))

                        elseif (teno) then
                            alpha0 = delta0*dR(0, j)
                            alpha1 = delta1*dR(1, j)
                            alpha2 = delta2*dR(2, j)
                        end if

                        sum_alpha = alpha0 + alpha1 + alpha2
                        omega0 = alpha0/sum_alpha
                        omega1 = alpha1/sum_alpha
                        omega2 = alpha2/sum_alpha
                        vR(j, k, l, i) = omega0*poly0 + omega1*poly1 + omega2*poly2
                    end do
                end do
            end do
        end do
        !$omp end target teams distribute parallel do

    end subroutine s_weno5_kernel

    subroutine s_weno5_kernel_x(vL, vR)

        real(wp), intent(inout), target :: vL(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR(is1_beg:, is2_beg:, is3_beg:, 1:)

        real(wp) :: dvd_m2, dvd_m1, dvd_0, dvd_p1
        real(wp) :: poly0, poly1, poly2
        real(wp) :: beta0, beta1, beta2
        real(wp) :: alpha0, alpha1, alpha2
        real(wp) :: omega0, omega1, omega2
        real(wp) :: delta0, delta1, delta2
        real(wp) :: tau, sum_alpha, vj, vjm2, vjm1, vjp1, vjp2
        integer :: i, j, k, l

        !$omp target teams distribute parallel do collapse(3) &
        !$omp& private(dvd_m2, dvd_m1, dvd_0, dvd_p1, poly0, poly1, poly2, beta0, beta1, beta2, &
        !$omp&         alpha0, alpha1, alpha2, omega0, omega1, omega2, tau, delta0, delta1, delta2, &
        !$omp&         sum_alpha, vj, vjm2, vjm1, vjp1, vjp2)
        do l = is3_beg, is3_end
            do k = is2_beg, is2_end
                do j = core1_beg, core1_end
                    do i = 1, v_size
                        alpha0 = 0._wp
                        alpha1 = 0._wp
                        alpha2 = 0._wp
                        omega0 = 0._wp
                        omega1 = 0._wp
                        omega2 = 0._wp
                        delta0 = 0._wp
                        delta1 = 0._wp
                        delta2 = 0._wp

                        vj = v_rs_ws_x(j, k, l, i)
                        vjp2 = v_rs_ws_x(j + 2, k, l, i)
                        vjp1 = v_rs_ws_x(j + 1, k, l, i)
                        vjm1 = v_rs_ws_x(j - 1, k, l, i)
                        vjm2 = v_rs_ws_x(j - 2, k, l, i)

                        dvd_p1 = vjp2 - vjp1
                        dvd_0  = vjp1 - vj
                        dvd_m1 = vj - vjm1
                        dvd_m2 = vjm1 - vjm2

                        poly0 = vj + poly_coef_cbL_x(j, 0, 0)*dvd_p1 + poly_coef_cbL_x(j, 0, 1)*dvd_0
                        poly1 = vj + poly_coef_cbL_x(j, 1, 0)*dvd_0 + poly_coef_cbL_x(j, 1, 1)*dvd_m1
                        poly2 = vj + poly_coef_cbL_x(j, 2, 0)*dvd_m1 + poly_coef_cbL_x(j, 2, 1)*dvd_m2

                        beta0 = beta_coef_x(j, 0, 0)*dvd_p1*dvd_p1 + beta_coef_x(j, 0, 1)*dvd_p1*dvd_0 &
                              + beta_coef_x(j, 0, 2)*dvd_0*dvd_0 + weno_eps
                        beta1 = beta_coef_x(j, 1, 0)*dvd_0*dvd_0 + beta_coef_x(j, 1, 1)*dvd_0*dvd_m1 &
                              + beta_coef_x(j, 1, 2)*dvd_m1*dvd_m1 + weno_eps
                        beta2 = beta_coef_x(j, 2, 0)*dvd_m1*dvd_m1 + beta_coef_x(j, 2, 1)*dvd_m1*dvd_m2 &
                              + beta_coef_x(j, 2, 2)*dvd_m2*dvd_m2 + weno_eps

                        if (wenojs) then
                            alpha0 = d_cbL_x(0, j)/(beta0**2._wp)
                            alpha1 = d_cbL_x(1, j)/(beta1**2._wp)
                            alpha2 = d_cbL_x(2, j)/(beta2**2._wp)
                        elseif (mapped_weno) then
                            alpha0 = d_cbL_x(0, j)/(beta0**2._wp)
                            alpha1 = d_cbL_x(1, j)/(beta1**2._wp)
                            alpha2 = d_cbL_x(2, j)/(beta2**2._wp)
                            sum_alpha = alpha0 + alpha1 + alpha2
                            omega0 = alpha0/sum_alpha
                            omega1 = alpha1/sum_alpha
                            omega2 = alpha2/sum_alpha
                            alpha0 = (d_cbL_x(0, j)*(1._wp + d_cbL_x(0, j) - 3._wp*omega0) + omega0**2._wp) &
                                   *(omega0/(d_cbL_x(0, j)**2._wp + omega0*(1._wp - 2._wp*d_cbL_x(0, j))))
                            alpha1 = (d_cbL_x(1, j)*(1._wp + d_cbL_x(1, j) - 3._wp*omega1) + omega1**2._wp) &
                                   *(omega1/(d_cbL_x(1, j)**2._wp + omega1*(1._wp - 2._wp*d_cbL_x(1, j))))
                            alpha2 = (d_cbL_x(2, j)*(1._wp + d_cbL_x(2, j) - 3._wp*omega2) + omega2**2._wp) &
                                   *(omega2/(d_cbL_x(2, j)**2._wp + omega2*(1._wp - 2._wp*d_cbL_x(2, j))))
                        elseif (wenoz) then
                            tau = abs(beta2 - beta0)
                            alpha0 = d_cbL_x(0, j)*(1._wp + (tau/beta0))
                            alpha1 = d_cbL_x(1, j)*(1._wp + (tau/beta1))
                            alpha2 = d_cbL_x(2, j)*(1._wp + (tau/beta2))
                        elseif (teno) then
                            tau = abs(beta2 - beta0)
                            alpha0 = (1._wp + tau/beta0)**6._wp
                            alpha1 = (1._wp + tau/beta1)**6._wp
                            alpha2 = (1._wp + tau/beta2)**6._wp
                            sum_alpha = alpha0 + alpha1 + alpha2
                            omega0 = alpha0/sum_alpha
                            omega1 = alpha1/sum_alpha
                            omega2 = alpha2/sum_alpha

                            if (omega0 < teno_CT) then
                                delta0 = 0._wp
                            else
                                delta0 = 1._wp
                            end if
                            if (omega1 < teno_CT) then
                                delta1 = 0._wp
                            else
                                delta1 = 1._wp
                            end if
                            if (omega2 < teno_CT) then
                                delta2 = 0._wp
                            else
                                delta2 = 1._wp
                            end if
                            alpha0 = delta0*d_cbL_x(0, j)
                            alpha1 = delta1*d_cbL_x(1, j)
                            alpha2 = delta2*d_cbL_x(2, j)
                        end if

                        sum_alpha = alpha0 + alpha1 + alpha2
                        omega0 = alpha0/sum_alpha
                        omega1 = alpha1/sum_alpha
                        omega2 = alpha2/sum_alpha
                        vL(j, k, l, i) = omega0*poly0 + omega1*poly1 + omega2*poly2

                        poly0 = vj + poly_coef_cbR_x(j, 0, 0)*dvd_p1 + poly_coef_cbR_x(j, 0, 1)*dvd_0
                        poly1 = vj + poly_coef_cbR_x(j, 1, 0)*dvd_0 + poly_coef_cbR_x(j, 1, 1)*dvd_m1
                        poly2 = vj + poly_coef_cbR_x(j, 2, 0)*dvd_m1 + poly_coef_cbR_x(j, 2, 1)*dvd_m2

                        if (wenojs) then
                            alpha0 = d_cbR_x(0, j)/(beta0**2._wp)
                            alpha1 = d_cbR_x(1, j)/(beta1**2._wp)
                            alpha2 = d_cbR_x(2, j)/(beta2**2._wp)
                        elseif (mapped_weno) then
                            alpha0 = d_cbR_x(0, j)/(beta0**2._wp)
                            alpha1 = d_cbR_x(1, j)/(beta1**2._wp)
                            alpha2 = d_cbR_x(2, j)/(beta2**2._wp)
                            sum_alpha = alpha0 + alpha1 + alpha2
                            omega0 = alpha0/sum_alpha
                            omega1 = alpha1/sum_alpha
                            omega2 = alpha2/sum_alpha
                            alpha0 = (d_cbR_x(0, j)*(1._wp + d_cbR_x(0, j) - 3._wp*omega0) + omega0**2._wp) &
                                   *(omega0/(d_cbR_x(0, j)**2._wp + omega0*(1._wp - 2._wp*d_cbR_x(0, j))))
                            alpha1 = (d_cbR_x(1, j)*(1._wp + d_cbR_x(1, j) - 3._wp*omega1) + omega1**2._wp) &
                                   *(omega1/(d_cbR_x(1, j)**2._wp + omega1*(1._wp - 2._wp*d_cbR_x(1, j))))
                            alpha2 = (d_cbR_x(2, j)*(1._wp + d_cbR_x(2, j) - 3._wp*omega2) + omega2**2._wp) &
                                   *(omega2/(d_cbR_x(2, j)**2._wp + omega2*(1._wp - 2._wp*d_cbR_x(2, j))))
                        elseif (wenoz) then
                            alpha0 = d_cbR_x(0, j)*(1._wp + (tau/beta0))
                            alpha1 = d_cbR_x(1, j)*(1._wp + (tau/beta1))
                            alpha2 = d_cbR_x(2, j)*(1._wp + (tau/beta2))
                        elseif (teno) then
                            alpha0 = delta0*d_cbR_x(0, j)
                            alpha1 = delta1*d_cbR_x(1, j)
                            alpha2 = delta2*d_cbR_x(2, j)
                        end if

                        sum_alpha = alpha0 + alpha1 + alpha2
                        omega0 = alpha0/sum_alpha
                        omega1 = alpha1/sum_alpha
                        omega2 = alpha2/sum_alpha
                        vR(j, k, l, i) = omega0*poly0 + omega1*poly1 + omega2*poly2
                    end do
                end do
            end do
        end do
        !$omp end target teams distribute parallel do

    end subroutine s_weno5_kernel_x

    subroutine s_weno5_kernel_y(vL, vR)

        real(wp), intent(inout), target :: vL(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR(is1_beg:, is2_beg:, is3_beg:, 1:)

        real(wp) :: dvd_m2, dvd_m1, dvd_0, dvd_p1
        real(wp) :: poly0, poly1, poly2
        real(wp) :: beta0, beta1, beta2
        real(wp) :: alpha0, alpha1, alpha2
        real(wp) :: omega0, omega1, omega2
        real(wp) :: delta0, delta1, delta2
        real(wp) :: tau, sum_alpha, vj, vjm2, vjm1, vjp1, vjp2
        integer :: i, j, k, l

        !$omp target teams distribute parallel do collapse(3) &
        !$omp& private(dvd_m2, dvd_m1, dvd_0, dvd_p1, poly0, poly1, poly2, beta0, beta1, beta2, &
        !$omp&         alpha0, alpha1, alpha2, omega0, omega1, omega2, tau, delta0, delta1, delta2, &
        !$omp&         sum_alpha, vj, vjm2, vjm1, vjp1, vjp2)
        do l = is3_beg, is3_end
            do k = is2_beg, is2_end
                do j = core1_beg, core1_end
                    do i = 1, v_size
                        alpha0 = 0._wp
                        alpha1 = 0._wp
                        alpha2 = 0._wp
                        omega0 = 0._wp
                        omega1 = 0._wp
                        omega2 = 0._wp
                        delta0 = 0._wp
                        delta1 = 0._wp
                        delta2 = 0._wp

                        vj = v_rs_ws_y(j, k, l, i)
                        vjp2 = v_rs_ws_y(j + 2, k, l, i)
                        vjp1 = v_rs_ws_y(j + 1, k, l, i)
                        vjm1 = v_rs_ws_y(j - 1, k, l, i)
                        vjm2 = v_rs_ws_y(j - 2, k, l, i)

                        dvd_p1 = vjp2 - vjp1
                        dvd_0  = vjp1 - vj
                        dvd_m1 = vj - vjm1
                        dvd_m2 = vjm1 - vjm2

                        poly0 = vj + poly_coef_cbL_y(j, 0, 0)*dvd_p1 + poly_coef_cbL_y(j, 0, 1)*dvd_0
                        poly1 = vj + poly_coef_cbL_y(j, 1, 0)*dvd_0 + poly_coef_cbL_y(j, 1, 1)*dvd_m1
                        poly2 = vj + poly_coef_cbL_y(j, 2, 0)*dvd_m1 + poly_coef_cbL_y(j, 2, 1)*dvd_m2

                        beta0 = beta_coef_y(j, 0, 0)*dvd_p1*dvd_p1 + beta_coef_y(j, 0, 1)*dvd_p1*dvd_0 &
                              + beta_coef_y(j, 0, 2)*dvd_0*dvd_0 + weno_eps
                        beta1 = beta_coef_y(j, 1, 0)*dvd_0*dvd_0 + beta_coef_y(j, 1, 1)*dvd_0*dvd_m1 &
                              + beta_coef_y(j, 1, 2)*dvd_m1*dvd_m1 + weno_eps
                        beta2 = beta_coef_y(j, 2, 0)*dvd_m1*dvd_m1 + beta_coef_y(j, 2, 1)*dvd_m1*dvd_m2 &
                              + beta_coef_y(j, 2, 2)*dvd_m2*dvd_m2 + weno_eps

                        if (wenojs) then
                            alpha0 = d_cbL_y(0, j)/(beta0**2._wp)
                            alpha1 = d_cbL_y(1, j)/(beta1**2._wp)
                            alpha2 = d_cbL_y(2, j)/(beta2**2._wp)
                        elseif (mapped_weno) then
                            alpha0 = d_cbL_y(0, j)/(beta0**2._wp)
                            alpha1 = d_cbL_y(1, j)/(beta1**2._wp)
                            alpha2 = d_cbL_y(2, j)/(beta2**2._wp)
                            sum_alpha = alpha0 + alpha1 + alpha2
                            omega0 = alpha0/sum_alpha
                            omega1 = alpha1/sum_alpha
                            omega2 = alpha2/sum_alpha
                            alpha0 = (d_cbL_y(0, j)*(1._wp + d_cbL_y(0, j) - 3._wp*omega0) + omega0**2._wp) &
                                   *(omega0/(d_cbL_y(0, j)**2._wp + omega0*(1._wp - 2._wp*d_cbL_y(0, j))))
                            alpha1 = (d_cbL_y(1, j)*(1._wp + d_cbL_y(1, j) - 3._wp*omega1) + omega1**2._wp) &
                                   *(omega1/(d_cbL_y(1, j)**2._wp + omega1*(1._wp - 2._wp*d_cbL_y(1, j))))
                            alpha2 = (d_cbL_y(2, j)*(1._wp + d_cbL_y(2, j) - 3._wp*omega2) + omega2**2._wp) &
                                   *(omega2/(d_cbL_y(2, j)**2._wp + omega2*(1._wp - 2._wp*d_cbL_y(2, j))))
                        elseif (wenoz) then
                            tau = abs(beta2 - beta0)
                            alpha0 = d_cbL_y(0, j)*(1._wp + (tau/beta0))
                            alpha1 = d_cbL_y(1, j)*(1._wp + (tau/beta1))
                            alpha2 = d_cbL_y(2, j)*(1._wp + (tau/beta2))
                        elseif (teno) then
                            tau = abs(beta2 - beta0)
                            alpha0 = (1._wp + tau/beta0)**6._wp
                            alpha1 = (1._wp + tau/beta1)**6._wp
                            alpha2 = (1._wp + tau/beta2)**6._wp
                            sum_alpha = alpha0 + alpha1 + alpha2
                            omega0 = alpha0/sum_alpha
                            omega1 = alpha1/sum_alpha
                            omega2 = alpha2/sum_alpha

                            if (omega0 < teno_CT) then
                                delta0 = 0._wp
                            else
                                delta0 = 1._wp
                            end if
                            if (omega1 < teno_CT) then
                                delta1 = 0._wp
                            else
                                delta1 = 1._wp
                            end if
                            if (omega2 < teno_CT) then
                                delta2 = 0._wp
                            else
                                delta2 = 1._wp
                            end if
                            alpha0 = delta0*d_cbL_y(0, j)
                            alpha1 = delta1*d_cbL_y(1, j)
                            alpha2 = delta2*d_cbL_y(2, j)
                        end if

                        sum_alpha = alpha0 + alpha1 + alpha2
                        omega0 = alpha0/sum_alpha
                        omega1 = alpha1/sum_alpha
                        omega2 = alpha2/sum_alpha
                        vL(j, k, l, i) = omega0*poly0 + omega1*poly1 + omega2*poly2

                        poly0 = vj + poly_coef_cbR_y(j, 0, 0)*dvd_p1 + poly_coef_cbR_y(j, 0, 1)*dvd_0
                        poly1 = vj + poly_coef_cbR_y(j, 1, 0)*dvd_0 + poly_coef_cbR_y(j, 1, 1)*dvd_m1
                        poly2 = vj + poly_coef_cbR_y(j, 2, 0)*dvd_m1 + poly_coef_cbR_y(j, 2, 1)*dvd_m2

                        if (wenojs) then
                            alpha0 = d_cbR_y(0, j)/(beta0**2._wp)
                            alpha1 = d_cbR_y(1, j)/(beta1**2._wp)
                            alpha2 = d_cbR_y(2, j)/(beta2**2._wp)
                        elseif (mapped_weno) then
                            alpha0 = d_cbR_y(0, j)/(beta0**2._wp)
                            alpha1 = d_cbR_y(1, j)/(beta1**2._wp)
                            alpha2 = d_cbR_y(2, j)/(beta2**2._wp)
                            sum_alpha = alpha0 + alpha1 + alpha2
                            omega0 = alpha0/sum_alpha
                            omega1 = alpha1/sum_alpha
                            omega2 = alpha2/sum_alpha
                            alpha0 = (d_cbR_y(0, j)*(1._wp + d_cbR_y(0, j) - 3._wp*omega0) + omega0**2._wp) &
                                   *(omega0/(d_cbR_y(0, j)**2._wp + omega0*(1._wp - 2._wp*d_cbR_y(0, j))))
                            alpha1 = (d_cbR_y(1, j)*(1._wp + d_cbR_y(1, j) - 3._wp*omega1) + omega1**2._wp) &
                                   *(omega1/(d_cbR_y(1, j)**2._wp + omega1*(1._wp - 2._wp*d_cbR_y(1, j))))
                            alpha2 = (d_cbR_y(2, j)*(1._wp + d_cbR_y(2, j) - 3._wp*omega2) + omega2**2._wp) &
                                   *(omega2/(d_cbR_y(2, j)**2._wp + omega2*(1._wp - 2._wp*d_cbR_y(2, j))))
                        elseif (wenoz) then
                            alpha0 = d_cbR_y(0, j)*(1._wp + (tau/beta0))
                            alpha1 = d_cbR_y(1, j)*(1._wp + (tau/beta1))
                            alpha2 = d_cbR_y(2, j)*(1._wp + (tau/beta2))
                        elseif (teno) then
                            alpha0 = delta0*d_cbR_y(0, j)
                            alpha1 = delta1*d_cbR_y(1, j)
                            alpha2 = delta2*d_cbR_y(2, j)
                        end if

                        sum_alpha = alpha0 + alpha1 + alpha2
                        omega0 = alpha0/sum_alpha
                        omega1 = alpha1/sum_alpha
                        omega2 = alpha2/sum_alpha
                        vR(j, k, l, i) = omega0*poly0 + omega1*poly1 + omega2*poly2
                    end do
                end do
            end do
        end do
        !$omp end target teams distribute parallel do

    end subroutine s_weno5_kernel_y

    subroutine s_weno5_kernel_z(vL, vR)

        real(wp), intent(inout), target :: vL(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR(is1_beg:, is2_beg:, is3_beg:, 1:)

        real(wp) :: dvd_m2, dvd_m1, dvd_0, dvd_p1
        real(wp) :: poly0, poly1, poly2
        real(wp) :: beta0, beta1, beta2
        real(wp) :: alpha0, alpha1, alpha2
        real(wp) :: omega0, omega1, omega2
        real(wp) :: delta0, delta1, delta2
        real(wp) :: tau, sum_alpha, vj, vjm2, vjm1, vjp1, vjp2
        integer :: i, j, k, l

        !$omp target teams distribute parallel do collapse(3) &
        !$omp& private(dvd_m2, dvd_m1, dvd_0, dvd_p1, poly0, poly1, poly2, beta0, beta1, beta2, &
        !$omp&         alpha0, alpha1, alpha2, omega0, omega1, omega2, tau, delta0, delta1, delta2, &
        !$omp&         sum_alpha, vj, vjm2, vjm1, vjp1, vjp2)
        do l = is3_beg, is3_end
            do k = is2_beg, is2_end
                do j = core1_beg, core1_end
                    do i = 1, v_size
                        alpha0 = 0._wp
                        alpha1 = 0._wp
                        alpha2 = 0._wp
                        omega0 = 0._wp
                        omega1 = 0._wp
                        omega2 = 0._wp
                        delta0 = 0._wp
                        delta1 = 0._wp
                        delta2 = 0._wp

                        vj = v_rs_ws_z(j, k, l, i)
                        vjp2 = v_rs_ws_z(j + 2, k, l, i)
                        vjp1 = v_rs_ws_z(j + 1, k, l, i)
                        vjm1 = v_rs_ws_z(j - 1, k, l, i)
                        vjm2 = v_rs_ws_z(j - 2, k, l, i)

                        dvd_p1 = vjp2 - vjp1
                        dvd_0  = vjp1 - vj
                        dvd_m1 = vj - vjm1
                        dvd_m2 = vjm1 - vjm2

                        poly0 = vj + poly_coef_cbL_z(j, 0, 0)*dvd_p1 + poly_coef_cbL_z(j, 0, 1)*dvd_0
                        poly1 = vj + poly_coef_cbL_z(j, 1, 0)*dvd_0 + poly_coef_cbL_z(j, 1, 1)*dvd_m1
                        poly2 = vj + poly_coef_cbL_z(j, 2, 0)*dvd_m1 + poly_coef_cbL_z(j, 2, 1)*dvd_m2

                        beta0 = beta_coef_z(j, 0, 0)*dvd_p1*dvd_p1 + beta_coef_z(j, 0, 1)*dvd_p1*dvd_0 &
                              + beta_coef_z(j, 0, 2)*dvd_0*dvd_0 + weno_eps
                        beta1 = beta_coef_z(j, 1, 0)*dvd_0*dvd_0 + beta_coef_z(j, 1, 1)*dvd_0*dvd_m1 &
                              + beta_coef_z(j, 1, 2)*dvd_m1*dvd_m1 + weno_eps
                        beta2 = beta_coef_z(j, 2, 0)*dvd_m1*dvd_m1 + beta_coef_z(j, 2, 1)*dvd_m1*dvd_m2 &
                              + beta_coef_z(j, 2, 2)*dvd_m2*dvd_m2 + weno_eps

                        if (wenojs) then
                            alpha0 = d_cbL_z(0, j)/(beta0**2._wp)
                            alpha1 = d_cbL_z(1, j)/(beta1**2._wp)
                            alpha2 = d_cbL_z(2, j)/(beta2**2._wp)
                        elseif (mapped_weno) then
                            alpha0 = d_cbL_z(0, j)/(beta0**2._wp)
                            alpha1 = d_cbL_z(1, j)/(beta1**2._wp)
                            alpha2 = d_cbL_z(2, j)/(beta2**2._wp)
                            sum_alpha = alpha0 + alpha1 + alpha2
                            omega0 = alpha0/sum_alpha
                            omega1 = alpha1/sum_alpha
                            omega2 = alpha2/sum_alpha
                            alpha0 = (d_cbL_z(0, j)*(1._wp + d_cbL_z(0, j) - 3._wp*omega0) + omega0**2._wp) &
                                   *(omega0/(d_cbL_z(0, j)**2._wp + omega0*(1._wp - 2._wp*d_cbL_z(0, j))))
                            alpha1 = (d_cbL_z(1, j)*(1._wp + d_cbL_z(1, j) - 3._wp*omega1) + omega1**2._wp) &
                                   *(omega1/(d_cbL_z(1, j)**2._wp + omega1*(1._wp - 2._wp*d_cbL_z(1, j))))
                            alpha2 = (d_cbL_z(2, j)*(1._wp + d_cbL_z(2, j) - 3._wp*omega2) + omega2**2._wp) &
                                   *(omega2/(d_cbL_z(2, j)**2._wp + omega2*(1._wp - 2._wp*d_cbL_z(2, j))))
                        elseif (wenoz) then
                            tau = abs(beta2 - beta0)
                            alpha0 = d_cbL_z(0, j)*(1._wp + (tau/beta0))
                            alpha1 = d_cbL_z(1, j)*(1._wp + (tau/beta1))
                            alpha2 = d_cbL_z(2, j)*(1._wp + (tau/beta2))
                        elseif (teno) then
                            tau = abs(beta2 - beta0)
                            alpha0 = (1._wp + tau/beta0)**6._wp
                            alpha1 = (1._wp + tau/beta1)**6._wp
                            alpha2 = (1._wp + tau/beta2)**6._wp
                            sum_alpha = alpha0 + alpha1 + alpha2
                            omega0 = alpha0/sum_alpha
                            omega1 = alpha1/sum_alpha
                            omega2 = alpha2/sum_alpha

                            if (omega0 < teno_CT) then
                                delta0 = 0._wp
                            else
                                delta0 = 1._wp
                            end if
                            if (omega1 < teno_CT) then
                                delta1 = 0._wp
                            else
                                delta1 = 1._wp
                            end if
                            if (omega2 < teno_CT) then
                                delta2 = 0._wp
                            else
                                delta2 = 1._wp
                            end if
                            alpha0 = delta0*d_cbL_z(0, j)
                            alpha1 = delta1*d_cbL_z(1, j)
                            alpha2 = delta2*d_cbL_z(2, j)
                        end if

                        sum_alpha = alpha0 + alpha1 + alpha2
                        omega0 = alpha0/sum_alpha
                        omega1 = alpha1/sum_alpha
                        omega2 = alpha2/sum_alpha
                        vL(j, k, l, i) = omega0*poly0 + omega1*poly1 + omega2*poly2

                        poly0 = vj + poly_coef_cbR_z(j, 0, 0)*dvd_p1 + poly_coef_cbR_z(j, 0, 1)*dvd_0
                        poly1 = vj + poly_coef_cbR_z(j, 1, 0)*dvd_0 + poly_coef_cbR_z(j, 1, 1)*dvd_m1
                        poly2 = vj + poly_coef_cbR_z(j, 2, 0)*dvd_m1 + poly_coef_cbR_z(j, 2, 1)*dvd_m2

                        if (wenojs) then
                            alpha0 = d_cbR_z(0, j)/(beta0**2._wp)
                            alpha1 = d_cbR_z(1, j)/(beta1**2._wp)
                            alpha2 = d_cbR_z(2, j)/(beta2**2._wp)
                        elseif (mapped_weno) then
                            alpha0 = d_cbR_z(0, j)/(beta0**2._wp)
                            alpha1 = d_cbR_z(1, j)/(beta1**2._wp)
                            alpha2 = d_cbR_z(2, j)/(beta2**2._wp)
                            sum_alpha = alpha0 + alpha1 + alpha2
                            omega0 = alpha0/sum_alpha
                            omega1 = alpha1/sum_alpha
                            omega2 = alpha2/sum_alpha
                            alpha0 = (d_cbR_z(0, j)*(1._wp + d_cbR_z(0, j) - 3._wp*omega0) + omega0**2._wp) &
                                   *(omega0/(d_cbR_z(0, j)**2._wp + omega0*(1._wp - 2._wp*d_cbR_z(0, j))))
                            alpha1 = (d_cbR_z(1, j)*(1._wp + d_cbR_z(1, j) - 3._wp*omega1) + omega1**2._wp) &
                                   *(omega1/(d_cbR_z(1, j)**2._wp + omega1*(1._wp - 2._wp*d_cbR_z(1, j))))
                            alpha2 = (d_cbR_z(2, j)*(1._wp + d_cbR_z(2, j) - 3._wp*omega2) + omega2**2._wp) &
                                   *(omega2/(d_cbR_z(2, j)**2._wp + omega2*(1._wp - 2._wp*d_cbR_z(2, j))))
                        elseif (wenoz) then
                            alpha0 = d_cbR_z(0, j)*(1._wp + (tau/beta0))
                            alpha1 = d_cbR_z(1, j)*(1._wp + (tau/beta1))
                            alpha2 = d_cbR_z(2, j)*(1._wp + (tau/beta2))
                        elseif (teno) then
                            alpha0 = delta0*d_cbR_z(0, j)
                            alpha1 = delta1*d_cbR_z(1, j)
                            alpha2 = delta2*d_cbR_z(2, j)
                        end if

                        sum_alpha = alpha0 + alpha1 + alpha2
                        omega0 = alpha0/sum_alpha
                        omega1 = alpha1/sum_alpha
                        omega2 = alpha2/sum_alpha
                        vR(j, k, l, i) = omega0*poly0 + omega1*poly1 + omega2*poly2
                    end do
                end do
            end do
        end do
        !$omp end target teams distribute parallel do

    end subroutine s_weno5_kernel_z

    subroutine s_weno5_left_kernel(v_rs_ws, polyL, dL, beta_coef, vL)

        real(wp), intent(in), target :: v_rs_ws(is1_beg-weno_polyn_max:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(in), target :: polyL(core1_beg:, 0:, 0:)
        real(wp), intent(in), target :: dL(0:, core1_beg:)
        real(wp), intent(in), target :: beta_coef(core1_beg:, 0:, 0:)
        real(wp), intent(inout), target :: vL(is1_beg:, is2_beg:, is3_beg:, 1:)

        real(wp) :: dvd_m2, dvd_m1, dvd_0, dvd_p1
        real(wp) :: poly0, poly1, poly2
        real(wp) :: beta0, beta1, beta2
        real(wp) :: alpha0, alpha1, alpha2
        real(wp) :: omega0, omega1, omega2
        real(wp) :: delta0, delta1, delta2
        real(wp) :: tau, sum_alpha, vj, vjm2, vjm1, vjp1, vjp2
        integer, parameter :: weno_num_stencils = 2
        integer :: i, j, k, l

        !$omp target teams distribute parallel do collapse(3) &
        !$omp& private(dvd_m2, dvd_m1, dvd_0, dvd_p1, poly0, poly1, poly2, beta0, beta1, beta2, &
        !$omp&         alpha0, alpha1, alpha2, omega0, omega1, omega2, tau, delta0, delta1, delta2, &
        !$omp&         sum_alpha, vj, vjm2, vjm1, vjp1, vjp2)
        do l = is3_beg, is3_end
            do k = is2_beg, is2_end
                do j = core1_beg, core1_end
                    do i = 1, v_size
                        alpha0 = 0._wp
                        alpha1 = 0._wp
                        alpha2 = 0._wp
                        omega0 = 0._wp
                        omega1 = 0._wp
                        omega2 = 0._wp
                        delta0 = 0._wp
                        delta1 = 0._wp
                        delta2 = 0._wp

                        vj = v_rs_ws(j, k, l, i)
                        vjp2 = v_rs_ws(j + 2, k, l, i)
                        vjp1 = v_rs_ws(j + 1, k, l, i)
                        vjm1 = v_rs_ws(j - 1, k, l, i)
                        vjm2 = v_rs_ws(j - 2, k, l, i)

                        dvd_p1 = vjp2 - vjp1
                        dvd_0  = vjp1 - vj
                        dvd_m1 = vj - vjm1
                        dvd_m2 = vjm1 - vjm2

                        poly0 = vj + polyL(j, 0, 0)*dvd_p1 + polyL(j, 0, 1)*dvd_0
                        poly1 = vj + polyL(j, 1, 0)*dvd_0 + polyL(j, 1, 1)*dvd_m1
                        poly2 = vj + polyL(j, 2, 0)*dvd_m1 + polyL(j, 2, 1)*dvd_m2

                        beta0 = beta_coef(j, 0, 0)*dvd_p1*dvd_p1 + beta_coef(j, 0, 1)*dvd_p1*dvd_0 &
                              + beta_coef(j, 0, 2)*dvd_0*dvd_0 + weno_eps
                        beta1 = beta_coef(j, 1, 0)*dvd_0*dvd_0 + beta_coef(j, 1, 1)*dvd_0*dvd_m1 &
                              + beta_coef(j, 1, 2)*dvd_m1*dvd_m1 + weno_eps
                        beta2 = beta_coef(j, 2, 0)*dvd_m1*dvd_m1 + beta_coef(j, 2, 1)*dvd_m1*dvd_m2 &
                              + beta_coef(j, 2, 2)*dvd_m2*dvd_m2 + weno_eps

                        if (wenojs) then
                            alpha0 = dL(0, j)/(beta0**2._wp)
                            alpha1 = dL(1, j)/(beta1**2._wp)
                            alpha2 = dL(2, j)/(beta2**2._wp)

                        elseif (mapped_weno) then
                            alpha0 = dL(0, j)/(beta0**2._wp)
                            alpha1 = dL(1, j)/(beta1**2._wp)
                            alpha2 = dL(2, j)/(beta2**2._wp)
                            sum_alpha = alpha0 + alpha1 + alpha2
                            omega0 = alpha0/sum_alpha
                            omega1 = alpha1/sum_alpha
                            omega2 = alpha2/sum_alpha
                            alpha0 = (dL(0, j)*(1._wp + dL(0, j) - 3._wp*omega0) + omega0**2._wp) &
                                   *(omega0/(dL(0, j)**2._wp + omega0*(1._wp - 2._wp*dL(0, j))))
                            alpha1 = (dL(1, j)*(1._wp + dL(1, j) - 3._wp*omega1) + omega1**2._wp) &
                                   *(omega1/(dL(1, j)**2._wp + omega1*(1._wp - 2._wp*dL(1, j))))
                            alpha2 = (dL(2, j)*(1._wp + dL(2, j) - 3._wp*omega2) + omega2**2._wp) &
                                   *(omega2/(dL(2, j)**2._wp + omega2*(1._wp - 2._wp*dL(2, j))))

                        elseif (wenoz) then
                            tau = abs(beta2 - beta0)
                            alpha0 = dL(0, j)*(1._wp + (tau/beta0))
                            alpha1 = dL(1, j)*(1._wp + (tau/beta1))
                            alpha2 = dL(2, j)*(1._wp + (tau/beta2))

                        elseif (teno) then
                            tau = abs(beta2 - beta0)
                            alpha0 = (1._wp + tau/beta0)**6._wp
                            alpha1 = (1._wp + tau/beta1)**6._wp
                            alpha2 = (1._wp + tau/beta2)**6._wp
                            sum_alpha = alpha0 + alpha1 + alpha2
                            omega0 = alpha0/sum_alpha
                            omega1 = alpha1/sum_alpha
                            omega2 = alpha2/sum_alpha

                            if (omega0 < teno_CT) then
                                delta0 = 0._wp
                            else
                                delta0 = 1._wp
                            end if
                            if (omega1 < teno_CT) then
                                delta1 = 0._wp
                            else
                                delta1 = 1._wp
                            end if
                            if (omega2 < teno_CT) then
                                delta2 = 0._wp
                            else
                                delta2 = 1._wp
                            end if
                            alpha0 = delta0*dL(0, j)
                            alpha1 = delta1*dL(1, j)
                            alpha2 = delta2*dL(2, j)
                        end if

                        sum_alpha = alpha0 + alpha1 + alpha2
                        omega0 = alpha0/sum_alpha
                        omega1 = alpha1/sum_alpha
                        omega2 = alpha2/sum_alpha
                        vL(j, k, l, i) = omega0*poly0 + omega1*poly1 + omega2*poly2
                    end do
                end do
            end do
        end do
        !$omp end target teams distribute parallel do

    end subroutine s_weno5_left_kernel

    subroutine s_weno5_right_kernel(v_rs_ws, polyR, dR, beta_coef, vR)

        real(wp), intent(in), target :: v_rs_ws(is1_beg-weno_polyn_max:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(in), target :: polyR(core1_beg:, 0:, 0:)
        real(wp), intent(in), target :: dR(0:, core1_beg:)
        real(wp), intent(in), target :: beta_coef(core1_beg:, 0:, 0:)
        real(wp), intent(inout), target :: vR(is1_beg:, is2_beg:, is3_beg:, 1:)

        real(wp) :: dvd_m2, dvd_m1, dvd_0, dvd_p1
        real(wp) :: poly0, poly1, poly2
        real(wp) :: beta0, beta1, beta2
        real(wp) :: alpha0, alpha1, alpha2
        real(wp) :: omega0, omega1, omega2
        real(wp) :: delta0, delta1, delta2
        real(wp) :: tau, sum_alpha, vj, vjm2, vjm1, vjp1, vjp2
        integer, parameter :: weno_num_stencils = 2
        integer :: i, j, k, l

        !$omp target teams distribute parallel do collapse(3) &
        !$omp& private(dvd_m2, dvd_m1, dvd_0, dvd_p1, poly0, poly1, poly2, beta0, beta1, beta2, &
        !$omp&         alpha0, alpha1, alpha2, omega0, omega1, omega2, tau, delta0, delta1, delta2, &
        !$omp&         sum_alpha, vj, vjm2, vjm1, vjp1, vjp2)
        do l = is3_beg, is3_end
            do k = is2_beg, is2_end
                do j = core1_beg, core1_end
                    do i = 1, v_size
                        alpha0 = 0._wp
                        alpha1 = 0._wp
                        alpha2 = 0._wp
                        omega0 = 0._wp
                        omega1 = 0._wp
                        omega2 = 0._wp
                        delta0 = 0._wp
                        delta1 = 0._wp
                        delta2 = 0._wp

                        vj = v_rs_ws(j, k, l, i)
                        vjp2 = v_rs_ws(j + 2, k, l, i)
                        vjp1 = v_rs_ws(j + 1, k, l, i)
                        vjm1 = v_rs_ws(j - 1, k, l, i)
                        vjm2 = v_rs_ws(j - 2, k, l, i)

                        dvd_p1 = vjp2 - vjp1
                        dvd_0  = vjp1 - vj
                        dvd_m1 = vj - vjm1
                        dvd_m2 = vjm1 - vjm2

                        poly0 = vj + polyR(j, 0, 0)*dvd_p1 + polyR(j, 0, 1)*dvd_0
                        poly1 = vj + polyR(j, 1, 0)*dvd_0 + polyR(j, 1, 1)*dvd_m1
                        poly2 = vj + polyR(j, 2, 0)*dvd_m1 + polyR(j, 2, 1)*dvd_m2

                        beta0 = beta_coef(j, 0, 0)*dvd_p1*dvd_p1 + beta_coef(j, 0, 1)*dvd_p1*dvd_0 &
                              + beta_coef(j, 0, 2)*dvd_0*dvd_0 + weno_eps
                        beta1 = beta_coef(j, 1, 0)*dvd_0*dvd_0 + beta_coef(j, 1, 1)*dvd_0*dvd_m1 &
                              + beta_coef(j, 1, 2)*dvd_m1*dvd_m1 + weno_eps
                        beta2 = beta_coef(j, 2, 0)*dvd_m1*dvd_m1 + beta_coef(j, 2, 1)*dvd_m1*dvd_m2 &
                              + beta_coef(j, 2, 2)*dvd_m2*dvd_m2 + weno_eps

                        if (wenojs) then
                            alpha0 = dR(0, j)/(beta0**2._wp)
                            alpha1 = dR(1, j)/(beta1**2._wp)
                            alpha2 = dR(2, j)/(beta2**2._wp)

                        elseif (mapped_weno) then
                            alpha0 = dR(0, j)/(beta0**2._wp)
                            alpha1 = dR(1, j)/(beta1**2._wp)
                            alpha2 = dR(2, j)/(beta2**2._wp)
                            sum_alpha = alpha0 + alpha1 + alpha2
                            omega0 = alpha0/sum_alpha
                            omega1 = alpha1/sum_alpha
                            omega2 = alpha2/sum_alpha
                            alpha0 = (dR(0, j)*(1._wp + dR(0, j) - 3._wp*omega0) + omega0**2._wp) &
                                   *(omega0/(dR(0, j)**2._wp + omega0*(1._wp - 2._wp*dR(0, j))))
                            alpha1 = (dR(1, j)*(1._wp + dR(1, j) - 3._wp*omega1) + omega1**2._wp) &
                                   *(omega1/(dR(1, j)**2._wp + omega1*(1._wp - 2._wp*dR(1, j))))
                            alpha2 = (dR(2, j)*(1._wp + dR(2, j) - 3._wp*omega2) + omega2**2._wp) &
                                   *(omega2/(dR(2, j)**2._wp + omega2*(1._wp - 2._wp*dR(2, j))))

                        elseif (wenoz) then
                            tau = abs(beta2 - beta0)
                            alpha0 = dR(0, j)*(1._wp + (tau/beta0))
                            alpha1 = dR(1, j)*(1._wp + (tau/beta1))
                            alpha2 = dR(2, j)*(1._wp + (tau/beta2))

                        elseif (teno) then
                            tau = abs(beta2 - beta0)
                            alpha0 = (1._wp + tau/beta0)**6._wp
                            alpha1 = (1._wp + tau/beta1)**6._wp
                            alpha2 = (1._wp + tau/beta2)**6._wp
                            sum_alpha = alpha0 + alpha1 + alpha2
                            omega0 = alpha0/sum_alpha
                            omega1 = alpha1/sum_alpha
                            omega2 = alpha2/sum_alpha

                            if (omega0 < teno_CT) then
                                delta0 = 0._wp
                            else
                                delta0 = 1._wp
                            end if
                            if (omega1 < teno_CT) then
                                delta1 = 0._wp
                            else
                                delta1 = 1._wp
                            end if
                            if (omega2 < teno_CT) then
                                delta2 = 0._wp
                            else
                                delta2 = 1._wp
                            end if
                            alpha0 = delta0*dR(0, j)
                            alpha1 = delta1*dR(1, j)
                            alpha2 = delta2*dR(2, j)
                        end if

                        sum_alpha = alpha0 + alpha1 + alpha2
                        omega0 = alpha0/sum_alpha
                        omega1 = alpha1/sum_alpha
                        omega2 = alpha2/sum_alpha
                        vR(j, k, l, i) = omega0*poly0 + omega1*poly1 + omega2*poly2
                    end do
                end do
            end do
        end do
        !$omp end target teams distribute parallel do

    end subroutine s_weno5_right_kernel

    subroutine s_weno5_kernel_dispatch(v_rs_ws, polyL, polyR, dL, dR, beta_coef, vL, vR)

        real(wp), intent(in), target :: v_rs_ws(is1_beg-weno_polyn_max:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(in), target :: polyL(core1_beg:, 0:, 0:)
        real(wp), intent(in), target :: polyR(core1_beg:, 0:, 0:)
        real(wp), intent(in), target :: dL(0:, core1_beg:)
        real(wp), intent(in), target :: dR(0:, core1_beg:)
        real(wp), intent(in), target :: beta_coef(core1_beg:, 0:, 0:)
        real(wp), intent(inout), target :: vL(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR(is1_beg:, is2_beg:, is3_beg:, 1:)

        if (weno5_split_kernels) then
            call s_weno5_left_kernel(v_rs_ws, polyL, dL, beta_coef, vL)
            call s_weno5_right_kernel(v_rs_ws, polyR, dR, beta_coef, vR)
        else
            call s_weno5_kernel(v_rs_ws, polyL, polyR, dL, dR, beta_coef, vL, vR)
        end if

    end subroutine s_weno5_kernel_dispatch

    subroutine s_weno7_kernel(v_rs_ws, polyL, polyR, dL, dR, beta_coef, vL, vR)

        real(wp), intent(in), target :: v_rs_ws(is1_beg-weno_polyn_max:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(in), target :: polyL(core1_beg:, 0:, 0:)
        real(wp), intent(in), target :: polyR(core1_beg:, 0:, 0:)
        real(wp), intent(in), target :: dL(0:, core1_beg:)
        real(wp), intent(in), target :: dR(0:, core1_beg:)
        real(wp), intent(in), target :: beta_coef(core1_beg:, 0:, 0:)
        real(wp), intent(inout), target :: vL(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR(is1_beg:, is2_beg:, is3_beg:, 1:)

        real(wp) :: dvd(-3:2), poly(0:4), alpha(0:4), omega(0:4), beta(0:4), delta(0:4)
        real(wp) :: v(-3:3), tau
        integer, parameter :: weno_num_stencils = 3
        integer :: i, j, k, l, q

        ! Source: src/simulation/m_weno.fpp:947 ($:GPU_PARALLEL_LOOP)
        ! Inner sequential GPU_LOOP sites in source: 951, 1075, 1087, 1147, 1153
        !$omp target teams distribute parallel do collapse(3) &
        !$omp& private(poly, beta, alpha, omega, tau, delta, dvd, v, q)
        do l = is3_beg, is3_end
            do k = is2_beg, is2_end
                do j = core1_beg, core1_end
                    do i = 1, v_size
                        alpha(:) = 0._wp
                        omega(:) = 0._wp
                        delta(:) = 0._wp
                        beta(:) = weno_eps

                        if (teno) v = v_rs_ws(j - 3:j + 3, k, l, i)

                        dvd(2) = v_rs_ws(j + 3, k, l, i) - v_rs_ws(j + 2, k, l, i)
                        dvd(1) = v_rs_ws(j + 2, k, l, i) - v_rs_ws(j + 1, k, l, i)
                        dvd(0) = v_rs_ws(j + 1, k, l, i) - v_rs_ws(j, k, l, i)
                        dvd(-1) = v_rs_ws(j, k, l, i) - v_rs_ws(j - 1, k, l, i)
                        dvd(-2) = v_rs_ws(j - 1, k, l, i) - v_rs_ws(j - 2, k, l, i)
                        dvd(-3) = v_rs_ws(j - 2, k, l, i) - v_rs_ws(j - 3, k, l, i)

                        poly(3) = v_rs_ws(j, k, l, i) + polyL(j, 0, 0)*dvd(2) + polyL(j, 0, 1)*dvd(1) + polyL(j, 0, 2)*dvd(0)
                        poly(2) = v_rs_ws(j, k, l, i) + polyL(j, 1, 0)*dvd(1) + polyL(j, 1, 1)*dvd(0) + polyL(j, 1, 2)*dvd(-1)
                        poly(1) = v_rs_ws(j, k, l, i) + polyL(j, 2, 0)*dvd(0) + polyL(j, 2, 1)*dvd(-1) + polyL(j, 2, 2)*dvd(-2)
                        poly(0) = v_rs_ws(j, k, l, i) + polyL(j, 3, 0)*dvd(-1) + polyL(j, 3, 1)*dvd(-2) + polyL(j, 3, 2)*dvd(-3)

                        beta(3) = beta_coef(j, 0, 0)*dvd(2)*dvd(2) + beta_coef(j, 0, 1)*dvd(2)*dvd(1) &
                                  + beta_coef(j, 0, 2)*dvd(2)*dvd(0) + beta_coef(j, 0, 3)*dvd(1)*dvd(1) &
                                  + beta_coef(j, 0, 4)*dvd(1)*dvd(0) + beta_coef(j, 0, 5)*dvd(0)*dvd(0) + weno_eps
                        beta(2) = beta_coef(j, 1, 0)*dvd(1)*dvd(1) + beta_coef(j, 1, 1)*dvd(1)*dvd(0) &
                                  + beta_coef(j, 1, 2)*dvd(1)*dvd(-1) + beta_coef(j, 1, 3)*dvd(0)*dvd(0) &
                                  + beta_coef(j, 1, 4)*dvd(0)*dvd(-1) + beta_coef(j, 1, 5)*dvd(-1)*dvd(-1) + weno_eps
                        beta(1) = beta_coef(j, 2, 0)*dvd(0)*dvd(0) + beta_coef(j, 2, 1)*dvd(0)*dvd(-1) &
                                  + beta_coef(j, 2, 2)*dvd(0)*dvd(-2) + beta_coef(j, 2, 3)*dvd(-1)*dvd(-1) &
                                  + beta_coef(j, 2, 4)*dvd(-1)*dvd(-2) + beta_coef(j, 2, 5)*dvd(-2)*dvd(-2) + weno_eps
                        beta(0) = beta_coef(j, 3, 0)*dvd(-1)*dvd(-1) + beta_coef(j, 3, 1)*dvd(-1)*dvd(-2) &
                                  + beta_coef(j, 3, 2)*dvd(-1)*dvd(-3) + beta_coef(j, 3, 3)*dvd(-2)*dvd(-2) &
                                  + beta_coef(j, 3, 4)*dvd(-2)*dvd(-3) + beta_coef(j, 3, 5)*dvd(-3)*dvd(-3) + weno_eps

                        if (wenojs) then
                            alpha(0:weno_num_stencils) = dL(0:weno_num_stencils, j)/(beta(0:weno_num_stencils)**2._wp)

                        elseif (mapped_weno) then
                            alpha(0:weno_num_stencils) = dL(0:weno_num_stencils, j)/(beta(0:weno_num_stencils)**2._wp)
                            omega = alpha/sum(alpha)
                            alpha(0:weno_num_stencils) = (dL(0:weno_num_stencils, j)*(1._wp + dL(0:weno_num_stencils, j) - 3._wp*omega(0:weno_num_stencils)) + omega(0:weno_num_stencils)**2._wp) &
                                                         *(omega(0:weno_num_stencils)/(dL(0:weno_num_stencils, j)**2._wp + omega(0:weno_num_stencils)*(1._wp - 2._wp*dL(0:weno_num_stencils, j))))

                        elseif (wenoz) then
                            tau = abs(beta(3) - beta(0))
                            do q = 0, weno_num_stencils
                                alpha(q) = dL(q, j)*(1._wp + (tau/beta(q))**wenoz_q)
                            end do

                        elseif (teno) then
                            tau = abs(beta(4) - beta(3))
                            alpha = 1._wp + tau/beta
                            alpha = (alpha**3._wp)**2._wp
                            omega = alpha/sum(alpha)

                            do q = 0, weno_num_stencils
                                if (omega(q) < teno_CT) then
                                    delta(q) = 0._wp
                                else
                                    delta(q) = 1._wp
                                end if
                                alpha(q) = delta(q)*dL(q, j)
                            end do
                        end if

                        omega = alpha/sum(alpha)
                        vL(j, k, l, i) = omega(0)*poly(0) + omega(1)*poly(1) + omega(2)*poly(2) + omega(3)*poly(3)

                        poly(3) = v_rs_ws(j, k, l, i) + polyR(j, 0, 0)*dvd(2) + polyR(j, 0, 1)*dvd(1) + polyR(j, 0, 2)*dvd(0)
                        poly(2) = v_rs_ws(j, k, l, i) + polyR(j, 1, 0)*dvd(1) + polyR(j, 1, 1)*dvd(0) + polyR(j, 1, 2)*dvd(-1)
                        poly(1) = v_rs_ws(j, k, l, i) + polyR(j, 2, 0)*dvd(0) + polyR(j, 2, 1)*dvd(-1) + polyR(j, 2, 2)*dvd(-2)
                        poly(0) = v_rs_ws(j, k, l, i) + polyR(j, 3, 0)*dvd(-1) + polyR(j, 3, 1)*dvd(-2) + polyR(j, 3, 2)*dvd(-3)

                        if (wenojs) then
                            alpha(0:weno_num_stencils) = dR(0:weno_num_stencils, j)/(beta(0:weno_num_stencils)**2._wp)

                        elseif (mapped_weno) then
                            alpha(0:weno_num_stencils) = dR(0:weno_num_stencils, j)/(beta(0:weno_num_stencils)**2._wp)
                            omega = alpha/sum(alpha)
                            alpha(0:weno_num_stencils) = (dR(0:weno_num_stencils, j)*(1._wp + dR(0:weno_num_stencils, j) - 3._wp*omega(0:weno_num_stencils)) + omega(0:weno_num_stencils)**2._wp) &
                                                         *(omega(0:weno_num_stencils)/(dR(0:weno_num_stencils, j)**2._wp + omega(0:weno_num_stencils)*(1._wp - 2._wp*dR(0:weno_num_stencils, j))))

                        elseif (wenoz) then
                            do q = 0, weno_num_stencils
                                alpha(q) = dR(q, j)*(1._wp + (tau/beta(q))**wenoz_q)
                            end do

                        elseif (teno) then
                            do q = 0, weno_num_stencils
                                alpha(q) = delta(q)*dR(q, j)
                            end do
                        end if

                        omega = alpha/sum(alpha)
                        vR(j, k, l, i) = omega(0)*poly(0) + omega(1)*poly(1) + omega(2)*poly(2) + omega(3)*poly(3)
                    end do
                end do
            end do
        end do
        !$omp end target teams distribute parallel do

    end subroutine s_weno7_kernel

    subroutine s_preserve_monotonicity(v_rs_ws, vL, vR)

        real(wp), intent(in), target :: v_rs_ws(is1_beg-weno_polyn_max:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vL(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR(is1_beg:, is2_beg:, is3_beg:, 1:)

        integer :: i, j, k, l
        real(wp) :: d(-1:1)
        real(wp) :: d_MD, d_LC
        real(wp) :: vL_UL, vR_UL, vL_MD, vR_MD, vL_LC, vR_LC
        real(wp) :: vL_min, vR_min, vL_max, vR_max
        real(wp), parameter :: alpha_mp = 2.0_wp
        real(wp), parameter :: beta_mp  = 4.0_wp/3.0_wp

        ! Source: src/simulation/m_weno.fpp:1304 ($:GPU_PARALLEL_LOOP, private='[d]')
        !$omp target teams distribute parallel do collapse(4) private(d)
        do l = is3_beg, is3_end
            do k = is2_beg, is2_end
                do j = core1_beg, core1_end
                    do i = 1, v_size
                        d(-1) = v_rs_ws(j, k, l, i) + v_rs_ws(j - 2, k, l, i) - 2.0_wp*v_rs_ws(j - 1, k, l, i)
                        d(0)  = v_rs_ws(j + 1, k, l, i) + v_rs_ws(j - 1, k, l, i) - 2.0_wp*v_rs_ws(j, k, l, i)
                        d(1)  = v_rs_ws(j + 2, k, l, i) + v_rs_ws(j, k, l, i) - 2.0_wp*v_rs_ws(j + 1, k, l, i)

                        d_MD = (sign(1._wp, 4._wp*d(-1) - d(0)) + sign(1._wp, 4._wp*d(0) - d(-1))) &
                               *abs((sign(1._wp, 4._wp*d(-1) - d(0)) + sign(1._wp, d(-1))) &
                                    *(sign(1._wp, 4._wp*d(-1) - d(0)) + sign(1._wp, d(0)))) &
                               *min(abs(4._wp*d(-1) - d(0)), abs(d(-1)), abs(4._wp*d(0) - d(-1)), abs(d(0)))/8._wp

                        d_LC = (sign(1._wp, 4._wp*d(0) - d(1)) + sign(1._wp, 4._wp*d(1) - d(0))) &
                               *abs((sign(1._wp, 4._wp*d(0) - d(1)) + sign(1._wp, d(0))) &
                                    *(sign(1._wp, 4._wp*d(0) - d(1)) + sign(1._wp, d(1)))) &
                               *min(abs(4._wp*d(0) - d(1)), abs(d(0)), abs(4._wp*d(1) - d(0)), abs(d(1)))/8._wp

                        vL_UL = v_rs_ws(j, k, l, i) - alpha_mp*(v_rs_ws(j + 1, k, l, i) - v_rs_ws(j, k, l, i))
                        vL_MD = (v_rs_ws(j, k, l, i) + v_rs_ws(j - 1, k, l, i) - d_MD)*5.e-1_wp
                        vL_LC = v_rs_ws(j, k, l, i) - (v_rs_ws(j + 1, k, l, i) - v_rs_ws(j, k, l, i))*5.e-1_wp + beta_mp*d_LC

                        vL_min = max(min(v_rs_ws(j, k, l, i), v_rs_ws(j - 1, k, l, i), vL_MD), &
                                     min(v_rs_ws(j, k, l, i), vL_UL, vL_LC))
                        vL_max = min(max(v_rs_ws(j, k, l, i), v_rs_ws(j - 1, k, l, i), vL_MD), &
                                     max(v_rs_ws(j, k, l, i), vL_UL, vL_LC))
                        vL(j, k, l, i) = vL(j, k, l, i) + (sign(5.e-1_wp, vL_min - vL(j, k, l, i)) + sign(5.e-1_wp, vL_max - vL(j, k, l, i))) &
                                         *min(abs(vL_min - vL(j, k, l, i)), abs(vL_max - vL(j, k, l, i)))

                        vR_UL = v_rs_ws(j, k, l, i) + alpha_mp*(v_rs_ws(j, k, l, i) - v_rs_ws(j - 1, k, l, i))
                        d_MD = (sign(1._wp, 4._wp*d(0) - d(1)) + sign(1._wp, 4._wp*d(1) - d(0))) &
                               *abs((sign(1._wp, 4._wp*d(0) - d(1)) + sign(1._wp, d(0))) &
                                    *(sign(1._wp, 4._wp*d(0) - d(1)) + sign(1._wp, d(1)))) &
                               *min(abs(4._wp*d(0) - d(1)), abs(d(0)), abs(4._wp*d(1) - d(0)), abs(d(1)))/8._wp

                        d_LC = (sign(1._wp, 4._wp*d(-1) - d(0)) + sign(1._wp, 4._wp*d(0) - d(-1))) &
                               *abs((sign(1._wp, 4._wp*d(-1) - d(0)) + sign(1._wp, d(-1))) &
                                    *(sign(1._wp, 4._wp*d(-1) - d(0)) + sign(1._wp, d(0)))) &
                               *min(abs(4._wp*d(-1) - d(0)), abs(d(-1)), abs(4._wp*d(0) - d(-1)), abs(d(0)))/8._wp

                        vR_MD = (v_rs_ws(j, k, l, i) + v_rs_ws(j + 1, k, l, i) - d_MD)*5.e-1_wp
                        vR_LC = v_rs_ws(j, k, l, i) + (v_rs_ws(j, k, l, i) - v_rs_ws(j - 1, k, l, i))*5.e-1_wp + beta_mp*d_LC

                        vR_min = max(min(v_rs_ws(j, k, l, i), v_rs_ws(j + 1, k, l, i), vR_MD), &
                                     min(v_rs_ws(j, k, l, i), vR_UL, vR_LC))
                        vR_max = min(max(v_rs_ws(j, k, l, i), v_rs_ws(j + 1, k, l, i), vR_MD), &
                                     max(v_rs_ws(j, k, l, i), vR_UL, vR_LC))
                        vR(j, k, l, i) = vR(j, k, l, i) + (sign(5.e-1_wp, vR_min - vR(j, k, l, i)) + sign(5.e-1_wp, vR_max - vR(j, k, l, i))) &
                                         *min(abs(vR_min - vR(j, k, l, i)), abs(vR_max - vR(j, k, l, i)))
                    end do
                end do
            end do
        end do
        !$omp end target teams distribute parallel do

    end subroutine s_preserve_monotonicity

    real(wp) function f_checksum(arr) result(chk)
        real(wp), intent(in) :: arr(:,:,:,:)
        chk = sum(abs(arr))
    end function f_checksum

    integer function f_get_env_int(name, default_value) result(value)
        character(len=*), intent(in) :: name
        integer, intent(in) :: default_value

        character(len=64) :: buffer
        integer :: length, status, parsed_value

        value = default_value
        buffer = ''

        call get_environment_variable(name, buffer, length=length, status=status)
        if (status /= 0 .or. length <= 0) return

        read(buffer(1:length), *, iostat=status) parsed_value
        if (status == 0) value = parsed_value
    end function f_get_env_int

    subroutine s_get_env_string(name, default_value, value)
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: default_value
        character(len=*), intent(out) :: value

        integer :: length, status

        value = default_value
        call get_environment_variable(name, value, length=length, status=status)
        if (status /= 0 .or. length <= 0) value = default_value
    end subroutine s_get_env_string

    subroutine s_run_all_kernels(vL_x, vR_x, vL_y, vR_y, vL_z, vR_z)
        real(wp), intent(inout), target :: vL_x(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR_x(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vL_y(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR_y(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vL_z(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR_z(is1_beg:, is2_beg:, is3_beg:, 1:)

        call s_initialize_weno_x()
        call s_initialize_weno_y()
        call s_initialize_weno_z()

        call s_order1_copy_x(vL_x, vR_x)
        call s_order1_copy_y(vL_y, vR_y)
        call s_order1_copy_z(vL_z, vR_z)

        call s_weno3_kernel(v_rs_ws_x, poly_coef_cbL_x, poly_coef_cbR_x, d_cbL_x, d_cbR_x, beta_coef_x, vL_x, vR_x)
        call s_weno3_kernel(v_rs_ws_y, poly_coef_cbL_y, poly_coef_cbR_y, d_cbL_y, d_cbR_y, beta_coef_y, vL_y, vR_y)
        call s_weno3_kernel(v_rs_ws_z, poly_coef_cbL_z, poly_coef_cbR_z, d_cbL_z, d_cbR_z, beta_coef_z, vL_z, vR_z)

        call s_weno5_kernel_dispatch(v_rs_ws_x, poly_coef_cbL_x, poly_coef_cbR_x, d_cbL_x, d_cbR_x, beta_coef_x, vL_x, vR_x)
        call s_weno5_kernel_dispatch(v_rs_ws_y, poly_coef_cbL_y, poly_coef_cbR_y, d_cbL_y, d_cbR_y, beta_coef_y, vL_y, vR_y)
        call s_weno5_kernel_dispatch(v_rs_ws_z, poly_coef_cbL_z, poly_coef_cbR_z, d_cbL_z, d_cbR_z, beta_coef_z, vL_z, vR_z)

        call s_weno7_kernel(v_rs_ws_x, poly_coef_cbL_x, poly_coef_cbR_x, d_cbL_x, d_cbR_x, beta_coef_x, vL_x, vR_x)
        call s_weno7_kernel(v_rs_ws_y, poly_coef_cbL_y, poly_coef_cbR_y, d_cbL_y, d_cbR_y, beta_coef_y, vL_y, vR_y)
        call s_weno7_kernel(v_rs_ws_z, poly_coef_cbL_z, poly_coef_cbR_z, d_cbL_z, d_cbR_z, beta_coef_z, vL_z, vR_z)

        call s_preserve_monotonicity(v_rs_ws_x, vL_x, vR_x)
        call s_preserve_monotonicity(v_rs_ws_y, vL_y, vR_y)
        call s_preserve_monotonicity(v_rs_ws_z, vL_z, vR_z)
    end subroutine s_run_all_kernels

    subroutine s_run_weno5_kernels(vL_x, vR_x, vL_y, vR_y, vL_z, vR_z)
        real(wp), intent(inout), target :: vL_x(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR_x(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vL_y(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR_y(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vL_z(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR_z(is1_beg:, is2_beg:, is3_beg:, 1:)

        call s_prepare_weno5_inputs(vL_x, vR_x, vL_y, vR_y, vL_z, vR_z)
        call s_run_weno5_reconstruction(vL_x, vR_x, vL_y, vR_y, vL_z, vR_z)
    end subroutine s_run_weno5_kernels

    subroutine s_prepare_weno5_inputs(vL_x, vR_x, vL_y, vR_y, vL_z, vR_z)
        real(wp), intent(inout), target :: vL_x(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR_x(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vL_y(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR_y(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vL_z(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR_z(is1_beg:, is2_beg:, is3_beg:, 1:)

        call s_initialize_weno_x()
        call s_initialize_weno_y()
        call s_initialize_weno_z()

        call s_order1_copy_x(vL_x, vR_x)
        call s_order1_copy_y(vL_y, vR_y)
        call s_order1_copy_z(vL_z, vR_z)
    end subroutine s_prepare_weno5_inputs

    subroutine s_run_weno5_reconstruction(vL_x, vR_x, vL_y, vR_y, vL_z, vR_z)
        real(wp), intent(inout), target :: vL_x(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR_x(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vL_y(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR_y(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vL_z(is1_beg:, is2_beg:, is3_beg:, 1:)
        real(wp), intent(inout), target :: vR_z(is1_beg:, is2_beg:, is3_beg:, 1:)

        if (weno5_specialized_combined .and. .not. weno5_split_kernels) then
            call s_weno5_kernel_x(vL_x, vR_x)
            call s_weno5_kernel_y(vL_y, vR_y)
            call s_weno5_kernel_z(vL_z, vR_z)
        else
            call s_weno5_kernel_dispatch(v_rs_ws_x, poly_coef_cbL_x, poly_coef_cbR_x, d_cbL_x, d_cbR_x, beta_coef_x, vL_x, vR_x)
            call s_weno5_kernel_dispatch(v_rs_ws_y, poly_coef_cbL_y, poly_coef_cbR_y, d_cbL_y, d_cbR_y, beta_coef_y, vL_y, vR_y)
            call s_weno5_kernel_dispatch(v_rs_ws_z, poly_coef_cbL_z, poly_coef_cbR_z, d_cbL_z, d_cbR_z, beta_coef_z, vL_z, vR_z)
        end if
    end subroutine s_run_weno5_reconstruction

end module m_weno_standalone

program weno_test_driver

    use omp_lib
    use m_weno_standalone
    implicit none

    real(wp), allocatable, target :: vL_x(:,:,:,:), vR_x(:,:,:,:)
    real(wp), allocatable, target :: vL_y(:,:,:,:), vR_y(:,:,:,:)
    real(wp), allocatable, target :: vL_z(:,:,:,:), vR_z(:,:,:,:)
    real(wp) :: chk_x, chk_y, chk_z, chk_total
    real(wp) :: t_start, t_end, elapsed_total, elapsed_average
    integer :: num_devices, dev_id, warmup_iters, benchmark_iters, iter
    character(len=16) :: kernel_mode
    character(len=16) :: weno5_bench_scope
    character(len=256) :: input_file

    input_file = 'inputs/weno_input.nml'
    if (command_argument_count() >= 1) then
        call get_command_argument(1, input_file)
    end if

    num_devices = omp_get_num_devices()
    dev_id = omp_get_default_device()
    warmup_iters = max(0, f_get_env_int('WENO_WARMUP_ITERS', 1))
    benchmark_iters = max(1, f_get_env_int('WENO_BENCH_ITERS', 5))
    call s_get_env_string('WENO_KERNEL_MODE', 'all', kernel_mode)
    call s_get_env_string('WENO5_BENCH_SCOPE', 'full', weno5_bench_scope)
    weno5_split_kernels = f_get_env_int('WENO5_SPLIT_KERNELS', 0) /= 0
    weno5_specialized_combined = f_get_env_int('WENO5_SPECIALIZED_COMBINED', 0) /= 0

    write(*,'(a)') '==============================================='
    write(*,'(a)') '  Standalone MFC WENO OpenMP Offload Harness'
    write(*,'(a)') '==============================================='
    write(*,'(a,i0)') '  OpenMP devices available : ', num_devices
    write(*,'(a,i0)') '  Default device           : ', dev_id
    write(*,'(a,a)')  '  Kernel mode              : ', trim(kernel_mode)
    write(*,'(a,a)')  '  WENO5 benchmark scope    : ', trim(weno5_bench_scope)
    write(*,'(a,l1)') '  WENO5 split kernels      : ', weno5_split_kernels
    write(*,'(a,l1)') '  WENO5 specialized comb.  : ', weno5_specialized_combined
    write(*,'(a,i0)') '  Warmup iterations        : ', warmup_iters
    write(*,'(a,i0)') '  Benchmark iterations     : ', benchmark_iters

    call s_read_input(trim(input_file))
    call s_initialize_data()

    allocate(vL_x(is1_beg:is1_end, is2_beg:is2_end, is3_beg:is3_end, 1:v_size))
    allocate(vR_x(is1_beg:is1_end, is2_beg:is2_end, is3_beg:is3_end, 1:v_size))
    allocate(vL_y(is1_beg:is1_end, is2_beg:is2_end, is3_beg:is3_end, 1:v_size))
    allocate(vR_y(is1_beg:is1_end, is2_beg:is2_end, is3_beg:is3_end, 1:v_size))
    allocate(vL_z(is1_beg:is1_end, is2_beg:is2_end, is3_beg:is3_end, 1:v_size))
    allocate(vR_z(is1_beg:is1_end, is2_beg:is2_end, is3_beg:is3_end, 1:v_size))

    vL_x = 0.0_wp; vR_x = 0.0_wp
    vL_y = 0.0_wp; vR_y = 0.0_wp
    vL_z = 0.0_wp; vR_z = 0.0_wp

    !$omp target enter data map(alloc:vL_x, vR_x, vL_y, vR_y, vL_z, vR_z)
    !$omp target update to(vL_x, vR_x, vL_y, vR_y, vL_z, vR_z)

    select case (trim(kernel_mode))
    case ('weno5')
        if (trim(weno5_bench_scope) == 'kernel') then
            call s_prepare_weno5_inputs(vL_x, vR_x, vL_y, vR_y, vL_z, vR_z)
            do iter = 1, warmup_iters
                call s_run_weno5_reconstruction(vL_x, vR_x, vL_y, vR_y, vL_z, vR_z)
            end do
        else
            do iter = 1, warmup_iters
                call s_run_weno5_kernels(vL_x, vR_x, vL_y, vR_y, vL_z, vR_z)
            end do
        end if
    case default
        do iter = 1, warmup_iters
            call s_run_all_kernels(vL_x, vR_x, vL_y, vR_y, vL_z, vR_z)
        end do
    end select

    t_start = omp_get_wtime()
    select case (trim(kernel_mode))
    case ('weno5')
        if (trim(weno5_bench_scope) == 'kernel') then
            do iter = 1, benchmark_iters
                call s_run_weno5_reconstruction(vL_x, vR_x, vL_y, vR_y, vL_z, vR_z)
            end do
        else
            do iter = 1, benchmark_iters
                call s_run_weno5_kernels(vL_x, vR_x, vL_y, vR_y, vL_z, vR_z)
            end do
        end if
    case default
        do iter = 1, benchmark_iters
            call s_run_all_kernels(vL_x, vR_x, vL_y, vR_y, vL_z, vR_z)
        end do
    end select
    t_end = omp_get_wtime()

    !$omp taskwait
    !$omp target update from(vL_x, vR_x, vL_y, vR_y, vL_z, vR_z)

    elapsed_total = t_end - t_start
    elapsed_average = elapsed_total/real(benchmark_iters, wp)

    chk_x = f_checksum(vL_x) + f_checksum(vR_x)
    chk_y = f_checksum(vL_y) + f_checksum(vR_y)
    chk_z = f_checksum(vL_z) + f_checksum(vR_z)
    chk_total = chk_x + chk_y + chk_z

    write(*,'(a,es18.10)') '  Checksum X : ', chk_x
    write(*,'(a,es18.10)') '  Checksum Y : ', chk_y
    write(*,'(a,es18.10)') '  Checksum Z : ', chk_z
    write(*,'(a,es18.10)') '  Checksum T : ', chk_total
    write(*,'(a,f12.6,a)') '  Kernel time total        : ', elapsed_total, ' s'
    write(*,'(a,f12.6,a)') '  Kernel time per iter     : ', elapsed_average, ' s'

    if (.not. ieee_is_finite(chk_total)) then
        write(*,'(a)') '  WARNING: non-finite values detected in standalone output.'
    else if (chk_total == 0.0_wp) then
        write(*,'(a)') '  WARNING: all-zero output; offload may not have executed.'
    else
        write(*,'(a)') '  All kernel families executed with finite non-zero output.'
    end if

    !$omp target exit data map(release:vL_x, vR_x, vL_y, vR_y, vL_z, vR_z)
    deallocate(vL_x, vR_x, vL_y, vR_y, vL_z, vR_z)

    call s_finalize_data()

end program weno_test_driver
