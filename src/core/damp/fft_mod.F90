module fft_mod

  use flogger
  use const_mod
  use block_mod
  use parallel_mod

  implicit none

  private

  public fft_init
  public fft_final
  public fft_array_on_vtx

  real, allocatable :: wave_array(:)
  real, allocatable :: work_array(:)

  real, allocatable :: full_fft_factor(:,:)
  real, allocatable :: half_fft_factor(:,:)

  logical, allocatable :: fft_on_full_lat(:)
  logical, allocatable :: fft_on_half_lat(:)

  integer :: fft_cutoff_wavenumber(100) = 0

contains

  subroutine fft_init()

    integer i, j, jr, n, ierr

    call fft_final()

    fft_cutoff_wavenumber(1:8) = [1, 1, 2, 3, 3, 4, 4, 4]

    allocate(fft_on_full_lat(global_mesh%num_full_lat))
    allocate(fft_on_half_lat(global_mesh%num_half_lat))

    fft_on_full_lat(:) = .false.
    fft_on_half_lat(:) = .false.

    do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
      if (global_mesh%full_lat(j) <= 0) then
        jr = j - global_mesh%full_lat_ibeg_no_pole + 1
      else
        jr = global_mesh%full_lat_iend_no_pole - j + 1
      end if
      if (fft_cutoff_wavenumber(jr) > 0) then
        fft_on_full_lat(j) = .true.
      end if
    end do
    do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
      if (global_mesh%half_lat(j) <= 0) then
        jr = j - global_mesh%half_lat_ibeg_no_pole + 1
      else
        jr = global_mesh%half_lat_iend_no_pole - j + 1
      end if
#ifdef V_POLE
      if (fft_cutoff_wavenumber(jr) > 0 .or. fft_cutoff_wavenumber(jr-1) > 0) then
#else
      if (fft_cutoff_wavenumber(jr) > 0 .or. fft_cutoff_wavenumber(jr+1) > 0) then
#endif
        fft_on_half_lat(j) = .true.
      end if
    end do

    ! N + INT(LOG(REAL(N))) + 4
    allocate(wave_array(global_mesh%num_full_lon + int(log(real(global_mesh%num_full_lon)) / log(2.0)) + 4))
    allocate(work_array(global_mesh%num_full_lon))

    allocate(full_fft_factor(global_mesh%num_full_lon,global_mesh%num_full_lat))
    allocate(half_fft_factor(global_mesh%num_full_lon,global_mesh%num_half_lat))

    call rfft1i(global_mesh%num_full_lon, wave_array, size(wave_array), ierr)
    if (ierr /= 0) then
      call log_error('Failed to initialize FFTPACK!')
    end if

    full_fft_factor = 0.0
    half_fft_factor = 0.0
    n = global_mesh%num_full_lon / 2
    do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
      if (global_mesh%full_lat(j) <= 0) then
        jr = j - global_mesh%full_lat_ibeg_no_pole + 1
      else
        jr = global_mesh%full_lat_iend_no_pole - j + 1
      end if
      if (fft_cutoff_wavenumber(jr) > 0) then
        do i = 1, fft_cutoff_wavenumber(jr) + 1
          full_fft_factor(2*i-1,j) = 1.0
          full_fft_factor(2*i,  j) = 1.0
        end do
      end if
    end do
    do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
      if (global_mesh%half_lat(j) <= 0) then
        jr = j - global_mesh%half_lat_ibeg_no_pole + 1
      else
        jr = global_mesh%half_lat_iend_no_pole - j + 1
      end if
#ifdef V_POLE
      if (fft_cutoff_wavenumber(jr) > 0 .or. fft_cutoff_wavenumber(jr-1) > 0) then
        do i = 1, min(fft_cutoff_wavenumber(jr), fft_cutoff_wavenumber(jr-1)) + 1
#else
      if (fft_cutoff_wavenumber(jr) > 0 .or. fft_cutoff_wavenumber(jr+1) > 0) then
        do i = 1, min(fft_cutoff_wavenumber(jr), fft_cutoff_wavenumber(jr+1)) + 1
#endif
          half_fft_factor(2*i-1,j) = 1.0
          half_fft_factor(2*i,  j) = 1.0
        end do
      end if
    end do

  end subroutine fft_init

  subroutine fft_final()

    if (allocated(fft_on_full_lat)) deallocate(fft_on_full_lat)
    if (allocated(fft_on_half_lat)) deallocate(fft_on_half_lat)
    if (allocated(wave_array     )) deallocate(wave_array)
    if (allocated(work_array     )) deallocate(work_array)
    if (allocated(full_fft_factor)) deallocate(full_fft_factor)
    if (allocated(half_fft_factor)) deallocate(half_fft_factor)

  end subroutine fft_final

  subroutine fft_array_on_vtx(block, f)

    type(block_type), intent(in), target :: block
    real(r8), intent(inout) :: f(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    real local_f(global_mesh%num_full_lon)
    integer i, j, k, ierr

    mesh => block%mesh

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        if (fft_on_half_lat(j)) then
          local_f = f(mesh%half_lon_ibeg:mesh%half_lon_iend,j,k)
          call rfft1f(global_mesh%num_full_lon, 1, local_f, global_mesh%num_full_lon, wave_array, size(wave_array), work_array, size(work_array), ierr)
          if (ierr /= 0) then
            call log_error('Failed to do forward FFT!')
          end if

          local_f = local_f * half_fft_factor(:,j)
          call rfft1b(global_mesh%num_full_lon, 1, local_f, global_mesh%num_full_lon, wave_array, size(wave_array), work_array, size(work_array), ierr)
          if (ierr /= 0) then
            call log_error('Failed to do backward FFT!')
          end if

          f(mesh%half_lon_ibeg:mesh%half_lon_iend,j,k) = local_f

          call fill_zonal_halo(block, mesh%lon_halo_width, f(:,j,k))
        end if
      end do
    end do

  end subroutine fft_array_on_vtx

end module fft_mod
