module allocator_mod

  use flogger
  use mesh_mod

  implicit none

  private

  public allocate_array

  interface allocate_array
    module procedure allocate_array_1d_r4
    module procedure allocate_array_2d_r4
    module procedure allocate_array_3d_r4
    module procedure allocate_array_1d_r8
    module procedure allocate_array_2d_r8
    module procedure allocate_array_3d_r8
    module procedure allocate_array_3d_i4
    module procedure allocate_array_pointer_3d_r4
    module procedure allocate_array_pointer_3d_r8
    module procedure allocate_array_1d_r16
    module procedure allocate_array_2d_r16
    module procedure allocate_array_3d_r16
  end interface allocate_array

contains

  subroutine allocate_array_1d_r4(mesh, array, full_lon, half_lon, full_lat, half_lat)

    type(mesh_type), intent(in )              :: mesh
    real(4)        , intent(out), allocatable :: array(:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat

    if (allocated(array)) deallocate(array)

    if (present(full_lon)) then
      if (full_lon) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub))
    else if (present(half_lon)) then
      if (half_lon) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub))
    else if (present(full_lat)) then
      if (full_lat) allocate(array(mesh%full_lat_lb:mesh%full_lat_ub))
    else if (present(half_lat)) then
      if (half_lat) allocate(array(mesh%half_lat_lb:mesh%half_lat_ub))
    else
      call log_error('allocate_array_1d_r4: Missing full_lon, half_lon, full_lat or half_lat argument!')
    end if

    array = 0.0d0

  end subroutine allocate_array_1d_r4

  subroutine allocate_array_2d_r4(mesh, array, full_lon, half_lon, full_lat, half_lat)

    type(mesh_type), intent(in )              :: mesh
    real(4)        , intent(out), allocatable :: array(:,:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat

    if (present(full_lon) .and. present(full_lat)) then
      if (full_lon .and. full_lat) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub))
    else if (present(full_lon) .and. present(half_lat)) then
      if (full_lon .and. half_lat) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub))
    else if (present(half_lon) .and. present(full_lat)) then
      if (half_lon .and. full_lat) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub))
    else if (present(half_lon) .and. present(half_lat)) then
      if (half_lon .and. half_lat) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub))
    else
      allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub))
    end if

    array = 0.0d0

  end subroutine allocate_array_2d_r4

  subroutine allocate_array_3d_r4(mesh, array, full_lon, half_lon, full_lat, half_lat, full_lev, half_lev)

    type(mesh_type), intent(in )              :: mesh
    real(4)        , intent(out), allocatable :: array(:,:,:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat
    logical        , intent(in ), optional    :: full_lev
    logical        , intent(in ), optional    :: half_lev

    if (present(full_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (full_lon .and. full_lat .and. full_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(full_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (full_lon .and. half_lat .and. full_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(full_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (full_lon .and. full_lat .and. half_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else if (present(full_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (full_lon .and. half_lat .and. half_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else if (present(half_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (half_lon .and. full_lat .and. full_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(half_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (half_lon .and. half_lat .and. full_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(half_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (half_lon .and. full_lat .and. half_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else if (present(half_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (half_lon .and. half_lat .and. half_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else
      allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    end if

    array = 0.0d0

  end subroutine allocate_array_3d_r4

  subroutine allocate_array_1d_r8(mesh, array, full_lon, half_lon, full_lat, half_lat)

    type(mesh_type), intent(in )              :: mesh
    real(8)        , intent(out), allocatable :: array(:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat

    if (allocated(array)) deallocate(array)

    if (present(full_lon)) then
      if (full_lon) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub))
    else if (present(half_lon)) then
      if (half_lon) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub))
    else if (present(full_lat)) then
      if (full_lat) allocate(array(mesh%full_lat_lb:mesh%full_lat_ub))
    else if (present(half_lat)) then
      if (half_lat) allocate(array(mesh%half_lat_lb:mesh%half_lat_ub))
    else
      call log_error('allocate_array_1d_r8: Missing full_lon, half_lon, full_lat or half_lat argument!')
    end if

    array = 0.0d0

  end subroutine allocate_array_1d_r8

  subroutine allocate_array_2d_r8(mesh, array, full_lon, half_lon, full_lat, half_lat)

    type(mesh_type), intent(in )              :: mesh
    real(8)        , intent(out), allocatable :: array(:,:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat

    if (present(full_lon) .and. present(full_lat)) then
      if (full_lon .and. full_lat) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub))
    else if (present(full_lon) .and. present(half_lat)) then
      if (full_lon .and. half_lat) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub))
    else if (present(half_lon) .and. present(full_lat)) then
      if (half_lon .and. full_lat) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub))
    else if (present(half_lon) .and. present(half_lat)) then
      if (half_lon .and. half_lat) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub))
    else
      allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub))
    end if

    array = 0.0d0

  end subroutine allocate_array_2d_r8

  subroutine allocate_array_3d_r8(mesh, array, full_lon, half_lon, full_lat, half_lat, full_lev, half_lev)

    type(mesh_type), intent(in )              :: mesh
    real(8)        , intent(out), allocatable :: array(:,:,:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat
    logical        , intent(in ), optional    :: full_lev
    logical        , intent(in ), optional    :: half_lev

    if (present(full_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (full_lon .and. full_lat .and. full_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(full_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (full_lon .and. half_lat .and. full_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(full_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (full_lon .and. full_lat .and. half_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else if (present(full_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (full_lon .and. half_lat .and. half_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else if (present(half_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (half_lon .and. full_lat .and. full_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(half_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (half_lon .and. half_lat .and. full_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(half_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (half_lon .and. full_lat .and. half_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else if (present(half_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (half_lon .and. half_lat .and. half_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else
      allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    end if

    array = 0.0d0

  end subroutine allocate_array_3d_r8

  subroutine allocate_array_3d_i4(mesh, array, full_lon, half_lon, full_lat, half_lat, full_lev, half_lev)

    type(mesh_type), intent(in )              :: mesh
    integer(4)     , intent(out), allocatable :: array(:,:,:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat
    logical        , intent(in ), optional    :: full_lev
    logical        , intent(in ), optional    :: half_lev

    if (present(full_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (full_lon .and. full_lat .and. full_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(full_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (full_lon .and. half_lat .and. full_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(full_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (full_lon .and. full_lat .and. half_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else if (present(full_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (full_lon .and. half_lat .and. half_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else if (present(half_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (half_lon .and. full_lat .and. full_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(half_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (half_lon .and. half_lat .and. full_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(half_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (half_lon .and. full_lat .and. half_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else if (present(half_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (half_lon .and. half_lat .and. half_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else
      allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    end if

    array = 0

  end subroutine allocate_array_3d_i4

  subroutine allocate_array_pointer_3d_r4(mesh, array, full_lon, half_lon, full_lat, half_lat, full_lev, half_lev)

    type(mesh_type), intent(in )              :: mesh
    real(4)        , intent(out), pointer     :: array(:,:,:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat
    logical        , intent(in ), optional    :: full_lev
    logical        , intent(in ), optional    :: half_lev

    if (present(full_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (full_lon .and. full_lat .and. full_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(full_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (full_lon .and. half_lat .and. full_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(full_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (full_lon .and. full_lat .and. half_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else if (present(full_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (full_lon .and. half_lat .and. half_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else if (present(half_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (half_lon .and. full_lat .and. full_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(half_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (half_lon .and. half_lat .and. full_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(half_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (half_lon .and. full_lat .and. half_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else if (present(half_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (half_lon .and. half_lat .and. half_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else
      allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    end if

    array = 0.0d0

  end subroutine allocate_array_pointer_3d_r4

  subroutine allocate_array_pointer_3d_r8(mesh, array, full_lon, half_lon, full_lat, half_lat, full_lev, half_lev)

    type(mesh_type), intent(in )              :: mesh
    real(8)        , intent(out), pointer     :: array(:,:,:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat
    logical        , intent(in ), optional    :: full_lev
    logical        , intent(in ), optional    :: half_lev

    if (present(full_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (full_lon .and. full_lat .and. full_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(full_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (full_lon .and. half_lat .and. full_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(full_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (full_lon .and. full_lat .and. half_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else if (present(full_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (full_lon .and. half_lat .and. half_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else if (present(half_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (half_lon .and. full_lat .and. full_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(half_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (half_lon .and. half_lat .and. full_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(half_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (half_lon .and. full_lat .and. half_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else if (present(half_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (half_lon .and. half_lat .and. half_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else
      allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    end if

    array = 0.0d0

  end subroutine allocate_array_pointer_3d_r8

  subroutine allocate_array_1d_r16(mesh, array, full_lon, half_lon, full_lat, half_lat)

    type(mesh_type), intent(in )              :: mesh
    real(16)        , intent(out), allocatable :: array(:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat

    if (allocated(array)) deallocate(array)

    if (present(full_lon)) then
      if (full_lon) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub))
    else if (present(half_lon)) then
      if (half_lon) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub))
    else if (present(full_lat)) then
      if (full_lat) allocate(array(mesh%full_lat_lb:mesh%full_lat_ub))
    else if (present(half_lat)) then
      if (half_lat) allocate(array(mesh%half_lat_lb:mesh%half_lat_ub))
    else
      call log_error('allocate_array_1d_r16: Missing full_lon, half_lon, full_lat or half_lat argument!')
    end if

    array = 0.0d0

  end subroutine allocate_array_1d_r16

  subroutine allocate_array_2d_r16(mesh, array, full_lon, half_lon, full_lat, half_lat)

    type(mesh_type), intent(in )              :: mesh
    real(16)        , intent(out), allocatable :: array(:,:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat

    if (present(full_lon) .and. present(full_lat)) then
      if (full_lon .and. full_lat) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub))
    else if (present(full_lon) .and. present(half_lat)) then
      if (full_lon .and. half_lat) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub))
    else if (present(half_lon) .and. present(full_lat)) then
      if (half_lon .and. full_lat) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub))
    else if (present(half_lon) .and. present(half_lat)) then
      if (half_lon .and. half_lat) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub))
    else
      allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub))
    end if

    array = 0.0d0

  end subroutine allocate_array_2d_r16

  subroutine allocate_array_3d_r16(mesh, array, full_lon, half_lon, full_lat, half_lat, full_lev, half_lev)

    type(mesh_type), intent(in )              :: mesh
    real(16)        , intent(out), allocatable :: array(:,:,:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat
    logical        , intent(in ), optional    :: full_lev
    logical        , intent(in ), optional    :: half_lev

    if (present(full_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (full_lon .and. full_lat .and. full_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(full_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (full_lon .and. half_lat .and. full_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(full_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (full_lon .and. full_lat .and. half_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else if (present(full_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (full_lon .and. half_lat .and. half_lev) allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else if (present(half_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (half_lon .and. full_lat .and. full_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(half_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (half_lon .and. half_lat .and. full_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    else if (present(half_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (half_lon .and. full_lat .and. half_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else if (present(half_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (half_lon .and. half_lat .and. half_lev) allocate(array(mesh%half_lon_lb:mesh%half_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
    else
      allocate(array(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
    end if

    array = 0.0d0

  end subroutine allocate_array_3d_r16

end module allocator_mod
