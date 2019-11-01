module nest_mod

  use const_mod
  use namelist_mod
  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod

  implicit none

  private

  public nest_init
  public nest_downward
  public nest_upward
  public nest_max_dom
  public nested_meshes
  public nested_statics
  public nested_states
  public nested_tends

  integer :: num_nested_dom = 0
  type(mesh_type  ), allocatable, dimension(:  ), target :: nested_meshes
  type(static_type), allocatable, dimension(:  ), target :: nested_statics
  type(state_type ), allocatable, dimension(:,:), target :: nested_states
  type(tend_type  ), allocatable, dimension(:,:), target :: nested_tends

contains

  subroutine nest_init()

    integer i

    do i = 1, nest_max_dom
      call nest_add_domain(nest_start_lon(i), nest_end_lon(i), nest_start_lat(i), nest_end_lat(i), nest_parent_id(i))
    end do

  end subroutine nest_init

  subroutine nest_add_domain(start_lon, end_lon, start_lat, end_lat, parent_id)

    real(r8), intent(in) :: start_lon
    real(r8), intent(in) :: end_lon
    real(r8), intent(in) :: start_lat
    real(r8), intent(in) :: end_lat
    integer , intent(in) :: parent_id

    type(mesh_type), pointer :: parent_mesh
    type(mesh_type), pointer :: nested_mesh
    type(state_type), pointer :: parent_states(:)
    real(r8) dlon, dlat
    integer num_lon, num_lat
    integer i, j

    if (.not. allocated(nested_meshes)) allocate(nested_meshes(nest_max_dom))

    if (parent_id == 0) then
      parent_mesh => mesh
      parent_states => states
    else
      do i = 1, size(nested_meshes)
        if (nested_meshes(i)%id == parent_id) then
          parent_mesh => nested_meshes(i)
          parent_states => nested_states(:,i)
          exit
        end if
      end do
    end if

    num_nested_dom = num_nested_dom + 1
    nested_mesh => nested_meshes(num_nested_dom)
    nested_mesh%id = num_nested_dom
    nested_mesh%parent => parent_mesh

    ! Find the parent grid indices of the four corners.
    do i = parent_mesh%full_lon_start_idx, parent_mesh%full_lon_end_idx - 1
      if (parent_mesh%full_lon_deg(i) <= start_lon .and. start_lon < parent_mesh%full_lon_deg(i+1)) then
        nested_mesh%parent_lon_start_idx = i
        exit
      end if
    end do

    do i = parent_mesh%full_lon_start_idx + 1, parent_mesh%full_lon_end_idx
      if (parent_mesh%full_lon_deg(i) < end_lon .and. end_lon <= parent_mesh%full_lon_deg(i+1)) then
        nested_mesh%parent_lon_end_idx = i
        exit
      end if
    end do

    do j = parent_mesh%full_lat_start_idx, parent_mesh%full_lat_end_idx - 1
      if (parent_mesh%full_lat_deg(j) <= start_lat .and. start_lat < parent_mesh%full_lat_deg(j+1)) then
        nested_mesh%parent_lat_start_idx = j
        exit
      end if
    end do

    do j = parent_mesh%full_lat_start_idx + 1, parent_mesh%full_lat_end_idx
      if (parent_mesh%full_lat_deg(j) < end_lat .and. end_lat <= parent_mesh%full_lat_deg(j+1)) then
        nested_mesh%parent_lat_end_idx = j
        exit
      end if
    end do

    dlon = parent_mesh%dlon / 3.0_r8
    dlat = parent_mesh%dlat / 3.0_r8

    num_lon = (                                              &
      parent_mesh%full_lon(nested_mesh%parent_lon_end_idx) - &
      parent_mesh%full_lon(nested_mesh%parent_lon_start_idx) &
    ) / dlon
    num_lat = (                                              &
      parent_mesh%full_lat(nested_mesh%parent_lat_end_idx) - &
      parent_mesh%full_lat(nested_mesh%parent_lat_start_idx) &
    ) / dlat

    call nested_mesh%init(num_lon, num_lat, num_nested_dom, parent_mesh%halo_width, &
                          parent_mesh%full_lon(nested_mesh%parent_lon_start_idx)  , &
                          parent_mesh%full_lon(nested_mesh%parent_lon_end_idx  )  , &
                          parent_mesh%full_lat(nested_mesh%parent_lat_start_idx)  , &
                          parent_mesh%full_lat(nested_mesh%parent_lat_end_idx  ))

    if (.not. allocated(nested_states)) then
      select case (nest_time_scheme)
      case ('pc2', 'rk2')
        allocate(nested_states(0:2,nest_max_dom))
        allocate(nested_tends (0:2,nest_max_dom))
      case ('pc3', 'rk3')
        allocate(nested_states(0:3,nest_max_dom))
        allocate(nested_tends (0:3,nest_max_dom))
      end select
      allocate(nested_statics(nest_max_dom))
    end if

    do i = lbound(nested_states, 1), ubound(nested_states, 1)
      call nested_states(i,num_nested_dom)%init(nested_mesh)
      nested_states(i,num_nested_dom)%id = nest_max_dom
      nested_states(i,num_nested_dom)%parent => parent_states(i)
      call nested_tends(i,num_nested_dom)%init(nested_mesh)
    end do
    call nested_statics(num_nested_dom)%init(nested_mesh)

  end subroutine nest_add_domain

  subroutine nest_downward(time_idx)

    integer, intent(in) :: time_idx

    integer i

    do i = 1, nest_max_dom
      call nest_coarse_to_fine(nested_states(time_idx,i)%parent, nested_states(time_idx,i))
    end do

  end subroutine nest_downward

  subroutine nest_coarse_to_fine(parent_state, nested_state)

    type(state_type), intent(in   ) :: parent_state
    type(state_type), intent(inout) :: nested_state

    ! Interpolate from coarse grids to halo of fine grids.

  end subroutine nest_coarse_to_fine

  subroutine nest_upward(time_idx)

    integer, intent(in) :: time_idx

  end subroutine nest_upward

  subroutine nest_fine_to_coarse(parent_state, nested_state)

    type(state_type), intent(inout) :: parent_state
    type(state_type), intent(in   ) :: nested_state

    ! Update coarse grids from fine grids.

  end subroutine nest_fine_to_coarse

end module nest_mod