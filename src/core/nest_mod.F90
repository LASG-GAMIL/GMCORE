module nest_mod

  use const_mod
  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod

  implicit none

  private

  public nest_add_domain
  public nest_from_coarse_to_fine
  public nest_from_fine_to_coarse

  type(mesh_type  ), allocatable, dimension(:  ) :: nested_mesh
  type(static_type), allocatable, dimension(:  ) :: nested_static
  type(state_type ), allocatable, dimension(:,:) :: nested_state
  type(tend_type  ), allocatable, dimension(:,:) :: nested_tend

contains

  subroutine nest_add_domain(start_lon, end_lon, start_lat, end_lat, root_mesh)

    real(r8)       , intent(in) :: start_lon
    real(r8)       , intent(in) :: end_lon
    real(r8)       , intent(in) :: start_lat
    real(r8)       , intent(in) :: end_lat
    type(mesh_type), intent(in) :: root_mesh
    
    integer :: time_substep_num
    
    select case(nest_time_scheme)
      case('pc2'.or.'rk2')
        time_substep_num = 2
      case('pc3'.or.'rk3')
        time_substep_num = 3
    end select
    
    allocate(nested_mesh  (nest_max_dom                    ))
    allocate(nested_static(nest_max_dom                    ))
    allocate(nested_state (nest_max_dom, 0:time_substep_num))
    allocate(nested_tend  (nest_max_dom, 0:time_substep_num))
    
  end subroutine nest_add_domain

  subroutine nest_init_domain(static, state)

    type(static_type), intent(in) :: static
    type(state_type) , intent(in) :: state

    ! static -> nested_static
    ! state -> nested_state

  end subroutine nest_init_domain

  subroutine nest_from_coarse_to_fine(root_mesh, root_state, nested_mesh, nested_state)

    type(mesh_type) , intent(in   ) :: root_mesh
    type(state_type), intent(in   ) :: root_state
    type(mesh_type) , intent(inout) :: nested_mesh
    type(state_type), intent(inout) :: nested_state

    ! Interpolate from coarse grids to halo of fine grids.

  end subroutine nest_from_coarse_to_fine

  subroutine nest_from_fine_to_coarse(root_mesh, root_state, nested_mesh, nested_state)

    type(mesh_type) , intent(inout) :: root_mesh
    type(state_type), intent(inout) :: root_state
    type(mesh_type) , intent(inout) :: nested_mesh
    type(state_type), intent(inout) :: nested_state

    ! Update coarse grids from fine grids.

  end subroutine nest_from_fine_to_coarse

end module nest_mod