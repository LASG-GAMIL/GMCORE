module gmcore_mod

  use flogger
  use parallel_mod
  use time_mod, old => old_time_idx, new => new_time_idx
  use history_mod
  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod

  implicit none

  private

  public gmcore_init
  public gmcore_run
  public gmcore_final

  ! Diagnostics
  real(real_kind) total_mass
  real(real_kind) total_energy
  real(real_kind) total_enstrophy

contains

  subroutine gmcore_init()

    call log_init()
    call parallel_init()
    call time_init()
    call create_meshes()
    call create_states()
    call create_static()
    call create_tends()
    call history_init()

  end subroutine gmcore_init

  subroutine gmcore_run()

    call diagnose(old)
    call output(old)
    call log_add_diag('total_mass', total_mass)
    call log_add_diag('total_energy', total_energy)
    call log_print_diag(curr_time%isoformat())

  end subroutine gmcore_run

  subroutine gmcore_final()

    call parallel_final()

  end subroutine gmcore_final

  subroutine output(time_idx)

    integer, intent(in) :: time_idx

    if (time_is_alerted('history_write')) call history_write(state(time_idx), static)

  end subroutine output

  subroutine diagnose(time_idx)

    integer, intent(in) :: time_idx

    type(mesh_type), pointer :: mesh
    real(real_kind), pointer :: u(:,:)
    real(real_kind), pointer :: v(:,:)
    real(real_kind), pointer :: gd(:,:)
    real(real_kind), pointer :: pv(:,:)
    real(real_kind), pointer :: hat_gd(:,:)
    real(real_kind), pointer :: ghs(:,:)
    integer i, j

    mesh => state(time_idx)%mesh

    u   => state(time_idx)%u
    v   => state(time_idx)%v
    gd  => state(time_idx)%gd
    pv  => state(time_idx)%pv
    ghs => static%ghs

    total_mass = 0.0
    do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        total_mass = total_mass + gd(i,j) * mesh%cell_area(j)
      end do
    end do

    total_energy = 0.0
    do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        total_energy = total_energy + (gd(i,j) + ghs(i,j))**2 * mesh%cell_area(j)
      end do
    end do
    do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        total_energy = total_energy + gd(i,j) * u(i,j)**2 * 2.0d0 * mesh%lon_edge_area(j)
      end do
    end do
    do j = mesh%half_lat_start_idx, mesh%half_lat_end_idx
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        total_energy = total_energy + gd(i,j) * v(i,j)**2 * 2.0d0 * mesh%lat_edge_area(j)
      end do
    end do

    do j = mesh%half_lat_start_idx, mesh%half_lat_end_idx
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        ! hat_gd(i,j) = 
      end do
    end do

    do j = mesh%half_lat_start_idx, mesh%half_lat_end_idx
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        pv(i,j) = (u(i,j  ) * mesh%cell_lon_distance(j  )  - &
                   u(i,j+1) * mesh%cell_lon_distance(j+1)  + &
                   v(i,j  ) * mesh%cell_lat_distance(j  )  - &
                   v(i-1,j) * mesh%cell_lat_distance(j  )) / &
                  mesh%vertex_area(j)
      end do
    end do

  end subroutine diagnose

end module gmcore_mod
