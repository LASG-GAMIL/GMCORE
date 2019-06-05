module namelist_mod

  use mesh_mod

  implicit none

  namelist /gmcore_swm_control/ &
    num_total_lon,              &
    num_total_lat

contains

  subroutine parse_namelist(file_path)

    character(*), intent(in) :: file_path

    open(10, file=file_path, status='old')
    read(10, nml=gmcore_swm_control)
    close(10)

  end subroutine parse_namelist

end module namelist_mod