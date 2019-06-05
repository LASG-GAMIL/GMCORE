program gmcore_swm_driver

  use namelist_mod
  use gmcore_mod

  implicit none

  call parse_namelist('namelist.gmcore_swm')
  call gmcore_init()
    
end program gmcore_swm_driver
