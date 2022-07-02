module hybrid_coord_mars_mod

  use flogger
  use const_mod
  use process_mod
  use mesh_mod

  implicit none

  private

  public hybrid_coord_mars_emars28

contains

  subroutine hybrid_coord_mars_emars28(p0, ptop, hyai, hybi)

    real(r8), intent(out) :: p0
    real(r8), intent(out) :: ptop
    real(r8), intent(out) :: hyai(29)
    real(r8), intent(out) :: hybi(29)

    if (global_mesh%num_full_lev /= 28 .and. is_root_proc()) then
      call log_error('num_lev should be 28 in namelist!')
    end if

    hyai = [       &
      0.02       , & !  1
      0.05738127 , & !  2
      0.1958398  , & !  3
      0.5922958  , & !  4
      1.566023   , & !  5
      2.445497   , & !  6
      2.768375   , & !  7
      2.885169   , & !  8
      2.917223   , & !  9
      2.908704   , & ! 10
      2.859894   , & ! 11
      2.768765   , & ! 12
      2.632705   , & ! 13
      2.450922   , & ! 14
      2.226681   , & ! 15
      1.968468   , & ! 16
      1.689483   , & ! 17
      1.405581   , & ! 18
      1.132426   , & ! 19
      0.8828918  , & ! 20
      0.6654847  , & ! 21
      0.4840102  , & ! 22
      0.3382412  , & ! 23
      0.225107   , & ! 24
      0.1399572  , & ! 25
      0.07761155 , & ! 26
      0.0330855  , & ! 27
      0.002      , & ! 28
      0.0          & ! 29
    ]
    hybi = [       &
      0.0        , & !  1
      0.0        , & !  2
      0.0        , & !  3
      0.0        , & !  4
      0.0        , & !  5
      0.001936639, & !  6
      0.007441913, & !  7
      0.01622727 , & !  8
      0.02707519 , & !  9
      0.043641   , & ! 10
      0.0681068  , & ! 11
      0.1028024  , & ! 12
      0.1497195  , & ! 13
      0.2098713  , & ! 14
      0.2827023  , & ! 15
      0.3658161  , & ! 16
      0.4552023  , & ! 17
      0.545936   , & ! 18
      0.6331097  , & ! 19
      0.7126763  , & ! 20
      0.7819615  , & ! 21
      0.8397753  , & ! 22
      0.8862035  , & ! 23
      0.9222317  , & ! 24
      0.9493454  , & ! 25
      0.9691962  , & ! 26
      0.9833726  , & ! 27
      0.9932694  , & ! 28
      1.0          & ! 29
    ]

    p0 = 1
    ptop = p0 * hyai(1)

  end subroutine hybrid_coord_mars_emars28

end module hybrid_coord_mars_mod
