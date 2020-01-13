module state_mod

  use const_mod
  use mesh_mod
  use namelist_mod
  use allocator_mod
  use parallel_mod

  implicit none

  private

  public state_type

  type state_type
    type(mesh_type), pointer :: mesh => null()
    ! For nesting
    integer :: id = 0
    type(state_type), pointer :: parent => null()
    ! Prognostic variables
    real(r8), allocatable, dimension(:,:) :: u
    real(r8), allocatable, dimension(:,:) :: v
    real(r8), allocatable, dimension(:,:) :: gd
    ! Diagnostic variables
    real(r8), allocatable, dimension(:,:) :: pv
    real(r8), allocatable, dimension(:,:) :: m_vtx
    real(r8), allocatable, dimension(:,:) :: m_lon
    real(r8), allocatable, dimension(:,:) :: m_lat
    real(r8), allocatable, dimension(:,:) :: mf_lon_n
    real(r8), allocatable, dimension(:,:) :: mf_lat_n
    real(r8), allocatable, dimension(:,:) :: mf_lat_t
    real(r8), allocatable, dimension(:,:) :: mf_lon_t
    real(r8), allocatable, dimension(:,:) :: pv_lon
    real(r8), allocatable, dimension(:,:) :: pv_lat
    real(r8), allocatable, dimension(:,:) :: dpv_lon_t
    real(r8), allocatable, dimension(:,:) :: dpv_lon_n
    real(r8), allocatable, dimension(:,:) :: dpv_lat_t
    real(r8), allocatable, dimension(:,:) :: dpv_lat_n
    real(r8), allocatable, dimension(:,:) :: ke
    real(r8) total_m
    real(r8) total_ke
    real(r8) total_e
    real(r8) total_av
    real(r8) total_pe
  contains
    procedure :: init => state_init
    procedure :: update_m_lon_m_lat => state_update_m_lon_m_lat
    procedure :: update_m_vtx => state_update_m_vtx
    procedure :: update_mf_lon_n_mf_lat_n => state_update_mf_lon_n_mf_lat_n
    procedure :: update_mf_lon_t_mf_lat_t => state_update_mf_lon_t_mf_lat_t
    procedure :: update_pv_vtx => state_update_pv_vtx
    procedure :: update_ke => state_update_ke
    procedure :: update_all => state_update_all
    procedure :: clear => state_clear
    final :: state_final
  end type state_type

contains

  subroutine state_init(this, mesh)

    class(state_type), intent(inout)         :: this
    type(mesh_type  ), intent(in   ), target :: mesh

    call this%clear()

    this%mesh => mesh

    call allocate_array(mesh, this%u        , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%v        , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%gd       , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%pv       , half_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%m_vtx    , half_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%m_lon    , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%m_lat    , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%mf_lon_n , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%mf_lon_t , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%mf_lat_n , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%mf_lat_t , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%pv_lon   , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%pv_lat   , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%dpv_lon_t, half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dpv_lon_n, half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dpv_lat_t, full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%dpv_lat_n, full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%ke       , full_lon=.true., full_lat=.true.)

  end subroutine state_init

  subroutine state_update_m_lon_m_lat(this)

    class(state_type), intent(inout) :: this

    integer i, j

    do j = this%mesh%full_lat_ibeg_no_pole, this%mesh%full_lat_iend_no_pole
      do i = this%mesh%half_lon_ibeg, this%mesh%half_lon_iend
        this%m_lon(i,j) = (this%mesh%lon_edge_left_area (j) * this%gd(i,  j) + &
                            this%mesh%lon_edge_right_area(j) * this%gd(i+1,j)   &
                           ) / this%mesh%lon_edge_area(j) / g
      end do
    end do

    do j = this%mesh%half_lat_ibeg_no_pole, this%mesh%half_lat_iend_no_pole
      do i = this%mesh%full_lon_ibeg, this%mesh%full_lon_iend
#ifdef V_POLE
        this%m_lat(i,j) = (this%mesh%lat_edge_up_area  (j) * this%gd(i,j  ) + &
                            this%mesh%lat_edge_down_area(j) * this%gd(i,j-1)   &
                           ) / this%mesh%lat_edge_area(j) / g
#else
        this%m_lat(i,j) = (this%mesh%lat_edge_up_area  (j) * this%gd(i,j+1) + &
                            this%mesh%lat_edge_down_area(j) * this%gd(i,j  )   &
                           ) / this%mesh%lat_edge_area(j) / g
#endif
      end do
    end do

  end subroutine state_update_m_lon_m_lat

  subroutine state_update_m_vtx(this)

    class(state_type), intent(inout) :: this

    integer i, j
    real(r8) pole

    do j = this%mesh%half_lat_ibeg_no_pole, this%mesh%half_lat_iend_no_pole
      do i = this%mesh%half_lon_ibeg, this%mesh%half_lon_iend
#ifdef V_POLE
        this%m_vtx(i,j) = (                                                       &
          (this%gd(i,j-1) + this%gd(i+1,j-1)) * this%mesh%subcell_area(2,j-1) + &
          (this%gd(i,j  ) + this%gd(i+1,j  )) * this%mesh%subcell_area(1,j  )   &
        ) / this%mesh%vertex_area(j) / g
#else
        this%m_vtx(i,j) = (                                                       &
          (this%gd(i,j  ) + this%gd(i+1,j  )) * this%mesh%subcell_area(2,j  ) + &
          (this%gd(i,j+1) + this%gd(i+1,j+1)) * this%mesh%subcell_area(1,j+1)   &
        ) / this%mesh%vertex_area(j) / g
#endif
      end do
    end do
#ifdef V_POLE
    if (this%mesh%has_south_pole()) then
      j = this%mesh%half_lat_ibeg
      pole = 0.0_r8
      do i = this%mesh%full_lon_ibeg, this%mesh%full_lon_iend
        pole = pole + this%gd(i,j)
      end do
      call zonal_sum(pole)
      pole = pole / this%mesh%num_half_lon / g
      do i = this%mesh%half_lon_ibeg, this%mesh%half_lon_iend
        this%m_vtx(i,j) = pole
      end do
    end if
    if (this%mesh%has_north_pole()) then
      j = this%mesh%half_lat_iend
      pole = 0.0_r8
      do i = this%mesh%full_lon_ibeg, this%mesh%full_lon_iend
        pole = pole + this%gd(i,j-1)
      end do
      call zonal_sum(pole)
      pole = pole / this%mesh%num_half_lon / g
      do i = this%mesh%half_lon_ibeg, this%mesh%half_lon_iend
        this%m_vtx(i,j) = pole
      end do
    end if
#endif

  end subroutine state_update_m_vtx

  subroutine state_update_mf_lon_n_mf_lat_n(this)

    class(state_type), intent(inout) :: this

    integer i, j

    do j = this%mesh%full_lat_ibeg_no_pole, this%mesh%full_lat_iend_no_pole
      do i = this%mesh%half_lon_ibeg, this%mesh%half_lon_iend
        this%mf_lon_n(i,j) = this%m_lon(i,j) * this%u(i,j)
      end do
    end do
    call fill_halo(this%mesh, this%mf_lon_n)

    do j = this%mesh%half_lat_ibeg_no_pole, this%mesh%half_lat_iend_no_pole
      do i = this%mesh%full_lon_ibeg, this%mesh%full_lon_iend
        this%mf_lat_n(i,j) = this%m_lat(i,j) * this%v(i,j)
      end do
    end do
    call fill_halo(this%mesh, this%mf_lat_n)

  end subroutine state_update_mf_lon_n_mf_lat_n

  subroutine state_update_mf_lon_t_mf_lat_t(this)

    class(state_type), intent(inout) :: this

    integer i, j

    do j = this%mesh%full_lat_ibeg_no_pole, this%mesh%full_lat_iend_no_pole
      do i = this%mesh%half_lon_ibeg, this%mesh%half_lon_iend
#ifdef V_POLE
        this%mf_lon_t(i,j) = this%mesh%full_tangent_wgt(1,j) * (this%mf_lat_n(i,j  ) + this%mf_lat_n(i+1,j  )) + &
                             this%mesh%full_tangent_wgt(2,j) * (this%mf_lat_n(i,j+1) + this%mf_lat_n(i+1,j+1))
#else
        this%mf_lon_t(i,j) = this%mesh%full_tangent_wgt(1,j) * (this%mf_lat_n(i,j-1) + this%mf_lat_n(i+1,j-1)) + &
                             this%mesh%full_tangent_wgt(2,j) * (this%mf_lat_n(i,j  ) + this%mf_lat_n(i+1,j  ))
#endif
      end do
    end do
    call fill_halo(this%mesh, this%mf_lon_t)

    do j = this%mesh%half_lat_ibeg_no_pole, this%mesh%half_lat_iend_no_pole
      do i = this%mesh%full_lon_ibeg, this%mesh%full_lon_iend
#ifdef V_POLE
        this%mf_lat_t(i,j) = this%mesh%half_tangent_wgt(1,j) * (this%mf_lon_n(i-1,j-1) + this%mf_lon_n(i,j-1)) + &
                             this%mesh%half_tangent_wgt(2,j) * (this%mf_lon_n(i-1,j  ) + this%mf_lon_n(i,j  ))
#else
        this%mf_lat_t(i,j) = this%mesh%half_tangent_wgt(1,j) * (this%mf_lon_n(i-1,j  ) + this%mf_lon_n(i,j  )) + &
                             this%mesh%half_tangent_wgt(2,j) * (this%mf_lon_n(i-1,j+1) + this%mf_lon_n(i,j+1))
#endif
      end do
    end do
    call fill_halo(this%mesh, this%mf_lat_t)

  end subroutine state_update_mf_lon_t_mf_lat_t

  subroutine state_update_pv_vtx(this)

    class(state_type), intent(inout) :: this

    type(mesh_type), pointer :: mesh
    real(r8) pole
    integer i, j

    mesh => this%mesh

    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
        this%pv(i,j) = (                                                              &
          (                                                                           &
            this%u(i  ,j-1) * mesh%de_lon(j-1) - this%u(i  ,j  ) * mesh%de_lon(j  ) + &
            this%v(i+1,j  ) * mesh%de_lat(j  ) - this%v(i  ,j  ) * mesh%de_lat(j  )   &
          ) / mesh%vertex_area(j) + mesh%half_f(j)                                    &
        ) / this%m_vtx(i,j)
#else
        this%pv(i,j) = (                                                              &
          (                                                                           &
            this%u(i  ,j  ) * mesh%de_lon(j  ) - this%u(i  ,j+1) * mesh%de_lon(j+1) + &
            this%v(i+1,j  ) * mesh%de_lat(j  ) - this%v(i  ,j  ) * mesh%de_lat(j  )   &
          ) / mesh%vertex_area(j) + mesh%half_f(j)                                    &
        ) / this%m_vtx(i,j)
#endif
      end do
    end do
#ifdef V_POLE
    if (mesh%has_south_pole()) then
      j = mesh%half_lat_ibeg
      pole = 0.0_r8
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        pole = pole - this%u(i,j) * mesh%de_lon(j)
      end do
      call zonal_sum(pole)
      pole = pole / mesh%num_half_lon / mesh%vertex_area(j)
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        this%pv(i,j) = (pole + mesh%half_f(j)) / this%m_vtx(i,j)
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%half_lat_iend
      pole = 0.0_r8
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        pole = pole + this%u(i,j-1) * mesh%de_lon(j-1)
      end do
      call zonal_sum(pole)
      pole = pole / mesh%num_half_lon / mesh%vertex_area(j)
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        this%pv(i,j) = (pole + mesh%half_f(j)) / this%m_vtx(i,j)
      end do
    end if
#else
    if (pv_pole_stokes) then
      ! Special treatment of vorticity around Poles
      if (mesh%has_south_pole()) then
        j = mesh%half_lat_ibeg
        pole = 0.0_r8
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          pole = pole - this%u(i,j+1) * mesh%de_lon(j+1)
        end do
        call zonal_sum(pole)
        pole = pole / mesh%num_half_lon / mesh%vertex_area(j)
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          this%pv(i,j) = (pole + mesh%half_f(j)) / this%m_vtx(i,j)
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%half_lat_iend
        pole = 0.0_r8
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          pole = pole + this%u(i,j) * mesh%de_lon(j)
        end do
        call zonal_sum(pole)
        pole = pole / mesh%num_half_lon / mesh%vertex_area(j)
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          this%pv(i,j) = (pole + mesh%half_f(j)) / this%m_vtx(i,j)
        end do
      end if
    end if
#endif
    call fill_halo(mesh, this%pv)

  end subroutine state_update_pv_vtx

  subroutine state_update_ke(this)

    class(state_type), intent(inout) :: this

    type(mesh_type), pointer :: mesh
    integer i, j
    real(r8) pole

    mesh => this%mesh

    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        this%ke(i,j) = (mesh%lon_edge_right_area(j  ) * this%u(i-1,j  )**2 + &
                         mesh%lon_edge_left_area (j  ) * this%u(i  ,j  )**2 + &
#ifdef V_POLE
                         mesh%lat_edge_up_area   (j  ) * this%v(i  ,j  )**2 + &
                         mesh%lat_edge_down_area (j+1) * this%v(i  ,j+1)**2   &
#else
                         mesh%lat_edge_up_area   (j-1) * this%v(i  ,j-1)**2 + &
                         mesh%lat_edge_down_area (j  ) * this%v(i  ,j  )**2   &
#endif
                        ) / mesh%cell_area(j)
      end do
    end do
#ifndef V_POLE
    ! Note: lat_edge_down_area and lat_edge_up_area at the Poles is the same as cell_area.
    if (mesh%has_south_pole()) then
      j = mesh%full_lat_ibeg
      pole = 0.0d0
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        pole = pole + this%v(i,j)**2
      end do
      call zonal_sum(pole)
      pole = pole / mesh%num_full_lon * 0.5_r8
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        this%ke(i,j) = pole
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_iend
      pole = 0.0d0
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        pole = pole + this%v(i,j-1)**2
      end do
      call zonal_sum(pole)
      pole = pole / mesh%num_full_lon * 0.5_r8
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        this%ke(i,j) = pole
      end do
    end if
#endif
    call fill_halo(mesh, this%ke)

  end subroutine state_update_ke

  subroutine state_update_all(this)

    class(state_type), intent(inout) :: this

    call this%update_m_lon_m_lat()
    call this%update_m_vtx()
    call this%update_mf_lon_n_mf_lat_n()
    call this%update_mf_lon_t_mf_lat_t()
    call this%update_pv_vtx()
    call this%update_ke()

  end subroutine state_update_all

  subroutine state_clear(this)

    class(state_type), intent(inout) :: this

    if (allocated(this%u        )) deallocate(this%u        )
    if (allocated(this%v        )) deallocate(this%v        )
    if (allocated(this%gd       )) deallocate(this%gd       )
    if (allocated(this%pv       )) deallocate(this%pv       )
    if (allocated(this%m_vtx    )) deallocate(this%m_vtx    )
    if (allocated(this%m_lon    )) deallocate(this%m_lon    )
    if (allocated(this%m_lat    )) deallocate(this%m_lat    )
    if (allocated(this%mf_lon_n )) deallocate(this%mf_lon_n )
    if (allocated(this%mf_lat_n )) deallocate(this%mf_lat_n )
    if (allocated(this%mf_lat_t )) deallocate(this%mf_lat_t )
    if (allocated(this%mf_lon_t )) deallocate(this%mf_lon_t )
    if (allocated(this%pv_lon   )) deallocate(this%pv_lon   )
    if (allocated(this%pv_lat   )) deallocate(this%pv_lat   )
    if (allocated(this%dpv_lon_t)) deallocate(this%dpv_lon_t)
    if (allocated(this%dpv_lon_n)) deallocate(this%dpv_lon_n)
    if (allocated(this%dpv_lat_t)) deallocate(this%dpv_lat_t)
    if (allocated(this%dpv_lat_n)) deallocate(this%dpv_lat_n)
    if (allocated(this%ke       )) deallocate(this%ke       )

  end subroutine state_clear

  subroutine state_final(this)

    type(state_type), intent(inout) :: this

    call this%clear()

  end subroutine state_final

end module state_mod
