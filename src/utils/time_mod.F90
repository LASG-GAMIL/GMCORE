module time_mod

  use datetime_mod
  use timedelta_mod
  use hash_table_mod

  implicit none

  private

  public time_init
  public time_reset_start_time
  public time_swap_indices
  public time_advance
  public time_elapsed_seconds
  public time_is_finished
  public time_add_alert
  public time_is_alerted

  public curr_time
  public start_time_format
  public curr_time_format
  public time_step
  public old_time_idx
  public new_time_idx

  type alert_type
    type(timedelta_type) period
    type(datetime_type) last_time
    logical :: ring = .false.
  end type alert_type

  type(datetime_type) start_time
  type(datetime_type) end_time
  type(datetime_type) curr_time
  type(timedelta_type) time_step_size
  real(8) elapsed_seconds
  type(hash_table_type) alerts
  integer time_step
  integer old_time_idx
  integer new_time_idx
  character(30) start_time_format
  character(30) curr_time_format

contains

  subroutine time_init(start_time_in, end_time_in, time_step_size_in)

    integer, intent(in) :: start_time_in(:)
    integer, intent(in) :: end_time_in(:)
    real, intent(in) :: time_step_size_in

    start_time = create_datetime(year=start_time_in(1), month=start_time_in(2), day=start_time_in(3), &
      hour=start_time_in(4), minute=start_time_in(5))
    end_time = create_datetime(year=end_time_in(1), month=end_time_in(2), day=end_time_in(3), &
      hour=end_time_in(4), minute=end_time_in(5))

    time_step = 0
    elapsed_seconds = 0
    old_time_idx = 1
    new_time_idx = 2
    time_step_size = timedelta(seconds=time_step_size_in)

    curr_time = start_time

    start_time_format = start_time%format('%Y-%m-%dT%H_%M_%S')
    curr_time_format  = curr_time %format('%Y-%m-%dT%H_%M_%S')

    alerts = hash_table()

  end subroutine time_init

  subroutine time_reset_start_time(time)

    type(datetime_type), intent(in) :: time

    start_time = time
    curr_time = start_time

    start_time_format = start_time%format('%Y-%m-%dT%H_%M_%S')
    curr_time_format  = curr_time %format('%Y-%m-%dT%H_%M_%S')

  end subroutine time_reset_start_time

  subroutine time_swap_indices(i, j)

    integer, intent(inout) :: i
    integer, intent(inout) :: j

    integer tmp

    tmp = i
    i = j
    j = tmp

  end subroutine time_swap_indices

  subroutine time_advance()
    
    type(hash_table_iterator_type) iter

    ! Update alerts.
    iter = hash_table_iterator(alerts)
    do while (.not. iter%ended())
      select type (alert => iter%value)
      type is (alert_type)
        if (alert%ring) then
          alert%last_time = curr_time
          alert%ring = .false.
        end if
      end select
      call iter%next()
    end do

    call time_swap_indices(old_time_idx, new_time_idx)

    time_step = time_step + 1
    elapsed_seconds = elapsed_seconds + time_step_size%total_seconds()
    curr_time = curr_time + time_step_size
    curr_time_format = curr_time%isoformat()

  end subroutine time_advance

  real(8) function time_elapsed_seconds() result(res)

    res = elapsed_seconds

  end function time_elapsed_seconds

  logical function time_is_finished() result(res)

    res = curr_time >= end_time

  end function time_is_finished

  subroutine time_add_alert(name, months, days, hours, minutes, seconds)

    character(*), intent(in) :: name
    real, intent(in), optional :: months
    real, intent(in), optional :: days
    real, intent(in), optional :: hours
    real, intent(in), optional :: minutes
    real, intent(in), optional :: seconds

    real months_
    real days_
    real hours_
    real minutes_
    real seconds_
    type(alert_type) alert

    if (present(months)) then
      months_ = months
    else
      months_ = 0.0
    end if
    if (present(days)) then
      days_ = days
    else
      days_ = 0.0
    end if
    if (present(hours)) then
      hours_ = hours
    else
      hours_ = 0.0
    end if
    if (present(minutes)) then
      minutes_ = minutes
    else
      minutes_ = 0.0
    end if
    if (present(seconds)) then
      seconds_ = seconds
    else
      seconds_ = 0.0
    end if

    alert%period = timedelta(months=months_, days=days_, hours=hours_, minutes=minutes_, seconds=seconds_)
    alert%last_time = start_time
    call alerts%insert(trim(name), alert)

  end subroutine time_add_alert

  function time_is_alerted(name) result(res)

    character(*), intent(in) :: name
    logical res

    type(alert_type), pointer :: alert => null()
    type(datetime_type) time

    alert => get_alert(name)
    if (associated(alert)) then
      time = alert%last_time + alert%period
      if (time <= curr_time) then
        alert%ring = .true.
        res = .true.
      else
        res = .false.
      end if
    else
      res = .false.
    end if

  end function time_is_alerted

  function get_alert(name) result(res)

    character(*), intent(in) :: name
    type(alert_type), pointer :: res

    class(*), pointer :: value

    if (alerts%hashed(name)) then
      value => alerts%value(name)
      select type (value)
      type is (alert_type)
        res => value
      end select
    else
      nullify(res)
    end if

  end function get_alert

end module time_mod
