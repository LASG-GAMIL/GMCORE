module time_mod

  use datetime
  use string
  use container
  use flogger
  use const_mod

  implicit none

  private

  public time_init
  public time_final
  public time_reset_start_time
  public time_swap_indices
  public time_advance
  public time_fast_forward
  public time_elapsed_seconds
  public time_is_first_step
  public time_is_finished
  public time_add_alert
  public time_has_alert
  public time_is_alerted

  public curr_time
  public start_time_str
  public curr_time_str
  public elapsed_seconds
  public time_step
  public old_time_idx
  public new_time_idx

  type alert_type
    type(timedelta_type) period
    type(datetime_type) last_time
    logical :: ring = .false.
  end type alert_type

  ! Namelist parameters
  integer, public :: start_time_array(5)  = 0
  integer, public :: end_time_array(5)    = 0
  real(8), public :: run_hours            = 0
  real(8), public :: run_days             = 0
  real(8), public :: run_years            = 0

  type(datetime_type) start_time
  type(datetime_type) end_time
  type(datetime_type) curr_time
  type(timedelta_type) dt
  real(8) elapsed_seconds
  type(hash_table_type) alerts
  integer time_step
  integer old_time_idx
  integer new_time_idx
  character(30) start_time_str
  character(30) curr_time_str

contains

  subroutine time_init(dt_in_seconds)

    real(8), intent(in) :: dt_in_seconds

    if (sum(start_time_array) > 0) then
      call start_time%init(year=start_time_array(1),  &
                           month=start_time_array(2), &
                           day=start_time_array(3),   &
                           hour=start_time_array(4),  &
                           minute=start_time_array(5))
    else
      call start_time%init(year=1, month=1, day=1, hour=0, minute=0)
    end if
    end_time = start_time
    if (run_years > 0 .or. run_days > 0 .or. run_hours > 0) then
      call end_time%add(years=run_years, days=run_days, hours=run_hours)
    else if (sum(end_time_array) > 0) then
      call end_time%init(year=end_time_array(1),  &
                         month=end_time_array(2), &
                         day=end_time_array(3),   &
                         hour=end_time_array(4),  &
                         minute=end_time_array(5))
    end if

    time_step = 0
    elapsed_seconds = 0
    old_time_idx = 1
    new_time_idx = 2
    call dt%init(seconds=dt_in_seconds)

    curr_time = start_time

    start_time_str = start_time%format('%Y-%m-%dT%H_%M_%S')
    curr_time_str = curr_time%format('%Y-%m-%dT%H_%M_%S')

    call alerts%init()

  end subroutine time_init

  subroutine time_final()

    call alerts%clear()

  end subroutine time_final

  subroutine time_reset_start_time(time)

    type(datetime_type), intent(in) :: time

    start_time = time
    curr_time = start_time

  end subroutine time_reset_start_time

  subroutine time_swap_indices(i, j)

    integer, intent(inout) :: i
    integer, intent(inout) :: j

    integer tmp

    tmp = i
    i = j
    j = tmp

  end subroutine time_swap_indices

  subroutine time_advance(dt_in_seconds)

    real(8), intent(in), optional :: dt_in_seconds

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
    if (present(dt_in_seconds)) then
      elapsed_seconds = elapsed_seconds + dt_in_seconds
      call curr_time%add(seconds=dt_in_seconds)
    else
      elapsed_seconds = elapsed_seconds + dt%total_seconds()
      curr_time = curr_time + dt
    end if
    curr_time_str = curr_time%format('%Y-%m-%dT%H_%M_%S')

  end subroutine time_advance

  subroutine time_fast_forward(time_value, time_units)

    real(8), intent(in) :: time_value
    character(*), intent(in) :: time_units

    type(timedelta_type) skipped_time
    character(30) tmp1, tmp2

    tmp1 = split_string(time_units, ' ', 1)
    tmp2 = split_string(time_units, ' ', 3)

    call curr_time%init(tmp2)
    select case (tmp1)
    case ('minutes')
      call curr_time%add_minutes(time_value)
    case ('hours')
      call curr_time%add_hours(time_value)
    case ('days')
      call curr_time%add_days(time_value)
    case ('seconds')
      call curr_time%add_seconds(time_value)
    case default
      call log_error('Unsupported time units ' // trim(time_units) // '!')
    end select

    skipped_time = curr_time - start_time
    elapsed_seconds = skipped_time%total_seconds()
    curr_time_str = curr_time%format('%Y-%m-%dT%H_%M_%S')

  end subroutine time_fast_forward

  real(8) function time_elapsed_seconds() result(res)

    res = elapsed_seconds

  end function time_elapsed_seconds

  logical function time_is_first_step() result(res)

    res = elapsed_seconds == 0

  end function time_is_first_step

  logical function time_is_finished() result(res)

    res = curr_time >= end_time

  end function time_is_finished

  subroutine time_add_alert(name, months, days, hours, minutes, seconds)

    character(*), intent(in)           :: name
    real(8)     , intent(in), optional :: months
    real(8)     , intent(in), optional :: days
    real(8)     , intent(in), optional :: hours
    real(8)     , intent(in), optional :: minutes
    real(8)     , intent(in), optional :: seconds

    real(8) months_
    real(8) days_
    real(8) hours_
    real(8) minutes_
    real(8) seconds_
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

    call alert%period%init(months=months_, days=days_, hours=hours_, minutes=minutes_, seconds=seconds_)
    alert%last_time = start_time - alert%period
    call alerts%insert(trim(name), alert)

  end subroutine time_add_alert

  function time_has_alert(name) result(res)

    character(*), intent(in) :: name
    logical res

    type(alert_type), pointer :: alert => null()

    alert => get_alert(name)
    res = associated(alert)

  end function time_has_alert

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
