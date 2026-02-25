module controller_m
    use hoch_m
    ! use atmosphere_m
    use jsonx_m
    use linalg_mod
    use connection_m

    implicit none
    
    type pid_type 
        character(len=:), allocatable :: name
        character(len=:), allocatable :: units
        real :: KP, KI, KD
        real :: error, prev_error, error_int, error_deriv
        real :: prev_time, prev_ans, update_rate
        real, allocatable :: limit(:)
        real :: display_units = 1.0
    end type pid_type

    type controller_type
        type(pid_type), allocatable :: pid_controllers(:)
        integer :: num_pid
        type(connection) :: pilot_conn
        logical :: running = .false.
    end type controller_type


contains

    subroutine controller_init(this, j_main)
        implicit none
        type(controller_type), intent(inout) :: this
        type(json_value), pointer :: j_main, j_connections, j_pid, j_temp
        integer :: i 
        logical :: found

        write(*,*) "     Initializing controller..."
        this%running = .true.
        
        ! initialize controller connections
        call jsonx_get(j_main, "connections", j_connections)
        call jsonx_get(j_connections, "receive_pilot_commands", j_temp)
        call this%pilot_conn%init(j_temp, time=0.0)
        
        ! initialize PID controllers
        call jsonx_get(j_main, "PID", j_pid)
        this%num_pid = json_value_count(j_pid)
        allocate(this%pid_controllers(this%num_pid))
        
        do i=1,this%num_pid
            call json_value_get(j_pid, i, j_temp)
            call pid_init(this%pid_controllers(i), j_temp)
        end do
        
        write(*,*) "     controller initialization complete."
    end subroutine controller_init


    function controller_update(this, states, time) result(ans)
        implicit none
        type(controller_type), intent(inout) :: this
        real, intent(in) :: states(21), time
        real :: ans(4)
        real :: omega_command(3), omega_actual(3)
        real :: Z, Temp, P, rho, a, mu, dynamic_pressure

        ! hard coded trim solution
        ans = [0.0, -9.086019165449*PI/180.0, 0.0, 0.070182002357]

        omega_actual = states(4:6)

        ! calculate rho for dynamic pressure
        call std_atm_English(-states(9), Z, Temp, P, rho, a, mu)
        dynamic_pressure = 0.5*rho*(states(1)**2 + states(2)**2 + states(3)**2)
        
        omega_command = this%pilot_conn%recv([time], time)
        ans(1) = pid_get_command(this%pid_controllers(1), omega_command(1), omega_actual(1), time, dynamic_pressure)
        ! ans(2) = pid_get_command(this%pid_controllers(2), omega_command(2), omega_actual(2), time, dynamic_pressure)
        ! ans(3) = pid_get_command(this%pid_controllers(3), omega_command(3), omega_actual(3), time, dynamic_pressure)

    end function controller_update


    subroutine pid_init(this, j_pid)
        implicit none
        type(pid_type), intent(inout) :: this
        type(json_value), pointer :: j_pid

        this%name = j_pid%name
        write(*,*) "          Initializing PID for ", this%name
        call jsonx_get(j_pid, "update_rate[hz]", this%update_rate)
        call jsonx_get(j_pid, "KP", this%KP)
        call jsonx_get(j_pid, "KI", this%KI)
        call jsonx_get(j_pid, "KD", this%KD)
        
        this%prev_error = 0.0
        this%error_int = 0.0
        this%prev_ans = 0.0
        this%prev_time = -1.0


        call jsonx_get(j_pid, "units", this%units, "none")
        if (this%units == "deg") then
            this%display_units = 180.0/PI
        else
            this%display_units = 1.0
        end if
        call jsonx_get(j_pid, "limits", this%limit, 0.0, 2)
        this%limit(:) = this%limit(:)/this%display_units

    end subroutine pid_init

    function pid_get_command(this, commanded, actual, time, dynamic_pressure) result(ans)
        implicit none
        type(pid_type), intent(inout) :: this
        real, intent(in) :: commanded, actual, time, dynamic_pressure
        real :: dt, ans

        ! if at first time step
        if (this%prev_time < 0.0) then
            this%prev_time = time
            this%prev_error = commanded - actual
            this%error_int = 0.0
            this%error_deriv = 0.0
            ans = this%KP*this%prev_error/dynamic_pressure
            this%prev_ans = ans
            return
        end if

        if (time - this%prev_time >= 1.0/this%update_rate - TOLERANCE) then 
            ! update controller

            dt = time - this%prev_time
            
            ! proportional 
            this%error = commanded - actual

            ! integral
            if ((this%prev_ans > this%limit(1)) .and. (this%prev_ans < this%limit(2))) then
                this%error_int = this%error_int + 0.5*(this%prev_error + this%error)*dt
            else
                write(*,*) this%name, " PID controller saturated at ", this%prev_ans*this%display_units,", using integrator clamping"
            end if

            ! derivataive
            if (dt> TOLERANCE) this%error_deriv = (this%error - this%prev_error)/dt

            ans = (this%KP*this%error + this%KI*this%error_int + this%KD*this%error_deriv)/dynamic_pressure

            this%prev_error = this%error
            this%prev_time = time
            this%prev_ans = ans
        else
            ! do nothing
            ans = this%prev_ans
        end if
    end function pid_get_command

end module controller_m