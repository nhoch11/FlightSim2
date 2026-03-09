module controller_m
    use hoch_m
    use jsonx_m
    use linalg_mod
    use connection_m
    use atmosphere_m

    implicit none
    
    type pid_type 
        character(len=:), allocatable :: name
        character(len=:), allocatable :: units
        real :: KP, KI, KD
        real :: error, prev_error, error_int, error_deriv
        real :: prev_time, prev_ans, update_rate
        real :: bias
        real, allocatable :: limit(:)
        real :: display_units = 1.0
        logical :: dynamic_pressure_schedule
    end type pid_type

    type controller_type
        ! type(pid_type), allocatable :: pid_controllers(:)
        type(pid_type) :: p_da, q_de, r_dr, bank_p, gamma_q, V_tau
        integer :: num_pid
        type(connection) :: pilot_conn
        logical :: running = .false.
        logical :: set_trim_bias
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
        call jsonx_get(j_main, "set_trim_bias", this%set_trim_bias, .false.)
        call jsonx_get(j_main, "connections", j_connections)
        call jsonx_get(j_connections, "receive_pilot_commands", j_temp)
        call this%pilot_conn%init(j_temp, time=0.0)
        
        ! initialize PID controllers
        write(*,*) "     Initializing PID Controllers"
        call jsonx_get(j_main, "PID", j_pid)
        ! this%num_pid = json_value_count(j_pid)
        ! allocate(this%pid_controllers(this%num_pid))
        
        ! old, just rates
        ! do i=1,this%num_pid
        !     call json_value_get(j_pid, i, j_temp)
        !     call pid_init(this%pid_controllers(i), j_temp)
        ! end do
        
        call jsonx_get(j_pid, "p->aileron", j_temp)
        call pid_init(this%p_da, j_temp)

        call jsonx_get(j_pid, "q->elevator", j_temp)
        call pid_init(this%q_de, j_temp)
        
        call jsonx_get(j_pid, "r->rudder", j_temp)
        call pid_init(this%r_dr, j_temp)
        
        call jsonx_get(j_pid, "bank->p", j_temp)
        call pid_init(this%bank_p, j_temp)
        
        call jsonx_get(j_pid, "gamma->q", j_temp)
        call pid_init(this%gamma_q, j_temp)
        
        call jsonx_get(j_pid, "V->throttle", j_temp)
        call pid_init(this%V_tau, j_temp)
        
        write(*,*) "     controller initialization complete."

    end subroutine controller_init


    ! function controller_update(this, states, time) result(ans)
    !     implicit none
    !     type(controller_type), intent(inout) :: this
    !     real, intent(in) :: states(21), time
    !     real :: ans(4)
    !     real :: omega_command(3), omega_actual(3)
    !     real :: Z, Temp, P, rho, a, mu, dynamic_pressure

    !     ! hard coded trim solution
    !     ! ans = [0.0, -9.086019165449*PI/180.0, 0.0, 0.070182002357] ! 12.6.1
    !     ! ans = [0.0, 0.0, 0.0, 0.0760529852986603] ! 12.6.2
        
    !     omega_actual = states(4:6)
        
    !     ! calculate rho for dynamic pressure
    !     call std_atm_English(-states(9), Z, Temp, P, rho, a, mu)
    !     dynamic_pressure = 0.5*rho*(states(1)**2 + states(2)**2 + states(3)**2)

    !     omega_command = this%pilot_conn%recv([time], time)

    !     ans(1) = pid_get_command(this%pid_controllers(1), omega_command(1), omega_actual(1), time, dynamic_pressure)
    !     ans(2) = pid_get_command(this%pid_controllers(2), omega_command(2), omega_actual(2), time, dynamic_pressure)
    !     ans(3) = pid_get_command(this%pid_controllers(3), omega_command(3), omega_actual(3), time, dynamic_pressure)

    ! end function controller_update
    
    function controller_update(this, states, time) result(ans)
        implicit none
        type(controller_type), intent(inout) :: this
        real, intent(in) :: states(21), time
        real :: ans(4)
        real :: pilot_command(3)
        real :: u, v, w, p, q, r, eul(3), sp, cp, st, ct, V_mag, g, gamma
        real :: p_sp, q_sp, r_sp
        real :: bank_sp, gamma_sp, V_sp
        real :: Z, Temp, Pres, rho, a, mu, dynamic_pressure
        
        ! hard coded trim solution
        ! ans = [0.0, -9.086019165449*PI/180.0, 0.0, 0.070182002357] ! 12.6.1
        ! ans = [0.0, 0.0, 0.0, 0.0760529852986603] ! 12.6.2
        
        u = states(1)
        v = states(2)
        w = states(3)
        p = states(4)
        q = states(5)
        r = states(6)
        eul = quat_to_euler(states(10:13))
        sp = sin(eul(1))
        cp = cos(eul(1))
        st = sin(eul(2))
        ct = cos(eul(2))
        V_mag = sqrt(u**2 + v**2 + w**2)
        gamma = asin((u*st - (v*sp + w*cp)*ct)/V_mag)
        
        ! calculate rho for dynamic pressure
        call std_atm_English(-states(9), Z, Temp, Pres, rho, a, mu)
        dynamic_pressure = 0.5*rho*V_mag**2
        g = gravity_English(-states(9))

        ! ! 12.6.1 ------------------------------------------------------------------------
        ! pilot_command = this%pilot_conn%recv([time], time)
        
        ! ! Outer loop controllers to get p,q,r set points
        ! p_sp = pilot_command(1)
        
        ! ! inner loop to get control deflection set points
        ! ans(1) = pid_get_command(this%p_da, p_sp, p, time, dynamic_pressure)
        ! ans(2) = -9.086019165449*PI/180.0
        ! ans(3) = 0.0
        ! ans(4) = 0.070182002357   
        

        ! ! 12.6.2 ------------------------------------------------------------------------
        ! pilot_command = this%pilot_conn%recv([time], time)
        
        ! ! Outer loop controllers to get p,q,r set points
        ! p_sp = pilot_command(1)
        ! q_sp = pilot_command(2)
        ! r_sp = pilot_command(3)
        
        ! ! inner loop to get control deflection set points
        ! ans(1) = pid_get_command(this%p_da, p_sp, p, time, dynamic_pressure)
        ! ans(2) = pid_get_command(this%q_de, q_sp, q, time, dynamic_pressure)
        ! ans(3) = pid_get_command(this%r_dr, r_sp, r, time, dynamic_pressure)
        ! ans(4) = 0.076052985299   
        
        
        ! 12.6.3 (normal operation)-------------------------------------------------------
        pilot_command = this%pilot_conn%recv([time], time)
        bank_sp = pilot_command(1)*PI/180.0
        gamma_sp = pilot_command(2)*PI/180.0
        V_sp = pilot_command(3)

        ! Outer loop controllers to get p,q,r set points
        p_sp = pid_get_command(this%bank_p, bank_sp, eul(1), time, dynamic_pressure)
        q_sp = pid_get_command(this%gamma_q, gamma_sp, gamma, time, dynamic_pressure)
        r_sp = (g*sp*ct + p*w)/u  ! neglects gravity relief
        
        ! inner loop to get control deflection set points
        ans(1) = pid_get_command(this%p_da, p_sp, p, time, dynamic_pressure)
        ans(2) = pid_get_command(this%q_de, q_sp, q, time, dynamic_pressure)
        ans(3) = pid_get_command(this%r_dr, r_sp, r, time, dynamic_pressure)
        ans(4) = pid_get_command(this%V_tau, V_sp, V_mag, time, dynamic_pressure)
        
    end function controller_update

        
    subroutine pid_init(this, j_pid)
        implicit none
        type(pid_type), intent(inout) :: this
        type(json_value), pointer :: j_pid

        this%name = j_pid%name
        write(*,*) "          Initializing PID for ", this%name
        call jsonx_get(j_pid, "update_rate[hz]", this%update_rate)
        call jsonx_get(j_pid, "kp", this%kp)
        call jsonx_get(j_pid, "ki", this%ki)
        call jsonx_get(j_pid, "kd", this%kd)
        call jsonx_get(j_pid, "dynamic_pressure_schedule", this%dynamic_pressure_schedule, .false.)

        this%bias = 0.0
        
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
            ans = this%kp*this%prev_error
            if (this%dynamic_pressure_schedule) then
                ans = ans/dynamic_pressure
            end if
            ans = ans + this%bias
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

            ans = (this%kp*this%error + this%ki*this%error_int + this%kd*this%error_deriv)
            if (this%dynamic_pressure_schedule) ans = ans/dynamic_pressure
            ans = ans + this%bias

            this%prev_error = this%error
            this%prev_time = time
            this%prev_ans = ans
        else
            ! do nothing
            ans = this%prev_ans
        end if
    end function pid_get_command

end module controller_m