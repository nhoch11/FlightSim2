module sim_m

    use vehicle_m
    use jsonx_m

    implicit none
    
    ! old time settings
    real :: t0, dt, tf
    real :: integrated_time, actual_time
    logical :: real_time
    real :: cpu_start_time, cpu_end_time
    real :: count_rate
    integer :: num_vehicles
    type(vehicle_t),dimension(:),allocatable :: vehicles

contains


    subroutine init(input_file)
        implicit none
        character(len=*), intent(in) :: input_file

        ! read json file and get all simulation parameters like this
        call jsonx_load(input_file,j_main)
        
        call jsonx_get(j_main, "simulation.time_step[s]", dt, 0.0)
        if (abs(dt) <= TOLERANCE) then
            real_time = .true.
        else 
            real_time = .false.
        end if 
        call jsonx_get(j_main, "simulation.total_time[s]", tf)
        call jsonx_get(j_main, "simulation.verbose", verbose, .false.)
        call jsonx_get(j_main, "simulation.verbose2", verbose2, .false.)

        ! initialize each vehicle
        do i=1,num_vehicles
            vehicles(i)%vehicle_init(vehicle dictionary)
        end do

    end subroutine


    subroutine run()
        
        implicit none
        
        real :: t
        real, dimension(13) :: y, y1, y_init
        real, dimension(14) :: s
        real, dimension(3) :: eul
        real, dimension(5) :: junk
        integer :: ios
        real :: time1, time2

        call mass_inertia()
        call std_atm_English(0.0, junk(1), junk(2), junk(3), rho0, junk(4), junk(5))

        if (init_type == "state") then

            y_init(1)  = V0*cos(alpha0)*cos(beta0)   ! u
            y_init(2)  = V0*sin(beta0)               ! v
            y_init(3)  = V0*sin(alpha0)*cos(beta0)   ! w
            y_init(4)  = p0              ! p
            y_init(5)  = q0              ! q
            y_init(6)  = r0              ! r
            y_init(7)  = 0.0             ! xf
            y_init(8)  = 0.0             ! yf
            y_init(9)  = zf0             ! zf
            y_init(10) = phi0            ! phi
            y_init(11) = theta0          ! theta
            y_init(12) = psi0            ! psi

            !  convert phi, theta, and psi to a quat
            y_init(10:13) = euler_to_quat(y_init(10:12))

            ! populate controls with initial values
            controls(1) = da0
            controls(2) = de0
            controls(3) = dr0
            controls(4) = tau0

        else ! call trim solver
            y_init = trim_solver()
            trimming = .false.
        end if
        ! stop ! just for debugging trim
        
        !  normalize the quat
        call quat_norm(y_init(10:13))
        
        if (verbose) then
            write(*,*) ""
            write(*,*) "y_init = ", y_init
        end if
        ! stop



        t = t0
        y = y_init
        if (real_time) then
            ! call cpu_time(time1)
            ! call system_clock(time1, count_rate)
            ! write(*,*) "count rate = ", count_rate
            time1 = get_time()
            write(*,*) "     time1 = ", time1
            y1 = rk4(t,y, dt)
            call quat_norm(y1(10:13))
            ! call system_clock(time2)
            time2 = get_time()
            write(*,*) "     time2 = ", time2
            ! dt = real(time2 - time2)/count_rate
            dt = real(time2 - time1)
            ! if (dt <= TOLERANCE) dt = 0.00001 / real(count_rate)
            t = t0
            y = y_init

        end if 


        ! Open output file (unit number 10 here, adjust as needed)
        open(unit=10, file=trim(output_file), status="replace", action="write", iostat=ios)
        if (ios /= 0) then
            write(*,*) "Error opening file!"
            stop
        end if

        write(10,'(15(1x,A20))') "time[s]", "dt[s]","u[ft/s]", "v[ft/s]", "w[ft/s]", "p[rad/s]", "q[rad/s]", "r[rad/s]", &
        "xf[ft]", "yf[ft]", "zf[ft]", "e0", "ex", "ey", "ez"

        write(10, '(15(1x,ES20.12))') t, dt, y(1:13)

        if (real_time) then
            ! call system_clock(cpu_start_time)
            cpu_start_time = get_time()
            time1 = cpu_start_time
            integrated_time = 0.0
        end if
        
        do while (t < tf)

            do i=1,num_vehicles
                vehicles(i)%vehicle_tick_state(t, time, dt)
            end do
        end do 
        
        close(10)
        
        if (real_time) then
            ! call system_clock(cpu_end_time)
            cpu_end_time = get_time()
            ! actual_time = real(cpu_end_time - cpu_start_time)/count_rate
            actual_time = cpu_end_time - cpu_start_time
            write(*,*) "       Integrated time[s] = ", integrated_time
            write(*,*) "           Actual time[s] = ", actual_time
            write(*,*) "  Total error in time [s] = ", integrated_time - actual_time
        end if

    end subroutine


end module sim_m