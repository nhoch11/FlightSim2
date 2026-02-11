module sim_m

    use vehicle_m
    use jsonx_m

    implicit none
    
    type(json_value), pointer :: j_main
    type(vehicle_type), dimension(:), allocatable :: vehicles
    integer :: num_vehicles
    character(len=:), allocatable, dimension(:) :: headers
    integer :: print_rate

    ! logical :: save_states, verbose, rk4_verbose
    
    contains
    
    
    subroutine init(input_file)
        implicit none
        character(len=*), intent(in) :: input_file
        type(json_value), pointer :: j_vehicles, j_temp
        integer :: i

        
        headers = [character(len=20) :: &
        " time[s]", " dt[s]", " u[ft/s]", " v[ft/s]", " w[ft/s]", &
        " p[rad/s]", " q[rad/s]", " r[rad/s]", " x[ft]", " y[ft]", " z[ft]",  &
        " e0", " ex", " ey", " ez"]
        
        write(*,*) "Initializing Simulation..."       
        ! read json file and get all simulation parameters 
        call jsonx_load(input_file, j_main)
        
        ! global sim settings
        call jsonx_get(j_main, "simulation.save_states", save_states, .true.)
        call jsonx_get(j_main, "simulation.verbose", verbose, .false.)
        call jsonx_get(j_main, "simulation.verbose2", verbose2, .false.)
        call jsonx_get(j_main, "simulation.rk4_verbose", rk4_verbose, .false.)
        call jsonx_get(j_main, "simulation.gravity_relief", gravity_relief, .true.)
        if (gravity_relief) then
            gravity_relief_factor = 1.0
        else
            gravity_relief_factor = 0.0
        end if
        
        call jsonx_get(j_main, "simulation.geographic_model", geographic_model, "none")
        geographic_model_ID = 0
        if (geographic_model == "sphere") geographic_model_ID = 1
        if (geographic_model == "ellipse") geographic_model_ID = 2

        if (geographic_model_ID > 0) then
            headers = [character(len=20) :: &
             headers, &
             " latitude[deg]", &
             " longitude[deg]", &
             " azimuth[deg]"]
        end if

        write(*,*) "Initializing Vehicles..."
        call jsonx_get(j_main, "vehicles", j_vehicles)
        num_vehicles = json_value_count(j_vehicles)
        allocate(vehicles(num_vehicles))
        
        ! initialize each vehicle
        do i=1,num_vehicles
            call json_value_get(j_vehicles, i, j_temp)
            call vehicle_init(vehicles(i), j_temp)
        end do
        write(*,*)"" 
        write(*,*) "Initialization Complete"

    end subroutine
    
    
    subroutine run()
        
        implicit none
        
        real :: time, dt, tf, stop_time
        real :: time1, time2
        real :: integrated_time, actual_time
        real :: cpu_start_time, cpu_end_time
        real :: count_rate
        real :: eul(3)
        logical :: real_time = .false.
        integer :: i, istep


        call jsonx_get(j_main, "simulation.time_step[s]", dt, 0.0)
        call jsonx_get(j_main, "simulation.end_time[s]", tf)
        call jsonx_get(j_main, "simulation.print_rate", print_rate, 1)

        time = 0.0
        istep = 0

        if (abs(dt) <= TOLERANCE) then
            real_time = .true.

            ! figure out how long it will take to tick all vehicles
            ! write(*,*) "count rate = ", count_rate ! this is for mac
            time1 = get_time()
            do i=1,num_vehicles
                if (vehicles(i)%run_physics) then
                    call vehicle_tick_state(vehicles(i), time, dt)
                end if
            end do
            
            ! set dt
            ! dt = real(time2 - time2)/count_rate  ! this is for mac
            ! if (dt <= TOLERANCE) dt = 0.00001 / real(count_rate) ! this is for mac

            time2 = get_time()
            dt = real(time2 - time1)
            
            ! reset init_state for all vehicles
            do i=1,num_vehicles
                if (vehicles(i)%run_physics) then
                    vehicles(i)%state = vehicles(i)%init_state
                end if
            end do
        end if 
        
        ! write to screen
        if (verbose) then
            ! write(*, '(*(A16))') headers
            write(*,*) ""
            write(*,'(*(A20))') "vehicle name", headers
            do i=1,num_vehicles
                if (vehicles(i)%run_physics) then
                    write(*,'(A19, *(ES20.12))', advance='no') vehicles(i)%name, time, dt, vehicles(i)%state(:)
                    if (geographic_model_ID>0) then
                        vehicles(i)%azimuth_deg = vehicles(i)%init_eul(3)*180.0/PI
                        write(*,'(*(ES20.12))', advance='no') vehicles(i)%latitude_deg, vehicles(i)%longitude_deg, vehicles(i)%azimuth_deg
                    end if
                    write(*,*) "" 
                end if
            end do
        end if
            
        ! get start time
        cpu_start_time = get_time()
        time1 = cpu_start_time
        integrated_time = 0.0

        do while (time <= tf)
            ! tick each vehicle forward one step
            do i=1,num_vehicles
                if (vehicles(i)%run_physics) then
                    call vehicle_tick_state(vehicles(i), time, dt)
                end if 
            end do
            
            time = time + dt
            istep = istep + 1

            ! write to terminal
            if (verbose .and. (modulo(istep, print_rate) == 0)) then
                do i=1,num_vehicles
                    if (vehicles(i)%run_physics) then
                        write(*,'(A19, *(ES20.12))', advance='no') vehicles(i)%name, time, dt, vehicles(i)%state(:)
                        if (geographic_model_ID>0) then
                            eul = quat_to_euler(vehicles(i)%state(10:13))
                            vehicles(i)%azimuth_deg = eul(3)*180.0/PI                            
                            write(*,'(*(ES20.12))', advance='no') vehicles(i)%latitude_deg, vehicles(i)%longitude_deg, vehicles(i)%azimuth_deg
                        end if
                        write(*,*) "" 
                    end if
                end do
            end if

            if (real_time) then
                time2 = get_time()
                dt = time2 - time1
                integrated_time = integrated_time + dt
                time1 = time2
            end if

            ! send data to connection
            ! s(1) = t
            ! s(2:14) = y(1:13)
            ! call graphics%send(s)

        end do 

        cpu_end_time = get_time()
        actual_time = cpu_end_time - cpu_start_time

        if (verbose) then
            write(*,*) ""
            write(*,*) "           Actual time[s] = ", actual_time
            if (real_time) then
                write(*,*) "       Integrated time[s] = ", integrated_time
                write(*,*) "  Total error in time [s] = ", integrated_time - actual_time
            end if 
        end if

    end subroutine


end module sim_m