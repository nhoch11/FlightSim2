program main
    
    use sim_m
    use udp_m, only: udp_initialize, udp_finalize
    ! use f16_566_m

    implicit none
    ! character(100) :: filename
    integer :: start_count, end_count, i_unit, ios
    real :: count_rate_main, runtime
    character(len=100) :: filename

    ! set up
    call udp_initialize()

    ! Start timer
    call system_clock(start_count, count_rate_main)

    ! call init("9.5.2.json")
    call get_command_argument(1,filename)
    call init(filename)
    call run()
    

    ! Figure out how long this took
    call system_clock(end_count)
    runtime = real(end_count - start_count)/count_rate_main
    
        ! Goodbye
    write(*,*)
    write(*,'(a, f10.4, a)') " sim exited successfully. Execution time: ", runtime, " s"
    

    ! set down
    call udp_finalize()

    call exit(0)
    
end program main