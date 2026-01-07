program main
    
    use sim_m
    use udp_m, only: udp_initialize, udp_finalize
    ! use f16_566_m

    implicit none
    ! character(100) :: filename
    integer :: start_count, end_count, i_unit, ios
    real :: count_rate_main, runtime
    
    ! set up
    call udp_initialize()

    ! Start timer
    call system_clock(start_count, count_rate)

    ! call init("7.2.2.B.json")
    call init("test_connection.json")
    call run()
    

    ! Figure out how long this took
    call system_clock(end_count)
    runtime = real(end_count - start_count)/count_rate_main
    
        ! Goodbye
    write(*,*)
    write(*,'(a, f10.4, a)') " f16 exited successfully. Execution time: ", runtime, " s"
    

    ! set down
    call udp_finalize()

    call exit(0)
    
end program main