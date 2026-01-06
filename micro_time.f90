module micro_time_m
    use, intrinsic :: iso_c_binding
    implicit none

    private
    public :: get_time

    interface
        subroutine QueryPerformanceCounter(lpPerformanceCount) bind(C, name="QueryPerformanceCounter")
            use iso_c_binding
            integer(c_int64_t), intent(out) :: lpPerformanceCount
        end subroutine

        subroutine QueryPerformanceFrequency(lpFrequency) bind(C, name="QueryPerformanceFrequency")
            use iso_c_binding
            integer(c_int64_t), intent(out) :: lpFrequency
        end subroutine

    end interface

contains

    function get_time() result(seconds)
        ! Returns current high-resolution time in seconds (double precision)
        real(c_double) :: seconds
        integer(c_int64_t) :: count, freq
        call QueryPerformanceCounter(count)
        call QueryPerformanceFrequency(freq)
        seconds = real(count, c_double) / real(freq, c_double)
    end function get_time

end module micro_time_m
