module connection_m
    use udp_m
    use iso_c_binding
    use jsonx_m
    use database_m, only: database, new_database, char_length, int2str
    implicit none
    
    
    type, abstract :: channel
        integer :: n
        logical :: only_send
        real, allocatable :: vals(:)
    contains
        procedure(interface_ch_init), deferred :: init
        procedure(interface_ch_send), deferred :: send
        procedure(interface_ch_recv), deferred :: recv
    end type channel
    
    abstract interface
        subroutine interface_ch_init(t, p, n)
            import :: channel, json_value
            class(channel), intent(out) :: t
            type(json_value), pointer, intent(in) :: p
            integer, intent(in), optional :: n
        end subroutine interface_ch_init
    end interface
    
    abstract interface
        subroutine interface_ch_send(t, vals)
            import :: channel
            class(channel), intent(in) :: t
            real, intent(in) :: vals(t%n)
        end subroutine interface_ch_send
    end interface
    
    abstract interface
        function   interface_ch_recv(t, x) result(y)
            import :: channel
            class(channel), intent(inout) :: t
            real, intent(in), optional :: x(:)      !! TODO: change the length maybe
            real :: y(t%n)
        end function interface_ch_recv
    end interface
    
    !!==================================================================
    
    type, extends(channel) :: udp_channel
        integer(c_int) :: sock
        character(len=:), allocatable :: ip
        integer :: port
        logical :: double_precision, blocking
    contains
        procedure :: init => udp_channel_init
        procedure :: send => udp_channel_send
        procedure :: recv => udp_channel_recv
    end type udp_channel
    
    type, extends(channel) :: file_channel
        integer :: fid
        class(database), allocatable :: db
        character(len=:), allocatable :: fn
    contains
        procedure :: init => file_channel_init
        procedure :: send => file_channel_send
        procedure :: recv => file_channel_recv
    end type file_channel
    
    !!==================================================================
    
    type :: connection
        ! real, allocatable :: vals(:)
        class(channel), allocatable :: ch
        ! type(zTimer) :: refresh
        real :: refresh_time, time
        ! logical :: only_send
    contains
        procedure :: init          => connection_init
        procedure :: send          => connection_send
        procedure :: recv          => connection_recv
        procedure :: check_refresh => connection_check_refresh
    end type connection
    
contains
    
    subroutine channel_init(t, p, n)
        class(channel), intent(out) :: t
        type(json_value), pointer, intent(in) :: p
        integer, intent(in), optional :: n
        
        character(len=:), allocatable :: char_temp
        
        if (present(n)) then
            t%n = n
        else
            call jsonx_get(p, 'number_of_values', t%n)
        end if
        
        call jsonx_get(p, 'type', char_temp)
        if     (char_temp == 'send') then
            t%only_send = .true.
        elseif (char_temp == 'receive') then
            t%only_send = .false.
            allocate(t%vals(t%n))
            t%vals = 0.
        else
            write(*,*) 'Error initializing channel '//trim(p%name)//'. Type not recognized as either send or receive! Quitting...'
            stop
        end if
        deallocate(char_temp)
        
    end subroutine channel_init
    
    !! udp channel =====================================================
    
    subroutine udp_channel_init(t, p, n)
        class(udp_channel), intent(out) :: t
        type(json_value), pointer, intent(in) :: p
        integer, intent(in), optional :: n
        
        character(len=:), allocatable :: char_temp
        
        call channel_init(t, p, n)
        
        call jsonx_get(p, 'channel_type', char_temp)
        if (char_temp /= 'udp' .and. char_temp /= 'UDP') then
            write(*,*) 'Error initializing channel '//trim(p%name)//'. Initializing as UDP channel but channel_type does not match! Quitting...'
            stop
        end if
        deallocate(char_temp)
        
        t%sock = udp_open_socket()
        call jsonx_get(p, 'port_ID', t%port)
        call jsonx_get(p, 'double_precision', t%double_precision, .false.)
        
        if (t%only_send) then
            call jsonx_get(p, 'IP_address', t%ip, '127.0.0.1')
        else
            call jsonx_get(p, 'wait_for_input', t%blocking, .false.)
            call udp_bind_socket(t%sock, t%port)
            if (.not. t%blocking) call udp_set_nonblocking(t%sock)
        end if
        
    end subroutine udp_channel_init
    
    subroutine udp_channel_send(t, vals)
        class(udp_channel), intent(in) :: t
        real, intent(in) :: vals(t%n)
        
        if (.not. t%only_send) then
            write(*,*) 'Error with UDP channel port '//int2str(t%port)//'. Attempting to send data, but channel setup as receive only! Quitting...'
            stop
        end if
        
        if (t%double_precision) then
            call udp_send_real8(t%sock, trim(t%ip), t%port, vals)
        else
            call udp_send_real4(t%sock, trim(t%ip), t%port, vals)
        end if
        
    end subroutine udp_channel_send
    
    function   udp_channel_recv(t, x) result(vals)
        class(udp_channel), intent(inout) :: t
        real, intent(in), optional :: x(:)
        
        real :: vals(t%n), temp(t%n)
        integer :: got_data
        logical :: flag
        
        if (t%only_send) then
            write(*,*) 'Error with UDP channel port '//int2str(t%port)//'. Attempting to receive data, but channel setup as send only! Quitting...'
            stop
        end if
        
        if (t%blocking) then
            if (t%double_precision) then
                call udp_recv_real8(t%sock, vals)
            else
                call udp_recv_real4(t%sock, vals)
            end if
            t%vals = vals
        else
            flag = .false.
            if (t%double_precision) then
                do
                    got_data = udp_recv_real8_nb(t%sock, temp)
                    if (got_data /= 1) then
                        exit
                    else
                        flag = .true.
                    end if
                end do
            else
                do
                    got_data = udp_recv_real4_nb(t%sock, temp)
                    if (got_data /= 1) then
                        exit
                    else
                        flag = .true.
                    end if
                end do
            end if
            if (flag) then
                vals = temp
                t%vals = vals
            else
                vals = t%vals
            end if
        end if
        
    end function udp_channel_recv
    
    !! file channel ====================================================
    
    subroutine file_channel_init(t, p, n)
        class(file_channel), intent(out) :: t
        type(json_value), pointer, intent(in) :: p
        integer, intent(in), optional :: n
        
        character(len=char_length), allocatable :: var_names(:), names(:)
        character(len=:), allocatable :: char_temp, fn, pn, input_file
        integer :: i
        
        call channel_init(t, p, n)
        
        call jsonx_get(p, 'channel_type', char_temp)
        if (char_temp /= 'file') then
            write(*,*) 'Error initializing channel '//trim(p%name)//'. Initializing as file channel but channel_type does not match! Quitting...'
            stop
        end if
        deallocate(char_temp)
        
        call jsonx_get(p, 'filename', fn)
        t%fn = fn
        call jsonx_get(p, 'pathname', pn, '')
        if (pn == '') deallocate(pn)
        
        if (t%only_send) then
            if (allocated(pn)) then
                i = len_trim(pn)
                if (pn(i:i) == '/' .or. pn(i:i) == '\') then
                    allocate(character(len=i+len_trim(fn)) :: input_file)
                    input_file = pn // fn
                else
                    allocate(character(len=i+1+len_trim(fn)) :: input_file)
                    input_file = pn // '/' // fn
                end if
            else
                allocate(character(len=len_trim(fn)) :: input_file)
                input_file = fn
            end if
            open(newunit=t%fid, file=trim(input_file), status='replace', action='write')
        else
            call jsonx_get(p, 'database_type', char_temp)
            if (allocated(pn)) then
                t%db = new_database(char_temp, t%fn, pn)
            else
                t%db = new_database(char_temp, t%fn)
            end if
            if (t%n /= t%db%n_dv) then
                write(*,*) 'Error initializing file receive channel '//trim(p%name)//'. Database output does not match expected amount of values! Quitting...'
                stop
            end if
        end if
        
    end subroutine file_channel_init
    
    subroutine file_channel_send(t, vals)
        class(file_channel), intent(in) :: t
        real, intent(in) :: vals(t%n)
        
        if (.not. t%only_send) then
            write(*,*) 'Error with file channel to file '//trim(t%fn)//'. Attempting to send data, but channel setup as receive only! Quitting...'
            stop
        end if
        
        ! call io_csvLineWrite(t%fid, vals)
        write(t%fid,'(*(E26.16E3, :","))') vals
        
    end subroutine file_channel_send
    
    function   file_channel_recv(t, x) result(vals)
        class(file_channel), intent(inout) :: t
        real, intent(in), optional :: x(:)
        real :: vals(t%n)
        
        real :: dummy
        
        if (t%only_send) then
            write(*,*) 'Error with file channel to file '//trim(t%fn)//'. Attempting to receive data, but channel setup as send only! Quitting...'
            stop
        end if
        if (.not. present(x)) then
            write(*,*) 'Error with file channel to file '//trim(t%fn)//'. Attempting to interpolate database, but independent variable values not given! Quitting...'
            stop
        end if
        if (size(x) /= t%db%n_iv) then
            write(*,*) 'Error with file channel to file '//trim(t%fn)//'. Attempting to interpolate database, but only received '//int2str(size(x))//' independent variable values but expecting '//int2str(t%db%n_iv)//'! Quitting...'
            stop
        end if
        
        vals = t%db%interpolate(x)
        t%vals = vals
        
    end function file_channel_recv
    
    !! channel =========================================================
    
    function   new_channel(p, n) result(obj)
        type(json_value), pointer, intent(in) :: p
        integer, intent(in), optional :: n
        
        class(channel), allocatable           :: obj
        type(udp_channel), allocatable        :: temp_udp
        type(file_channel), allocatable       :: temp_file
        
        character(len=:), allocatable :: ch_type
        
        call jsonx_get(p, 'channel_type', ch_type)
        
        select case (ch_type)
            case('udp', 'UDP')
                allocate(temp_udp)
                call temp_udp%init(p, n)
                obj = temp_udp
            case('file')
                allocate(temp_file)
                call temp_file%init(p, n)
                obj = temp_file
            case default
                write(*,*) 'Error initializing channel for '//trim(p%name)//'. Channel type '//trim(ch_type)//' not recognized! Quitting...'
                stop
        end select
    end function new_channel
    
    !! connection ======================================================
    
    subroutine connection_init(t, p, n)
        class(connection), intent(out) :: t
        type(json_value), pointer, intent(in) :: p
        integer, intent(in), optional :: n
        
        real :: refresh_rate
        
        t%ch = new_channel(p, n)
        
        ! t%only_send = t%ch%only_send
        
        call jsonx_get(p, 'refresh_rate', refresh_rate, 0.)
        if (refresh_rate <= 0.) then
            t%refresh_time = 1.0e-12
        else
            t%refresh_time = 1. / refresh_rate
        end if
        
        ! call t%refresh%init()
        ! t%refresh%lapTime = t%refresh%startTime-t%refresh_time    !! guarantees the first call will trigger a send/recv
        
        call cpu_time(t%time)
        t%time = t%time - t%refresh_time
        
    end subroutine connection_init
    
    subroutine connection_send(t, vals)
        class(connection), intent(inout) :: t
        real, intent(in) :: vals(t%ch%n)
        
        if (t%check_refresh()) then
            call t%ch%send(vals)
            ! t%refresh%lapTime = t%refresh%lapTime + t%refresh_time
            t%time = t%time + t%refresh_time
        end if
        
    end subroutine connection_send
    
    function   connection_recv(t, x) result(vals)
        class(connection), intent(inout) :: t
        real, intent(in), optional :: x(:)
        real :: vals(t%ch%n)
        
        if (t%check_refresh()) then
            vals = t%ch%recv(x)
            ! t%refresh%lapTime = t%refresh%lapTime + t%refresh_time
            t%time = t%time + t%refresh_time
        else
            vals = t%ch%vals
        end if
        
    end function connection_recv
    
    function   connection_check_refresh(t) result(flag)
        class(connection), intent(in) :: t
        logical :: flag
        real :: current
        ! flag = t%refresh%getLapTimeSeconds() >= t%refresh_time
        call cpu_time(current)
        flag = current - t%time >= t%refresh_time
    end function connection_check_refresh
    
end module connection_m
