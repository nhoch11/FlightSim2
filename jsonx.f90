! Extends the original JSON module for usability
module jsonx_m
    use json_m
    implicit none
    
    ! Everything is public. Not great coding practice, but allows module re-export downstream to minimize "use" statements
    ! private
    ! public :: jsonx_get, jsonx_load
    
    logical :: json_found
    
    interface jsonx_get
        module procedure :: jsonx_value_get_real, jsonx_file_get_real
        module procedure :: jsonx_value_get_integer, jsonx_file_get_integer
        module procedure :: jsonx_value_get_string, jsonx_file_get_string
        module procedure :: jsonx_value_get_logical, jsonx_file_get_logical
        module procedure :: jsonx_value_get_real_array
        module procedure :: jsonx_value_get_char_array
        module procedure :: jsonx_value_get_real_spanwise_array
        module procedure :: jsonx_value_get_char_spanwise_array
        module procedure :: jsonx_value_get_logical_array
        module procedure :: jsonx_value_get_json
    end interface jsonx_get
    
contains

    
    subroutine jsonx_load(fn, json)
        character(len=*), intent(in) :: fn
        type(json_value), pointer, intent(out) :: json
        
        call json_parse(fn, json)
        call json_check()
    end subroutine jsonx_load
    
    subroutine jsonx_value_get_json(json_in, name, value)
        implicit none
        type(json_value), pointer, intent(in) :: json_in
        character(len=*), intent(in) :: name
        type(json_value), pointer, intent(out) :: value
        
        type(json_value), pointer :: json
        
        call json_check_4_new_file(json_in, json)
        
        call json_value_get(json, name, value)
        
        if (json_failed()) then
            write(*,*) 'Error: Unable to read required value: '//trim(name)
            stop
        end if
    end subroutine jsonx_value_get_json
    
    subroutine jsonx_value_get_real(json_in, name, value, default_value)
        implicit none
        type(json_value),intent(in),pointer :: json_in
        character(len=*), intent(in) :: name
        real, intent(out) :: value
        real, intent(in), optional :: default_value
        type(json_value), pointer :: json
        
        call json_check_4_new_file(json_in, json)
        
        ! integer :: i
        ! type(json_value), pointer :: p
        
        call json_get(json, name, value, json_found)
        if(json_failed() .or. (.not. json_found)) then
            ! if (index(name, '[') > 0) then
                ! call json_clear_exceptions()
                ! do i=1, json_value_count(json)
                    ! call json_value_get(json, i, p)
                    ! if (trim(p%name) == trim(name)) then
                        ! call json_get(p, value=value, found=json_found)
                        ! if (.not. json_failed() .and. json_found) return
                    ! end if
                ! end do
            ! end if
            
            if (present(default_value)) then
                value = default_value
                call json_clear_exceptions()
            else
                write(*,*) 'Error: Unable to read required value: ',name
                STOP
            end if
        end if
    end subroutine jsonx_value_get_real
    
    
    subroutine jsonx_value_get_integer(json_in, name, value, default_value)
        implicit none
        type(json_value),intent(in),pointer :: json_in
        character(len=*), intent(in) :: name
        integer, intent(out) :: value
        integer, intent(in), optional :: default_value
        type(json_value), pointer :: json
        
        call json_check_4_new_file(json_in, json)
        
        ! integer :: i
        ! type(json_value), pointer :: p
        
        call json_get(json, name, value, json_found)
        if((.not.json_found) .or. json_failed()) then
            ! if (index(name, '[') > 0) then
                ! call json_clear_exceptions()
                ! do i=1, json_value_count(json)
                    ! call json_value_get(json, i, p)
                    ! if (trim(p%name) == trim(name)) then
                        ! call json_get(p, value=value, found=json_found)
                        ! if (.not. json_failed() .and. json_found) return
                    ! end if
                ! end do
            ! end if
            
            if (present(default_value)) then
                value = default_value
                call json_clear_exceptions()
            else
                write(*,*) 'Error: Unable to read required value: ', name
                STOP
            end if
        end if
        
    end subroutine jsonx_value_get_integer
    
    
    subroutine jsonx_value_get_string(json_in, name, value, default_value)
        implicit none
        type(json_value), intent(in), pointer :: json_in
        character(len=*), intent(in) :: name
        character(:), allocatable, intent(out) :: value
        character(len=*), intent(in), optional :: default_value
        type(json_value), pointer :: json
        
        call json_check_4_new_file(json_in, json)
        
        ! integer :: i
        ! type(json_value), pointer :: p
        
        call json_get(json, name, value, json_found)
        if((.not.json_found) .or. json_failed()) then
            ! if (index(name, '[') > 0) then
                ! call json_clear_exceptions()
                ! do i=1, json_value_count(json)
                    ! call json_value_get(json, i, p)
                    ! if (trim(p%name) == trim(name)) then
                        ! call json_get(p, value=value, found=json_found)
                        ! if (.not. json_failed() .and. json_found) return
                    ! end if
                ! end do
            ! end if
            
            if (present(default_value)) then
                value = default_value
                call json_clear_exceptions()
            else
                write(*,*) 'Error: Unable to read required value: ', name
                STOP
            end if
        end if
        
    end subroutine jsonx_value_get_string
    
    
    subroutine jsonx_value_get_logical(json_in, name, value, default_value)
        implicit none
        type(json_value), intent(in), pointer :: json_in
        character(len=*), intent(in) :: name
        logical, intent(out) :: value
        logical, intent(in), optional :: default_value
        type(json_value), pointer :: json
        
        call json_check_4_new_file(json_in, json)
        
        ! integer :: i
        ! type(json_value), pointer :: p
        
        call json_get(json, name, value, json_found)
        if((.not.json_found) .or. json_failed()) then
            ! if (index(name, '[') > 0) then
                ! call json_clear_exceptions()
                ! do i=1, json_value_count(json)
                    ! call json_value_get(json, i, p)
                    ! if (trim(p%name) == trim(name)) then
                        ! call json_get(p, value=value, found=json_found)
                        ! if (.not. json_failed() .and. json_found) return
                    ! end if
                ! end do
            ! end if
            
            if (present(default_value)) then
                value = default_value
                call json_clear_exceptions()
            else
                write(*,*) 'Error: Unable to read required value: ', name
                STOP
            end if
        end if
        
    end subroutine jsonx_value_get_logical
    
    
    subroutine jsonx_file_get_real(json, name, value, default_value)
        implicit none
        type(json_file) :: json
        character(len=*), intent(in) :: name
        real, intent(out) :: value
        real, intent(in), optional :: default_value
        
        call json%get(name, value)
        if(json_failed()) then
            if (present(default_value)) then
                value = default_value
                call json_clear_exceptions()
            else
                write(*,*) 'Error: Unable to read required value: ', trim(name)
                STOP
            end if
        end if
        
    end subroutine jsonx_file_get_real
    
    
    subroutine jsonx_file_get_integer(json, name, value, default_value)
        implicit none
        type(json_file) :: json
        character(len=*), intent(in) :: name
        integer, intent(out) :: value
        integer, intent(in), optional :: default_value
        
        call json%get(name, value)
        if(json_failed()) then
            if (present(default_value)) then
                value = default_value
                call json_clear_exceptions()
            else
                write(*,*) 'Error: Unable to read required value: ',name
                STOP
            end if
        end if
    end subroutine jsonx_file_get_integer
    
    
    subroutine jsonx_file_get_string(json, name, value)
        implicit none
        type(json_file) :: json
        character(len=*), intent(in) :: name
        character(:), allocatable, intent(out) :: value
        
        call json%get(name, value)
        if(json_failed()) then
            write(*,*) 'Error: Unable to read required value: ',name
            STOP
        end if
        
        value = trim(value)
    end subroutine jsonx_file_get_string
    
    
    subroutine jsonx_file_get_logical(json, name, value)
        implicit none
        type(json_file) :: json
        character(len=*), intent(in) :: name
        logical, intent(out) :: value
        
        call json%get(name, value)
        if(json_failed()) then
            write(*,*) 'Error: Unable to read required value: ',name
            STOP
        end if
        
    end subroutine jsonx_file_get_logical
    
    
    subroutine json_check()
        if(json_failed()) then
            call print_json_error_message()
            STOP
        end if
    end subroutine json_check
        
    subroutine print_json_error_message()
        implicit none
        character(len=:),allocatable :: error_msg
        logical :: status_ok
        
        !get error message:
        call json_check_for_errors(status_ok, error_msg)
        
        !print it if there is one:
        if (.not. status_ok) then
            write(*,'(A)') error_msg
            deallocate(error_msg)
            call json_clear_exceptions()
        end if
        
    end subroutine print_json_error_message
    
    subroutine init_spanArray_constant(a, c)
        
        real, allocatable, dimension(:,:), intent(out) :: a
        real, intent(in) :: c
        
        allocate(a(2,2))
        a(:,1) = [0.0, 1.0]
        a(:,2) = c
        
    end subroutine init_spanArray_constant
    
    subroutine jsonx_value_get_real_array(json_in, name, value, default_value, l)
        implicit none
        type(json_value), intent(in), pointer :: json_in
        character(len=*), intent(in) :: name
        real, intent(out), allocatable, dimension(:) :: value
        real, intent(in), optional :: default_value
        integer, intent(in), optional :: l
        type(json_value), pointer :: json
        
        call json_check_4_new_file(json_in, json)
        
        ! integer :: i
        ! type(json_value), pointer :: p
        
        call json_get(json, name, value, json_found)
        if (json_failed() .or. (.not. json_found)) then
            if (present(default_value)) then
                if (present(l)) then
                    allocate(value(l))
                else
                    allocate(value(2))
                end if
                value = default_value
                call json_clear_exceptions()
            else
                write(*,*) 'Error: Unable to read required value: ', name
                stop
            end if
        end if
    end subroutine jsonx_value_get_real_array
    
    subroutine jsonx_value_get_real_spanwise_array(json_in, name, value, default_value)
        implicit none
        type(json_value), intent(in), pointer :: json_in
        type(json_value), pointer :: temp
        character(len=*), intent(in) :: name
        real, intent(out), allocatable, dimension(:,:) :: value
        real, intent(in), optional :: default_value
        real, allocatable, dimension(:) :: span, v
        real :: singleValue
        integer :: l, i
        !~ logical :: found
        type(json_value), pointer :: json
        
        call json_check_4_new_file(json_in, json)
        
        call json_get(json, name, singleValue, json_found)
        if(json_failed() .or. (.not. json_found)) then
            call json_clear_exceptions()
            
            call json_get(json, name, temp, json_found)
            if(json_failed() .or. (.not. json_found)) then
                if (present(default_value)) then
                    !~ allocate(value(2,2))
                    !~ value(1,1) = 0.0
                    !~ value(2,1) = 1.0
                    !~ value(1,2) = default_value
                    !~ value(2,2) = default_value
                    call init_spanArray_constant(value, default_value)
                    call json_clear_exceptions()
                    return
                else
                    write(*,*) 'Error: Unable to read required value: ',name
                    STOP
                end if
            end if
            
            call json_get(temp, 'value', v, json_found)
            if (json_failed() .or. (.not. json_found)) then
                if (present(default_value)) then
                    !~ allocate(value(2,2))
                    !~ value(1,1) = 0.0
                    !~ value(2,1) = 1.0
                    !~ value(1,2) = default_value
                    !~ value(2,2) = default_value
                    call init_spanArray_constant(value, default_value)
                    call json_clear_exceptions()
                    return
                else
                    write(*,*) 'Error: Unable to read required value: ', name
                    stop
                end if
            end if
            
            call json_get(temp, 'span', span, json_found)
            if (json_failed() .or. (.not. json_found)) then
                if (present(default_value)) then
                    !~ allocate(value(2,2))
                    !~ value(1,1) = 0.0
                    !~ value(2,1) = 1.0
                    !~ value(1,2) = default_value
                    !~ value(2,2) = default_value
                    call init_spanArray_constant(value, default_value)
                    call json_clear_exceptions()
                    return
                else
                    write(*,*) 'Error: Unable to read required value: ', name
                    stop
                end if
            end if
            
            l = size(v)
            if (l == size(span)) then
                allocate(value(l, 2))
                do i=1,l
                    value(i,1) = span(i)
                    value(i,2) = v(i)
                end do
            else
                write(*,*) 'Error: Number of spanwise points and value points are not equivalent for value: ', name
                stop
            end if
        else
            !~ allocate(value(2,2))
            !~ value(1,1) = 0.0
            !~ value(2,1) = 1.0
            !~ value(1,2) = singleValue
            !~ value(2,2) = singleValue
            call init_spanArray_constant(value, singleValue)
        end if
        
    end subroutine jsonx_value_get_real_spanwise_array
    
    subroutine jsonx_value_get_char_array(json_in, name, value, default_value)
        implicit none
        type(json_value), intent(in), pointer :: json_in
        character(len=*), intent(in) :: name
        character(len=*), intent(out), allocatable, dimension(:) :: value
        character(len=*), intent(in), optional :: default_value
        
        character(len=:), allocatable :: singleValue
        type(json_value), pointer :: json
        
        call json_check_4_new_file(json_in, json)
        
        call json_get(json, name, singleValue, json_found)
        if (json_failed() .or. (.not. json_found)) then
            call json_clear_exceptions()
            call json_get(json, name, value, json_found)
            if (json_failed() .or. (.not. json_found)) then
                if (present(default_value)) then
                    allocate(value(1))
                    value(1) = default_value
                    call json_clear_exceptions()
                else
                    write(*,*) 'Error: Unable to read required value: ', name
                    stop
                end if
            end if
        else
            allocate(value(1))
            value(1) = trim(singleValue)
            call json_clear_exceptions()
        end if
    end subroutine jsonx_value_get_char_array
    
    subroutine jsonx_value_get_char_spanwise_array(json_in, name, span, value, default_value)
        implicit none
        type(json_value), intent(in), pointer :: json_in
        type(json_value), pointer :: temp
        character(len=*), intent(in) :: name
        real, intent(out), allocatable, dimension(:) :: span
        character(len=*), intent(out), allocatable, dimension(:) :: value
        character(len=*), intent(in), optional :: default_value
        !~ character(len=*), allocatable, dimension(:) :: v
        real, allocatable, dimension(:) :: s
        character(len=:), allocatable :: singleValue
        integer :: l, i
        !~ logical :: found
        type(json_value), pointer :: json
        
        call json_check_4_new_file(json_in, json)
        
        call json_get(json, name, singleValue, json_found)
        if(json_failed() .or. (.not. json_found)) then
            call json_clear_exceptions()
            
            call json_get(json, name, temp, json_found)
            if(json_failed() .or. (.not. json_found)) then
                if (present(default_value)) then
                    allocate(value(2), span(2))
                    span(1) = 0.0
                    span(2) = 1.0
                    value(1) = default_value
                    value(2) = default_value
                    call json_clear_exceptions()
                    return
                else
                    write(*,*) 'Error: Unable to read required value: ',name
                    STOP
                end if
            end if
            
            call json_get(temp, 'value', value, json_found)
            if (json_failed() .or. (.not. json_found)) then
                if (present(default_value)) then
                    allocate(value(2), span(2))
                    span(1) = 0.0
                    span(2) = 1.0
                    value(1) = default_value
                    value(2) = default_value
                    call json_clear_exceptions()
                    return
                else
                    write(*,*) 'Error: Unable to read required value: ', name
                    stop
                end if
            end if
            
            call json_get(temp, 'span', s, json_found)
            if (json_failed() .or. (.not. json_found)) then
                if (present(default_value)) then
                    allocate(value(2), span(2))
                    span(1) = 0.0
                    span(2) = 1.0
                    value(1) = default_value
                    value(2) = default_value
                    call json_clear_exceptions()
                    return
                else
                    write(*,*) 'Error: Unable to read required value: ', name
                    stop
                end if
            end if
            
            l = size(value)
            if (l == size(s)) then
                allocate(span(l))
                do i=1,l
                    span(i)  = s(i)
                end do
            else
                write(*,*) 'Error: Number of spanwise points and value points are not equivalent for value: ', name
                stop
            end if
        else
            allocate(value(2), span(2))
            span(1) = 0.0
            span(2) = 1.0
            value(1) = trim(singleValue)
            value(2) = trim(singleValue)
            call json_clear_exceptions()
            return
        end if
        
    end subroutine jsonx_value_get_char_spanwise_array
    
    subroutine jsonx_value_get_logical_array(json_in, name, value, default_value, l)
        implicit none
        type(json_value), intent(in), pointer :: json_in
        character(len=*), intent(in) :: name
        logical, intent(out), allocatable, dimension(:) :: value
        logical, intent(in), optional :: default_value
        integer, intent(in), optional :: l
        type(json_value), pointer :: json
        
        call json_check_4_new_file(json_in, json)
        
        call json_get(json, name, value, json_found)
        if (json_failed() .or. (.not. json_found)) then
            call json_clear_exceptions()
            if (present(default_value)) then
                if (present(l)) then
                    allocate(value(l))
                else
                    allocate(value(2))
                end if
                value = default_value
            else
                write(*,*) 'Error: Unable to read required value: ', name
                stop
            end if
        end if
    end subroutine jsonx_value_get_logical_array
    
    subroutine json_check_4_new_file(s_in, s_out)
        type(json_value), pointer, intent(in) :: s_in
        type(json_value), pointer, intent(out) :: s_out
        
        ! type(json_value), pointer :: temp
        character(len=:), allocatable :: fn
        logical :: exists
        
        call json_get(s_in, "filepath", fn, json_found)
        
        if (json_failed() .or. (.not. json_found)) then
            call json_clear_exceptions()
            s_out => s_in
        else
            !! Check if file exists
            inquire(file=fn, exist=exists)
            if (.not. exists) then
                write(*,*) "!!! The file ", trim(fn), " does not exist. Quitting..."
                stop
            end if
            !! Load settings from input file
            call json_parse(fn, s_out)
            call json_check()
        end if
    end subroutine json_check_4_new_file
    
end module jsonx_m
