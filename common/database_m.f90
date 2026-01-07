module database_m
    implicit none
    
    integer, parameter :: char_length = 200
    
    type, abstract :: database
        integer :: n_iv, n_dv, n_pts
        real, allocatable :: x(:,:), y(:,:)
        character(len=char_length), allocatable, dimension(:) :: ind_vars, dep_vars
    contains
        procedure(interface_interpolate), deferred :: interpolate
    end type
    
    abstract interface
        function   interface_interpolate(t, x) result(y)
            import :: database
            class(database), intent(in) :: t
            real, intent(in) :: x(t%n_iv)
            real :: y(t%n_dv)
        end function interface_interpolate
    end interface
    
    !!==================================================================
    
    type, extends(database) :: db_rect
        integer, allocatable :: n_pts_ind_vars(:)
        real, allocatable :: unique_ind_vars(:,:)
    contains
        procedure :: init        => db_rect_init
        procedure :: interpolate => db_rect_interpolate
    end type
    
contains
    
    function   int2str(i) result(s)
        integer, intent(in) :: i
        integer :: j
        character(:), allocatable :: s
        character(len=100) :: t
        
        write(t, '(I100)') i
        j = len(trim(adjustl(t)))
        allocate(character(j) :: s)
        s = trim(adjustl(t))
        
    end function int2str
    
    !! database helper functions =======================================
    
    subroutine append_unique_real_values(arr, new_vals)
        real, allocatable, intent(inout), dimension(:) :: arr
        real, intent(in), dimension(:) :: new_vals
        real, allocatable, dimension(:) :: temp
        integer :: i, cnt
        if (.not. allocated(arr)) then
            allocate(arr(1))
            arr(1) = new_vals(1)
        end if
        cnt = size(arr)
        do i=1, size(new_vals)
            if (.not. any(arr == new_vals(i))) then
                cnt = cnt + 1
                allocate(temp(cnt))
                temp(:cnt-1) = arr
                temp(cnt) = new_vals(i)
                deallocate(arr)
                allocate(arr(cnt))
                arr = temp
                deallocate(temp)
            end if
        end do
    end subroutine append_unique_real_values
    
    subroutine sort_array(a)
        
        real, dimension(:), intent(inout) :: a
        real :: temp
        integer :: k, i, j, m
        
        k = size(a)
        
        do m=1,k-1
            do j=k,m+1,-1
                i = j-1
                if (a(j) < a(i)) then
                    temp = a(j)
                    a(j) = a(i)
                    a(i) = temp
                end if
            end do
        end do
        
    end subroutine sort_array
    
    subroutine get_interp_indices_and_fraction(a, s, il, ih, f)
        real, intent(in), dimension(:) :: a
        real, intent(in) :: s
        integer, intent(out) :: il, ih
        real, intent(out) :: f
        
        integer :: i
        
        do i=1, size(a)-1
            if (a(i+1) < a(i)) then
                write(*,*) 'Spanwise array not sorted: ', a
                stop
            else if (a(i) <= s .and. s <= a(i+1)) then
                il = i
                ih = i+1
                f = (s - a(i)) / (a(i+1) - a(i))
                return
            end if
        end do
        write(*,*) 'Value not found in spanwise array: ', a(1) , s, a(size(a))
        stop
    end subroutine get_interp_indices_and_fraction
    
    function   decompose_j(j, Nvec) result(n)
        integer, intent(in) :: j
        integer, dimension(:), intent(in) :: Nvec
        integer :: VV
        integer, dimension(:), allocatable :: n
        
        integer :: v, denom, summ, prod, u, s, w
        
        if (j < 1 .or. j > calc_numb_coefs(Nvec)) then
            write(*,*) 'Error decomposing j '//int2str(j)//' into n for Nvec with '//int2str(calc_numb_coefs(Nvec))//' allowable coefficients! Quitting...'
            stop
        end if
        
        VV = size(Nvec)
        allocate(n(VV))
        
        n = -1
        
        do v=VV, 1, -1
            denom = 1
            do w=v+1, VV
                denom = denom * Nvec(w)
            end do
            summ = 0
            do u=v+1, VV
                prod = 1
                do s=u+1, VV
                    prod = prod * Nvec(s)
                end do
                summ = summ + n(u) * prod
            end do
            n(v) = mod((j - 1 - summ) / denom, Nvec(v))
        end do
        n = n + 1
        
    end function decompose_j
    
    function   compose_j(n, Nvec) result(j)
        integer, intent(in), dimension(:) :: n, Nvec
        integer :: j
        integer :: VV
        
        integer :: v, prod, w
        
        VV = size(Nvec)
        
        if (size(n) /= VV) then
            write(*,*) 'Error composing j for n and Nvec. Size of Nvec '//int2str(VV)//' does not equal size of n '//int2str(size(n))//'! Quitting...'
            stop
        end if
        
        j = 1
        do v=1, VV
            prod = 1
            do w=v+1, VV
                prod = prod * Nvec(w)
            end do
            j = j + (n(v)-1) * prod
        end do
    end function compose_j
    
    function   calc_numb_coefs(Nvec) result(J)
        integer, intent(in), dimension(:) :: Nvec
        integer :: J
        
        integer :: V, i
        
        V = size(Nvec)
        
        J = 1
        do i=1, V
            J = J * Nvec(i)
        end do
    end function calc_numb_coefs
    
    subroutine split_line(line, tokens, sep)
        character(len=*), intent(in) :: line
        character(len=:), allocatable, dimension(:), intent(out) :: tokens
        character(len=1), optional :: sep
        
        integer :: i, start_pos, len_line, token_count, cnt
        integer :: token_len, max_token_len
        character(len=1) :: sepp
        
        if (.not. present(sep)) then
            sepp = ','
        else
            sepp = sep
        end if
        
        len_line = len_trim(line)
        
        !! find the number of tokens and max token length
        start_pos = 1
        token_count = 0
        max_token_len = 0
        
        do i = 1, len_line + 1
            if (i > len_line .or. line(i:i) == sepp) then
                token_len = i - start_pos
                if (token_len > 0) then
                    token_count = token_count + 1
                    if (token_len > max_token_len) max_token_len = token_len
                end if
                start_pos = i + 1
            end if
        end do
        
        !! Allocate token
        allocate(character(len=max_token_len) :: tokens(token_count))
        
        !! fill the tokens
        start_pos = 1
        cnt = 0
        do i = 1, len_line + 1
            if (i > len_line .or. line(i:i) == sepp) then
                token_len = i - start_pos
                if (token_len > 0) then
                    cnt = cnt + 1
                    tokens(cnt) = line(start_pos:start_pos + token_len - 1)
                end if
                start_pos = i + 1
            end if
        end do
        
    end subroutine split_line
    
    !!==================================================================
    
    subroutine database_init(t, fn, pn)
        class(database), intent(inout) :: t
        character(len=*), intent(in) :: fn
        character(len=*), intent(in), optional :: pn
        
        character(len=:), allocatable :: input_file
        logical :: flag
        integer :: fid, ios
        character(len=char_length*100) :: line                          !! TODO: change to allocatable length?
        integer :: n_comments, n_parameters, n_rows, n_pts, headings_row
        character(len=:), allocatable, dimension(:) :: cols
        integer :: i, j, k1, k2, k3, k4, ii
        character(len=:), allocatable :: temp
        real, allocatable :: vals(:)
        
        !! and on the pathname if given
        if (present(pn)) then
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
        
        !! Check if it exists
        inquire(file=input_file, exist=flag)
        if (.not. flag) then
            write(*,*) 'Error initializing database. CSV file: '//trim(input_file)//' does not exist. Quitting...'
            stop
        end if
        
        !! open the file and determine the number of
        !!      dependent variables
        !!      independent variables
        !!      data points
        !!      length of the file / rows
        open(newunit=fid, file=input_file, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Error opening file '//trim(input_file)//'! Quitting...'
            stop
        end if
        n_rows = 0
        n_comments = 0
        n_parameters = 0
        n_pts = 0
        flag = .false.
        outer: do
            
            read(fid, '(A)', iostat=ios) line
            if (ios /= 0) exit
            
            n_rows = n_rows + 1
            
            if (line(1:1) == '#') then
                n_comments = n_comments + 1
                cycle outer
            end if
            
            if (flag) then
                n_pts = n_pts + 1
                call split_line(line, cols)
                if (size(cols) /= t%n_iv + t%n_dv) then
                    write(*,*) 'Error reading database csv '//trim(input_file)//'! Row '//&
                        int2str(n_rows)//' has a different number of columns! Quitting...'
                    stop
                end if
                deallocate(cols)
                cycle outer
            end if
            
            call split_line(line, cols)
            
            if (trim(cols(1)) == 'parameter') then
                n_parameters = n_parameters + 1
                if (trim(cols(2)) == 'independent_variables') then
                    read(cols(3), *) t%n_iv
                end if
            else
                flag = .true.           !! we've reached the column headings row
                headings_row = n_rows
                if (t%n_iv == 0) then
                    write(*,*) 'Error reading database csv '//trim(input_file)//'! Number of independent variables not given! Quitting...'
                    close(fid)
                    stop
                end if
                t%n_dv = size(cols) - t%n_iv
                allocate(t%ind_vars(t%n_iv), t%dep_vars(t%n_dv))
                do i=1, t%n_iv
                    t%ind_vars(i) = cols(i)
                end do
                do i=1, t%n_dv
                    j = i + t%n_iv
                    t%dep_vars(i) = cols(j)
                end do
            end if
            deallocate(cols)
        end do outer
        
        !! check results
        if (n_rows /= n_comments + n_parameters + 1 + n_pts) then
            write(*,*) 'Error reading database csv '//trim(input_file)//'! Number of rows, comments, parameters, and data points do not agree!'
            write(*,*) '    rows:        '//int2str(n_rows)
            write(*,*) '    comments:    '//int2str(n_comments)
            write(*,*) '    parameters:  '//int2str(n_parameters)
            write(*,*) '    data points: '//int2str(n_pts)
            close(fid)
            stop
        end if
        
        !! get data
        t%n_pts = n_pts
        allocate(t%x(n_pts, t%n_iv), t%y(n_pts, t%n_dv), vals(t%n_iv + t%n_dv))
        rewind(fid)
        do i=1, headings_row
            read(fid, *)
        end do
        j = 0
        do i=1, n_rows - headings_row
            read(fid, '(A)', iostat=ios) line
            if (ios /= 0) then
                write(*,*) 'Error reading csv file: '//trim(input_file) //'! Quitting...'
                stop
            end if
            ! call io_split_line(line, cols)
            ! do j=1, t%n_iv
                ! read(cols(j), *) t%x(i,j)
            ! end do
            ! do j=1, t%n_dv
                ! read(cols(t%n_iv + j), *) t%y(i,j)
            ! end do
            ! deallocate(cols)
            if (line(1:1) == '#') cycle
            j = j + 1
            read(line, *) vals
            t%x(j,:) = vals(:t%n_iv)
            t%y(j,:) = vals(t%n_iv+1:)
        end do
        close(fid)
    end subroutine database_init
    
    !!==================================================================
    
    subroutine db_rect_init(t, fn, pn)  !! potential edit: add verbosity?
        class(db_rect), intent(out) :: t
        character(len=*), intent(in) :: fn
        character(len=*), intent(in), optional :: pn
        
        real, dimension(:,:), allocatable :: x, y
        ! character(len=char_length), dimension(:), intent(in) :: names_x, names_y
        
        integer :: i, j, k, l, ii
        real, allocatable, dimension(:) :: unique_temp
        integer, allocatable, dimension(:) :: n
        
        call database_init(t, fn, pn)
        
        k = t%n_pts
        allocate(x(k,t%n_iv), y(k,t%n_dv))
        x = t%x
        y = t%y
        t%x = 0.
        t%y = 0.
        
        !!==============================================================
        !! determine the unique ind var values
        !!==============================================================
        if (.true.) then
        allocate(t%n_pts_ind_vars(t%n_iv), n(t%n_iv))
        !! loop thru the independent variables
        do i=1, t%n_iv
            !! append the value onto the unique temp array
            call append_unique_real_values(unique_temp, x(:,i))
            !! find the number of unique values
            t%n_pts_ind_vars(i) = size(unique_temp)
            !! reset the temp array
            deallocate(unique_temp)
        end do
        !! allocate the unique array
        allocate(t%unique_ind_vars(maxval(t%n_pts_ind_vars), t%n_iv))
        t%unique_ind_vars = -1.0e9
        !! loop thru the independent variables
        do i=1, t%n_iv
            !! append the values onto the unique temp array
            call append_unique_real_values(unique_temp, x(:,i))
            !! sort the temp array
            call sort_array(unique_temp)
            !! store the values
            t%unique_ind_vars(:t%n_pts_ind_vars(i),i) = unique_temp
            !! reset the temp array
            deallocate(unique_temp)
        end do
        end if
        
        !!==============================================================
        !! set the data arrays sorted by j indexing
        !!==============================================================
        if (.true.) then
        !! loop thru j indices
        do j=1, calc_numb_coefs(t%n_pts_ind_vars)
            !! decompose j
            n = decompose_j(j, t%n_pts_ind_vars)
            !! loop thru points
            ii = 0
            outer: do i=1, k
                !! check if the point matches each ind var value
                !! loop thru ind var values
                do l=1, t%n_iv
                    if (x(i,l) /= t%unique_ind_vars(n(l), l)) cycle outer
                end do
                ii = i
                exit outer
            end do outer
            !! store values
            if (ii /= 0) then
                t%x(j,:) = x(ii,:)
                t%y(j,:) = y(ii,:)
            !! failure
            else
                write(*,*) 'Failed to find data in a rectilinear grid format! Quitting...'
                write(*,*)
                write(*,*) 'unique variable values'
                do i=1, t%n_iv
                    write(*,*) trim(t%ind_vars(i))
                    write(*,*) t%unique_ind_vars(:t%n_pts_ind_vars(i),i)
                    write(*,*)
                end do
                write(*,*) 'Values not found:'
                do l=1, t%n_iv
                    write(*,*) t%unique_ind_vars(n(l), l)
                end do
                stop
            end if
        end do
        end if
        
    end subroutine db_rect_init
    
    function   db_rect_interpolate(t, x) result(y)
        class(db_rect), intent(in) :: t
        real, intent(in) :: x(t%n_iv)
        real :: y(t%n_dv)
        
        integer, dimension(t%n_iv) :: il, ih
        real, dimension(t%n_iv) :: f
        integer :: i, l, j, jh, V, l_new, jl
        real, allocatable, dimension(:,:) :: d, d_new
        integer, allocatable, dimension(:) :: Nvec, nh, n, nl
        
        !!==============================================================
        !! down select the data
        !!==============================================================
        !! this function needs to be able to handle an arbitrary number
        !! of independent variables. However, fortran cannot dynamically
        !! change the number of dimensions on an array. This complicates
        !! the process somewhat. To overcome this, the would-be multi
        !! dimensional arrays are flattened using j-indexing controlled by
        !! the decompose_j and compose_j functions, along with the func
        !! calc_numb_coefs.
        
        !! perform 1D interpolations
        do i=1, t%n_iv
            call get_interp_indices_and_fraction(t%unique_ind_vars(:t%n_pts_ind_vars(i), i), x(i), il(i), ih(i), f(i))
        end do
        
        !! get the bounding rectilinear box
        l = 2 ** t%n_iv
        allocate(d(l,t%n_dv), Nvec(t%n_iv), nh(t%n_iv), n(t%n_iv))
        Nvec = 2
        do j=1, calc_numb_coefs(Nvec)
            !! get local n vec
            n = decompose_j(j, Nvec)
            !! convert to global values
            nh = 0
            do i=1, t%n_iv
                if (n(i) == 1) then
                    nh(i) = il(i)
                else
                    nh(i) = ih(i)
                end if
            end do
            jh = compose_j(nh, t%n_pts_ind_vars)
            !! store results
            d(j,:) = t%y(jh,:)
        end do
        deallocate(nh, n)
        
        !!==============================================================
        !! down select the data
        !!==============================================================
        !! d has 2^V points
        !! interpolate on the last dimension to down select the data
        !! the number of points will half, i.e. 2^V -> 2^(V-1)
        !! repeat this until only one point remains
        
        !! loop thru the number of dimensions - 1
        do V=t%n_iv-1, 1, -1
            !! determine the size of the new down select array and allocate
            l_new = 2**V
            allocate(d_new(l_new, t%n_dv), nh(V+1), nl(V+1), n(V))
            !! loop thru the j indices of the new array
            do j=1, l_new
                !! decompose j to get the n values
                n = decompose_j(j, Nvec(:V))
                !! create nn values for hi and lo points
                nh(:V) = n
                nl(:V) = n
                nh(V+1) = 2
                nl(V+1) = 1
                !! convert these to j hi and lo points
                jh = compose_j(nh, Nvec)
                jl = compose_j(nl, Nvec)
                !! set the new d value
                d_new(j,:) = (d(jh,:) - d(jl,:)) * f(V+1) + d(jl,:)
                !! reset arrays as needed
                
            end do
            !! reset arrays as needed
            deallocate(Nvec, nh, nl, d, n)
            allocate(Nvec(V), d(l_new, t%n_dv))
            Nvec = 2
            d = d_new
            deallocate(d_new)
        end do
        
        !! perform last interpolation on first dimension
        y = (d(2,:) - d(1,:)) * f(1) + d(1,:)
        deallocate(Nvec, d)
        
    end function db_rect_interpolate
    
    !!==================================================================
    
    function   new_database(db_type, fn, pn) result(obj)
        character(len=*), intent(in) :: db_type, fn
        character(len=*), intent(in), optional :: pn
        class(database), allocatable           :: obj
        
        type(db_rect), allocatable             :: temp_rect
        ! type(db_recu), allocatable             :: temp_recu
        ! type(db_scat), allocatable             :: temp_scat
        
        select case (db_type)
            case('rectilinear', 'Rectilinear')
                allocate(temp_rect)
                call temp_rect%init(fn, pn)
                obj = temp_rect
            case('recursive', 'Recursive')
                ! allocate(temp_recu)
                ! call temp_recu%init(fn, pn)
                ! obj = temp_recu
                write(*,*) 'Error initializing database from file '//trim(fn)//'. Recursive database type currently not available! Quitting...'
                stop
            case('scatter', 'Scatter')
                ! allocate(temp_scat)
                ! call temp_scat%init(fn, pn)
                ! obj = temp_scat
                write(*,*) 'Error initializing database from file '//trim(fn)//'. Scatter database type currently not available! Quitting...'
                stop
            case default
                write(*,*) 'Error initializing database from file '//trim(fn)//'. Database type '//trim(db_type)//' not recognized!'
                write(*,*) 'Must be either rectilinear, recursive, or scatter (currently only rectilinear is available). Quitting...'
                stop
        end select
    end function new_database
    
end module database_m
