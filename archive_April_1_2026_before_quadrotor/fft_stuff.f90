subroutine psd(x, dt, psd_norm, filename)
    implicit none
    real, intent(in) :: x(:), dt
    real, intent(out), optional :: psd_norm(:,:)
    character(len=*), intent(in), optional :: filename
    integer :: n, n_pad, k, i, iunit        ! <-- added n_pad
    real :: fs, f, mean, stdev
    real, allocatable :: Pxx(:), loc_psd_norm(:,:)
    complex, allocatable :: cx(:)           ! <-- replaced eye/temp with cx

    n = size(x)
    if (mod(n,2) /= 0) error stop "periodogram: N must be even for N/2 bin."

    ! <-- find next power of 2 and allocate padded complex array
    n_pad = 1
    do while (n_pad < n)
        n_pad = n_pad * 2
    end do

    allocate(cx(n_pad))
    allocate(Pxx(n/2+1))
    allocate(loc_psd_norm(n/2+1,2))

    ! compute mean and stdev (unchanged)
    mean = 0.0
    stdev = 0.0
    do i=1,n
        mean = mean + x(i)/n
    end do
    do i=1,n
        stdev = stdev + (x(i)-mean)**2
    end do
    stdev = sqrt(stdev/(n-1))

    if(present(filename)) then
        open(newunit=iunit, file=filename, status="replace", action="write")
        write(iunit,'(A)') 'fs[Hz],PSD[ft^2/Hz],Normalized_PSD,Mean,STDEV'
        write(*,*) trim(filename),'   Mean = ',mean,'   Std. Dev = ',stdev
    end if

    ! <-- build zero-padded complex input and run FFT
    do i = 1, n
        cx(i) = cmplx(x(i) - mean, 0.0)
    end do
    cx(n+1:n_pad) = cmplx(0.0, 0.0)
    call fft(cx)

    ! <-- replaced DFT k-loop body with FFT output lookup
    Pxx(:) = 0.0
    fs = 1.0/dt
    do k=0,n/2
        f = real(k)*fs/real(n)
        Pxx(k+1) = 1/fs/real(n)*abs(cx(k+1))**2
        if (k /= 0 .and. k /= n/2) Pxx(k+1) = 2.0*Pxx(k+1)
        loc_psd_norm(k+1,1) = f
        loc_psd_norm(k+1,2) = Pxx(k+1)/stdev**2
        if(present(filename)) then
            if(k==0) then
                write(iunit,'(ES20.12,",",ES20.12,",",ES20.12,",",ES20.12,",",ES20.12)') f, Pxx(k+1), loc_psd_norm(k+1,2),mean,stdev
            else
                write(iunit,'(ES20.12,",",ES20.12,",",ES20.12)') f, Pxx(k+1), loc_psd_norm(k+1,2)
            end if
        end if
    end do

    if(present(filename)) close(iunit)
    if(present(psd_norm)) psd_norm = loc_psd_norm
end subroutine psd

function next_power_of_2(n) result(m)
        implicit none
        integer, intent(in) :: n
        integer :: m
        m = 1
        do while (m < n)
            m = m * 2
        end do
    end function next_power_of_2

    subroutine fft(x)
    complex, intent(inout) :: x(:)
    integer :: n, m, j, k, half, i
    complex :: t, w, wn

    n = size(x)

    ! Bit-reversal permutation
    j = 1
    do i = 1, n-1
        if (i < j) then
            t = x(j); x(j) = x(i); x(i) = t
        end if
        m = n/2
        do while (m >= 1 .and. j > m)
            j = j - m; m = m/2
        end do
        j = j + m
    end do

    ! Cooley-Tukey butterfly stages
    half = 1
    do while (half < n)
        wn = exp(cmplx(0.0, -2.0*pi/real(2*half)))
        w  = cmplx(1.0, 0.0)
        do k = 1, half
            do j = k, n, 2*half
                t          = w * x(j + half)
                x(j + half) = x(j) - t
                x(j)       = x(j) + t
            end do
            w = w * wn
        end do
        half = half * 2
    end do
end subroutine fft