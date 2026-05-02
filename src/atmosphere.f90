module atmosphere_m
    use hoch_m
    use jsonx_m
    implicit none

    type atmosphere_type
        real, allocatable :: wind(:)
        character(len=:), allocatable :: turb_model, turb_intensity
        real :: wingspan, hstab_dist, vstab_dist
        logical :: turb_repeatable, turb_on, turb_model_allowed
        integer :: hist_len

        real :: light_height_ag(3) = (/2000., 8000., 17000./)
        real :: light_sigma(3) = (/5., 5., 3./)
        real :: moderate_height_ag(3) = (/2000., 11000., 45000./)
        real :: moderate_sigma(3) = (/10., 10., 3./)
        real :: severe_height_ag(4) = (/2000., 4000., 20000., 80000./)
        real :: severe_sigma(4) = (/15., 21., 21., 3./)

        real, allocatable :: turb_height_ag(:), turb_sigma(:)
        real :: prev_turb(6), prev_xyz(3), prev_f, prev_g
        real :: Lu, Lv, Lw, Lb, xff
        real, allocatable :: w_tilde_hist(:), v_tilde_hist(:), xff_hist(:)
    end type atmosphere_type

contains

    subroutine atmosphere_init(this, j_atmosphere)
        implicit none
        type(atmosphere_type) :: this
        type(json_value), pointer :: j_atmosphere, j_turb, j_sample
        logical :: found
        integer :: i, n, seed
        integer, allocatable :: seed_array(:)

        write(*,*) "Initializing Atmospheric Model..."
        call jsonx_get(j_atmosphere, "constant_wind[ft/s]", this%wind, 0.0, 3)
        write(*,*) " Constant Wind [ft/s] = ",this%wind(:)

        call json_get(j_atmosphere, "turbulence", j_turb, found)
        if (found) then
            call jsonx_get(j_turb, "model", this%turb_model, 'none')
            if (this%turb_model .ne. 'none') then
                call jsonx_get(j_turb, "wingspan[ft]", this%wingspan)
                call jsonx_get(j_turb, "hstab_distance[ft]", this%hstab_dist)
                call jsonx_get(j_turb, "vstab_distance[ft]", this%vstab_dist)
                call jsonx_get(j_turb, "intensity", this%turb_intensity)
                call jsonx_get(j_turb, "history_array_length", this%hist_len)
                call jsonx_get(j_turb, "repeatable", this%turb_repeatable)
                allocate(this%w_tilde_hist(this%hist_len), source=0.0)
                allocate(this%v_tilde_hist(this%hist_len), source=0.0)
                allocate(this%xff_hist(this%hist_len), source=0.0)
                
                write(*,*) "    Turbulence Intesity = ", this%turb_intensity
                write(*,*) "       Turbulence Model = ", this%turb_model
                
                ! set up random number generator
                if (this%turb_repeatable) then
                    seed = 12345
                else
                    call system_clock(count=seed)
                end if 
                call random_seed(size=n)
                allocate(seed_array(n))
                seed_array = seed + 7 + [(i-1, i=1,n)]
                call random_seed(put=seed_array)
                deallocate(seed_array)
                
                select case(trim(this%turb_intensity))
                case("light")
                    allocate(this%turb_height_ag(3))
                    allocate(this%turb_sigma(3))
                    this%turb_height_ag = this%light_height_ag
                    this%turb_sigma = this%light_sigma
                case("moderate")
                    allocate(this%turb_height_ag(3))
                    allocate(this%turb_sigma(3))
                    this%turb_height_ag = this%moderate_height_ag
                    this%turb_sigma = this%moderate_sigma
                case("severe")
                    allocate(this%turb_height_ag(4))
                    allocate(this%turb_sigma(4))
                    this%turb_height_ag = this%severe_height_ag
                    this%turb_sigma = this%severe_sigma
                end select
                    
                select case(trim(this%turb_model))
                case("dryden_beal")
                    this%Lu = 1750.0
                    this%Lv = 875.0
                    this%Lw = 875.0
                    this%Lb = 4.0*this%wingspan/PI
                    this%turb_model_allowed = .true.
                    this%turb_on = .true.
                case default
                    write(*,*) ' Error in input turbulence model. Turning turbulence off'
                    this%turb_model_allowed = .false.
                    this%turb_on = .false.
                end select

                ! initialize values
                this%prev_turb(:) = 0.0
                this%prev_xyz(:) = 0.0
                this%prev_f = 0.0
                this%prev_g = 0.0
                this%xff = 0.0
                this%xff_hist = 0.0

                call json_get(j_turb, "sample", j_sample, found)
                if (found) call turbulence_sample(this, j_sample)
            else
                ! turb model was specified as 'none'
                this%turb_model_allowed = .false.
                this%turb_on = .false.
            end if
        end if
    end subroutine atmosphere_init


    function atmosphere_get_turbulence(this, states) result(ans)
        implicit none
        type(atmosphere_type) :: this
        real :: states(21)
        real :: ans(6)
        real :: dx, sigma

        if (this%turb_on) then
            dx = sqrt((states(7) - this%prev_xyz(1))**2 + (states(8) - this%prev_xyz(2))**2 + (states(9) - this%prev_xyz(3))**2)
            if (dx > TOLERANCE) then
                sigma = interpolate_1D(this%turb_height_ag, this%turb_sigma, -states(9))
            
                ans(:) = get_turbulence(this, dx, sigma, sigma, sigma)
                this%prev_xyz(:) = states(7:9)
            else
                write(*,*) "  WARNING: dx in atmosphere_get_turbulence is near zero, using previous turb value"
                ans(:) = this%prev_turb(:)
            end if
        else
            ans(:) = 0.0
        end if 
        
    end function atmosphere_get_turbulence


    subroutine turbulence_sample(this, j_sample)
        implicit none
        type(atmosphere_type) :: this
        type(json_value), pointer :: j_sample
        character(len=:), allocatable :: fn, analysis_variable
        character(len=100) :: filename, filename_psd_variance
        integer :: i, j, n, n_psd, iunit, psd_mean_unit, psd_variance_unit, psd_var
        real, allocatable :: vals(:,:), psd_mean(:,:), psd_temp(:,:)
        real :: height_ag, sigma, dx, turb(6)
        real :: psd_stdev
        logical :: found
        
        ! test random number generator function
        call test_rand_normal()

        write(*,*) "    Sampling Atmospheric Turbulence..."
        call jsonx_get(j_sample, "save_filename", fn)
        open(newunit=iunit, file=fn, status='REPLACE')
        write(*,*) "    saving sample to ", fn
        call jsonx_get(j_sample, "height_above_ground[ft]", height_ag)
        call jsonx_get(j_sample, "dx[ft]", dx)
        call jsonx_get(j_sample, "number_of_points", n)
        allocate(vals(n,n))
        
        sigma = interpolate_1D(this%turb_height_ag, this%turb_sigma, height_ag)

        write(*,*) "    Altitude = ", height_ag
        write(*,*) "    Turbulence Standard Deviation, sigma = ", sigma

        write(iunit,*)"distance[ft], uprime[ft/s], vprime[ft/s], wprime[ft/s], pprime[rad/s], qprime[rad/s], rprime[rad/s]"
        do i=1,n
            turb(:) = get_turbulence(this, dx, sigma, sigma, sigma)
            write(iunit,*) dx*real(i-1),",", turb(1), ",", turb(2), ",", turb(3), ",", turb(4), ",", turb(5), ",", turb(6)
            vals(i,:) = turb(:)
        end do 
        close(iunit)

        call json_get(j_sample, "psd_analyses", n_psd, found)
        if (found) then
            call jsonx_get(j_sample, "psd_analyses", n_psd)
            call jsonx_get(j_sample, "analysis_variable", analysis_variable)
            allocate(psd_mean(n/2+1,2)) 
            allocate(psd_temp(n/2+1,2)) 
            psd_mean = 0.0
            psd_temp = 0.0

            write(*,*) "    Turbulence Normalized Mean PSD Analysis"
            write(filename,'("psd_mean_analysis_",A,"_",I0,".csv")') trim(analysis_variable), n_psd
            write(filename_psd_variance,'("psd_variance_",A,"_",I0,".csv")') trim(analysis_variable), n_psd
            open(newunit=psd_mean_unit, file=trim(filename), status="REPLACE")
            open(newunit=psd_variance_unit, file=trim(filename_psd_variance), status="REPLACE")
            write(*,*) "    - saving normalized PSD analysis to ", filename
            do j=1,n_psd
                write(*,*) "PSD ", j, " of ", n_psd
                do i=1,n
                    turb(:) = get_turbulence(this, dx, sigma, sigma, sigma)
                    vals(i,:) = turb(:)
                end do
                select case(trim(analysis_variable))
                case("u")
                    psd_var = 1
                case("v")
                    psd_var = 2
                case("w")
                    psd_var = 3
                case("p")
                    psd_var = 4
                case("q")
                    psd_var = 5
                case("r")
                    psd_var = 6
                case default
                    write(*,*) " JSON INPUT ERROR! PSD analysis variable must be u,v,w,p,q, or r.... Quitting"
                    stop
                end select
                call psd(vals(:,psd_var), dx, psd_norm=psd_temp, stdev=psd_stdev) ! 1=u, 2=v, 3=w, 4=p, 5=q, 6=r
                psd_mean(:,2) = psd_mean(:,2) + psd_temp(:,2)/n_psd
                write(psd_variance_unit,*) j,",", psd_stdev**2
            end do

            do i=1,n/2+1
                write(psd_mean_unit,*) psd_temp(i,1)*2*PI,",",psd_mean(i,2)/2/PI
            end do
            close(psd_mean_unit)
        end if

        deallocate(vals)
    end subroutine turbulence_sample


    function get_turbulence(this, dx, sigma_u, sigma_v, sigma_w) result(ans)
        implicit none
        type(atmosphere_type) :: this
        real :: dx, sigma_u, sigma_v, sigma_w
        real :: ans(6)
        
        select case(trim(this%turb_model))
        case("dryden_beal")
            this%xff = this%xff + dx
            ans(:) =  dryden_beal(this, dx, sigma_u, sigma_v, sigma_w)
        end select
    end function get_turbulence
    

    function dryden_beal(this, dx, sigma_u, sigma_v, sigma_w) result(ans)
        implicit none
        type(atmosphere_type) :: this
        real :: dx, sigma_u, sigma_v, sigma_w
        real :: ans(6)
        real :: Au, Av, Aw, Ap, Aq, Ar
        real :: etau, etav, etaw, etap, etaq, etar, f, g
        real :: w_tilde_at_hstab, v_tilde_at_vstab

        Au = 0.5*dx/this%Lu
        Av = 0.25*dx/this%Lv
        Aw = 0.25*dx/this%Lw
        Ap = 0.5*dx/this%Lb

        etau = rand_normal()*sigma_u*sqrt(2.0*this%Lu/dx)
        etav = rand_normal()*sigma_v*sqrt(2.0*this%Lv/dx)
        etaw = rand_normal()*sigma_w*sqrt(2.0*this%Lw/dx)
        etap = rand_normal()*sigma_w*sqrt(0.8*PI*(this%Lw/this%Lb)**(1.0/3.0)/(this%Lw*dx))

        f = ((1.0 - Av)*this%prev_f + 2.0*Av*etav)/(1.0 + Av)
        g = ((1.0 - Aw)*this%prev_g + 2.0*Aw*etaw)/(1.0 + Aw)

        
        ans(1) = ((1.0 - Au)*this%prev_turb(1) + 2.0*Au*etau)/(1.0 + Au)
        ans(2) = ((1.0 - Av)*this%prev_turb(2) + Av*(f + this%prev_f) + sqrt(3.0)*(f - this%prev_f))/(1.0 + Av)
        ans(3) = ((1.0 - Aw)*this%prev_turb(3) + Aw*(g + this%prev_g) + sqrt(3.0)*(g - this%prev_g))/(1.0 + Aw)
        ans(4) = ((1.0 - Ap)*this%prev_turb(4) + 2.0*Ap*etaP)/(1.0 + Ap)
        
        this%xff_hist(1:this%hist_len - 1) = this%xff_hist(2:this%hist_len)
        this%xff_hist(this%hist_len) = this%xff

        this%w_tilde_hist(1:this%hist_len - 1) = this%w_tilde_hist(2:this%hist_len)
        this%w_tilde_hist(this%hist_len) = ans(3)
        
        this%v_tilde_hist(1:this%hist_len - 1) = this%v_tilde_hist(2:this%hist_len)
        this%v_tilde_hist(this%hist_len) = ans(2)

        w_tilde_at_hstab = interpolate_1D(this%xff_hist, this%w_tilde_hist, this%xff - this%hstab_dist)
        v_tilde_at_vstab = interpolate_1D(this%xff_hist, this%v_tilde_hist, this%xff - this%vstab_dist)
        ! write(*,*) ""
        ! write(*,*) "     xff  history", this%xff_hist
        ! write(*,*) "  w tilde history", this%w_tilde_hist
        ! write(*,*) " w tilde at hstab", w_tilde_at_hstab
        
        ans(5) =  (ans(3) - w_tilde_at_hstab)/this%hstab_dist  ! eq. 9.5.18
        ans(6) = -(ans(2) - v_tilde_at_vstab)/this%vstab_dist  ! eq. 9.5.19

        this%prev_f = f
        this%prev_g = g
        this%prev_turb(:) = ans(:)

    end function dryden_beal     



end module atmosphere_m