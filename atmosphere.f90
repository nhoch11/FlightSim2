module atmosphere_m
    use hoch_m
    use jsonx_m
    implicit none

    type atmosphere_type
        real, allocatable :: wind(:)
        character(len=:), allocatable :: turb_model, turb_intensity
        real :: wingspan, hstab_dist, vstab_dist
        logical :: turb_repeatable

        real :: light_height_ag(3) = (/2000., 8000., 17000./)
        real :: light_sigma(3) = (/5., 5., 3./)
        real :: moderate_height_ag(3) = (/2000., 11000., 45000./)
        real :: moderate_sigma(3) = (/10., 10., 3./)
        real :: severe_height_ag(4) = (/2000., 4000., 20000., 80000./)
        real :: severe_sigma(4) = (/15., 21., 21., 3./)

        real, allocatable :: turb_height_ag(:), turb_sigma(:)
        real :: prev_turb(3), prev_xyz(3), prev_f, prev_g
        real :: Lu, Lv, Lw
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
        call jsonx_get(j_atmosphere, "constant_wind[ft]", this%wind, 0.0, 3)
        write(*,*) " Constant Wind [ft/s] = ",this%wind(:)

        call json_get(j_atmosphere, "turbulence", j_turb, found)
        if (found) then
            call jsonx_get(j_turb, "model", this%turb_model, 'none')
            if (this%turb_model .ne. 'none') then
                call jsonx_get(j_turb, "wingspan[ft]", this%wingspan)
                call jsonx_get(j_turb, "hstab_distance[ft]", this%hstab_dist)
                call jsonx_get(j_turb, "vstab_distance[ft]", this%vstab_dist)
                call jsonx_get(j_turb, "intensity", this%turb_intensity)
                call jsonx_get(j_turb, "repeatable", this%turb_repeatable)

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
                end select

                ! initialize values
                this%prev_turb(:) = 0.0
                this%prev_xyz(:) = 0.0
                this%prev_f = 0.0
                this%prev_g = 0.0

                call json_get(j_turb, "sample", j_sample, found)
                if (found) call turbulence_sample(this, j_sample)
            end if
        end if
    end subroutine atmosphere_init


    function atmosphere_get_turbulence(this, states) result(ans)
        implicit none
        type(atmosphere_type) :: this
        real :: states(21)
        real :: ans(3)
        real :: dx, sigma

        dx = sqrt((states(7) - this%prev_xyz(1))**2 + (states(8) - this%prev_xyz(2))**2 + (states(9) - this%prev_xyz(3))**2)
        sigma = interpolate_1D(this%turb_height_ag, this%turb_sigma, -states(9))

        ans(:) = get_turbulence(this, dx, sigma, sigma, sigma)
        
        this%prev_xyz(:) = states(7:9)
    end function atmosphere_get_turbulence


    subroutine turbulence_sample(this, j_sample)
        implicit none
        type(atmosphere_type) :: this
        type(json_value), pointer :: j_sample
        character(len=:), allocatable :: fn 
        integer :: i, j, n, n_psd, iunit, psd_mean_unit
        real, allocatable :: vals(:,:), psd_mean(:,:), psd_temp(:,:)
        real :: height_ag, sigma, dx, turb(3)
        real :: mean, stdev
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

        write(iunit,*)"distance[ft], uprime[ft/s], vprime[ft/s], wprime[ft/s]"
        do i=1,n
            turb(:) = get_turbulence(this, dx, sigma, sigma, sigma)
            write(iunit,*) dx*real(i-1),",", turb(1), ",", turb(2), ",", turb(3)
            vals(i,:) = turb(:)
        end do 
        close(iunit)

        call json_get(j_sample, "psd_analyses", n_psd, found)
        if (found) then
            call jsonx_get(j_sample, "psd_analyses", n_psd)
            allocate(psd_mean(n/2+1,2)) 
            allocate(psd_temp(n/2+1,2)) 
            psd_mean = 0.0
            psd_temp = 0.0

            write(*,*) "    Turbulence Normalized Mean PSD Analysis"
            open(newunit=psd_mean_unit, file="psd_mean_analysis.csv", status="REPLACE")
            write(*,*) "    - saving normalized PSD analysis to psd_mean_analysis.csv"
            do j=1,n_psd
                write(*,*) "PSD ", j, " of ", n_psd
                do i=1,n
                    turb(:) = get_turbulence(this, dx, sigma, sigma, sigma)
                    vals(i,:) = turb(:)
                end do
                call psd(vals(:,1), dx, psd_norm=psd_temp) ! 1=u, 2=v, 3=w
                psd_mean(:,2) = psd_mean(:,2) + psd_temp(:,2)/n_psd
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
        real :: ans(3)
        
        select case(trim(this%turb_model))
        case("dryden_beal")
            ans(:) =  dryden_beal(this, dx, sigma_u, sigma_v, sigma_w)
        end select
    end function get_turbulence
    

    function dryden_beal(this, dx, sigma_u, sigma_v, sigma_w) result(ans)
        implicit none
        type(atmosphere_type) :: this
        real :: dx, sigma_u, sigma_v, sigma_w
        real :: ans(3)
        real :: Au, Av, Aw, etau, etav, etaw, f, g

        Au = 0.5*dx/this%Lu
        Av = 0.25*dx/this%Lv
        Aw = 0.25*dx/this%Lw

        etau = rand_normal()*sigma_u*sqrt(2.0*this%Lu/dx)
        etav = rand_normal()*sigma_v*sqrt(2.0*this%Lv/dx)
        etaw = rand_normal()*sigma_w*sqrt(2.0*this%Lw/dx)

        f = ((1.0 - Av)*this%prev_f + 2.0*Av*etav)/(1.0 + Av)
        g = ((1.0 - Aw)*this%prev_g + 2.0*Aw*etaw)/(1.0 + Aw)
        
        ans(1) = ((1.0 - Au)*this%prev_turb(1) + 2.0*Au*etau)/(1.0 + Au)
        ans(2) = ((1.0 - Av)*this%prev_turb(2) + Av*(f + this%prev_f) + sqrt(3.0)*(f - this%prev_f))/(1.0 + Av)
        ans(3) = ((1.0 - Aw)*this%prev_turb(3) + Aw*(g + this%prev_g) + sqrt(3.0)*(g - this%prev_g))/(1.0 + Aw)
        
        this%prev_f = f
        this%prev_g = g
        this%prev_turb(:) = ans(:)

    end function dryden_beal     



end module atmosphere_m