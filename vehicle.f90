module vehicle_m
    use hoch_m
    use jsonx_m
    use micro_time_m
    use linalg_mod
    use connection_m
    ! use json_xtnsn_mod
    
    implicit none

    type stall_settings_type
        real :: alpha_0, alpha_s, lambda_b, min_val 
    end type stall_settings_type

    type trim_settings_type
        ! old trim settings
        logical :: has_sideslip, has_gamma, verbose
        character(:), allocatable :: type
        real :: fd_step, r_factor, tol
        integer :: max_iter
        real :: sideslip0, gamma0
        real :: u, v, w, p, q, r
        real :: phi0, theta0, psi0
    end type trim_settings_type
        
    type vehicle_type
        
        type(json_value), pointer :: j_vehicle
        
        character(len=:), allocatable :: name
        character(len=:), allocatable :: type
        character(100) :: states_filename, rk4_filename
        
        integer :: iunit_states, iunit_rk4, iunit_trim
        logical :: run_physics
        logical :: trimming
        
        ! mass constants
        real :: weight, mass
        real, dimension(3,3) :: I, I_inv, h
        real, dimension(6) :: FM
        real :: Ixx0, Iyy0, Izz0, Ixz0, Ixy0, Iyz0, hx0, hy0, hz0
        
        ! aero constants
        real, dimension(:), allocatable :: aero_ref_loc
        real :: ref_area, ref_len, ref_long, ref_lat
        real :: CL0, CLa, CLahat, CLqbar, CLde
        real :: CD0, CD1, CD2, CDS2, CDqbar, CDaqbar, CDde, CDade, CDde2
        real :: CSb, CSpbar, CSapbar, CSrbar, CSda, CSdr
        real :: Cll0, Clb, Clpbar, Clrbar, Clarbar, Clda, Cldr
        real :: Cm0, Cma, Cmqbar, Cmahat, Cmde
        real :: Cnb, Cnpbar, Cnapbar, Cnrbar, Cnda, Cnada, Cndr
        
        
        ! stall model constants
        logical :: include_stall
        type(stall_settings_type) :: CL_stall, CD_stall, Cm_stall

        ! thrust
        real :: thrust_T0, thrust_Ta
        real, dimension(:), allocatable :: thrust_loc
        real, dimension(4) :: thrust_quat

        ! intitialization constants
        character(len=:), allocatable :: init_type
        real, dimension(:), allocatable :: init_eul
        real, dimension(13) :: init_state
        real :: init_V, init_alt
        real :: rho0

        ! variables
        real, dimension(13) :: state
        real, dimension(4) :: controls
        real :: lat, long

        type(trim_settings_type) :: trim

        
        ! old misc settings
        type(connection) :: graphics, control_udp
        real, dimension(3) :: Vwf

    end type vehicle_type
        

contains 

    subroutine vehicle_init(this, vehicle_json)
        implicit none
        type(vehicle_type) :: this 
        type(json_value), pointer :: vehicle_json
        real, dimension(:), allocatable :: thrust_orientation, hxyz
        real :: det
        real :: junk1, junk2, junk3, junk4, junk5

        this%j_vehicle => vehicle_json
        this%name = this%j_vehicle%name

        call jsonx_get(this%j_vehicle, "type", this%type)
        call jsonx_get(this%j_vehicle, "run_physics", this%run_physics)
        
        write(*,*) "Initializing ", this%name
        write(*,*) "    - type = ", this%type
        write(*,*) "    - run physics = ", this%run_physics
        
        if (this%run_physics) then 
            if (save_states) then
                this%states_filename = trim(this%name)//"_states.csv"
                open(newunit=this%iunit_states, file=this%states_filename, status="REPLACE")
                write(this%iunit_states,*) "time[s],u[ft/s],v[ft/s],w[ft/s],p[rad/s],q[rad/s],r[rad/s],x[ft],y[ft],z[ft],e0,ex,ey,ez"
                write(*,*) "    - states will be written to ", this%states_filename
            end if 
            
            if (rk4_verbose) then
                this%rk4_filename = trim(this%name)//"_RK4.csv"
                open(newunit=this%iunit_rk4, file=this%rk4_filename, status="REPLACE")
                write(*,*) "    - RK4 results will be written to ", this%states_filename
            end if
            
            write(*,*) "    - mass"
            call jsonx_get(this%j_vehicle, "mass.weight[lbf]",  this%weight)
            call jsonx_get(this%j_vehicle, "mass.Ixx[slug-ft^2]",  this%I(1,1))
            call jsonx_get(this%j_vehicle, "mass.Iyy[slug-ft^2]",  this%I(2,2))
            call jsonx_get(this%j_vehicle, "mass.Izz[slug-ft^2]",  this%I(3,3))
            call jsonx_get(this%j_vehicle, "mass.Ixy[slug-ft^2]",  this%I(1,2), 0.0)
            call jsonx_get(this%j_vehicle, "mass.Ixz[slug-ft^2]",  this%I(1,3), 0.0)
            call jsonx_get(this%j_vehicle, "mass.Iyz[slug-ft^2]",  this%I(2,3), 0.0)
            
            ! calc mass
            this%mass = this%weight/gravity_English(0.0)
            write(*,*) "        - weight[lbf] = ", this%weight
            write(*,*) "        - mass[slug]  = ", this%mass
            
            this%I(1,2) = -this%I(1,2)
            this%I(1,3) = -this%I(1,3)
            this%I(2,1) =  this%I(1,2) ! this is negative at this point
            this%I(2,3) = -this%I(2,3)
            this%I(3,1) =  this%I(1,3) ! this is negative at this point
            this%I(3,2) =  this%I(2,3) ! this is negative at this point
            
            ! calculate I_inv by hand
            det = this%I(1,1)*(this%I(2,2)*this%I(3,3) - this%I(2,3)*this%I(3,2)) &
                - this%I(1,2)*(this%I(2,1)*this%I(3,3) - this%I(2,3)*this%I(3,1)) &
                + this%I(1,3)*(this%I(2,1)*this%I(3,2) - this%I(2,2)*this%I(3,1))
                
            if (det < 0) then
                ! The matrix is singular, cannot be inverted
                write(*,*) "Singular inertia matrix.... quitting"
                stop
            end if
            
            this%I_inv(1,1) = (this%I(2,2)*this%I(3,3) - this%I(2,3)*this%I(3,2)) / det;
            this%I_inv(1,2) = (this%I(1,3)*this%I(3,2) - this%I(1,2)*this%I(3,3)) / det;
            this%I_inv(1,3) = (this%I(1,2)*this%I(2,3) - this%I(1,3)*this%I(2,2)) / det;
            this%I_inv(2,1) = (this%I(2,3)*this%I(3,1) - this%I(2,1)*this%I(3,3)) / det;
            this%I_inv(2,2) = (this%I(1,1)*this%I(3,3) - this%I(1,3)*this%I(3,1)) / det;
            this%I_inv(2,3) = (this%I(1,1)*this%I(2,3) - this%I(1,3)*this%I(2,1)) / det;
            this%I_inv(3,1) = (this%I(2,1)*this%I(3,2) - this%I(2,2)*this%I(3,1)) / det;
            this%I_inv(3,2) = (this%I(1,2)*this%I(3,1) - this%I(1,1)*this%I(3,2)) / det;
            this%I_inv(3,3) = (this%I(1,1)*this%I(2,2) - this%I(1,2)*this%I(2,1)) / det;
            
            
            call jsonx_get(this%j_vehicle, "mass.h[slug-ft^2/s]",  hxyz, 0.0, 3)
            
            this%h = 0.
            this%h(1,2) = -hxyz(3)
            this%h(1,3) =  hxyz(2)
            this%h(2,1) =  hxyz(3)
            this%h(2,3) = -hxyz(1)
            this%h(3,1) = -hxyz(2)
            this%h(3,2) =  hxyz(1)
            
            if (verbose) then
                write(*,*) "        - I matrix = ", this%I
                write(*,*) "        - I_inv    = ", this%I_inv
                write(*,*) "        - h matrix = ", this%h
            end if
            
            
            write(*,*) "    - aerodynamics"
            call jsonx_get(this%j_vehicle, "aerodynamics.reference.area[ft^2]",  this%ref_area)
            call jsonx_get(this%j_vehicle, "aerodynamics.reference.longitudinal_length[ft]",  this%ref_long)
            call jsonx_get(this%j_vehicle, "aerodynamics.reference.lateral_length[ft]",  this%ref_lat)
            call jsonx_get(this%j_vehicle, "aerodynamics.reference.relative_location[ft]",  this%aero_ref_loc, 0.0, 3)
            
            if (this%type == "arrow") then
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CL.alpha", this%CLa)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.L0", this%CD0)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.CL1_CL1", this%CD2)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.0", this%Cll0)
                ! call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.1", this%Cl1)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.pbar", this%Clpbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cm.alpha", this%Cma)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cm.qbar", this%Cmqbar)
            end if
            
            if (this%type == "aircraft") then
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CL.0", this%CL0)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CL.alpha", this%CLa)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CL.alphahat", this%CLahat)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CL.qbar", this%CLqbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CL.elevator", this%CLde)
                
                ! write(*,*) ""
                ! write(*,*) "    CL 0 = ", this%CL0
                ! write(*,*) "    CL a = ", this%CLa
                ! write(*,*) "    CL ahat = ", this%CLahat
                ! write(*,*) "    CL qbar = ", this%CLqbar
                ! write(*,*) "    CL de = ", this%CLde
                
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CS.beta", this%CSb)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CS.pbar", this%CSpbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CS.alpha_pbar", this%CSapbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CS.rbar", this%CSrbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CS.aileron", this%CSda)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CS.rudder", this%CSdr)
                
                ! write(*,*) ""
                ! write(*,*) "    CS beta = ", this%CSb
                ! write(*,*) "    CS pbar = ", this%CSpbar
                ! write(*,*) "    CS apbar = ", this%CSapbar
                ! write(*,*) "    CS rbar = ", this%CSrbar
                ! write(*,*) "    CS da = ", this%CSda
                ! write(*,*) "    CS dr = ", this%CSdr
                
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.L0", this%CD0)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.CL1", this%CD1)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.CL1_CL1", this%CD2)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.CS_CS", this%CDS2)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.qbar", this%CDqbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.alpha_qbar", this%CDaqbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.elevator", this%CDde)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.alpha_elevator", this%CDade)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.elevator_elevator", this%CDde2)
                
                ! write(*,*) ""
                ! write(*,*) "    CD L0 = ", this%CD0
                ! write(*,*) "    CD CL1 = ", this%CD1
                ! write(*,*) "    CD CL2 = ", this%CD2
                ! write(*,*) "    CD S2 = ", this%CDS2
                ! write(*,*) "    CD qbar = ", this%CDqbar
                ! write(*,*) "    CD aqbar = ", this%CDaqbar
                ! write(*,*) "    CD de = ", this%CDde
                ! write(*,*) "    CD ade = ", this%CDade
                ! write(*,*) "    CD de2 = ", this%CDde2
                
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.beta", this%Clb)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.pbar", this%Clpbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.rbar", this%Clrbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.alpha_rbar", this%Clarbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.aileron", this%Clda)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.rudder", this%Cldr)
                
                ! write(*,*) ""
                ! write(*,*) "    Cl beta = ", this%Clb
                ! write(*,*) "    Cl pbar = ", this%Clpbar
                ! write(*,*) "    Cl rbar = ", this%Clrbar
                ! write(*,*) "    Cl arbar = ", this%Clarbar
                ! write(*,*) "    Cl da = ", this%Clda
                ! write(*,*) "    Cl dr = ", this%Cldr
                
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cm.0", this%Cm0)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cm.alpha", this%Cma)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cm.qbar", this%Cmqbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cm.alphahat", this%Cmahat)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cm.elevator", this%Cmde)
                
                ! write(*,*) ""
                ! write(*,*) "    Cm 0 = ", this%Cm0
                ! write(*,*) "    Cm a = ", this%Cma
                ! write(*,*) "    Cm qbar = ", this%Cmqbar
                ! write(*,*) "    Cm ahat = ", this%Cmahat
                ! write(*,*) "    Cm de = ", this%Cmde
                
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.beta", this%Cnb)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.pbar", this%Cnpbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.alpha_pbar", this%Cnapbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.rbar", this%Cnrbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.aileron", this%Cnda)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.alpha_aileron", this%Cnada)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.rudder", this%Cndr)
                
                ! write(*,*) ""
                ! write(*,*) "    Cn beta = ", this%Cnb
                ! write(*,*) "    Cn pbar = ", this%Cnpbar
                ! write(*,*) "    Cn apbar = ", this%Cnapbar
                ! write(*,*) "    Cn rbar = ", this%Cnrbar
                ! write(*,*) "    Cn da = ", this%Cnda
                ! write(*,*) "    Cn ada = ", this%Cnada
                ! write(*,*) "    Cn dr = ", this%Cndr

            end if

            ! stall
            if (this%type == "arrow" .or. this%type == "aircraft") then
                call jsonx_get(this%j_vehicle, "aerodynamics.stall.include_stall", this%include_stall)
                if (this%include_stall) then
                    ! CL stall
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.CL.alpha_0[deg]", this%CL_stall%alpha_0)
                    this%CL_stall%alpha_0 = this%CL_stall%alpha_0*PI/180.0
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.CL.alpha_s[deg]", this%CL_stall%alpha_s)
                    this%CL_stall%alpha_s = this%CL_stall%alpha_s*PI/180.0
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.CL.lambda_b", this%CL_stall%lambda_b)
                    
                    ! CD stall
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.CD.alpha_0[deg]", this%CD_stall%alpha_0)
                    this%CD_stall%alpha_0 = this%CD_stall%alpha_0*PI/180.0
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.CD.alpha_s[deg]", this%CD_stall%alpha_s)
                    this%CD_stall%alpha_s = this%CD_stall%alpha_s*PI/180.0
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.CD.lambda_b", this%CD_stall%lambda_b)
                    
                    ! Cm stall
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.Cm.min", this%Cm_stall%min_val)
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.Cm.alpha_0[deg]", this%Cm_stall%alpha_0)
                    this%Cm_stall%alpha_0 = this%Cm_stall%alpha_0*PI/180.0
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.Cm.alpha_s[deg]", this%Cm_stall%alpha_s)
                    this%Cm_stall%alpha_s = this%Cm_stall%alpha_s*PI/180.0
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.Cm.lambda_b", this%Cm_stall%lambda_b)
                end if
            end if
            
            call jsonx_get(this%j_vehicle, "thrust.T0[lbf]", this%thrust_T0, 0.0)
            call jsonx_get(this%j_vehicle, "thrust.Ta", this%thrust_Ta, 0.0)
            call jsonx_get(this%j_vehicle, "thrust.locationp[ft]", this%thrust_loc, 0.0, 3)
            call jsonx_get(this%j_vehicle, "thrust.locationp[ft]", thrust_orientation, 0.0, 3)
            thrust_orientation = thrust_orientation*PI/180.0
            this%thrust_quat = euler_to_quat(thrust_orientation)


            ! initial conditions
            this%init_state = 0.
            this%controls = 0.

            call jsonx_get(this%j_vehicle, "initial.airspeed[ft/s]",  this%init_V)
            call jsonx_get(this%j_vehicle, "initial.altitude[ft]",  this%init_alt)
            this%init_state(9) = -this%init_alt
            call jsonx_get(this%j_vehicle, "initial.latitude[deg]",  this%lat)
            this%lat = this%lat*PI/180.0
            call jsonx_get(this%j_vehicle, "initial.longitude[deg]",  this%long)
            this%long = this%long*PI/180.0
            call jsonx_get(this%j_vehicle, "initial.Euler_angles[deg]",  this%init_eul, 0.0, 3)
            this%init_eul = this%init_eul*PI/180.0
          

            call jsonx_get(this%j_vehicle, "initial.type", this%init_type)
            
            if (this%init_type == "state") then
                call init_to_state(this)
            else if (this%init_type == "trim") then
                call init_to_trim(this)
            else
                write(*,*) " Invalid initial type in json. must be trim or state. Quitting..."
                stop
            end if

            this%init_state(10:13) = euler_to_quat(this%init_eul)
            this%state = this%init_state

            if (save_states) then
                call vehicle_write_state(this, 0.0, this%state)
            end if

            ! get rho, mu, Reynolds
            call std_atm_English(0.0,junk1,junk2,junk3,this%rho0,junk4,junk5)

        end if ! run_physics

        ! connections
        ! call jsonx_get(j_main, "connections", j_connections)
        ! call jsonx_get(j_connections, "graphics", j_graphics)
        ! call graphics%init(j_graphics)

    end subroutine vehicle_init


    subroutine init_to_state(this)

        implicit none
        type(vehicle_type) :: this 
        real :: p0, q0, r0
        real :: alpha0, beta0, da0, de0, dr0, tau0, rho0

        this%trimming = .false.

        call jsonx_get(this%j_vehicle, "initial.state.angle_of_attack[deg]", alpha0)
        call jsonx_get(this%j_vehicle, "initial.state.sideslip_angle[deg]", beta0)
        alpha0 = alpha0*PI/180.
        beta0 = beta0*PI/180.

        this%init_state(1)  = this%init_V*cos(alpha0)*cos(beta0)   ! u
        this%init_state(2)  = this%init_V*sin(beta0)               ! v
        this%init_state(3)  = this%init_V*sin(alpha0)*cos(beta0)   ! w

        call jsonx_get(this%j_vehicle, "initial.state.p[deg/s]",  this%init_state(4))
        call jsonx_get(this%j_vehicle, "initial.state.q[deg/s]",  this%init_state(5))
        call jsonx_get(this%j_vehicle, "initial.state.r[deg/s]",  this%init_state(6))
        this%init_state(4) = this%init_state(4)*PI/180. ! p
        this%init_state(5) = this%init_state(5)*PI/180. ! q
        this%init_state(6) = this%init_state(6)*PI/180. ! r
        
        call jsonx_get(this%j_vehicle, "initial.state.aileron[deg]",  this%controls(1), 0.0)
        call jsonx_get(this%j_vehicle, "initial.state.elevator[deg]",  this%controls(2), 0.0)
        call jsonx_get(this%j_vehicle, "initial.state.rudder[deg]",  this%controls(3), 0.0)
        call jsonx_get(this%j_vehicle, "initial.state.throttle",  this%controls(4), 0.0)
        this%controls(1) = this%controls(1)*PI/180. ! da
        this%controls(2) = this%controls(2)*PI/180. ! de    
        this%controls(3) = this%controls(3)*PI/180. ! dr

    end subroutine init_to_state


    subroutine init_to_trim(this)

        implicit none
        type(vehicle_type) :: this 

        this%trimming = .true.

        ! read in trim settings
        call jsonx_get(this%j_vehicle, "initial.trim.type", this%trim%type)

        this%trim%has_sideslip = .false.
        call json_get(this%j_vehicle, "initial.trim.sideslip_angle[deg]", this%trim%sideslip0, this%trim%has_sideslip)
        if (this%trim%has_sideslip) then
            this%trim%sideslip0 = this%trim%sideslip0*PI/180.
        end if

        this%trim%has_gamma = .false.
        call json_get(this%j_vehicle, "initial.trim.ref_climb_angle[deg]", this%trim%gamma0, this%trim%has_gamma)
        write(*,*) " has gamma = ", this%trim%has_gamma
        if (this%trim%has_gamma) then
            this%trim%gamma0 = this%trim%gamma0*PI/180.
        end if
        
        this%trim%phi0 = this%init_eul(1)
        this%trim%theta0 = this%init_eul(2)
        this%trim%psi0 = this%init_eul(3)

        call jsonx_get(this%j_vehicle, "initial.trim.solver.finite_difference_step_size", this%trim%fd_step)
        call jsonx_get(this%j_vehicle, "initial.trim.solver.relaxation_factor", this%trim%r_factor)
        call jsonx_get(this%j_vehicle, "initial.trim.solver.tolerance", this%trim%tol)
        call jsonx_get(this%j_vehicle, "initial.trim.solver.max_iterations", this%trim%max_iter)
        call jsonx_get(this%j_vehicle, "initial.trim.solver.verbose", this%trim%verbose)
        
        ! jump into trim

        this%trimming = .false.
    
    end subroutine init_to_trim


    ! function calc_R(p,q,r, theta, G) result(R_)

    !     implicit none 

    !     real, intent(in) :: p, q, r, theta
    !     real, dimension(6), intent(in) :: G
    !     real, dimension(6) :: R_
    !     real, dimension(13) :: y, dydt 
    !     real :: alpha, beta, phi, psi

    !     alpha = G(1)
    !     ! theta = theta0
    !     psi = psi0
    !     controls(1) = G(3)
    !     controls(2) = G(4)
    !     controls(3) = G(5)
    !     controls(4) = G(6)

    !     if (trim_type == "shss" .and. has_sideslip) then
    !         phi = G(2)
    !         beta = sideslip0 

    !     else
    !         beta = G(2)
    !         phi = phi0
    !     end if
        
    !     y(1)  = V0*cos(alpha)*cos(beta)   ! u
    !     y(2)  = V0*sin(beta)               ! v
    !     y(3)  = V0*sin(alpha)*cos(beta)   ! w
    !     y(4)  = p             ! p
    !     y(5)  = q              ! q
    !     y(6)  = r              ! r
    !     y(7)  = 0.0             ! xf
    !     y(8)  = 0.0             ! yf
    !     y(9)  = zf0             ! zf
    !     y(10) = phi            ! phi
    !     y(11) = theta          ! theta
    !     y(12) = psi            ! psi

        
    !     !  convert phi, theta, and psi to a quat
    !     y(10:13) = euler_to_quat(y(10:12))
        
    !     !  normalize the quat
    !     call quat_norm(y(10:13))
        
    !     ! write(*,*) "before diff_eq"
    !     if (verbose2) then
    !         write(*,*) "y going into diff_eq in calc_R", y
    !     end if
        
    !     ! write(*,*)"y = ",y
    !     dydt = diff_eq(y)
    !     ! write(*,*) "after diff_eq"
    !     ! write(*,*)"dydt = ",dydt
    !     R_ = dydt(1:6)


    ! end function calc_R


    ! function trim_solver() result(ans)

    !     implicit none

    !     real :: alpha, beta, u, v, w, p, q, r, phi, theta, psi, da, de, dr, tau, eps, grav, c1, pw
    !     real :: a, b, c, theta_plus, theta_minus, err_theta_plus, err_theta_minus
    !     real, dimension(13) :: ans
    !     real, dimension(6) :: G, R_plus, R_minus, R_neg, R_
    !     real,dimension(:),allocatable :: delta_G
    !     real, dimension(6,6) :: jac
    !     integer :: i,j, iter

    !     p_trim = 0.0
    !     q_trim = 0.0
    !     r_trim = 0.0
    !     ! theta = theta0
    !     psi = psi0
    !     controls = 0.0
    !     G = 0.0

    !     alpha = G(1)
            
    !     if (trim_type == "shss" .and. has_sideslip) then
    !         phi = G(2)
    !         beta = sideslip0 
    !     else
    !         beta = G(2)
    !         phi = phi0
    !     end if

    !     eps = 1.0

    !     if (verbose) then
    !         write(*,*) ""
    !         if (trim_type == "shss") then

    !             if (has_sideslip) then
    !                 write(*,*) "Solving for trim condition: Steady Heading Sideslip"
    !                 if (has_gamma) then
    !                     write(*,*) "beta[deg] = ", sideslip0*180./pi, ",  climb_angle[deg] = ", gamma0*180./pi, &
    !                     ",  psi[deg] = ", psi*180./pi
    !                 else
    !                     write(*,*) "beta[deg] = ", sideslip0*180./pi, ",  theta[deg] = ", theta*180./pi, &
    !                     ",  psi[deg] = ", psi*180./pi
    !                 end if
    !             else
    !                 write(*,*) "Solving for trim condition: Steady Heading Sideslip"
    !                 if (has_gamma) then
    !                     write(*,*) "phi[deg] = ", phi*180./pi, ",  climb_angle[deg] = ", gamma0*180./pi, &
    !                     ",  psi[deg] = ", psi*180./pi
    !                 else
    !                     write(*,*) "phi[deg] = ", phi*180./pi, ",  theta[deg] = ", theta*180./pi, &
    !                     ",  psi[deg] = ", psi*180./pi
    !                 end if
    !             end if

    !         else if (trim_type == "sct") then
    !             write(*,*) "Solving for trim condition: Steady Coordinated Turn"
    !             if (has_gamma) then
    !                 write(*,*) "phi[deg] = ", phi*180./pi, ",  climb_angle[deg] = ", gamma0*180./pi, &
    !                 ",  psi[deg] = ", psi*180./pi
    !             else
    !                 write(*,*) "phi[deg] = ", phi*180./pi, ",  theta[deg] = ", theta*180./pi, &
    !                 ",  psi[deg] = ", psi*180./pi
    !             end if
    !         else 
    !             write(*,*) "Solving for trim condition: Vertical Barrel Roll"
    !             if (has_gamma) then
    !                 write(*,*) "phi[deg] = ", phi*180./pi, ",  climb_angle[deg] = ", gamma0*180./pi, &
    !                 ",  psi[deg] = ", psi*180./pi
    !             else
    !                 write(*,*) "phi[deg] = ", phi*180./pi, ",  theta[deg] = ", theta*180./pi, &
    !                 ",  psi[deg] = ", psi*180./pi
    !             end if
    !         end if
    !         write(*,*)""
    !         ! write(*,*) "initial eps = ", eps
    !     end if

    !     grav = gravity_English(-zf0)
    !     iter = 1
        
    !     do while(eps > trim_tol .and. iter < trim_max_iter)
    !         if (verbose2) then
    !             write(*,*)""
    !             write(*,*) "iter = ", iter
    !         end if
    !         ! limit throttle to 0 to 1
    !         if (G(6) < 0.0) then
    !             G(6) = 0.0
    !             if (verbose2) then
    !                 write(*,*)""
    !                 write(*,*) " Overwriting negative throttle"
    !             end if 
    !         else if (G(6) > 1.0) then
    !             G(6) = 1.0
    !             if (verbose2) then
    !                 write(*,*)""
    !                 write(*,*) " Limiting throttle to 100 percent"
    !             end if
    !         end if

    !         alpha = G(1)
            
    !         if (trim_type == "shss" .and. has_sideslip) then
    !             phi = G(2)
    !             beta = sideslip0 
    !         else
    !             beta = G(2)
    !             phi = phi0
    !         end if
            
    !         u  = V0*cos(alpha)*cos(beta)   ! u
    !         v  = V0*sin(beta)               ! v
    !         w  = V0*sin(alpha)*cos(beta)

    !         if (has_gamma) then
    !             a = u*V0*sin(gamma0)
    !             b = (v*sin(phi) + w*cos(phi))*sqrt(u*u + (v*sin(phi) + w*cos(phi))**2 - V0*V0*sin(gamma0)*sin(gamma0))
    !             c = u*u + (v*sin(phi) + w*cos(phi))**2
    !             theta_plus = asin((a + b)/c)
    !             theta_minus = asin((a - b)/c)
    !             err_theta_plus = abs(u*sin(theta_plus) - (v*sin(phi) + w*cos(phi))*cos(theta_plus) - V0*sin(gamma0))
    !             err_theta_minus = abs(u*sin(theta_minus) - (v*sin(phi) + w*cos(phi))*cos(theta_minus) - V0*sin(gamma0))
    !             ! write(*,*) "err plus = ", err_theta_plus
    !             ! write(*,*) "err minus = ", err_theta_minus
    !             if (err_theta_plus < err_theta_minus) then
    !                 theta = theta_plus
    !             else 
    !                 theta = theta_minus
    !             end if
    !             if (verbose2) then
    !                 write(*,*)""
    !                 write(*,'(A, 1(1x,ES20.12))') "      theta_1[deg] = ", theta_plus*180./pi
    !                 write(*,'(A, 1(1x,ES20.12))') "      theta_2[deg] = ", theta_minus*180./pi
    !                 write(*,'(A, 1(1x,ES20.12))') "correct theta[deg] = ", theta*180./pi
    !                 write(*,'(A, 1(1x,ES20.12))') "correct theta[rad] = ", theta
    !                 write(*,*)""
    !             end if 
    !         else
    !             theta = theta0
    !         end if
            
    !         c1 = grav*sin(phi)*cos(theta)/(u*cos(theta)*cos(phi) + w*sin(theta))
            
    !         if (trim_type == "sct") then ! eq 6.2.3
    !             ! write(*,*) "phi", phi
    !             ! write(*,*) "theta", theta
    !             ! write(*,*) "u", u
    !             ! write(*,*) "v", v
    !             ! write(*,*) "w", w
    !             ! write(*,*) "grav", grav
    !             ! write(*,*) "alt", -zf0
    !             p_trim = -c1*sin(theta)
    !             q_trim =  c1*sin(phi)*cos(theta)
    !             r_trim =  c1*cos(phi)*cos(theta)
    !         end if
            
    !         if (trim_type == "vbr") then
    !             pw = c1*(-u*sin(theta) + v*sin(phi)*cos(theta) + w*cos(phi)*cos(theta))/V0
    !             p_trim = pw*u/V0
    !             q_trim = pw*v/V0
    !             r_trim = pw*w/V0
    !         end if 


    !         if (verbose2) then
    !             write(*,*) ""
    !             write(*,*) "Updating p,q,r"
    !             write(*,'(A, 1(1x,ES20.12))') "p[deg/s] = ", p_trim*180./pi
    !             write(*,'(A, 1(1x,ES20.12))') "q[deg/s] = ", q_trim*180./pi
    !             write(*,'(A, 1(1x,ES20.12))') "r[deg/s] = ", r_trim*180./pi
    !         end if
            
    !         if (verbose2) then
    !             write(*,*) ""
    !             write(*,*) "G with Throttle forced to be 0-1"
    !             write(*,'(A, 6(1x,ES20.12))') "    G = ", G
    !         end if
    !         R_ = calc_R(p_trim, q_trim, r_trim, theta, G)
    !         if (verbose2) then
    !             write(*,'(A, 6(1x,ES20.12))') "   R_ = ", R_
    !         end if

    !         do j=1,6
    !             G(j) = G(j) + trim_fd_step
    !             if (verbose2) then
    !                 write(*,*) ""
    !                 write(*,*) "                 j = ", j
    !                 write(*,'(A, 6(1x,ES20.12))') "   G_up = ", G
    !             end if
    !             R_plus = calc_R(p_trim, q_trim, r_trim, theta, G)

    !             if (verbose2) then
    !                 write(*,'(A, 6(1x,ES20.12))') " R_plus = ", R_plus
    !             end if

    !             G(j) = G(j) - 2*trim_fd_step
    !             if (verbose2) then
    !                 write(*,*) ""
    !                 write(*,'(A, 6(1x,ES20.12))') "   G_dn = ", G
    !             end if
    !             R_minus = calc_R(p_trim, q_trim, r_trim, theta, G)

    !             if (verbose2) then
    !                 write(*,'(A, 6(1x,ES20.12))') "R_minus = ", R_minus
    !             end if

    !             do i = 1,6
    !                 jac(i,j) = (R_plus(i) - R_minus(i))/(2*trim_fd_step)
    !             end do
    !             G(j) = G(j) + trim_fd_step
                
    !         end do

    !         if (verbose2) then
    !             write(*,*) ""
    !             write(*,*) "Jacobian = "
    !             do i = 1, 6
    !                 write(*,'(6ES20.12)') jac(i, :)
    !             end do
    !         end if

    !         R_neg = -1*R_
    !         if (verbose2) then
    !             write(*,*) ""
    !             write(*,*) "G used to calculate R_neg which goes into lu_solve:"
    !             write(*,'(A, 6(1x,ES20.12))') "      G  = ", G
    !             write(*,'(A, 6(1x,ES20.12))') "  R_neg  = ", R_neg
    !         end if
    !         call lu_solve(6, jac, R_neg, delta_G) 
            


    !         G = G + trim_r_factor*delta_G
    !         if (verbose2) then
    !             write(*,*) ""
    !             write(*,'(A, 6(1x,ES20.12))') " Delta G = ", delta_G
    !             write(*,'(A, 6(1x,ES20.12))') "   G new = ", G
    !         end if
    !         R_ = calc_R(p_trim, q_trim, r_trim, theta, G)
            
    !         eps = norm2(R_)

            
    !         if (verbose2) then
    !             write(*,*) ""
    !             write(*,'(A, 6(1x,ES20.12))') "      R_ = ", R_
    !             write(*,'(A, 6(1x,ES20.12))') "     eps = ", eps
    !             write(*,*) "iter ", iter, " complete"
    !         end if
            
    !         iter = iter+1
    !     end do
        
        
        
    !     alpha = G(1)
        
    !     if (trim_type == "shss" .and. has_sideslip) then
    !         phi = G(2)
    !         beta = sideslip0 
    !     else
    !         beta = G(2)
    !         phi = phi0
    !     end if
            
    !     u  = V0*cos(alpha)*cos(beta)   ! u
    !     v  = V0*sin(beta)               ! v
    !     w  = V0*sin(alpha)*cos(beta)

    !     if (has_gamma) then
    !         a = u*V0*sin(gamma0)
    !         b = (v*sin(phi) + w*cos(phi))*sqrt(u*u + (v*sin(phi) + w*cos(phi))**2 - V0*V0*sin(gamma0)*sin(gamma0))
    !         c = u*u + (v*sin(phi) + w*cos(phi))**2
    !         theta_plus = asin((a + b)/c)
    !         theta_minus = asin((a - b)/c)
    !         err_theta_plus = abs(u*sin(theta_plus) - (v*sin(phi) + w*cos(phi))*cos(theta_plus) - V0*sin(gamma0))
    !         err_theta_minus = abs(u*sin(theta_minus) - (v*sin(phi) + w*cos(phi))*cos(theta_minus) - V0*sin(gamma0))
    !         ! write(*,*) "err plus = ", err_theta_plus
    !         ! write(*,*) "err minus = ", err_theta_minus
    !         if (err_theta_plus < err_theta_minus) then
    !             theta = theta_plus
    !         else 
    !             theta = theta_minus
    !         end if
    !         if (verbose2) then
    !             write(*,*)""
    !             write(*,'(A, 1(1x,ES20.12))') "      theta_1[deg] = ", theta_plus*180./pi
    !             write(*,'(A, 1(1x,ES20.12))') "      theta_2[deg] = ", theta_minus*180./pi
    !             write(*,'(A, 1(1x,ES20.12))') "correct theta[deg] = ", theta*180./pi
    !             write(*,'(A, 1(1x,ES20.12))') "correct theta[rad] = ", theta
    !             write(*,*)""
    !         end if 
    !         else
    !             theta = theta0
    !     end if  

    !     c1 = grav*sin(phi)*cos(theta)/(u*cos(theta)*cos(phi) + w*sin(theta))

    !     if (trim_type == "sct") then
    !         p_trim = -c1*sin(theta)
    !         q_trim =  c1*sin(phi)*cos(theta)
    !         r_trim =  c1*cos(phi)*cos(theta)
    !     end if

        
    !     if (trim_type == "vbr") then
    !         pw = c1*(-u*sin(theta) + v*sin(phi)*cos(theta) + w*cos(phi)*cos(theta))/V0
    !         p_trim = pw*u/V0
    !         q_trim = pw*v/V0
    !         r_trim = pw*w/V0
    !     end if 
        
        
    !     ans(1)  = V0*cos(alpha)*cos(beta)   ! u
    !     ans(2)  = V0*sin(beta)               ! v
    !     ans(3)  = V0*sin(alpha)*cos(beta)   ! w
    !     ans(4)  = p_trim              ! p
    !     ans(5)  = q_trim              ! q
    !     ans(6)  = r_trim              ! r
    !     ans(7)  = 0.0             ! xf
    !     ans(8)  = 0.0             ! yf
    !     ans(9)  = zf0             ! zf
    !     ans(10) = phi            ! phi
    !     ans(11) = theta          ! theta
    !     ans(12) = psi            ! psi
        
    !     !  convert phi, theta, and psi to a quat
    !     ans(10:13) = euler_to_quat(ans(10:12))
        
    !     controls(1) = G(3)
    !     controls(2) = G(4)
    !     controls(3) = G(5)
    !     controls(4) = G(6)
        
    !     if (verbose) then
    !         write(*,*) ""
    !         write(*,*) "Solved in ", iter-1, "iterations"
    !         write(*,'(A, 1(1x,ES20.12))') ""
    !         write(*,'(A, 1(1x,ES20.12))') "       Trim settings "
    !         write(*,'(A, 1(1x,ES20.12))') "      alpha[deg] = ", alpha*180./pi
    !         write(*,'(A, 1(1x,ES20.12))') "       beta[deg] = ", beta*180./pi
    !         write(*,'(A, 1(1x,ES20.12))') "        p[deg/s] = ", p_trim*180./pi
    !         write(*,'(A, 1(1x,ES20.12))') "        q[deg/s] = ", q_trim*180./pi
    !         write(*,'(A, 1(1x,ES20.12))') "        r[deg/s] = ", r_trim*180./pi
    !         write(*,'(A, 1(1x,ES20.12))') "         da[deg] = ", controls(1)*180./pi
    !         write(*,'(A, 1(1x,ES20.12))') "         de[deg] = ", controls(2)*180./pi
    !         write(*,'(A, 1(1x,ES20.12))') "         dr[deg] = ", controls(3)*180./pi
    !         write(*,'(A, 1(1x,ES20.12))') "        throttle = ", controls(4)
    !         write(*,*)"" 
    !         write(*,'(A, 1(1x,ES20.12))') "        phi[deg] = ", phi*180./pi
    !         if (has_gamma) then
    !             write(*,'(A, 1(1x,ES20.12))') "climb_angle[deg] = ", gamma0*180./pi
    !         end if              
    !         write(*,'(A, 1(1x,ES20.12))') "      theta[deg] = ", theta*180./pi
    !         write(*,'(A, 1(1x,ES20.12))') "      theta[deg] = ", theta*180./pi
    !         write(*,'(A, 1(1x,ES20.12))') "        psi[deg] = ", psi*180./pi
    !     end if

    ! end function trim_solver


    function pseudo_aero(this, y) result(FM)
        implicit none
        type(vehicle_type) :: this
        real, dimension(13) :: y
        real, dimension(6) :: FM
        real :: u, v, w, p, q, r, xf, yf, zf, phi, theta, psi, alpha, alpha_hat, beta, beta_f, V_mag
        real :: CL1, CL, CS, CD, C_l, Cm, Cn, Reynolds
        real :: CLnewt, CDnewt, Cmnewt, sigma, pos, neg
        real :: Z_pot,T,Press,rho,a, mu, ca, sa, cb, sb, sign_a, pbar, qbar, rbar
        real :: da, de, dr, tau
        real :: FTx, FTy, FTz, MTx, MTy, MTz

        if (verbose2 .and. .not. this%trimming) then
            write(*, '(A, 13(1x,ES20.12))')" state coming in = ", y
        end if
        
        ! make dummy variables for readibility
        u  = y(1)
        v  = y(2)
        w  = y(3)
        p  = y(4) 
        q  = y(5)
        r  = y(6)
        xf = y(7) 
        yf = y(8)
        zf = y(9)

        da = this%controls(1)
        de = this%controls(2)
        dr = this%controls(3)
        tau = this%controls(4)

        V_mag = sqrt(u*u + v*v + w*w)
        ! calculate alpha, V
        alpha  = atan2(w, u)
        sign_a = sign(1.0, alpha)
        alpha_hat = 0.
        beta = asin(v/V_mag)
        ca = cos(alpha)
        sa = sin(alpha)
        cb = cos(beta)
        sb = sin(beta)
        
        ! calc bars
        pbar = p*this%ref_lat/(2*V_mag)
        qbar = q*this%ref_long/(2*V_mag)
        rbar = r*this%ref_lat/(2*V_mag)
        
        ! get rho, mu, Reynolds
        call std_atm_English(-zf,Z_pot,T,Press,rho,a,mu)

        if (this%type == "sphere") then
            CL = 0.0
            CS = 0.0
            C_l = 0.0
            Cm = 0.0
            Cn = 0.0
            Reynolds = 2.*rho*V_mag*this%ref_long/mu

            if (Reynolds<0.01) then
                CD = 2405.0
            else if (Reynolds<450000.0) then 
                CD = 24.0/Reynolds + 6.0/(1.0 + sqrt(Reynolds)) + 0.4
            else if (Reynolds < 560000.0) then
                CD = 10.0**29. * Reynolds**(-5.211)
            else if (Reynolds < 1400000.0) then
                CD = -2.0*10.0**(-23) * Reynolds**3 - 1.0*10.0**(-16) * Reynolds**2 + 9.0*10.0**(-9)*Reynolds + 0.069
            else
                CD = 0.12
            end if
            ! write(*,*) "CD = ", CD
        end if
        
        
        if (this%type == "arrow") then
            beta_f = atan2(v,u)
            CL =  this%CLa*alpha
            CS = -this%CLa*beta_f
            CD =  this%CD0 + this%CD2*CL**2 + this%CD2*CS**2
            C_l = this%Cll0 + this%Clpbar*pbar
            Cm =  this%Cma*alpha + this%Cmqbar*qbar
            Cn = -this%Cma*beta_f + this%Cmqbar*rbar
        end if

        if (this%type == "aircraft") then           
            ! calc forces
            CL1 = this%CL0 + this%CLa*alpha
            CL = CL1 + this%CLqbar*qbar + this%CLahat*alpha_hat + this%CLde*de 
            CS = this%CSb*beta + (this%CSpbar + this%CSapbar*alpha)*pbar + this%CSrbar*rbar + this%CSda*da + this%CSdr*dr
            ! write(*,*) " alpha = ", alpha*180.0/PI
            ! write(*,*) " beta deg= ", beta*180.0/PI
            ! write(*,*) " pbar deg /s = ", pbar*180.0/PI
            ! write(*,*) " rbar deg /s = ", rbar*180.0/PI
            ! write(*,*) " da deg = ", da*180.0/PI
            ! write(*,*) " dr deg = ", dr*180.0/PI
            ! write(*,*) " CSb = ", this%CSb
            ! write(*,*) " CSpbar = ", this%CSpbar
            ! write(*,*) " CSapbar = ", this%CSapbar
            ! write(*,*) " CSrbar = ", this%CSrbar
            ! write(*,*) " CSda = ", this%CSda
            ! write(*,*) " CSdr = ", this%CSdr
            CD = this%CD0 + this%CD1*CL1 + this%CD2*CL1*CL1 + this%CDS2*CS**2 + (this%CDqbar + this%CDaqbar*alpha)*qbar &
            + (this%CDde + this%CDade*alpha)*de + this%CDde2*de*de
            C_l = this%Clb*beta + this%Clpbar*pbar + (this%Clrbar + this%Clarbar*alpha)*rbar + this%Clda*da + this%Cldr*dr
            Cm = this%Cm0 + this%Cma*alpha + this%Cmqbar*qbar + this%Cmahat*alpha_hat + this%Cmde*de
            Cn = this%Cnb*beta + (this%Cnpbar + this%Cnapbar*alpha)*pbar + this%Cnrbar*rbar + (this%Cnda + this%Cnada*alpha)*da + this%Cndr*dr
        
            if (this%include_stall) then
                ! stall CL
                CLnewt = 2.0*sign_a*sa*sa*ca
                pos = exp( this%CL_stall%lambda_b*(alpha - this%CL_stall%alpha_0 + this%CL_stall%alpha_s))
                neg = exp(-this%CL_stall%lambda_b*(alpha - this%CL_stall%alpha_0 - this%CL_stall%alpha_s))
                sigma = (1.0 + neg + pos)/((1.0 + neg)*(1.0 + pos))
                CL = (1.0 - sigma)*CL + sigma*CLnewt

                ! stall CD
                CDnewt = 2.0*sin(abs(alpha))**3
                pos = exp( this%CD_stall%lambda_b*(alpha - this%CD_stall%alpha_0 + this%CD_stall%alpha_s))
                neg = exp(-this%CD_stall%lambda_b*(alpha - this%CD_stall%alpha_0 - this%CD_stall%alpha_s))
                sigma = (1.0 + neg + pos)/((1.0 + neg)*(1.0 + pos))
                CD = (1.0 - sigma)*CD + sigma*CDnewt

                ! stall Cm
                Cmnewt = this%Cm_stall%min_val*sign_a*sa*sa
                pos = exp( this%Cm_stall%lambda_b*(alpha - this%Cm_stall%alpha_0 + this%Cm_stall%alpha_s))
                neg = exp(-this%Cm_stall%lambda_b*(alpha - this%Cm_stall%alpha_0 - this%Cm_stall%alpha_s))
                sigma = (1.0 + neg + pos)/((1.0 + neg)*(1.0 + pos))
                Cm = (1.0 - sigma)*Cm + sigma*Cmnewt
            end if
        end if
            
        ! limit throttle to 0 to 100
        if (tau < 0.0) then
            tau = 0.0
        else if (tau > 1.0) then
            tau = 1.0
        end if
        
        ! Calculate forces and moments due to throttle
        FTx = tau*this%thrust_T0*(rho/this%rho0)**this%thrust_Ta
        FTy = 0.
        FTz = 0.
        MTx = 0.
        MTy = 0.
        MTz = 0.
        ! trimming = .false. ! for debugging
        if (verbose2 .and. .not. this%trimming) then
            write(*,*) "Intermediate values in pseudo_aero"
            write(*,*) "u = ", u
            write(*,*) "v = ", v
            write(*,*) "w = ", w
            write(*,*) "p = ", p
            write(*,*) "q = ", q
            write(*,*) "r = ", r
            write(*,*) "xf = ", xf
            write(*,*) "yf = ", yf
            write(*,*) "zf = ", zf
            write(*,*) "V_mag = ", V_mag
            write(*,*) "pbar = ", pbar
            write(*,*) "qbar = ", qbar
            write(*,*) "rbar = ", rbar
            write(*,*) "sw = ", this%ref_lat
            write(*,*) "cw = ", this%ref_long
            write(*,*) "cos(alpha) = ", ca
            write(*,*) "sin(alpha) = ", sa
            write(*,*) "cos(beta) = ", cb
            write(*,*) "sin(beta) = ", sb
            write(*,*) "ref_area = ", this%ref_area
            write(*,*) "alpha = ", alpha
            write(*,*) "beta = ", beta
            write(*,*) "rho = ", rho
            write(*,*) "CL1 = ", CL1
            write(*,*) "CL = ", CL
            write(*,*) "CS = ", CS
            write(*,*) "CD = ", CD
            write(*,*) "C_l = ", C_l
            write(*,*) "Cm = ", Cm
            write(*,*) "Cn = ", Cn
            write(*,*) "rho0 = ", this%rho0
            write(*,*) "da [deg] = ", this%controls(1)*180./pi
            write(*,*) "de [deg] = ", this%controls(2)*180./pi
            write(*,*) "dr [deg] = ", this%controls(3)*180./pi
            write(*,*) "throttle = ", this%controls(4)

            write(*,*) ""
            write(*,*) "Forces and moments due to thrust"
            write(*,*) "FTx = ", FTx
            write(*,*) "FTy = ", FTy
            write(*,*) "FTz = ", FTz
            write(*,*) "MTx = ", MTx
            write(*,*) "MTy = ", MTy
            write(*,*) "MTz = ", MTz
        end if

        ! update forces and moments
        FM(1) =  FTx + -0.5*rho*(V_mag**2)*this%ref_area*(CD*ca*cb + CS*ca*sb - CL*sa) ! + tau*pow(rho/m_rho0, m_a)*(m_T0 + m_T1*V + m_T2*V*V) !F_xb
        FM(2) =  FTy +  0.5*rho*(V_mag**2)*this%ref_area*(CS*cb - CD*sb)           ! F_yb
        FM(3) =  FTz + -0.5*rho*(V_mag**2)*this%ref_area*(CD*sa*cb + CS*sa*sb + CL*ca) ! F_zb
        FM(4) =  MTx +  0.5*rho*(V_mag**2)*this%ref_area*this%ref_lat*C_l  ! M_xb
        FM(5) =  MTy +  0.5*rho*(V_mag**2)*this%ref_area*this%ref_long*Cm  ! M_yb
        FM(6) =  MTz +  0.5*rho*(V_mag**2)*this%ref_area*this%ref_lat*Cn  ! M_zb
        ! write(*,*) FM(5), rho, V_mag, ref_area, ref_long, Cm

        if (verbose2 .and. .not. this%trimming) then
            write(*,*) ""
            write(*,*) "Forces"
            write(*,*) "Fx = ", FM(1)
            write(*,*) "Fy = ", FM(2)
            write(*,*) "Fz = ", FM(3)
            write(*,*) ""
            write(*,*) "Moments before aero reference shift"
            write(*,*) "Mx = ", FM(4)
            write(*,*) "My = ", FM(5)
            write(*,*) "Mz = ", FM(6)
        end if

        ! shift aero ref location
        ! write(*,*) " forces", FM(1:3)
        ! write(*,*) " moments before cross", FM(4:6)
        FM(4:6) = FM(4:6) + cross3(this%aero_ref_loc, FM(1:3))

        if (verbose2 .and. .not. this%trimming) then
            write(*,*) ""
            write(*,*) "Moments after areo reference shift"
            write(*,*) "Mx = ", FM(4)
            write(*,*) "My = ", FM(5)
            write(*,*) "Mz = ", FM(6)
            write(*, '(A, 6(1x,ES20.12))')"      pseuo aero = ", FM
        end if
     
    end function pseudo_aero


    function diff_eq(this,y) result(dydt)
        
        implicit none

        type(vehicle_type) :: this 
        real, dimension(:), intent(in) :: y
        real, dimension(size(y)) :: dydt
        real, dimension(6) :: FM
        real :: Fxb, Fyb, Fzb, Mxb, Myb, Mzb
        real :: u, v, w, p, q, r, xf, yf, zf, e0, ex, ey, ez
        real :: phi, theta, psi
        real :: g, weight, ac
        real :: Ixx, Iyy, Izz, Ixz, Ixy, Iyz
        real, dimension(4) :: quat_A, quat_B, quat_AB, quat_xyz
        real, dimension(3) :: hpqr, pqr_dot_stuff, pqr_dot

        
        
        Ixx =  this%I(1,1)
        Iyy =  this%I(2,2)
        Izz =  this%I(3,3)
        Ixy = -this%I(1,2)
        Ixz = -this%I(1,3)
        Iyz = -this%I(2,3)

        ! if (verbose2 .and. .not. this%trimming) then
        !     write(*,*) "Ixy = ", Ixy
        !     write(*,*) "Ixz = ", Ixz
        !     write(*,*) "Iyz = ", Iyz
        ! end if
        
        FM = pseudo_aero(this, y)

        if (rk4_verbose .and. .not. this%trimming) then
            write(*, '(A, 13(ES20.12))')"     state vector coming in = ", y
            write(*, '(A,  6(ES20.12))')"          pseudo aero (F,M) = ", FM
        end if

        Fxb = FM(1)
        Fyb = FM(2)
        Fzb = FM(3)
        Mxb = FM(4)
        Myb = FM(5)
        Mzb = FM(6)

        ! declare variables to keep track of stuff
        u  = y(1)
        v  = y(2)
        w  = y(3)
        p  = y(4) 
        q  = y(5)
        r  = y(6)
        xf = y(7) 
        yf = y(8)
        zf = y(9)
        e0 = y(10)
        ex = y(11)
        ey = y(12)
        ez = y(13)
        
        g = gravity_English(-zf)
        ! write(*,*) "g = ", g

        ! weight = mass*g

        hpqr = matmul(this%h, y(4:6))
        ! write(*,*) " matmul [h][pqr] = ", hpqr
        ! write(*,*) " Moments ", Mxb, Myb, Mzb
        pqr_dot_stuff(1) = hpqr(1) + Mxb + (Iyy - Izz)*q*r + Iyz*(q*q - r*r) + Ixz*p*q - Ixy*p*r
        pqr_dot_stuff(2) = hpqr(2) + Myb + (Izz - Ixx)*p*r + Ixz*(r*r - p*p) + Ixy*q*r - Iyz*p*q
        pqr_dot_stuff(3) = hpqr(3) + Mzb + (Ixx - Iyy)*p*q + Ixy*(p*p - q*q) + Iyz*p*r - Ixz*q*r

        ! write(*,*) "I stuff = ", (Iyy - Izz)*q*r + Iyz*(q*q - r*r) + Ixz*p*q - Ixy*p*r
        ! write(*,*) "I stuff = ", (Izz - Ixx)*p*r + Ixz*(r*r - p*p) + Ixy*q*r - Iyz*p*q
        ! write(*,*) "I stuff = ", (Ixx - Iyy)*p*q + Ixy*(p*p - q*q) + Iyz*p*r - Ixz*q*r
        ! write(*,*) ""
        ! write(*,*) " I inv = ", I_inv
        ! write(*,*) ""
        pqr_dot = matmul(this%I_inv, pqr_dot_stuff)
        
        ! do this first to get xdot_f, ydot_f, zdot_f for gravity relief
        dydt(7:9) = quat_dependent_to_base(y(1:3), y(10:13)) ! + Vwf

        ac = gravity_relief_factor*(dydt(7)**2 + dydt(8)**2)/(RE/0.3048 - zf)
        ! write(*,*)" ac = ", ac

        dydt(1)  = (Fxb/this%mass) + (g-ac)*(2.0*(ex*ez - ey*e0)) + (r*v) - (q*w) ! udot
        dydt(2)  = (Fyb/this%mass) + (g-ac)*(2.0*(ey*ez + ex*e0)) + (p*w) - (r*u) ! vdot
        dydt(3)  = (Fzb/this%mass) + (g-ac)*(ez**2 + e0**2 - ex**2 - ey**2) + (q*u) - (p*v) ! wdot
        dydt(4)  = pqr_dot(1) ! pdot
        dydt(5)  = pqr_dot(2) ! qdot
        dydt(6)  = pqr_dot(3) ! rdot
        
        
        dydt(10) = 0.5*(-ex*p - ey*q - ez*r) ! e0 dot
        dydt(11) = 0.5*( e0*p - ez*q + ey*r) ! ex dot
        dydt(12) = 0.5*( ez*p + e0*q - ex*r) ! ey dot
        dydt(13) = 0.5*(-ey*p + ex*q + e0*r) ! ez dot

        if (rk4_verbose .and. .not. this%trimming) then
            write(*, '(A, 13(ES20.12))')"                    diff eq = ", dydt
            write(*,*)""
        end if
    end function diff_eq


    function rk4(this, t0, y, dt) result(y_next)
    
        implicit none
        
        type(vehicle_type) :: this
        real, intent(in) :: t0, dt
        real, dimension(:), intent(in) :: y
        real, dimension(size(y)) :: k1, k2, k3, k4, y_next
    
        ! k1:
        if (rk4_verbose) then
            write(*,*) ""
            write(*, *)"                   RK call =   1"
            write(*, '(A, (ES20.12))')"                       time = ", t0
        end if
        k1 = diff_eq(this, y)


        ! k2:  
        if (rk4_verbose) then
            write(*, *)"                   RK call =   2"
            write(*, '(A, (ES20.12))')"                       time = ", t0 + 0.5*dt
        end if
        k2 = diff_eq(this, y + 0.5*k1*dt)
        
        ! k3:
        if (rk4_verbose) then
            write(*, *)"                   RK call =   3"
            write(*, '(A, (ES20.12))')"                       time = ", t0 + 0.5*dt
        end if
        k3 = diff_eq(this, y + 0.5*k2*dt)
        
        ! k4:
        if (rk4_verbose) then
            write(*, *)"                   RK call =   4"
            write(*, '(A, (ES20.12))')"                       time = ", t0 + dt
        end if
        k4 = diff_eq(this, y + k3*dt)

        y_next = y + dt*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0  
   
   
    end function rk4


    subroutine vehicle_tick_state(this, time, dt)
        
        implicit none
        
        type(vehicle_type) :: this
        real, intent(in) :: time, dt
        real, dimension(13) :: y_next
            
        y_next = rk4(this, time, this%state, dt)

        ! write(*,'(A20, *(ES20.12))') this%name, time, dt, y_next(:)
        call quat_norm(y_next(10:13))
        ! write(*,'(A20, *(ES20.12))') this%name, time, dt, y_next

        this%state = y_next

        ! write to file
        if (save_states) then
            call vehicle_write_state(this, dt, y_next)
        end if
        
        if (verbose2) then
            write(*,*)
            write(*, '(A, (1x,ES20.12))') "           time = ", time
            write(*, '(A, 13(1x,ES20.12))')"          state  = ", y_next
            write(*,*) ""
        end if
            
    end subroutine vehicle_tick_state


    subroutine vehicle_write_state(this, time, y)
        implicit none
        type(vehicle_type) :: this
        real, intent(in) :: time
        real, dimension(13), intent(in) :: y
    end subroutine vehicle_write_state

end module vehicle_m