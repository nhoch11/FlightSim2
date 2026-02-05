module vehicle_m
    use hoch_m
    use jsonx_m
    use micro_time_m
    use linalg_mod
    use connection_m
    ! use json_xtnsn_mod
    
    implicit none

    logical :: save_states, verbose, verbose2, rk4_verbose, gravity_relief
    real :: gravity_relief_factor
    integer :: geographic_model_ID
    character(len=:), allocatable :: geographic_model

    type stall_settings_type
        real :: alpha_0, alpha_s, lambda_b, min_val 
    end type stall_settings_type

    type trim_settings_type
        ! old trim settings
        logical :: solve_for_sideslip, solve_for_fixed_climb_angle, solve_for_relative_climb_angle, solve_for_load_factor, verbose
        logical, allocatable, dimension(:) :: free_vars
        character(:), allocatable :: type
        real :: sideslip, climb_angle, load_factor
        real :: u, v, w, p, q, r
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
        real :: latitude, longitude, latitude_deg, longitude_deg, azimuth_deg

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
                
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CS.beta", this%CSb)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CS.pbar", this%CSpbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CS.alpha_pbar", this%CSapbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CS.rbar", this%CSrbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CS.aileron", this%CSda)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CS.rudder", this%CSdr)
                
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.L0", this%CD0)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.CL1", this%CD1)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.CL1_CL1", this%CD2)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.CS_CS", this%CDS2)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.qbar", this%CDqbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.alpha_qbar", this%CDaqbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.elevator", this%CDde)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.alpha_elevator", this%CDade)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.elevator_elevator", this%CDde2)
                
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.beta", this%Clb)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.pbar", this%Clpbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.rbar", this%Clrbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.alpha_rbar", this%Clarbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.aileron", this%Clda)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.rudder", this%Cldr)
                
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cm.0", this%Cm0)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cm.alpha", this%Cma)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cm.qbar", this%Cmqbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cm.alphahat", this%Cmahat)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cm.elevator", this%Cmde)
                
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.beta", this%Cnb)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.pbar", this%Cnpbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.alpha_pbar", this%Cnapbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.rbar", this%Cnrbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.aileron", this%Cnda)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.alpha_aileron", this%Cnada)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.rudder", this%Cndr)

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
            call jsonx_get(this%j_vehicle, "thrust.location[ft]", this%thrust_loc, 0.0, 3)
            call jsonx_get(this%j_vehicle, "thrust.orientation[ft]", thrust_orientation, 0.0, 3)
            thrust_orientation = thrust_orientation*PI/180.0
            this%thrust_quat = euler_to_quat(thrust_orientation)


            ! initial conditions
            this%init_state = 0.
            this%controls = 0.

            call jsonx_get(this%j_vehicle, "initial.airspeed[ft/s]",  this%init_V)
            call jsonx_get(this%j_vehicle, "initial.altitude[ft]",  this%init_alt)
            this%init_state(9) = -this%init_alt
            call jsonx_get(this%j_vehicle, "initial.latitude[deg]",  this%latitude_deg)
            call jsonx_get(this%j_vehicle, "initial.longitude[deg]",  this%longitude_deg)
            call jsonx_get(this%j_vehicle, "initial.Euler_angles[deg]",  this%init_eul, 0.0, 3)
            this%azimuth_deg = this%init_eul(3)
            this%init_eul = this%init_eul*PI/180.0
            this%latitude = this%latitude_deg*PI/180.
            this%longitude = this%longitude_deg*PI/180.

            
            ! get rho, mu, Reynolds
            call std_atm_English(0.0,junk1,junk2,junk3,this%rho0,junk4,junk5)
            
            call jsonx_get(this%j_vehicle, "initial.type", this%init_type)
            
            if (this%init_type == "state") then
                call init_to_state(this)
            else if (this%init_type == "trim") then
                call init_to_trim(this)
            else
                write(*,*) " Invalid initial type in json. must be trim or state. Quitting..."
                stop
            end if
            
            this%state = this%init_state
            
            if (save_states) then
                call vehicle_write_state(this, 0.0, this%state)
            end if
            
            
        end if ! run_physics
        
        ! connections
        ! call jsonx_get(j_main, "connections", j_connections)
        ! call jsonx_get(j_connections, "graphics", j_graphics)
        ! call graphics%init(j_graphics)
        
    end subroutine vehicle_init
    
    
    subroutine init_to_state(this)
        
        implicit none
        type(vehicle_type) :: this 
        real :: alpha0, beta0
        
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

        this%init_state(10:13) = euler_to_quat(this%init_eul)

    end subroutine init_to_state


    subroutine init_to_trim(this)
        
        implicit none
        type(vehicle_type) :: this 
        integer, parameter :: n_vars = 9
        integer :: max_iter, iter, n_free, i, j, k
        integer, dimension(:), allocatable :: idx_free
        logical :: free_vars(n_vars), found
        real :: fd_step, relaxation, trim_tol, max_R
        real :: alpha, beta, sa, sb, ca, cb
        real, allocatable :: R(:), R_neg(:), R_up(:), R_dn(:), jac(:,:), x(:), dx(:)
        
        this%trimming = .true.
        
        write(*,*) ""
        write(*,*) "---------- Initializing with a Trim State ------------ "
        
        ! read in trim settings
        call jsonx_get(this%j_vehicle, "initial.trim.type", this%trim%type)
        write(*,*) "     Trim type = ", this%trim%type
        
        call jsonx_get(this%j_vehicle, "initial.trim.solver.finite_difference_step_size", fd_step)
        call jsonx_get(this%j_vehicle, "initial.trim.solver.relaxation_factor", relaxation)
        call jsonx_get(this%j_vehicle, "initial.trim.solver.tolerance", trim_tol)
        call jsonx_get(this%j_vehicle, "initial.trim.solver.max_iterations", max_iter)
        call jsonx_get(this%j_vehicle, "initial.trim.solver.verbose", this%trim%verbose)
        
        
        allocate(x(n_vars))
        allocate(dx(n_vars))
        
        ! initialize stuff
        iter = 0
        x = 0.
        x(7:9) = this%init_eul
        free_vars(1:6) = .true.
        free_vars(7:9) = .false.
        this%trim%load_factor = 1.0
        
        
        this%trim%solve_for_sideslip = .false.
        this%trim%solve_for_load_factor = .false.
        this%trim%solve_for_fixed_climb_angle = .false.
        this%trim%solve_for_relative_climb_angle = .false.

        
        
        call json_get(this%j_vehicle, "initial.trim.fixed_climb_angle[deg]", this%trim%climb_angle, this%trim%solve_for_fixed_climb_angle)
        if (this%trim%solve_for_fixed_climb_angle) then
            write(*,*) ""
            write(*,*) "     Fixed Climb Angle Specified "
            write(*,*) "     climb_angle[deg] = ", this%trim%climb_angle
            write(*,*) "     Setting elevation angle to fixed climb angle as initial guess"
            this%trim%climb_angle = this%trim%climb_angle*PI/180
            x(8) = this%trim%climb_angle ! initial guess
            free_vars(8) = .true.
        end if
        
        call json_get(this%j_vehicle, "initial.trim.relative_climb_angle[deg]", this%trim%climb_angle, this%trim%solve_for_relative_climb_angle)
        if (this%trim%solve_for_relative_climb_angle) then
            write(*,*) ""
            write(*,*) "     Relative Climb Angle Specified "
            write(*,*) "     climb_angle[deg] = ", this%trim%climb_angle
            write(*,*) "     Setting elevation angle to relative climb angle as initial guess"
            this%trim%climb_angle = this%trim%climb_angle*PI/180
            x(8) = this%trim%climb_angle ! initial guess
            free_vars(8) = .true.
        end if

        if (this%trim%type == "sct") then
            call json_get(this%j_vehicle, "initial.trim.load_factor", this%trim%load_factor, this%trim%solve_for_load_factor)
            if (this%trim%solve_for_load_factor) then
                write(*,*) ""
                write(*,*) "     Load Factor Specified "
                write(*,*) "     load factor = ", this%trim%load_factor
                x(7) = acos(cos(x(8))/this%trim%load_factor)
                write(*,*) "     Setting load factor to initial guess based on approximation nL = cos(theta)/cos(phi)"
                free_vars(7) = .true.
            end if
        end if

        if (this%trim%type == "shss") then
            call json_get(this%j_vehicle, "initial.trim.sideslip_angle[deg]", this%trim%sideslip, this%trim%solve_for_sideslip)
            if (this%trim%solve_for_sideslip) then
                write(*,*) ""
                write(*,*) "     Sideslip Angle Specified "
                write(*,*) "     beta[deg] = ", this%trim%sideslip
                this%trim%sideslip = this%trim%sideslip*PI/180.
                x(2) = this%trim%sideslip
                free_vars(2) = .false.
                free_vars(7) = .true.
            end if 
        end if


        n_free = count(free_vars)
        write(*,*) "     n_free = ", n_free
        write(*,*) "     idx_free = ", idx_free
        
        allocate(R(n_free))
        allocate(R_up(n_free))
        allocate(R_dn(n_free))
        allocate(Jac(n_free, n_free))
        
        idx_free = pack([(i, i=1,n_vars)], free_vars)
        
        if (this%trim%verbose) then
            write(*,*) " "
            write(*,*) "     Initial guess "
            write(*,'(A, *(1x,ES24.16))') "      x = ", x
        end if 
        
        R = calc_R(this, x, n_free)
        max_R = maxval(abs(R))

        if (this%trim%verbose) then
            write(*,'(A, *(1x,ES24.16))') "      R = ", R
            write(*,'(A, (1x,ES24.16))') "  max_R = ", max_R
            write(*,*) ""
            write(*,*) ""
            write(*,*) "---------- Solving for Trim State ------------ "
            write(*,*) ""
        end if
        
        
        do while (max_R > trim_tol)
            if (iter + 1 > max_iter) exit
            iter = iter + 1
            
            if (this%trim%verbose) then
                write(*,*) ""
                write(*,*) "---------------------------------------------"
                write(*,*) "     Iteration ", iter
                write(*,*) "---------------------------------------------"
            end if
            
            do i=1,n_free
                k = idx_free(i)
                ! step up
                x(k) = x(k) + fd_step
                if (this%trim%verbose) then
                    write(*,*) ""
                    write(*,'(A, I0, A)') "     Perturbing up,   x[",k,"]"
                    write(*,'(A, *(1x,ES24.16))') "     x_up = ", x
                end if
                R_up = calc_R(this, x, n_free)
                    
                if (this%trim%verbose) then
                    write(*,'(A, *(1x,ES24.16))') "     R_up = ", R_up
                end if


                ! step dn
                x(k) = x(k) - 2*fd_step

                if (this%trim%verbose) then
                    write(*,*) ""
                    write(*,'(A, I0, A)') "     Perturbing dn,   x[",k,"]"
                    write(*,'(A, *(1x,ES24.16))') "     x_dn = ", x
                end if
                
                R_dn = calc_R(this, x, n_free)
                if (this%trim%verbose) then
                    write(*,'(A, *(1x,ES24.16))') "     R_dn = ", R_dn
                    write(*,*) ""
                    write(*,*) "     ---------------------------"
                end if
                
                ! do i=1,n_free
                !     jac(i,j) = (R_up(i) - R_dn(i))/(2*fd_step) 
                ! end do
                
                jac(:,i) = (R_up - R_dn)/(2.0*fd_step)

                x(k) = x(k) + fd_step
                
            end do 

            if (this%trim%verbose) then
                write(*,*) ""
                write(*,*) "     Jacobian = "
                do j = 1, n_free
                    write(*,'(*(ES24.16))') jac(j, :)
                end do
            end if
            
            call lu_solve(n_free, jac, -R, dx) 
            
            do i=1,n_free
                k = idx_free(i)
                x(k) = x(k) + relaxation*dx(i)
            end do
            
            if (this%trim%verbose) then
                write(*,*) ""
                write(*,'(A, *(1x,ES24.16))') "      dx = ", dx
                write(*,'(A, *(1x,ES24.16))') "   x new = ", x
            end if
            
            R = calc_R(this, x, n_free)
            max_R = maxval(abs(R))

            if (this%trim%verbose) then
                write(*,'(A, *(1x,ES24.16))') "   R new = ", R
                write(*,'(A, *(1x,ES24.16))') "     eps = ", max_R
                write(*,*) ""
                write(*,*) "     Iteration ", iter, " complete"
                write(*,*) "---------------------------------------------"
            end if
            
        end do
        
        this%controls = x(3:6)
        
        write(*,*) ""
        write(*,*) "     Solved in ", iter, "iterations"
        write(*,*) ""
        write(*,*) "---------------------------------------------"
        write(*,*) "           Trim Solution for ", this%trim%type
        write(*,*) "---------------------------------------------"
        write(*,'(A, 1(1x,ES22.14))') ""
        write(*,'(A, 1(1x,ES22.14))') "                      Trim settings "
        write(*,'(A, 1(1x,ES22.14))') "                 alpha[deg] = ", x(1)*180./pi
        write(*,'(A, 1(1x,ES22.14))') "                  beta[deg] = ", x(2)*180./pi
        write(*,*)"" 
        write(*,'(A, 1(1x,ES22.14))') "                   p[deg/s] = ", this%init_state(4)*180./pi
        write(*,'(A, 1(1x,ES22.14))') "                   q[deg/s] = ", this%init_state(5)*180./pi
        write(*,'(A, 1(1x,ES22.14))') "                   r[deg/s] = ", this%init_state(6)*180./pi
        write(*,*)"" 
        write(*,'(A, 1(1x,ES22.14))') "                    da[deg] = ", this%controls(1)*180./pi
        write(*,'(A, 1(1x,ES22.14))') "                    de[deg] = ", this%controls(2)*180./pi
        write(*,'(A, 1(1x,ES22.14))') "                    dr[deg] = ", this%controls(3)*180./pi
        write(*,'(A, 1(1x,ES22.14))') "                   throttle = ", this%controls(4)
        write(*,*)"" 
        write(*,'(A, 1(1x,ES22.14))') "                   phi[deg] = ", x(7)*180./pi
        write(*,'(A, 1(1x,ES22.14))') "                 theta[deg] = ", x(8)*180./pi
        write(*,'(A, 1(1x,ES22.14))') "                   psi[deg] = ", x(9)*180./pi
        if (this%trim%solve_for_fixed_climb_angle) then
            write(*,*)"" 
            write(*,'(A, 1(1x,ES22.14))') "     fixed climb angle[deg] = ", this%trim%climb_angle*180./pi
        end if
        if (this%trim%solve_for_relative_climb_angle) then
            write(*,*)"" 
            write(*,'(A, 1(1x,ES22.14))') "  relative climb angle[deg] = ", this%trim%climb_angle*180./pi
        end if
        if (this%trim%solve_for_load_factor) then
            write(*,*)"" 
            write(*,'(A, 1(1x,ES22.14))') "                load factor = ", this%trim%load_factor
        end if 
        
        
        this%trimming = .false.
    
    end subroutine init_to_trim


    function calc_R(this, x, n_free) result(R)
        implicit none
        type(vehicle_type) :: this
        real, intent(in) :: x(9)
        integer, intent(in) :: n_free
        real :: dydt(13), R(n_free), xyzdot(3), FM(6)
        real :: ca, sa, cb, sb, cp, sp, ct, st 
        real :: u, v, w, V0, grav, ac
        real :: c1, pw
        integer :: last

        ca = cos(x(1))
        sa = sin(x(1))
        cb = cos(x(2))
        sb = sin(x(2))
        cp = cos(x(7))
        sp = sin(x(7))
        ct = cos(x(8))
        st = sin(x(8))

        V0 = this%init_V
        u = V0*ca*cb
        v = V0*sb
        w = V0*sa*cb
        
        this%controls = x(3:6)

        ! build state for diff eq
        this%init_state(1) = u
        this%init_state(2) = v
        this%init_state(3) = w
        
        grav = gravity_English(-this%init_state(9))
        ! do this first to calc ac
        this%init_state(10:13) = euler_to_quat(x(7:9))
        xyzdot = quat_dependent_to_base(this%init_state(1:3), this%init_state(10:13)) ! + Vwf
        ac = gravity_relief_factor*(xyzdot(1)**2 + xyzdot(2)**2)/(RE/0.3048 - this%init_state(9))
        c1 = (grav-ac)*sp*ct/(u*ct*cp + w*st)
            
        if (this%trim%type == "shss") then
            this%init_state(4) = 0.0
            this%init_state(5) = 0.0
            this%init_state(6) = 0.0     
        else if (this%trim%type == "sct") then ! eq 6.2.3
            this%init_state(4) = -c1*st
            this%init_state(5) =  c1*sp*ct
            this%init_state(6) =  c1*cp*ct
        else if (this%trim%type == "br") then
            pw = c1*(-u*st + v*sp*ct + w*cp*ct)/V0
            this%init_state(4) = pw*u/V0
            this%init_state(5) = pw*v/V0
            this%init_state(6) = pw*w/V0
        end if 

        if (this%trim%verbose) then
            write(*,*) "     Updating p,q,r for ", this%trim%type
            write(*,'(A, 1(1x,ES24.16))') "        p[deg/s] = ", this%init_state(4)*180./pi
            write(*,'(A, 1(1x,ES24.16))') "        q[deg/s] = ", this%init_state(5)*180./pi
            write(*,'(A, 1(1x,ES24.16))') "        r[deg/s] = ", this%init_state(6)*180./pi
        end if


        dydt = diff_eq(this, this%init_state)
        R(1:6) = dydt(1:6) 
        last = 6

        if (this%trim%solve_for_load_factor) then
            last = last + 1
            FM = pseudo_aero(this, this%init_state)
            R(last) = this%trim%load_factor - (FM(1)*sa - FM(3)*ca)/(this%mass*(grav-ac))
        end if

        if (this%trim%solve_for_fixed_climb_angle) then
            last = last + 1
            R(last) = this%trim%climb_angle - calc_relative_climb_angle(this%init_state)
        end if

        if (this%trim%solve_for_relative_climb_angle) then
            last = last + 1
            R(last) = this%trim%climb_angle - calc_relative_climb_angle(this%init_state)
        end if

    end function calc_R


    ! function calc_load_factor(y) result(ans)
    !     implicit none
    !     real, dimension(13) :: y
    !     real :: ans 
    !     real, dimension(3) :: xyzdot
    !     real :: Vmag 

    !     xyzdot = quat_dependent_to_base(y(1:3), y(10:13))
    !     Vmag = sqrt(y(1)**2 + y(2)**2 + y(3)**2)
    !     ans = asin(-xyzdot(3)/Vmag)
    ! end function calc_load_factor


    function calc_fixed_climb_angle(y) result(ans)
        implicit none
        real, dimension(13) :: y
        real :: ans 
        real, dimension(3) :: xyzdot
        real :: Vmag 

        xyzdot = quat_dependent_to_base(y(1:3), y(10:13))
        ans = asin(-xyzdot(3)/sqrt(xyzdot(1)**2 + xyzdot(2)**2 + xyzdot(3)**2))
    end function calc_fixed_climb_angle


    function calc_relative_climb_angle(y) result(ans)
        implicit none
        real, dimension(13) :: y
        real :: ans 
        real, dimension(3) :: xyzdot
        real :: Vmag 

        xyzdot = quat_dependent_to_base(y(1:3), y(10:13))
        Vmag = sqrt(y(1)**2 + y(2)**2 + y(3)**2)
        ans = asin(-xyzdot(3)/Vmag)
    end function calc_relative_climb_angle



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
        ! if (tau < 0.0) then
        !     tau = 0.0
        ! else if (tau > 1.0) then
        !     tau = 1.0
        ! end if
        
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
        ! write(*,*)" Fxb/mass = ",Fxb/this%mass
        ! write(*,*)" other stuff = ",(Fxb/this%mass) + (g-ac)*(2.0*(ex*ez - ey*e0))
        
        
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
        real, dimension(13) :: y, y_next
        
        y = this%state

        y_next = rk4(this, time, y, dt)

        if (geographic_model_ID > 0) call update_geographic(this, y, y_next)

        call quat_norm(y_next(10:13))
        if (sqrt(y_next(4)**2 + y_next(5)**2 + y_next(6)**2)*dt/(2.0*PI) > 0.1) write(*,*) "WARNING: High Vehicle Rotation relative to integration timestep"
        
        this%state = y_next

        ! write to file
        if (save_states) then
            ! call vehicle_write_state(this, dt, y_next)
            write(this%iunit_states,*) "time[s],u[ft/s],v[ft/s],w[ft/s],p[rad/s],q[rad/s],r[rad/s],x[ft],y[ft],z[ft],e0,ex,ey,ez"
            write(*,*) "    - states will be written to ", this%states_filename
        end if
        
        if (verbose2) then
            write(*,*)
            write(*, '(A, (1x,ES20.12))') "           time = ", time
            write(*, '(A, 13(1x,ES20.12))')"          state  = ", y_next
            write(*,*) ""
        end if
            
    end subroutine vehicle_tick_state


    subroutine update_geographic(this, y1, y2)
        implicit none
        type(vehicle_type) :: this
        real, dimension(13), intent(in) :: y1, y2
        real :: dx, dy, dz, d
        real :: H, Phi1, Theta, Psi1, psi_g1
        real :: cPhi, sPhi, cTheta, sTheta, cpsi_g1, spsi_g1 
        real :: xhat, yhat, zhat, xhp, yhp, zhp
        real :: Chat, Shat, rhat, dg
        real :: quat(4), eul(3)

        dx = (y2(7) - y1(7))
        dy = (y2(8) - y1(8))
        d = sqrt(dx**2 + dy**2)
        
        if (d < TOLERANCE) then
            write(*,*) "WARNING: Geographic coordinates update- xf, yf displacement is too small, d = ", d
        else
            H = -y1(9)
            dz = (y2(9) - y1(9))
            Phi1 = this%latitude
            Psi1 = this%longitude
            cPhi = cos(Phi1)
            sPhi = sin(Phi1)
            if (geographic_model_ID == 1) then ! spherical earth model
                Theta = d/(RE + H - dz/2.0)
                cTheta = cos(Theta)
                sTheta = sin(Theta)
                
                psi_g1 = atan2(dy, dx)
                cpsi_g1 = cos(psi_g1)
                spsi_g1 = sin(psi_g1)

                xhat = cPhi*cTheta - sPhi*sTheta*cpsi_g1 
                yhat = sTheta*spsi_g1
                zhat = sPhi*cTheta + cPhi*sTheta*cpsi_g1

                xhp = -cPhi*sTheta - sPhi*cTheta*cpsi_g1
                yhp = cTheta*spsi_g1
                zhp = -sPhi*sTheta + cPhi*cTheta*cpsi_g1

                rhat = sqrt(xhat**2 + yhat**2)

                this%latitude = atan2(zhat, rhat)
                this%longitude = Psi1 + atan2(yhat, xhat)

                Chat = (xhat**2)*zhp
                Shat = (xhat*yhp - yhat*xhp)*cos(this%latitude)**2*cos(this%longitude - Psi1)**2
                dg = atan2(Shat, Chat) - psi_g1

            else if (geographic_model_ID == 2) then ! ellipse model
                ! implement ellipse model
            end if

            ! limit longitude
            if (this%longitude >  PI) this%longitude = this%longitude - 2.*PI
            if (this%longitude < -PI) this%longitude = this%longitude + 2.*PI

            quat(1) = -this%state(13)
            quat(2) = -this%state(12)
            quat(3) =  this%state(11)
            quat(4) =  this%state(10)

            this%state(10:13) = cos(dg/2.)*this%state(10:13) + sin(dg/2.)*quat(:)
            
            
            eul = quat_to_euler(this%state(10:13))
            this%azimuth_deg = eul(3)*180./PI
            this%latitude_deg = this%latitude*180./PI
            this%longitude_deg = this%longitude*180./PI
        end if

        


    end subroutine update_geographic


    subroutine vehicle_write_state(this, time, y)
        implicit none
        type(vehicle_type) :: this
        real, intent(in) :: time
        real, dimension(13), intent(in) :: y
    end subroutine vehicle_write_state

end module vehicle_m