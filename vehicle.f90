module vehicle_m
    use hoch_m
    use jsonx_m
    use micro_time_m
    use linalg_mod
    use connection_m
    ! use json_xtnsn_mod
    
    implicit none

    type stall_settings_t

    end type stall_settings_t

    type trim_settings_t
        ! old trim settings
        logical :: has_sideslip, has_gamma
        character(:), allocatable :: init_type, trim_type
        real :: trim_fd_step, trim_r_factor, trim_tol, u_trim, v_trim, w_trim, p_trim, q_trim, r_trim, sideslip0, gamma0
        integer :: trim_max_iter
        logical :: trimming
    end type trim_settings_t

    type vehicle_t

        type(json_value), pointer :: j_vehicle

        character(len=:), allocatable :: name
        character(len=:), allocatable :: type
        character(100) :: states_filename, rk4_filename

        logical :: run_physics
        logical :: save_states
        logical :: iunit_states, iunit_rk4, iunit_trim
        logical :: verbose, verbose2

        ! mass constants
        real :: weight, mass
        real, dimension(3,3) :: I, I_inv, h
        real, dimension(6) :: FM
        real :: Ixx0, Iyy0, Izz0, Ixz0, Ixy0, Iyz0, hx0, hy0, hz0
        
        ! aero constants
        real, dimension(:), allocatable :: aero_ref_loc
        real :: ref_area, ref_long, ref_lat
        real :: CL0, CLa, CLahat, CLqbar, CLde
        real :: CD0, CD1, CD2, CDS2, CDqbar, CDaqbar, CDde, CDade, CDde2
        real :: CSb, CSpbar, CSapbar, CSrbar, CSda, CSdr
        real :: Clb, Clpbar, Clrbar, Clarbar, Clda, Cldr
        real :: Cm0, Cma, Cmqbar, Cmahat, Cmde
        real :: Cnb, Cnpbar, Cnapbar, Cnrbar, Cnda, Cnada, Cndr
        real :: Throt_0, Throt_a
        
        ! stall model constants
        logical :: include stall
        type(stall_settings_t) :: CLstall, CDstall, Cmstall

        ! intitialization constants
        real, dimension(:), allocatable :: eul0
        real, dimension(13) :: init_state
        real :: V0, zf0

        ! old init constants
        real :: phi0, theta0, psi0, p0, q0, r0
        real :: alpha0, beta0, da0, de0, dr0, tau0, rho0


        ! variables
        real, dimension(13) :: states
        real, dimension(4) :: controls

        type(trim_settings_t) :: trim

        
        ! old misc settings
        type(connection) :: graphics, control_udp
        real, dimension(3) :: Vwf

    end type vehicle_t
        
        

contains 

    subroutine vehicle_init(input_file)
        implicit none
        character(len=*), intent(in) :: input_file
        type(json_value), pointer :: j_main
        type(json_value), pointer :: j_connections, j_graphics
        logical :: found
        integer :: dot_pos
        character(len=256) :: base_name

        ! Find position of ".json" (or last dot)
        dot_pos = index(input_file, ".json")

        if (dot_pos > 0) then
            base_name = input_file(1:dot_pos-1)
        else
            base_name = input_file   ! fallback if no .json found
        end if

        ! Build the output filename
        output_file = trim(base_name)//"_output.txt"
        
        call jsonx_load(input_file,j_main)
        
        call jsonx_get(j_main, "simulation.time_step[s]", dt, 0.0)
        if (abs(dt) <= TOLERANCE) then
            real_time = .true.
        else 
            real_time = .false.
        end if 
        call jsonx_get(j_main, "simulation.total_time[s]", tf)
        call jsonx_get(j_main, "simulation.verbose", verbose, .false.)
        call jsonx_get(j_main, "simulation.verbose2", verbose2, .false.)
        
        ! call json_xtnsn_get(mass_settings, "mass_slug",  mass)
        call jsonx_get(j_main, "vehicle.mass.weight[lbf]",  weight)
        call jsonx_get(j_main, "vehicle.mass.Ixx[slug-ft^2]",  Ixx0)
        call jsonx_get(j_main, "vehicle.mass.Iyy[slug-ft^2]",  Iyy0)
        call jsonx_get(j_main, "vehicle.mass.Izz[slug-ft^2]",  Izz0)
        call jsonx_get(j_main, "vehicle.mass.Ixy[slug-ft^2]",  Ixy0)
        call jsonx_get(j_main, "vehicle.mass.Ixz[slug-ft^2]",  Ixz0)
        call jsonx_get(j_main, "vehicle.mass.Iyz[slug-ft^2]",  Iyz0)
        call jsonx_get(j_main, "vehicle.mass.hx[slug-ft^2/s]",  hx0)
        call jsonx_get(j_main, "vehicle.mass.hy[slug-ft^2/s]",  hy0)
        call jsonx_get(j_main, "vehicle.mass.hz[slug-ft^2/s]",  hz0)

        call jsonx_get(j_main, "vehicle.thrust.T0[lbf]",  Throt_0)
        call jsonx_get(j_main, "vehicle.thrust.a",  Throt_a)

        call jsonx_get(j_main, "vehicle.aerodynamics.reference.area[ft^2]",  ref_area)
        call jsonx_get(j_main, "vehicle.aerodynamics.reference.longitudinal_length[ft]",  ref_long)
        call jsonx_get(j_main, "vehicle.aerodynamics.reference.lateral_length[ft]",  ref_lat)
        call jsonx_get(j_main, "vehicle.aerodynamics.reference.relative_location[ft]",  aero_ref_loc, 0.0, 3)
        
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CL.0", CL0)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CL.alpha", CLa)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CL.alphahat", CLahat)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CL.qbar", CLqbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CL.elevator", CLde)

        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CS.beta", CSb)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CS.pbar", CSpbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CS.alpha_pbar", CSapbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CS.rbar", CSrbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CS.aileron", CSda)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CS.rudder", CSdr)

        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.L0", CD0)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.CL1", CD1)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.CL1_CL1", CD2)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.CS_CS", CDS2)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.qbar", CDqbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.alpha_qbar", CDaqbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.elevator", CDde)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.alpha_elevator", CDade)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.elevator_elevator", CDde2)
        
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cl.beta", Clb)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cl.pbar", Clpbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cl.rbar", Clrbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cl.alpha_rbar", Clarbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cl.aileron", Clda)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cl.rudder", Cldr)
        
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cm.0", Cm0)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cm.alpha", Cma)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cm.qbar", Cmqbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cm.alphahat", Cmahat)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cm.elevator", Cmde)

        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.beta", Cnb)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.pbar", Cnpbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.alpha_pbar", Cnapbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.rbar", Cnrbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.aileron", Cnda)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.alpha_aileron", Cnada)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.rudder", Cndr)

        call jsonx_get(j_main, "initial.airspeed[ft/s]",  V0)
        call jsonx_get(j_main, "initial.altitude[ft]",  zf0)
        zf0 = -1.*zf0

        call jsonx_get(j_main, "initial.type", init_type)

        if (init_type == "state") then
            call init_to_state(t)
        else if (init_type == "trim") then
            call init_to_trim(t)
        else
            write(*,*) " Invalid initial type in json. must be trim or state. Quitting..."
            stop
        end if
     
        
        Vwf(1) = 0.
        Vwf(2) = 0.
        Vwf(3) = 0.
        !  print results in a file
        ! FILE* en_file = fopen("cannonball.txt", "w");
        ! fprintf(en_file, "  Time[s]              u[ft/s]            v[ft/s]                w[ft/s]              p[rad/s]             q[rad/s]             r[rad/s]             x[ft]                y[ft]                z[ft]                e0                   ex                   ey                   ez\n");
    

        ! connections
        call jsonx_get(j_main, "connections", j_connections)
        call jsonx_get(j_connections, "graphics", j_graphics)
        call graphics%init(j_graphics)
    end subroutine vehicle_init


    subroutine init_to_state(t)
        trimming = .false.
        call jsonx_get(j_main, "initial.state.alpha[deg]", alpha0)
        call jsonx_get(j_main, "initial.state.beta[deg]", beta0)
        alpha0 = alpha0*pi/180.
        beta0 = beta0*pi/180.
        write(*,*) "beta0 rad = ", beta0

        call jsonx_get(j_main, "initial.state.bank_angle[deg]", phi0)
        call jsonx_get(j_main, "initial.state.elevation_angle[deg]", theta0)
        call jsonx_get(j_main, "initial.state.heading_angle[deg]", psi0)
        phi0 = phi0*pi/180.
        theta0 = theta0*pi/180.
        psi0 = psi0*pi/180.
        
        call jsonx_get(j_main, "initial.state.p[deg/s]",  p0)
        call jsonx_get(j_main, "initial.state.q[deg/s]",  q0)
        call jsonx_get(j_main, "initial.state.r[deg/s]",  r0)
        p0 = p0*pi/180.
        q0 = q0*pi/180.
        r0 = r0*pi/180.
        
        call jsonx_get(j_main, "initial.state.aileron[deg]",  da0)
        call jsonx_get(j_main, "initial.state.elevator[deg]",  de0)
        call jsonx_get(j_main, "initial.state.rudder[deg]",  dr0)
        call jsonx_get(j_main, "initial.state.throttle",  tau0)
        da0 = da0*pi/180.
        de0 = de0*pi/180.
        dr0 = dr0*pi/180.
    end subroutine init_to_state


    subroutine init_to_trim(t)
        trimming = .true.
        call jsonx_get(j_main, "initial.trim.type", trim_type)

        has_sideslip = .false.
        call json_get(j_main, "initial.trim.sideslip[deg]", sideslip0, has_sideslip)
        if (has_sideslip) then
            sideslip0 = sideslip0*pi/180.
        end if

        call json_get(j_main, "initial.trim.bank_angle[deg]", phi0, found)
        if (found) then
            phi0 = phi0*pi/180.
        end if
        
        call jsonx_get(j_main, "initial.trim.heading_angle[deg]", psi0)
        psi0 = psi0*pi/180.

        call json_get(j_main, "initial.trim.climb_angle[deg]", gamma0, has_gamma)
        write(*,*) " has gamma = ", has_gamma
        if (has_gamma) then
            gamma0 = gamma0*pi/180.
        else
            call jsonx_get(j_main, "initial.trim.elevation_angle[deg]", theta0)
            theta0 = theta0*pi/180.
        end if

        call jsonx_get(j_main, "initial.trim.solver.finite_difference_step_size", trim_fd_step)
        call jsonx_get(j_main, "initial.trim.solver.relaxation_factor", trim_r_factor)
        call jsonx_get(j_main, "initial.trim.solver.tolerance", trim_tol)
        call jsonx_get(j_main, "initial.trim.solver.max_iterations", trim_max_iter)
    end subroutine init_to_trim


    subroutine mass_inertia()
        implicit none

        real :: det


        ! calc mass
        ! r = ref_length/2.
        mass = weight*.3048/GSSL
        if (verbose2 .and. .not. trimming) then
            write(*,*) "mass = ", mass
        end if

        ! Ixxyyzz = mass*2.*(r**5 - (r-0.00131)**5)/(5.*(r**3 - (r-0.00131)**3))
        h = 0.
        h(1,2) = -hz0
        h(1,3) =  hy0
        h(2,1) =  hz0
        h(2,3) = -hx0
        h(3,1) = -hy0
        h(3,2) =  hx0

        I(1,1) =  Ixx0
        I(2,2) =  Iyy0
        I(3,3) =  Izz0
        
        I(1,2) = -Ixy0
        I(1,3) = -Ixz0
        I(2,1) = -Ixy0
        I(2,3) = -Iyz0
        I(3,1) = -Ixz0
        I(3,2) = -Iyz0


        ! calculate I_inv by hand

        det = I(1,1)*(I(2,2)*I(3,3) - I(2,3)*I(3,2)) &
            - I(1,2)*(I(2,1)*I(3,3) - I(2,3)*I(3,1)) &
            + I(1,3)*(I(2,1)*I(3,2) - I(2,2)*I(3,1))
        write(*,*) "denom = ", det

        if (det < 0) then
            ! The matrix is singular, cannot be inverted
            write(*,*) "Singular inertia matrix.... quitting"
            stop
        end if

        I_inv(1,1) = (I(2,2)*I(3,3) - I(2,3)*I(3,2)) / det;
        I_inv(1,2) = (I(1,3)*I(3,2) - I(1,2)*I(3,3)) / det;
        I_inv(1,3) = (I(1,2)*I(2,3) - I(1,3)*I(2,2)) / det;
        I_inv(2,1) = (I(2,3)*I(3,1) - I(2,1)*I(3,3)) / det;
        I_inv(2,2) = (I(1,1)*I(3,3) - I(1,3)*I(3,1)) / det;
        I_inv(2,3) = (I(1,1)*I(2,3) - I(1,3)*I(2,1)) / det;
        I_inv(3,1) = (I(2,1)*I(3,2) - I(2,2)*I(3,1)) / det;
        I_inv(3,2) = (I(1,2)*I(3,1) - I(1,1)*I(3,2)) / det;
        I_inv(3,3) = (I(1,1)*I(2,2) - I(1,2)*I(2,1)) / det;

        if (verbose2 .and. .not. trimming) then
            write(*,*) "I matrix", I
            write(*,*) "I_inv = ", I_inv
            write(*,*) "h matrix = ", h
        end if
       
        return
    end subroutine mass_inertia


    function calc_R(p,q,r, theta, G) result(R_)

        implicit none 

        real, intent(in) :: p, q, r, theta
        real, dimension(6), intent(in) :: G
        real, dimension(6) :: R_
        real, dimension(13) :: y, dydt 
        real :: alpha, beta, phi, psi

        alpha = G(1)
        ! theta = theta0
        psi = psi0
        controls(1) = G(3)
        controls(2) = G(4)
        controls(3) = G(5)
        controls(4) = G(6)

        if (trim_type == "shss" .and. has_sideslip) then
            phi = G(2)
            beta = sideslip0 

        else
            beta = G(2)
            phi = phi0
        end if
        
        y(1)  = V0*cos(alpha)*cos(beta)   ! u
        y(2)  = V0*sin(beta)               ! v
        y(3)  = V0*sin(alpha)*cos(beta)   ! w
        y(4)  = p             ! p
        y(5)  = q              ! q
        y(6)  = r              ! r
        y(7)  = 0.0             ! xf
        y(8)  = 0.0             ! yf
        y(9)  = zf0             ! zf
        y(10) = phi            ! phi
        y(11) = theta          ! theta
        y(12) = psi            ! psi

        
        !  convert phi, theta, and psi to a quat
        y(10:13) = euler_to_quat(y(10:12))
        
        !  normalize the quat
        call quat_norm(y(10:13))
        
        ! write(*,*) "before diff_eq"
        if (verbose2) then
            write(*,*) "y going into diff_eq in calc_R", y
        end if
        
        ! write(*,*)"y = ",y
        dydt = diff_eq(y)
        ! write(*,*) "after diff_eq"
        ! write(*,*)"dydt = ",dydt
        R_ = dydt(1:6)


    end function calc_R


    function trim_solver() result(ans)

        implicit none

        real :: alpha, beta, u, v, w, p, q, r, phi, theta, psi, da, de, dr, tau, eps, grav, c1, pw
        real :: a, b, c, theta_plus, theta_minus, err_theta_plus, err_theta_minus
        real, dimension(13) :: ans
        real, dimension(6) :: G, R_plus, R_minus, R_neg, R_
        real,dimension(:),allocatable :: delta_G
        real, dimension(6,6) :: jac
        integer :: i,j, iter

        p_trim = 0.0
        q_trim = 0.0
        r_trim = 0.0
        ! theta = theta0
        psi = psi0
        controls = 0.0
        G = 0.0

        alpha = G(1)
            
        if (trim_type == "shss" .and. has_sideslip) then
            phi = G(2)
            beta = sideslip0 
        else
            beta = G(2)
            phi = phi0
        end if

        eps = 1.0

        if (verbose) then
            write(*,*) ""
            if (trim_type == "shss") then

                if (has_sideslip) then
                    write(*,*) "Solving for trim condition: Steady Heading Sideslip"
                    if (has_gamma) then
                        write(*,*) "beta[deg] = ", sideslip0*180./pi, ",  climb_angle[deg] = ", gamma0*180./pi, &
                        ",  psi[deg] = ", psi*180./pi
                    else
                        write(*,*) "beta[deg] = ", sideslip0*180./pi, ",  theta[deg] = ", theta*180./pi, &
                        ",  psi[deg] = ", psi*180./pi
                    end if
                else
                    write(*,*) "Solving for trim condition: Steady Heading Sideslip"
                    if (has_gamma) then
                        write(*,*) "phi[deg] = ", phi*180./pi, ",  climb_angle[deg] = ", gamma0*180./pi, &
                        ",  psi[deg] = ", psi*180./pi
                    else
                        write(*,*) "phi[deg] = ", phi*180./pi, ",  theta[deg] = ", theta*180./pi, &
                        ",  psi[deg] = ", psi*180./pi
                    end if
                end if

            else if (trim_type == "sct") then
                write(*,*) "Solving for trim condition: Steady Coordinated Turn"
                if (has_gamma) then
                    write(*,*) "phi[deg] = ", phi*180./pi, ",  climb_angle[deg] = ", gamma0*180./pi, &
                    ",  psi[deg] = ", psi*180./pi
                else
                    write(*,*) "phi[deg] = ", phi*180./pi, ",  theta[deg] = ", theta*180./pi, &
                    ",  psi[deg] = ", psi*180./pi
                end if
            else 
                write(*,*) "Solving for trim condition: Vertical Barrel Roll"
                if (has_gamma) then
                    write(*,*) "phi[deg] = ", phi*180./pi, ",  climb_angle[deg] = ", gamma0*180./pi, &
                    ",  psi[deg] = ", psi*180./pi
                else
                    write(*,*) "phi[deg] = ", phi*180./pi, ",  theta[deg] = ", theta*180./pi, &
                    ",  psi[deg] = ", psi*180./pi
                end if
            end if
            write(*,*)""
            ! write(*,*) "initial eps = ", eps
        end if

        grav = gravity_English(-zf0)
        iter = 1
        
        do while(eps > trim_tol .and. iter < trim_max_iter)
            if (verbose2) then
                write(*,*)""
                write(*,*) "iter = ", iter
            end if
            ! limit throttle to 0 to 1
            if (G(6) < 0.0) then
                G(6) = 0.0
                if (verbose2) then
                    write(*,*)""
                    write(*,*) " Overwriting negative throttle"
                end if 
            else if (G(6) > 1.0) then
                G(6) = 1.0
                if (verbose2) then
                    write(*,*)""
                    write(*,*) " Limiting throttle to 100 percent"
                end if
            end if

            alpha = G(1)
            
            if (trim_type == "shss" .and. has_sideslip) then
                phi = G(2)
                beta = sideslip0 
            else
                beta = G(2)
                phi = phi0
            end if
            
            u  = V0*cos(alpha)*cos(beta)   ! u
            v  = V0*sin(beta)               ! v
            w  = V0*sin(alpha)*cos(beta)

            if (has_gamma) then
                a = u*V0*sin(gamma0)
                b = (v*sin(phi) + w*cos(phi))*sqrt(u*u + (v*sin(phi) + w*cos(phi))**2 - V0*V0*sin(gamma0)*sin(gamma0))
                c = u*u + (v*sin(phi) + w*cos(phi))**2
                theta_plus = asin((a + b)/c)
                theta_minus = asin((a - b)/c)
                err_theta_plus = abs(u*sin(theta_plus) - (v*sin(phi) + w*cos(phi))*cos(theta_plus) - V0*sin(gamma0))
                err_theta_minus = abs(u*sin(theta_minus) - (v*sin(phi) + w*cos(phi))*cos(theta_minus) - V0*sin(gamma0))
                ! write(*,*) "err plus = ", err_theta_plus
                ! write(*,*) "err minus = ", err_theta_minus
                if (err_theta_plus < err_theta_minus) then
                    theta = theta_plus
                else 
                    theta = theta_minus
                end if
                if (verbose2) then
                    write(*,*)""
                    write(*,'(A, 1(1x,ES20.12))') "      theta_1[deg] = ", theta_plus*180./pi
                    write(*,'(A, 1(1x,ES20.12))') "      theta_2[deg] = ", theta_minus*180./pi
                    write(*,'(A, 1(1x,ES20.12))') "correct theta[deg] = ", theta*180./pi
                    write(*,'(A, 1(1x,ES20.12))') "correct theta[rad] = ", theta
                    write(*,*)""
                end if 
            else
                theta = theta0
            end if
            
            c1 = grav*sin(phi)*cos(theta)/(u*cos(theta)*cos(phi) + w*sin(theta))
            
            if (trim_type == "sct") then ! eq 6.2.3
                ! write(*,*) "phi", phi
                ! write(*,*) "theta", theta
                ! write(*,*) "u", u
                ! write(*,*) "v", v
                ! write(*,*) "w", w
                ! write(*,*) "grav", grav
                ! write(*,*) "alt", -zf0
                p_trim = -c1*sin(theta)
                q_trim =  c1*sin(phi)*cos(theta)
                r_trim =  c1*cos(phi)*cos(theta)
            end if
            
            if (trim_type == "vbr") then
                pw = c1*(-u*sin(theta) + v*sin(phi)*cos(theta) + w*cos(phi)*cos(theta))/V0
                p_trim = pw*u/V0
                q_trim = pw*v/V0
                r_trim = pw*w/V0
            end if 


            if (verbose2) then
                write(*,*) ""
                write(*,*) "Updating p,q,r"
                write(*,'(A, 1(1x,ES20.12))') "p[deg/s] = ", p_trim*180./pi
                write(*,'(A, 1(1x,ES20.12))') "q[deg/s] = ", q_trim*180./pi
                write(*,'(A, 1(1x,ES20.12))') "r[deg/s] = ", r_trim*180./pi
            end if
            
            if (verbose2) then
                write(*,*) ""
                write(*,*) "G with Throttle forced to be 0-1"
                write(*,'(A, 6(1x,ES20.12))') "    G = ", G
            end if
            R_ = calc_R(p_trim, q_trim, r_trim, theta, G)
            if (verbose2) then
                write(*,'(A, 6(1x,ES20.12))') "   R_ = ", R_
            end if

            do j=1,6
                G(j) = G(j) + trim_fd_step
                if (verbose2) then
                    write(*,*) ""
                    write(*,*) "                 j = ", j
                    write(*,'(A, 6(1x,ES20.12))') "   G_up = ", G
                end if
                R_plus = calc_R(p_trim, q_trim, r_trim, theta, G)

                if (verbose2) then
                    write(*,'(A, 6(1x,ES20.12))') " R_plus = ", R_plus
                end if

                G(j) = G(j) - 2*trim_fd_step
                if (verbose2) then
                    write(*,*) ""
                    write(*,'(A, 6(1x,ES20.12))') "   G_dn = ", G
                end if
                R_minus = calc_R(p_trim, q_trim, r_trim, theta, G)

                if (verbose2) then
                    write(*,'(A, 6(1x,ES20.12))') "R_minus = ", R_minus
                end if

                do i = 1,6
                    jac(i,j) = (R_plus(i) - R_minus(i))/(2*trim_fd_step)
                end do
                G(j) = G(j) + trim_fd_step
                
            end do

            if (verbose2) then
                write(*,*) ""
                write(*,*) "Jacobian = "
                do i = 1, 6
                    write(*,'(6ES20.12)') jac(i, :)
                end do
            end if

            R_neg = -1*R_
            if (verbose2) then
                write(*,*) ""
                write(*,*) "G used to calculate R_neg which goes into lu_solve:"
                write(*,'(A, 6(1x,ES20.12))') "      G  = ", G
                write(*,'(A, 6(1x,ES20.12))') "  R_neg  = ", R_neg
            end if
            call lu_solve(6, jac, R_neg, delta_G) 
            


            G = G + trim_r_factor*delta_G
            if (verbose2) then
                write(*,*) ""
                write(*,'(A, 6(1x,ES20.12))') " Delta G = ", delta_G
                write(*,'(A, 6(1x,ES20.12))') "   G new = ", G
            end if
            R_ = calc_R(p_trim, q_trim, r_trim, theta, G)
            
            eps = norm2(R_)

            
            if (verbose2) then
                write(*,*) ""
                write(*,'(A, 6(1x,ES20.12))') "      R_ = ", R_
                write(*,'(A, 6(1x,ES20.12))') "     eps = ", eps
                write(*,*) "iter ", iter, " complete"
            end if
            
            iter = iter+1
        end do
        
        
        
        alpha = G(1)
        
        if (trim_type == "shss" .and. has_sideslip) then
            phi = G(2)
            beta = sideslip0 
        else
            beta = G(2)
            phi = phi0
        end if
            
        u  = V0*cos(alpha)*cos(beta)   ! u
        v  = V0*sin(beta)               ! v
        w  = V0*sin(alpha)*cos(beta)

        if (has_gamma) then
            a = u*V0*sin(gamma0)
            b = (v*sin(phi) + w*cos(phi))*sqrt(u*u + (v*sin(phi) + w*cos(phi))**2 - V0*V0*sin(gamma0)*sin(gamma0))
            c = u*u + (v*sin(phi) + w*cos(phi))**2
            theta_plus = asin((a + b)/c)
            theta_minus = asin((a - b)/c)
            err_theta_plus = abs(u*sin(theta_plus) - (v*sin(phi) + w*cos(phi))*cos(theta_plus) - V0*sin(gamma0))
            err_theta_minus = abs(u*sin(theta_minus) - (v*sin(phi) + w*cos(phi))*cos(theta_minus) - V0*sin(gamma0))
            ! write(*,*) "err plus = ", err_theta_plus
            ! write(*,*) "err minus = ", err_theta_minus
            if (err_theta_plus < err_theta_minus) then
                theta = theta_plus
            else 
                theta = theta_minus
            end if
            if (verbose2) then
                write(*,*)""
                write(*,'(A, 1(1x,ES20.12))') "      theta_1[deg] = ", theta_plus*180./pi
                write(*,'(A, 1(1x,ES20.12))') "      theta_2[deg] = ", theta_minus*180./pi
                write(*,'(A, 1(1x,ES20.12))') "correct theta[deg] = ", theta*180./pi
                write(*,'(A, 1(1x,ES20.12))') "correct theta[rad] = ", theta
                write(*,*)""
            end if 
            else
                theta = theta0
        end if  

        c1 = grav*sin(phi)*cos(theta)/(u*cos(theta)*cos(phi) + w*sin(theta))

        if (trim_type == "sct") then
            p_trim = -c1*sin(theta)
            q_trim =  c1*sin(phi)*cos(theta)
            r_trim =  c1*cos(phi)*cos(theta)
        end if

        
        if (trim_type == "vbr") then
            pw = c1*(-u*sin(theta) + v*sin(phi)*cos(theta) + w*cos(phi)*cos(theta))/V0
            p_trim = pw*u/V0
            q_trim = pw*v/V0
            r_trim = pw*w/V0
        end if 
        
        
        ans(1)  = V0*cos(alpha)*cos(beta)   ! u
        ans(2)  = V0*sin(beta)               ! v
        ans(3)  = V0*sin(alpha)*cos(beta)   ! w
        ans(4)  = p_trim              ! p
        ans(5)  = q_trim              ! q
        ans(6)  = r_trim              ! r
        ans(7)  = 0.0             ! xf
        ans(8)  = 0.0             ! yf
        ans(9)  = zf0             ! zf
        ans(10) = phi            ! phi
        ans(11) = theta          ! theta
        ans(12) = psi            ! psi
        
        !  convert phi, theta, and psi to a quat
        ans(10:13) = euler_to_quat(ans(10:12))
        
        controls(1) = G(3)
        controls(2) = G(4)
        controls(3) = G(5)
        controls(4) = G(6)
        
        if (verbose) then
            write(*,*) ""
            write(*,*) "Solved in ", iter-1, "iterations"
            write(*,'(A, 1(1x,ES20.12))') ""
            write(*,'(A, 1(1x,ES20.12))') "       Trim settings "
            write(*,'(A, 1(1x,ES20.12))') "      alpha[deg] = ", alpha*180./pi
            write(*,'(A, 1(1x,ES20.12))') "       beta[deg] = ", beta*180./pi
            write(*,'(A, 1(1x,ES20.12))') "        p[deg/s] = ", p_trim*180./pi
            write(*,'(A, 1(1x,ES20.12))') "        q[deg/s] = ", q_trim*180./pi
            write(*,'(A, 1(1x,ES20.12))') "        r[deg/s] = ", r_trim*180./pi
            write(*,'(A, 1(1x,ES20.12))') "         da[deg] = ", controls(1)*180./pi
            write(*,'(A, 1(1x,ES20.12))') "         de[deg] = ", controls(2)*180./pi
            write(*,'(A, 1(1x,ES20.12))') "         dr[deg] = ", controls(3)*180./pi
            write(*,'(A, 1(1x,ES20.12))') "        throttle = ", controls(4)
            write(*,*)"" 
            write(*,'(A, 1(1x,ES20.12))') "        phi[deg] = ", phi*180./pi
            if (has_gamma) then
                write(*,'(A, 1(1x,ES20.12))') "climb_angle[deg] = ", gamma0*180./pi
            end if              
            write(*,'(A, 1(1x,ES20.12))') "      theta[deg] = ", theta*180./pi
            write(*,'(A, 1(1x,ES20.12))') "      theta[deg] = ", theta*180./pi
            write(*,'(A, 1(1x,ES20.12))') "        psi[deg] = ", psi*180./pi
        end if

    end function trim_solver


    subroutine pseudo_aero(y)
        implicit none

        real, dimension(13) :: y
        real :: u, v, w, p, q, r, xf, yf, zf, phi, theta, psi, alpha, alpha_hat, beta, beta_f, V_mag
        real :: CL1, CL, CS, CD, C_l, Cm, Cn, Re
        real :: Z_pot,T,Press,rho,a, mu, ca, sa, cb, sb, pbar, qbar, rbar
        real :: da, de, dr, tau
        real :: FTx, FTy, FTz, MTx, MTy, MTz

        if (verbose .and. .not. trimming) then
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
        ! phi = y(10)
        ! theta = y(11)
        ! psi = y(12)

        da = controls(1)
        de = controls(2)
        dr = controls(3)
        tau = controls(4)

        V_mag = sqrt(u*u + v*v + w*w)

        ! calculate alpha, V
        alpha  = atan2(w, u)
        alpha_hat = 0.
        beta = asin(v/V_mag)
        ! beta_f = atan2(v,u)
        beta_f = beta
        ca = cos(alpha)
        sa = sin(alpha)
        cb = cos(beta)
        sb = sin(beta)
        
        ! get rho, mu, Re
        call std_atm_English(-zf, Z_pot,T,Press,rho,a, mu)

        ! calc bars
        pbar = p*ref_lat/(2*V_mag)
        qbar = q*ref_long/(2*V_mag)
        rbar = r*ref_lat/(2*V_mag)
        
        ! calc forces
        CL1 = CL0 + CLa*alpha
        CL = CL1 + CLqbar*qbar + CLahat*alpha_hat + CLde*de 
        CS = CSb*beta + (CSpbar + CSapbar*alpha)*pbar + CSrbar*rbar + CSda*da + CSdr*dr
        CD = CD0 + CD1*CL1 + CD2*CL1**2 + CDS2*CS**2 &
        + (CDqbar + CDaqbar*alpha)*qbar + (CDde + CDade*alpha)*de + CDde2*de**2
        C_l = Clb*beta + Clpbar*pbar + (Clrbar + Clarbar*alpha)*rbar + Clda*da + Cldr*dr
        Cm = Cm0 + Cma*alpha + Cmqbar*qbar + Cmahat*alpha_hat + Cmde*de
        Cn = Cnb*beta + (Cnpbar + Cnapbar*alpha)*pbar + Cnrbar*rbar + (Cnda + Cnada*alpha)*da + Cndr*dr
        
        ! limit throttle to 0 to 100
        if (tau < 0.0) then
            tau = 0.0
        else if (tau > 1.0) then
            tau = 1.0
        end if
        
        ! Calculate forces and moments due to throttle
        FTx = tau*Throt_0*(rho/rho0)**Throt_a
        FTy = 0.
        FTz = 0.
        MTx = 0.
        MTy = 0.
        MTz = 0.
        ! trimming = .false. ! for debugging
        ! Re  = rho*V_mag*ref_length/mu
        if (verbose2 .and. .not. trimming) then
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
            write(*,*) "sw = ", ref_lat
            write(*,*) "cw = ", ref_long
            write(*,*) "cos(alpha) = ", ca
            write(*,*) "sin(alpha) = ", sa
            write(*,*) "cos(beta) = ", cb
            write(*,*) "sin(beta) = ", sb
            write(*,*) "ref_area = ", ref_area
            write(*,*) "alpha = ", alpha
            write(*,*) "beta = ", beta
            write(*,*) "beta_flank = ", beta_f
            write(*,*) "rho = ", rho
            write(*,*) "CL1 = ", CL1
            write(*,*) "CL = ", CL
            write(*,*) "CS = ", CS
            write(*,*) "CD = ", CD
            write(*,*) "C_l = ", C_l
            write(*,*) "Cm = ", Cm
            write(*,*) "Cn = ", Cn
            write(*,*) "rho0 = ", rho0
            write(*,*) "da [deg] = ", controls(1)*180./pi
            write(*,*) "de [deg] = ", controls(2)*180./pi
            write(*,*) "dr [deg] = ", controls(3)*180./pi
            write(*,*) "throttle = ", controls(4)

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
        FM(1)=  FTx + -0.5*rho*(V_mag**2)*ref_area*(CD*ca*cb + CS*ca*sb - CL*sa) ! + tau*pow(rho/m_rho0, m_a)*(m_T0 + m_T1*V + m_T2*V*V) !F_xb
        FM(2)=  FTy +  0.5*rho*(V_mag**2)*ref_area*(CS*cb - CD*sb)           ! F_yb
        FM(3)=  FTz + -0.5*rho*(V_mag**2)*ref_area*(CD*sa*cb + CS*sa*sb + CL*ca) ! F_zb
        FM(4)=  MTx +  0.5*rho*(V_mag**2)*ref_area*ref_lat*C_l  ! M_xb
        FM(5)=  MTy +  0.5*rho*(V_mag**2)*ref_area*ref_long*Cm  ! M_yb
        FM(6)=  MTz +  0.5*rho*(V_mag**2)*ref_area*ref_lat*Cn  ! M_zb
        ! write(*,*) FM(5), rho, V_mag, ref_area, ref_long, Cm
       
        if (verbose .and. .not. trimming) then
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
        FM(4:6) = FM(4:6) + cross3(aero_ref_loc, FM(1:3))

        if (verbose .and. .not. trimming) then
            write(*,*) ""
            write(*,*) "Moments after areo reference shift"
            write(*,*) "Mx = ", FM(4)
            write(*,*) "My = ", FM(5)
            write(*,*) "Mz = ", FM(6)
            write(*, '(A, 6(1x,ES20.12))')"      pseuo aero = ", FM
        end if
    
    end subroutine pseudo_aero


    function diff_eq(y) result(dydt)
        
        implicit none
        ! real, intent(in) :: t 
        real, dimension(:), intent(in) :: y
        real, dimension(size(y)) :: dydt
        real :: Fxb, Fyb, Fzb, Mxb, Myb, Mzb, u, v, w, p, q, r, xf, yf, zf, e0, ex, ey, ez, phi, theta, psi, &
        g, weight, Ixx, Iyy, Izz, Ixz, Ixy, Iyz
        real, dimension(4) :: quat_A, quat_B, quat_AB, quat_xyz
        real, dimension(3) :: hpqr, pqr_dot_stuff, pqr_dot

        
        
        Ixx =  I(1,1)
        Iyy =  I(2,2)
        Izz =  I(3,3)
        Ixy = -I(1,2)
        Ixz = -I(1,3)
        Iyz = -I(2,3)

        if (verbose2 .and. .not. trimming) then
            write(*,*) "Ixy = ", Ixy
            write(*,*) "Ixz = ", Ixz
            write(*,*) "Iyz = ", Iyz
        end if
        
        
        call pseudo_aero(y)
        
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

        ! weight = mass*g

        hpqr = matmul(h, y(4:6))
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
        pqr_dot = matmul(I_inv, pqr_dot_stuff)
        
        dydt(1)  = (Fxb/mass) + (g*2.0*(ex*ez - ey*e0)) + (r*v) - (q*w) ! udot
        dydt(2)  = (Fyb/mass) + (g*2.0*(ey*ez + ex*e0)) + (p*w) - (r*u) ! vdot
        dydt(3)  = (Fzb/mass) + (g*(ez**2 + e0**2 - ex**2 - ey**2)) + (q*u) - (p*v) ! wdot
        dydt(4)  = pqr_dot(1) ! pdot
        dydt(5)  = pqr_dot(2) ! qdot
        dydt(6)  = pqr_dot(3) ! rdot
        
        dydt(7:9) = quat_dependent_to_base(y(1:3), y(10:13)) + Vwf
        
        dydt(10) = 0.5*(-ex*p - ey*q - ez*r) ! e0 dot
        dydt(11) = 0.5*( e0*p - ez*q + ey*r) ! ex dot
        dydt(12) = 0.5*( ez*p + e0*q - ex*r) ! ey dot
        dydt(13) = 0.5*(-ey*p + ex*q + e0*r) ! ez dot
        if (verbose .and. .not. trimming) then
            write(*, '(A, 13(1x,ES20.12))')"         diff eq = ", dydt
            write(*,*)""
        end if
    end function diff_eq


    function rk4(t0, y0, dt) result(ans)
    
        implicit none
        
        real, intent(in) :: t0, dt
        real, dimension(:), intent(in) :: y0
        real, dimension(size(y0)) :: k1, k2, k3, k4, ans, y_temp
        integer :: i, n

        n = size(y0)
    
        ! k1:
        if (verbose) then
            write(*, *)"        RK call =  ", 1
            write(*, '(A, 1(1x,ES20.12))')"           time =  ", t0
        end if

        k1 = diff_eq(y0)

        ! multiply k1 by dt 
        ! for (int i = 0; i < size; i++)
        do i=1, n
            k1(i) = dt * k1(i)
            ! use k1 to update y_temporary
            y_temp(i) = y0(i) + 0.5*k1(i)
        end do

        ! k2:  
        if (verbose) then
            write(*, *)"        RK call =  ", 2
            write(*, '(A, 1(1x,ES20.12))')"           time =  ", t0 + 0.5*dt
        end if
        k2 = diff_eq(y_temp)
        
        ! multiply k2 by dt and then update y_temp
        do i = 1, n
        
            k2(i) = dt * k2(i)
            y_temp(i) = y0(i) + 0.5*k2(i)
        end do
        
        ! k3:
        if (verbose) then
            write(*, *)"        RK call =  ", 3
            write(*, '(A, 1(1x,ES20.12))')"           time =  ",t0 + 0.5*dt
        end if
        k3 = diff_eq(y_temp)
        
        !  multiply k2 by dt and then update y_temp
        do i = 1, n
          
            k3(i) = dt * k3(i)
            y_temp(i) = y0(i) + k3(i)
        end do
        
        ! k4:
        if (verbose) then
            write(*, *)"        RK call =  ", 4
            write(*, '(A, 1(1x,ES20.12))')"           time =  ", t0 + dt
        end if
        k4 = diff_eq(y_temp)

        ! multiply k4 by dt   
        do i = 1, n
        
            k4(i) = dt * k4(i)

            ! find dxdt, then dzdt
            ans(i) = y0(i) + (k1(i) + 2.0*k2(i) + 2.0*k3(i) + k4(i))/6.0  
        end do
   
   
    end function rk4


    function vehicle_tick_state(t,time,dt) result(y1)
        ! controls = controls_udp%recv()
            
        y1 = rk4(t, y, dt)
        call quat_norm(y1(10:13))
        y = y1
        t = t + dt
        integrated_time = integrated_time + dt
        
        if (verbose) then
            write(*,*)
            write(*, '(A, (1x,ES20.12))') "           time = ", t
            write(*, '(A, 13(1x,ES20.12))')"          state  = ", y
            write(*,*) ""
        end if
        
        ! write to file
        write(10, '(15(1x,ES20.12))') t, dt, y(1:13)

        ! send data to connection
        s(1) = t
        s(2:14) = y(1:13)
        call graphics%send(s)

        if (real_time) then
            ! call system_clock(time2)
            time2 = get_time()
            dt = time2 - time1
            ! dt = real(time2 - time1)/count_rate
            ! if (dt <= TOLERANCE) dt = 0.00001 / real(count_rate)
            time1 = time2
        end if

    end function vehicle_tick_state


end module vehicle_m