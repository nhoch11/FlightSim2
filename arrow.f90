module arrow_m

    use hoch_m
    use jsonx_m
    ! use json_xtnsn_mod
    
    implicit none
    
    real :: weight, mass
    real, dimension(3,3) :: I, I_inv, h
    real, dimension(6) :: FM
    real, dimension(3) :: Vwf

    real :: t0, dt, ref_area, ref_long, ref_lat
    real :: v_inf, zf0, theta0, p0, q0, r0, CLa, CD0, CD2, Cma, Cmqbar, Clpbar, C_l0
    real :: Ixx0, Iyy0, Izz0, Ixz0, Ixy0, Iyz0

    


contains 

    subroutine read_in(input_file)
        implicit none
        character(len=*), intent(in) :: input_file
        type(json_value), pointer :: j_main
        
        write(*,*) "made it here"
        write(*,*) input_file
        ! call get_command_argument(1,input_file)
        write(*,*) "after"
        
        call jsonx_load(input_file,j_main)
        
        call jsonx_get(j_main, "simulation.time_step[s]", dt)
        call jsonx_get(j_main, "reference.area[ft^2]",  ref_area)
        call jsonx_get(j_main, "reference.longitudinal[ft]",  ref_long)
        call jsonx_get(j_main, "reference.lateral[ft]",  ref_lat)
        
        call jsonx_get(j_main, "initial.airspeed[ft/s]",  v_inf)
        call jsonx_get(j_main, "initial.altitude[ft]",  zf0)
        zf0 = -1.*zf0
        call jsonx_get(j_main, "initial.elevation_angle[deg]",  theta0)
        theta0 = theta0*pi/180.
        ! call jsonx_get(j_main, "initial.azimuth[deg]]",  theta)
        call jsonx_get(j_main, "initial.initial_p[rad/s]",  p0)
        call jsonx_get(j_main, "initial.initial_q[rad/s]",  q0)
        call jsonx_get(j_main, "initial.initial_r[rad/s]",  r0)
        
        ! call json_xtnsn_get(mass_settings, "mass_slug",  mass)
        call jsonx_get(j_main, "mass.weight[lbf]",  weight)
        call jsonx_get(j_main, "mass.Ixx[slug*ft^2]",  Ixx0, 0.)
        call jsonx_get(j_main, "mass.Iyy[slug*ft^2]",  Iyy0, 0.)
        call jsonx_get(j_main, "mass.Izz[slug*ft^2]",  Izz0, 0.)
        call jsonx_get(j_main, "mass.Ixy[slug*ft^2]",  Ixy0, 0.)
        call jsonx_get(j_main, "mass.Iyz[slug*ft^2]",  Iyz0, 0.)
        
        call jsonx_get(j_main, "aerodynamics.CLa",  CLa)
        call jsonx_get(j_main, "aerodynamics.CD0",  CD0)
        call jsonx_get(j_main, "aerodynamics.CD2",  CD2)
        call jsonx_get(j_main, "aerodynamics.Cma",  Cma)
        call jsonx_get(j_main, "aerodynamics.Cmqbar",  Cmqbar)
        call jsonx_get(j_main, "aerodynamics.C_l0",  C_l0)
        call jsonx_get(j_main, "aerodynamics.Clpbar",  Clpbar)

        ! write(*,*)dt, ref_area, ref_long, ref_lat, v_inf, zf0, theta0, weight, &
        !  Ixx0, Iyy0, Izz0, CLa, CD0, CD2, Cma, Cmqbar, C_l0, Clpbar
        ! stop
        

        h(1,1) = 0.
        h(1,2) = 0.
        h(1,3) = 0.
        h(2,1) = 0.
        h(2,2) = 0.
        h(2,3) = 0.
        h(3,1) = 0.
        h(3,2) = 0.
        h(3,3) = 0.
        
        Vwf(1) = 0.
        Vwf(2) = 0.
        Vwf(3) = 0.
        !  print results in a file
        ! FILE* en_file = fopen("cannonball.txt", "w");
        ! fprintf(en_file, "  Time[s]              u[ft/s]            v[ft/s]                w[ft/s]              p[rad/s]             q[rad/s]             r[rad/s]             x[ft]                y[ft]                z[ft]                e0                   ex                   ey                   ez\n");
    
    end subroutine read_in

    function rk4(t0, y0, dt) result(ans)
    
        implicit none
        
        real, intent(in) :: t0, dt
        real, dimension(:), intent(in) :: y0
        real, dimension(size(y0)) :: k1, k2, k3, k4, ans, y_temp
        integer :: i, n

        n = size(y0)
    
        ! k1:
        k1 = diff_eq(t0, y0)

        ! multiply k1 by dt 
        ! for (int i = 0; i < size; i++)
        do i=1, n
            k1(i) = dt * k1(i)
            ! use k1 to update y_temporary
            y_temp(i) = y0(i) + 0.5*k1(i)
        end do

        ! k2:  
        k2 = diff_eq(t0 + 0.5*dt, y_temp)
        
        ! multiply k2 by dt and then update y_temp
        do i = 1, n
        
            k2(i) = dt * k2(i)
            y_temp(i) = y0(i) + 0.5*k2(i)
        end do
        
        ! k3:
        k3 = diff_eq(t0 + 0.5*dt, y_temp)
        
        !  multiply k2 by dt and then update y_temp
        do i = 1, n
          
            k3(i) = dt * k3(i)
            y_temp(i) = y0(i) + k3(i)
        end do
        
        ! k4:
        k4 = diff_eq(t0 + dt, y_temp)

        ! multiply k4 by dt   
        do i = 1, n
        
            k4(i) = dt * k4(i)

            ! find dxdt, then dzdt
            ans(i) = y0(i) + (k1(i) + 2.0*k2(i) + 2.0*k3(i) + k4(i))/6.0  
        end do
   
   
    end function rk4



    function diff_eq(t, y) result(dydt)
        
        implicit none
        real, intent(in) :: t 
        real, dimension(:), intent(in) :: y
        real, dimension(size(y)) :: dydt
        real :: Fxb, Fyb, Fzb, Mxb, Myb, Mzb, u, v, w, p, q, r, xf, yf, zf, e0, ex, ey, ez, phi, theta, psi, &
        g, weight, Ixx, Iyy, Izz, Ixz, Ixy, Iyz
        real, dimension(4) :: quat_A, quat_B, quat_AB, quat_xyz
        real, dimension(3) :: hpqr, pqr_dot_stuff, pqr_dot

        
        write(*, *)" time =  ", t
        Ixx =  I(1,1)
        Iyy =  I(2,2)
        Izz =  I(3,3)
        Ixy = -I(1,2)
        Ixz = -I(1,3)
        Iyz = -I(2,3)
        
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

        pqr_dot_stuff(1) = hpqr(1) + Mxb + (Iyy - Izz)*q*r + Iyz*(q*q - r*r) + Ixz*p*q - Ixy*p*r
        pqr_dot_stuff(2) = hpqr(2) + Myb + (Izz - Ixx)*p*r + Ixz*(r*r - p*p) + Ixy*q*r - Iyz*p*q
        pqr_dot_stuff(3) = hpqr(3) + Mzb + (Ixx - Iyy)*p*q + Ixy*(p*p - q*q) + Iyz*p*r - Ixz*q*r

        pqr_dot = matmul(I_inv, pqr_dot_stuff)
        dydt(1)  = (Fxb/mass) + (g*2.0*(ex*ez - ey*e0)) + (r*v) - (q*w) ! udot
        write(*,*) "Looj HEEREEEEE"
        write(*,*) (1./mass),Fxb, g, ex, ez, ey, e0, r, v, q, w
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
        write(*, '(A, 13(1x,ES20.12))')" diff eq = ", dydt
        write(*,*)""
        
    end function diff_eq


    subroutine mass_inertia()
        implicit none

        real :: det


        ! calc mass
        ! r = ref_length/2.
        mass = weight*.3048/GSSL
        write(*,*) "mass = ", mass
        
        ! Ixxyyzz = mass*2.*(r**5 - (r-0.00131)**5)/(5.*(r**3 - (r-0.00131)**3))

        I(1,1) =  Ixx0
        I(2,2) =  Iyy0
        I(3,3) =  Izz0
        
        I(1,3) = -0.0
        I(3,1) = -0.0
        I(1,2) = -0.0
        I(2,1) = -0.0
        I(2,3) = -0.0
        I(3,2) = -0.0


        ! calculate I_inv by hand

        det = I(1,1)*(I(2,2)*I(3,3) - I(2,3)*I(3,2)) &
            - I(1,2)*(I(2,1)*I(3,3) - I(2,3)*I(3,1)) &
            + I(1,3)*(I(2,1)*I(3,2) - I(2,2)*I(3,1))

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

        return
    end subroutine mass_inertia


    subroutine pseudo_aero(y)
        implicit none

        real, dimension(13) :: y
        real :: u, v, w, p, q, r, xf, yf, zf, phi, theta, psi, alpha, beta, beta_f, V_mag, CL, CS, CD, C_l, Cm, Cn, Re
        real :: Z_pot,T,Press,rho,a, mu, ca, sa, cb, sb, pbar, qbar, rbar

        write(*, '(A, 13(1x,ES20.12))')" state coming in = ", y
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
        
        V_mag = sqrt(u**2 + v**2 + w**2)

        ! calculate alpha, V
        alpha  = atan2(w, u)
        beta_f = atan2(v,u)
        beta = asin(v/V_mag)
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


        ! Re  = rho*V_mag*ref_length/mu
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
        
        
        ! calc forces
        CL = CLa*alpha
        CS = -CLa*beta_f 
        CD = CD0 + CD2*CL**2 + CD2*CS**2
        C_l = C_l0 + Clpbar*pbar
        Cm = Cma*alpha + Cmqbar*qbar
        Cn = -Cma*beta_f + Cmqbar*rbar
        
        write(*,*) "alpha = ", alpha
        write(*,*) "beta = ", beta
        write(*,*) "beta_flank = ", beta_f
        write(*,*) "rho = ", rho
        write(*,*) "CL = ", CL
        write(*,*) "CS = ", CS
        write(*,*) "CD = ", CD
        write(*,*) "C_l = ", C_l
        write(*,*) "Cm = ", Cm
        write(*,*) "Cn = ", Cn

        ! update forces and moments
        FM(1)=  -0.5*rho*(V_mag**2)*ref_area*(CD*ca*cb + CS*ca*sb - CL*sa) ! + tau*pow(rho/m_rho0, m_a)*(m_T0 + m_T1*V + m_T2*V*V) !F_xb
        FM(2)=   0.5*rho*(V_mag**2)*ref_area*(CS*cb - CD*sb)           ! F_yb
        FM(3)=  -0.5*rho*(V_mag**2)*ref_area*(CD*sa*cb + CS*sa*sb + CL*ca) ! F_zb
        FM(4)=   0.5*rho*(V_mag**2)*ref_area*ref_lat*C_l  ! M_xb
        FM(5)=   0.5*rho*(V_mag**2)*ref_area*ref_long*Cm  ! M_yb
        FM(6)=   0.5*rho*(V_mag**2)*ref_area*ref_lat*Cn  ! M_zb
        write(*, '(A, 6(1x,ES20.12))')" pseuo aero = ", FM

    
    end subroutine pseudo_aero


    subroutine run()

        implicit none

        real :: t, tf
        real, dimension(13) :: y, y1
        real, dimension(3) :: eul
        integer :: ios
        call mass_inertia()

        y(1)  = v_inf         ! u
        y(2)  = 0.0              ! v
        y(3)  = 0.0              ! w
        y(4)  = 0.0              ! p
        y(5)  = 0.0              ! q
        y(6)  = 0.0              ! r
        y(7)  = 0.0              ! xf
        y(8)  = 0.0              ! yf
        y(9)  = zf0 ! zf
        y(10)  = 0.0              ! phi
        y(11) = theta0   ! theta
        y(12) = 0.0              ! psi

        !  convert phi, theta, and psi to a quat
        y(10:13) = euler_to_quat(y(10:12))

        !  normalize the quat
        call quat_norm(y(10:13))

        t = t0
        
        write(*,*)t, dt

        ! Open output file (unit number 10 here, adjust as needed)
        open(unit=10, file="arrow_output.txt", status="replace", action="write", iostat=ios)
        if (ios /= 0) then
            write(*,*) "Error opening file!"
            stop
        end if

        write(10,'(14(1x,A20))') "time[s]", "u[ft/s]", "v[ft/s]", "w[ft/s]", "p[rad/s]", "q[rad/s]", "r[rad/s]", &
        "xf[ft]", "yf[ft]", "zf[ft]", "e0", "ex", "ey", "ez"

        write(10, '(14(1x,ES20.12))') t, y(1:13)
        
        do while (-y(9) > 0.0)
            ! write(*,*) -y(9)
            write(*,*) "time = ", t
            write(*, '(A, 13(1x,ES20.12))')" state  = ", y

            y1 = rk4(t, y, dt)
            
            ! write to file
            
            !  re-normalize the quat
            call quat_norm(y1(10:13))
            
            !  y = y0 ( copy y0 into y)
            y = y1
            
            !  add time step
            t = t + dt

            write(10, '(14(1x,ES20.12))') t, y(1:13)

        end do 
        ! while (t0<0.005);
        
        close(10)


    end subroutine


    


end module arrow_m