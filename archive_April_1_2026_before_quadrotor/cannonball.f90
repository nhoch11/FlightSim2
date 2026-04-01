module cannonball_m

    use hoch_m
    use json_mod
    use json_xtnsn_mod
    
    implicit none
    
    real :: weight, mass
    real, dimension(3,3) :: I, I_inv, h
    real, dimension(6) :: FM
    real, dimension(3) :: Vwf

    real :: t0, dt, ref_area, ref_length
    real :: v_inf, zf0, theta0, p0, q0, r0, CLa, CD0, CD2, Cma, Cmq, Clp, Cl_0
    real :: Ixx0, Iyy0, Izz0, Ixz0, Ixy0, Iyz0

    


contains 

    subroutine read_in(input_file)
        character(len=*), intent(in) :: input_file

        type(json_file) :: input_json
        type(json_value),pointer :: sim_settings, &
                                    ref_settings, &
                                    initial_settings, &
                                    mass_settings, &
                                    aero_settings
        logical :: exists, found

        real, dimension(4) :: quat

        ! Set up run
        call json_initialize()

        ! Get input file from command line
        ! call getarg(1, input_file)
        ! input_file = trim(input_file)

        ! input_file = "12_pounder.json"

        ! Get input file from user
        ! if (input_file == '') then
        !     write(*,*) "Please specify an input file:"
        !     read(*,'(a)') input_file
        !     input_file = trim(input_file)
        ! end if

        ! Check it exists
        inquire(file=input_file, exist=exists)
        if (.not. exists) then
            write(*,*) "!!! The file ", input_file, " does not exist. Quitting..."
            stop
        end if

        ! Load settings from input file
        call input_json%load_file(filename=input_file)
        call json_check()
        call input_json%get('simulation', sim_settings, found)
        call input_json%get('reference', ref_settings, found)
        call input_json%get('initial', initial_settings, found)
        call input_json%get('mass', mass_settings, found)
        call input_json%get('aerodynamics', aero_settings, found)

        call json_xtnsn_get(sim_settings, "time_step_s", dt, 1.0)
        write(*,*)dt
        call json_xtnsn_get(ref_settings, "area_ft2",  ref_area)
        call json_xtnsn_get(ref_settings, "length_ft",  ref_length)
        
        call json_xtnsn_get(initial_settings, "airspeed_ft_per_s",  v_inf)
        call json_xtnsn_get(initial_settings, "altitude_ft",  zf0)
        call json_xtnsn_get(initial_settings, "elevation_angle_deg",  theta0)
        ! call json_xtnsn_get(initial_settings, "azimuth[deg]]",  theta)
        call json_xtnsn_get(initial_settings, "initial_p_rad_per_s",  p0)
        call json_xtnsn_get(initial_settings, "initial_q_rad_per_s",  q0)
        call json_xtnsn_get(initial_settings, "initial_r_rad_per_s",  r0)
        
        ! call json_xtnsn_get(mass_settings, "mass_slug",  mass)
        call json_xtnsn_get(mass_settings, "weight_lbf",  weight)
        call json_xtnsn_get(mass_settings, "Ixx_slug_ft2",  Ixx0, 0.)
        call json_xtnsn_get(mass_settings, "Iyy_slug_ft2",  Iyy0, 0.)
        call json_xtnsn_get(mass_settings, "Izz_slug_ft2",  Izz0, 0.)
        call json_xtnsn_get(mass_settings, "Ixy_slug_ft2",  Ixy0, 0.)
        call json_xtnsn_get(mass_settings, "Iyz_slug_ft2",  Iyz0, 0.)
        
        call json_xtnsn_get(aero_settings, "CLa",  CLa)
        call json_xtnsn_get(aero_settings, "CD0",  CD0)
        call json_xtnsn_get(aero_settings, "CD2",  CD2)
        call json_xtnsn_get(aero_settings, "Cma",  Cma)
        call json_xtnsn_get(aero_settings, "Cmq",  Cmq)
        call json_xtnsn_get(aero_settings, "Clp",  Clp)
        call json_xtnsn_get(aero_settings, "Cl_0",  Cl_0)
        

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

        real :: r, det,Ixxyyzz


        ! calc mass
        r = ref_length/2.
        mass = weight*.3048/GSSL
        write(*,*) "mass = ", mass
        
        Ixxyyzz = mass*2.*(r**5 - (r-0.00131)**5)/(5.*(r**3 - (r-0.00131)**3))

        I(1,1) =  Ixxyyzz
        I(2,2) =  Ixxyyzz
        I(3,3) =  Ixxyyzz
        
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

    function sphere_CD(Re) result(ans)

        implicit none

        real, intent(in) :: Re
        real :: ans

        ! page 105

        ! check if i reynolds input is out of range
        if (0.< Re .and. Re <= 450000.) then
            ans =  24./Re + 6./(1 + sqrt(Re)) + 0.4

        else if (450000. < Re .and. Re <= 560000.) then
            ans = 1.0e29*Re**(-5.211) 

        else if (560000. < Re .and. Re <= 14000000.) then
            ans = -2.0e-23 * Re**3 -1.0e-6 * Re**2 + 9.0e-9*Re + 0.0069

        else if (Re > 14000000.) then
            ans = 0.12
        else
            write(*,*) "Error, invalid Reynolds number"
            stop
        end if
    
    end function sphere_CD



    subroutine pseudo_aero(y)
        implicit none

        real, dimension(13) :: y
        real :: u, v, w, p, q, r, xf, yf, zf, phi, theta, psi, alpha, beta, V_mag, CL, CS, CD, C_l, Cm, Cn, Re
        real :: Z_pot,T,Press,rho,a, mu

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
        phi = y(10)
        theta = y(11)
        psi = y(12)
    
        ! calculate alpha, V
        alpha  = atan2(w, u);
        beta = atan2(v,u);
        
        V_mag = sqrt(u**2 + v**2 + w**2)

        ! calculate CL, CD, and Cm
        ! CL =  m_CLa * alpha;
        ! CS =  m_CLa * beta;
        ! CD =  m_CD0 + (m_CD2 * (CL**2));
        ! C_l =  m_Cl0 + ((m_Clp * m_ref_length * p) / V_mag);
        ! Cm =  (m_Cma * alpha) + ((m_Cmq * m_ref_length * q) / V_mag);
        ! Cn = -(m_Cma * beta) +  ((m_Cmq * m_ref_length * r) / V_mag);

        ! get rho, mu, Re
        call std_atm_English(-zf, Z_pot,T,Press,rho,a, mu)
        Re  = rho*V_mag*ref_length/mu
        write(*,*) "u = ", u
        write(*,*) "v = ", v
        write(*,*) "w = ", w
        write(*,*) "rho = ", rho
        write(*,*) "mu = ", mu

        write(*,*) "Re = ", Re
        
        ! get drag of a sphere
        CD = sphere_CD(Re)

        write(*,*) "CD = ", CD

        ! update forces and moments
        ! FM(1)=  -0.5*rho*(V_mag**2)*pi*(ref_length/2.)**2*CD*cos(alpha)*cos(beta); ! F_xb
        ! FM(2)=  -0.5*rho*(V_mag**2)*pi*(ref_length/2.)**2*CD*sin(beta);            ! F_yb
        ! FM(3)=  -0.5*rho*(V_mag**2)*pi*(ref_length/2.)**2*CD*sin(alpha)*cos(beta); ! F_zb
        
        FM(1)=  -0.5*rho*(V_mag**2)*pi*(ref_length/2.)**2*CD*u/V_mag ! F_xb
        FM(2)=  -0.5*rho*(V_mag**2)*pi*(ref_length/2.)**2*CD*v/V_mag            ! F_yb
        FM(3)=  -0.5*rho*(V_mag**2)*pi*(ref_length/2.)**2*CD*w/V_mag ! F_zb
        FM(4)=  0.0; ! 0.5*rho*(V_mag**2)*m_ref_area*m_ref_length*Cl;  ! M_xb
        FM(5)=  0.0; ! 0.5*rho*(V_mag**2)*m_ref_area*m_ref_length*Cm;  ! M_yb
        FM(6)=  0.0; ! 0.5*rho*(V_mag**2)*m_ref_area*m_ref_length*Cn;  ! M_zb
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
        y(9)  = -zf0 ! zf
        y(10)  = 0.0              ! phi
        y(11) = theta0           ! theta
        y(12) = 0.0              ! psi

        !  convert phi, theta, and psi to a quat
        y(10:13) = euler_to_quat(y(10:12))

        !  normalize the quat
        call quat_norm(y(10:13))

        t = t0
        
        write(*,*)t, dt

        ! Open output file (unit number 10 here, adjust as needed)
        open(unit=10, file="cannonball_output.txt", status="replace", action="write", iostat=ios)
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


    


end module cannonball_m