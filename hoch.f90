module hoch_m

    implicit none
    
    real, parameter :: PI = 3.1415926535897932384626433832795
    real, parameter :: TOLERANCE = 1.0e-14
    real, parameter :: REZ = 6356766.0 ! [m] radius of Earth used to compute geopotential altitude
    real, parameter :: RE = 6366707.01949371 ! [m] mean sea-level radius of the Earth
    real, parameter :: GSSL = 9.80665 ! [m/s^2] gravity at sea level
    real, parameter :: RGAS = 287.0528 ! gas constant
    real, parameter :: PSSL = 101325.0 ! pressure at sea level
    real, parameter :: GAMMA = 1.4 ! heat capacity ratio for dry air
    real, parameter :: Z_const(8) = [   0.,11000.,20000.,32000.,47000.,52000.,61000.,79000.]
    real, parameter :: T_array(8) = [288.15,216.65,216.65,228.65,270.65,270.65,252.65,180.65]
    real, parameter :: T_prime(8) = [-0.0065, 0.0, 0.0010, 0.0028, 0.0, -0.0020, -0.0040, 0.0]
    ! real, parameter :: p_array(8) = [1.01325e5, 2.26320318222212e4, 5.47487352827083e3, 8.68014769086723e2, &
    !             1.10905588989225e2, 5.90005242789244e1, 1.82099249050177e1, 1.03770045489203]
    real, parameter :: p_array(8) = [1.0132500000000000E+5, 2.2632031822221168E+4 , &
    5.4748735282708267E+3 , 8.6801476908672271E+2 , 1.1090558898922531E+2 , &
    5.9000524278924367E+1 , 1.8209924905017658E+1 , 1.0377004548920223E+0]
    ! real, parameter ::  mu0 = 1.716e-5 ! dynamic viscosity constant for Sutherlands Formula
    ! real, parameter :: T0 = 273.15  ! Kelvin for Sutherlands formula
contains

    subroutine std_atm_SI(H, Z, T, P, rho, a, mu)
        implicit none
        ! input:
            ! H = geometric altitude [ft]
        ! output:
            ! Z = geopotential altitude [m]
            ! T = temperature [K]
            ! P = pressure [N/m**2]
            ! rho = density [kg/m**3]
            ! a = speed of sound [m/s]

        real, intent(in) :: H
        real, intent(out):: Z,T,P,rho,a,mu !lambda
        integer :: i
        ! real Z_const(0:8),T_array(0:7),T_prime(0:7),p_array(0:7)
        ! data Z_const/    0.,11000.,20000.,32000.,47000.,52000.,61000.,79000.,9.9e20/
        ! data T_array/ 288.15,216.65,216.65,228.65,270.65,270.65,252.65,180.65/ !,180.65/
        ! data T_prime/-0.0065, 0.0, 0.0010, 0.0028, 0.0, -0.0020, -0.0040, 0.0/
        ! data p_array/101325.0, 22632.0491189944,5474.88167406511, 868.016875642432, 110.905974487888,59.0007483456162, 18.210004975660, 1.03770045489203/ ! Pa
        real :: P0, T0, mu0, C 

        P0 = PSSL

        if (H < 0.0) then
            Z = 0.
        else
            Z=REZ*H/(REZ+H)
        end if
        ! write(*,*) H

        ! dynamic viscosity using sutherlands formula
        ! double mu; // kg/(m*s)
        ! double mu0 = 1.716e-5; // dynamic viscosity constant for Sutherlands Formula
        ! double T0 = 273.15; !Kelvin for Sutherlands formula

        ! const double R = 287.0528; ! Nm/kgK
        
        ! store pre calculated temperatures and pressures at each layer boundary
        ! T_array(7) = {288.150, 216.650, 216.650, 228.650, 270.650, 270.650, 252.650}; !K
        ! p_array(7) = {101325.0, 22632.0491189944,5474.88167406511, 868.016875642432, 110.905974487888,59.0007483456162, 18.2100049756603 }; ! Pa
        ! T_prime(7) = {-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.002, -0.004}; ! K/m
    
        !  first layer
        if (Z <= 11000.0) then
            T = T_array(1) + T_prime(1)*(Z-Z_const(1)); !  K
            p = p_array(1)*(T/T_array(1))**(-GSSL/(RGAS*T_prime(1)))

            ! write(*,*) Z, Z_const(1), T_array(1), T_prime(1), P_array(1), gssl, rgas, p, T
            
        
        !  second layer
        else if (Z > 11000.0 .and. Z <= 20000.0) then
        ! add the portion of the second altitude range to the previously integrated first layer
            T = T_array(2) + T_prime(2)*(Z-Z_const(2)) !  K
            ! use the Ti prime = 0 equation
            p = p_array(2)*exp(-GSSL*(Z-Z_const(2))/(RGAS*T_array(2)))

        
            
        
        ! third layer
        else if (32000.0 >= Z .and. Z > 20000.0) then
            ! add the portion of the third altitude range to the previously integrated 2 layers
            T = T_array(3) + T_prime(3)*(Z-Z_const(3)) !  K
            ! use the Ti prime /= 0 equation
            p = p_array(3)*(T/T_array(3))**(-GSSL/(RGAS*T_prime(3)))
                    
        
        !fourth layer
        else if (47000.0 >= Z .and. Z > 32000.0) then
            ! add the portion of the fourth altitude to the previously integrated 3 layers
            T = T_array(4) + T_prime(4)*(Z-Z_const(4)) !  K
            ! use the Ti prime /= 0 equation
            p = p_array(4)*(T/T_array(4))**(-GSSL/(RGAS*T_prime(4)))
                    
        
        ! fifth layer
        else if (52000.0 >= Z .and. Z > 47000.0) then
            ! add the portion of the fifth altitude range to the previously integrated 4 layers
            T = T_array(5) + T_prime(5)*(Z-Z_const(5)) !  K
            ! use the Ti prime /= 0 equation
            p = p_array(5)*exp((-GSSL*(Z-Z_const(5)))/(RGAS*T_array(5)))
        
        
        ! sixth layer
        else if (61000.0 >= Z .and. Z > 52000.0) then      
            ! add the portion of the sixth altitude range to the prevously integrated 5 layers
            T = T_array(6) + T_prime(6)*(Z-Z_const(6)) ! K
            ! use the Ti prime /= 0 equation
            p = p_array(6)*(T/T_array(6))**(-GSSL/(RGAS*T_prime(6)))
        
        ! seventh layer
        else if (79000.0 >= Z .and. Z > 61000.0) then
            ! add the portion of the seventh altitude range to the previously integrated 6 layers
            T = T_array(7) + T_prime(7)*(Z-Z_const(7)) ! K
            ! use the Ti prime /= 0 equation
            p = p_array(7)*(T/T_array(7))**(-GSSL/(RGAS*T_prime(7)))
            
        else if (90000.0 >= Z .and. Z > 79000.0) then
            ! add the portion of the seventh altitude range to the previously integrated 6 layers
            T = T_array(8) + T_prime(8)*(Z-Z_const(8)) ! K
            ! use the Ti prime /= 0 equation
            p = p_array(8)*exp((-GSSL*(Z-Z_const(8)))/(RGAS*T_array(8)))
            
        else
            write(*,*) "std_atm_si error, out of atmosphere"
        end if
        rho = p/(RGAS*T);

        ! write(*,*) "rho  , p, R, T"
        ! write(*,*) rho, p, RGAS, T
        a = sqrt(1.4*RGAS*T);
        T0 = 273.15; ! Kelvin
        mu0 = 0.00001716; ! kg/m-s
        C = 110.4; ! Sutherland's Constant for air
        mu = mu0*(T0+C)/(T+C)*(T/T0)**1.5;
        return
        
    end subroutine std_atm_SI


    subroutine std_atm_English(H,Z,T,P,rho,a, mu)
        ! input:
            ! H = geometric altitude [ft]
        ! output:
            ! Z = geopotential altitude [ft]
            ! T = temperature [R]
            ! P = pressure [lbf/ft^2]
            ! rho = density [slug/ft^3]
            ! a = speed of sound [ft/s]
        
        real :: H,Z,T,P,rho,a, mu
        real :: Z1, T1, P1, rho1, a1, mu1
        
        call std_atm_SI(H*0.3048, Z1, T1, P1, rho1, a1, mu1)
        Z = Z1/0.3048 ! meters to ft
        T = T1*1.8 ! Kelvin to Rankine
        P = P1/47.880258 ! Pa to lbf/ft^2
        rho = rho1/515.379 ! kg/m^3 to slugs/ ft^3
        a = a1/0.3048 ! m/s to ft/s
        ! kg/(m*s) to slug/(ft*s)
        ! mu = mu1*0.020885434225
        mu = mu1/47.880258

        return
    end subroutine std_atm_English


    function gravity_SI(H) result(g)
        implicit none
        real, intent(in) :: H 
        real :: g
        
        g = GSSL*(REZ/(REZ+H))**2
        
    end function gravity_SI


    function gravity_English(H) result(g)
    
        implicit none
        real, intent(in) :: H
        real :: g
        g = gravity_SI(0.3048*H)/0.3048

    end function gravity_English

    
    function quat_mult(A, B) result(quat_result)

        implicit none

        real, dimension(4), intent(in) :: A, B
        real, dimension(4) :: quat_result

        quat_result(1) = A(1)*B(1) - A(2)*B(2) - A(3)*B(3) - A(4)*B(4)
        quat_result(2) = A(1)*B(2) + A(2)*B(1) + A(3)*B(4) - A(4)*B(3)
        quat_result(3) = A(1)*B(3) - A(2)*B(4) + A(3)*B(1) + A(4)*B(2)
        quat_result(4) = A(1)*B(4) + A(2)*B(3) - A(3)*B(2) + A(4)*B(1)

    end function quat_mult

    
    function quat_base_to_dependent(base,quat) result(dependent)

        implicit none

        real, intent(in) :: base(3),quat(4)
        real :: T(4), dependent(3)
        
        T(1) = -base(1)*quat(2) - base(2)*quat(3) - base(3)*quat(4)
        T(2) =  base(1)*quat(1) + base(2)*quat(4) - base(3)*quat(3)
        T(3) = -base(1)*quat(4) + base(2)*quat(1) + base(3)*quat(2)
        T(4) =  base(1)*quat(3) - base(2)*quat(2) + base(3)*quat(1)
        
        dependent(1) = quat(1)*T(2) - quat(2)*T(1) - quat(3)*T(4) + quat(4)*T(3)
        dependent(2) = quat(1)*T(3) + quat(2)*T(4) - quat(3)*T(1) - quat(4)*T(2)
        dependent(3) = quat(1)*T(4) - quat(2)*T(3) + quat(3)*T(2) - quat(4)*T(1)

    end function quat_base_to_dependent


    function quat_dependent_to_base(dependent, quat) result(base)

        implicit none

        real, intent(in) :: dependent(3), quat(4)
        real :: T(4), base(3)
        
        T(1) =  dependent(1)*quat(2) + dependent(2)*quat(3) + dependent(3)*quat(4)
        T(2) =  dependent(1)*quat(1) - dependent(2)*quat(4) + dependent(3)*quat(3)
        T(3) =  dependent(1)*quat(4) + dependent(2)*quat(1) - dependent(3)*quat(2)
        T(4) = -dependent(1)*quat(3) + dependent(2)*quat(2) + dependent(3)*quat(1)
        
        base(1) = quat(1)*T(2) + quat(2)*T(1) + quat(3)*T(4) - quat(4)*T(3)
        base(2) = quat(1)*T(3) - quat(2)*T(4) + quat(3)*T(1) + quat(4)*T(2)
        base(3) = quat(1)*T(4) + quat(2)*T(3) - quat(3)*T(2) + quat(4)*T(1)

    end function quat_dependent_to_base


    subroutine quat_norm(quat)

        implicit none

        real, intent(inout) :: quat(4)
        real :: norm
        norm = sqrt(sum(quat**2))
        quat = quat/norm

    end subroutine quat_norm


    function euler_to_quat(euler) result(quat)
        ! euler angles in radians

        implicit none

        real, intent(in) :: euler(3)
        real :: C_half_phi, S_half_phi, C_half_theta, S_half_theta, C_half_psi, S_half_psi, quat(4)

        C_half_phi = cos(euler(1)/2.) 
        S_half_phi = sin(euler(1)/2.) 
        C_half_theta = cos(euler(2)/2.)
        S_half_theta = sin(euler(2)/2.)
        C_half_psi = cos(euler(3)/2.)
        S_half_psi = sin(euler(3)/2.)

        quat(1) = C_half_phi*C_half_theta*C_half_psi + S_half_phi*S_half_theta*S_half_psi
        quat(2) = S_half_phi*C_half_theta*C_half_psi - C_half_phi*S_half_theta*S_half_psi
        quat(3) = C_half_phi*S_half_theta*C_half_psi + S_half_phi*C_half_theta*S_half_psi
        quat(4) = C_half_phi*C_half_theta*S_half_psi - S_half_phi*S_half_theta*C_half_psi
        

    end function euler_to_quat

    
    function quat_to_euler(quat) result(euler)
    
        implicit none

        real, intent(in) :: quat(4)
        real :: a, euler(3)

        a = quat(1)*quat(3) - quat(2)*quat(4)
        
        if (abs(a - 0.5) < 1.e-12) then
            euler(1) = 2.0*asin(quat(2)/cos(PI/4.0))
            euler(2) = PI/2.0
            euler(3) = 0.
        else if (abs(a + 0.5) < 1.e-12) then
            euler(1) = 2.0*asin(quat(2)/cos(PI/4.0));
            euler(2) = -PI/2.0
            euler(3) = 0.0
        else
            ! write(*,*) "here"
            euler(1) = atan2( 2.0*(quat(1)*quat(2) + quat(3)*quat(4)), &
                quat(1)*quat(1) + quat(4)*quat(4) - quat(2)*quat(2) - quat(3)*quat(3))
            euler(2) = asin(2.0*(quat(1)*quat(3) - quat(2)*quat(4)))
            euler(3) = atan2( 2.0*(quat(1)*quat(4) + quat(2)*quat(3)), &
                quat(1)*quat(1) + quat(2)*quat(2) - quat(3)*quat(3) - quat(4)*quat(4))
        end if


    end function quat_to_euler


    function cross3(A, B) result(ans)

        implicit none

        real, dimension(3), intent(in) :: A, B
        real, dimension(3) :: ans

        ans(1) =   A(2)*B(3) - A(3)*B(2)
        ans(2) =   A(3)*B(1) - A(1)*B(3)
        ans(3) =   A(1)*B(2) - A(2)*B(1)


    end function


end module hoch_m 
