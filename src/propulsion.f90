module propulsion_m
    use hoch_m
    use atmosphere_m
    use jsonx_m

    implicit none

    type propulsion_type
        character(len=:), allocatable :: name, type, units
        real, allocatable :: location(:), orientation_eul(:)
        real :: orientation_quat(4), rho0
        integer :: control_ID ! ID of controls array associated with this propulsion element

        ! for type = T=f(V)
        real, allocatable :: T_coefs(:)
        real :: Ta

        ! for type = propeller_polynomial
        real :: diameter, Ixx
        integer :: rotation_delta
        real, allocatable :: CT_J(:), CP_J(:), CNa_J(:), Cnna_J(:)
        real :: motor_Kv, motor_no_load_amp, motor_ohm
        real :: batt_no_load_volt, batt_ohm
        real :: esc_ohm
        real :: last_Nr

    end type propulsion_type


contains

    subroutine propulsion_init(this, j_propulsion)
        implicit none
        type(propulsion_type), intent(inout) :: this
        type(json_value), pointer :: j_propulsion
        character(len=:), allocatable :: temp
        real :: junkZ, junkP, junkT, junka, junkmu

        this%name = j_propulsion%name
        write(*,*) "     Initializing Propulsion : ", trim(this%name)

        call jsonx_get(j_propulsion, "location[ft]", this%location, 0.0, 3)
        call jsonx_get(j_propulsion, "orientation[deg]", this%orientation_eul, 0.0, 3)
        this%orientation_eul = this%orientation_eul*PI/180.0
        this%orientation_quat = euler_to_quat(this%orientation_eul)

        call jsonx_get(j_propulsion, "type", this%type)
        select case (this%type)
        case("T=f(V)")
            call jsonx_get(j_propulsion, "T_coefficients[lbf]", this%T_coefs, 0.0)
            call jsonx_get(j_propulsion, "Ta", this%Ta)
            this%rotation_delta = 1
            this%Ixx = 0.0
        case("propeller_polynomial")
            call jsonx_get(j_propulsion, "diameter[ft]", this%diameter)
            call jsonx_get(j_propulsion, "Ixx[slug-ft^2]", this%Ixx)
            call jsonx_get(j_propulsion, "rotation", temp)
            this%rotation_delta = 1
            if(trim(temp) == "LH") this%rotation_delta = -1
            call jsonx_get(j_propulsion, "CT(J)", this%CT_J, 0.0)
            call jsonx_get(j_propulsion, "CPb(J)", this%CP_J, 0.0)
            call jsonx_get(j_propulsion, "CN,alpha(J)", this%CNa_J, 0.0)
            call jsonx_get(j_propulsion, "Cn,alpha(J)", this%Cnna_J, 0.0) 
            ! read in motor
            call jsonx_get(j_propulsion, "motor.Kv", this%motor_Kv, 0.0)
            call jsonx_get(j_propulsion, "motor.no_load_current[amp]", this%motor_no_load_amp, 0.0)
            call jsonx_get(j_propulsion, "motor.resistance[ohm]", this%motor_ohm, 0.0)
            ! read in battery
            call jsonx_get(j_propulsion, "battery.no_load_voltage[volt]", this%batt_no_load_volt, 0.0)
            call jsonx_get(j_propulsion, "battery.resistance[ohm]", this%batt_ohm, 0.0)
            ! read in esc
            call jsonx_get(j_propulsion, "esc.resistance[ohm]", this%esc_ohm, 0.0)
            this%last_Nr = 1000.0
        end select

        call std_atm_English(0.0, junkz, junkT, junkP, this%rho0, junka, junkmu)

    end subroutine propulsion_init

    
    function propulsion_get_FMh(this, states, tau) result(ans)
        implicit none
        type(propulsion_type), intent(inout) :: this
        real, intent(in) :: states(21), tau
        real :: ans(9) ! returns propulsion FMh
        real :: Vc(3), Vc_mag, uc(3), vN(3), vN_mag, uN(3), alpha_c
        real :: CT0, CT1, CT2, CP0, CP1, CP2
        real :: Fc(3), Mc(3) ! subscript c is for component coordinate system
        real :: thrust, normal, torque, yaw, hxx
        real :: rho, junkZ, junkT, junkP, junka, junkmu
        real :: Hz, omega, J, Nr

        
        ! velocity in component coord 
        Vc = quat_base_to_dependent(states(1:3) + cross3(states(4:6), this%location), this%orientation_quat) 
        Vc_mag = sqrt(Vc(1)**2 + Vc(2)**2 + Vc(3)**2)
        uc = Vc/Vc_mag
        if (Vc_mag < TOLERANCE) uc = [1.0, 0.0, 0.0]
        alpha_c = acos(uc(1)) 
        
        ! normal velocity component
        vN = -[0.0, uc(2), uc(3)]
        vN_mag = sqrt(vN(1)**2 + vN(2)**2 + vN(3)**2)
        uN = vN/Vn_mag
        if (vN_mag < TOLERANCE) uN = [0.0, 0.0, 1.0]

        call std_atm_English(-states(9), junkZ, junkT, junkP, rho, junka, junkmu) 

        select case (this%type)
            case("T=f(V)")
                thrust = tau*calc_polynomial(this%T_coefs, Vc_mag)*(rho/this%rho0)**this%Ta
                normal = 0.0
                torque = 0.0
                yaw = 0.0
                hxx = 0.0
            case("propeller_polynomial")
                select case (this%units)
                    case ("throttle")
                        ! tau is given as throttle, use electric propulsion function to solve for RPM
                        Nr = solve_electric_propulsion(this, tau, Vc_mag, rho)
                        Hz = Nr/60.0 ! Revolutions per second
                        omega = Hz*2*PI ! radians per second   Hz = omega/2pi = rotations per second
                    case ("RPM")
                        ! tau is given as RPM
                        Hz = tau/60.0 ! Revolutions per second
                        omega = Hz*2*PI ! radians per second   Hz = omega/2pi = rotations per second
                    case("thrust[lbf]")
                        ! tau is given as thrust[lbf]
                        CT0 = this%CT_J(1)
                        CT1 = this%CT_J(2)
                        CT2 = this%CT_J(3)
                        omega = PI/this%diameter/CT0*(-Vc_mag*CT1 + sqrt((Vc_mag*CT1)**2 - 4.*CT0*(Vc_mag**2*CT2 - tau/(rho*this%diameter**2))))
                        Hz = omega/(2*PI)
                    case("torque[ft-lbf/s]")
                        ! tau is given as torque[ft-lbf]
                        CP0 = this%CP_J(1)
                        CP1 = this%CP_J(2)
                        CP2 = this%CP_J(3)
                        omega = PI/this%diameter/CP0*(-Vc_mag*CP1 + sqrt((Vc_mag*CP1)**2 - 4.*CP0*(Vc_mag**2*CP2 - 2.*PI*tau/(rho*this%diameter**3))))
                        Hz = omega/(2*PI)
                end select

                J = 2.0*PI*Vc_mag/omega/this%diameter ! advance ratio

                thrust = rho*(Hz**2)*(this%diameter**4)*calc_polynomial(this%CT_J, J)
                torque = rho*(Hz**3)*(this%diameter**5)*calc_polynomial(this%CP_J, J) / omega
                normal = rho*(Hz**2)*(this%diameter**4)*calc_polynomial(this%CNa_J, J)*alpha_c
                yaw    = rho*(Hz**2)*(this%diameter**5)*calc_polynomial(this%Cnna_J, J)*alpha_c
                hxx    = this%rotation_delta*this%Ixx*omega 
                
                ! write(*,*) "rho    = ", rho
                ! write(*,*) "u,v,w  = ", states(1:3)
                ! write(*,*) "x,y,z  = ", states(4:6)
                ! write(*,*) "Vc     = ", Vc
                ! write(*,*) "Vc_mag = ", Vc_mag
                ! write(*,*) "uc     = ", uc
                ! write(*,*) "Vn_mag = ", Vn_mag
                ! write(*,*) "uN     = ", uN
                ! write(*,*) "Hz     = ", Hz
                ! write(*,*) "omega  = ", omega
                ! write(*,*) "J      = ", J
                ! write(*,*) ""
                ! write(*,*) "Ex 4.5.4"
                ! write(*,*) ""
                ! write(*,*) "thrust [lbf]      = ", thrust
                ! write(*,*) "torque [ft-lbf/s] = ", torque
                ! write(*,*) "normal [lbf]      = ", normal
                ! write(*,*) "yaw [lbf]         = ", yaw
                ! write(*,*) "hxx    = ", hxx
        end select

        Fc = [thrust, 0.0, 0.0] + normal*uN
        Mc = -real(this%rotation_delta)*([torque, 0.0, 0.0] + yaw*uN)

        ! convert to body coordinates
        ans(1:3) = quat_dependent_to_base(Fc, this%orientation_quat)
        ans(4:6) = quat_dependent_to_base(Mc, this%orientation_quat) + cross3(this%location, ans(1:3))
        ans(7:9) = quat_dependent_to_base([hxx, 0.0, 0.0], this%orientation_quat)
        ! write(*,*)"ans(7:9) = ", ans(7:9)
        ! write(*,*)"hxx = ", hxx
        ! write(*,*)"orientqtion_quat = ", this%orientation_quat
        ! write(*,*)"ans = ", ans

        ! write(*,*)""
        ! write(*,*)""
        ! write(*,*)"Rotor      Rotation                                               Fc[lbf]                                                     Mc[ft-lbf]"
        ! write(*,*)this%control_ID, ",", this%rotation_delta, ",", Fc(1), ",",  Fc(2), ",", Fc(3), ",", Mc(1), ",",Mc(2), ",",Mc(3)         
        ! write(*,*)"Rotor      Rotation                                               Fb[lbf]                                                     Mb[ft-lbf]"
        ! write(*,*)this%control_ID, ",", this%rotation_delta, ",", ans(1), ",",  ans(2), ",", ans(3), ",", ans(4), ",",ans(5), ",",ans(6)         
        ! write(*,*)"Rotor      Rotation                                               hc[slug-ft^2/s]                                            hb[slug-ft^2s] "
        ! write(*,*)this%control_ID, ",", this%rotation_delta, ",", hxx, ", 0.0, 0.0,", ans(7), ",",ans(8), ",",ans(9)         

    end function propulsion_get_FMh


    function solve_electric_propulsion(this, tau, Vc_mag, rho) result(Nr)
        implicit none
        type(propulsion_type), intent(inout) :: this
        real, intent(in) :: tau, Vc_mag, rho
        real :: Nr ! RPM of rotor
        real :: relax, error, Hz, omega, J
        real :: lr, ls, lm, Im, eta_C, B, C, Em, Nm, Ns, Eb, Ib
        integer :: iter

        ! initial guess for rotor RPM
        Nr = this%last_Nr
        error = 1.0
        relax = 0.9
        iter = 1
        if (this%name == "prop_1") then
            write(*,*) "      iter        Nr [RPM]                Ns [RPM]                  error [RPM]              "
        end if

        do while (abs(error) > TOLERANCE)
            ! compute rotor torque
            Hz = Nr/60.0 ! Revolutions per second
            ! if (this%name == "prop_1") then
            !     ! write(*,*) " First Iteration"
            !     ! write(*,*) " Nr    = ", Nr
            ! end if
            omega = Hz*2*PI
            J = 2.0*PI*Vc_mag/omega/this%diameter ! advance ratio
            lr = rho*(Hz**3)*(this%diameter**5)*calc_polynomial(this%CP_J, J) / omega
            ls = lr ! set shaft torque to rotor torque
            lm = ls ! motor torque lm = ls/eta_g/Gm if there is a gearbox
            Im = lm*this%motor_Kv/C_I + this%motor_no_load_amp ! motor current
            eta_C = 1.0 - 0.078*(1.0 - tau) ! esc efficiency
            B = 2.*Im*this%esc_ohm - tau*this%batt_no_load_volt + tau**2*Im*this%batt_ohm/eta_C
            C = Im**2*this%esc_ohm**2 - tau*this%batt_no_load_volt*Im*this%esc_ohm
            Em = 0.5*(-B + sqrt(B**2 - 4.0*C)) ! motor voltage
            Nm = this%motor_Kv*(Em - Im*this%motor_ohm) ! motor rotation rate
            Ns = Nm ! shaft rotation rate, modify if theres a gear box
            error = Ns - Nr
            if (this%name == "prop_1") then
                write(*,*)iter, Nr, Ns, error
                ! write(*,*) " lr    = ", lr
                ! write(*,*) " ls    = ", ls
                ! write(*,*) " lm    = ", lm
                ! write(*,*) " Im    = ", Im
                ! write(*,*) " eta_C = ", eta_C
                ! write(*,*) " B     = ", B
                ! write(*,*) " C     = ", C
                ! write(*,*) " Em    = ", Em
                ! write(*,*) " Nm    = ", Nm
                ! write(*,*) " Ns    = ", Ns
                ! write(*,*) " error = ", error
                ! write(*,*) " Nr    = ", Nr
                ! write(*,*) ""
                ! write(*,*) " Nr    = ", Nr
            end if
            
            iter = iter + 1
            Nr = Nr + relax*(Ns - Nr)
            
        end do
        Eb = (Nm/this%motor_Kv + Im*(this%motor_ohm + this%esc_ohm))/tau
        Ib = (this%batt_no_load_volt - Eb)/this%batt_ohm
        this%last_Nr = Nr
        if (this%name == "prop_1") then
            write(*,*) ""
            write(*,*) " Eb    = ", Eb
            write(*,*) " Ib    = ", Ib
            write(*,*) ""
        end if

    end function solve_electric_propulsion


end module propulsion_m