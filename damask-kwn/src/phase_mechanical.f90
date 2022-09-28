!----------------------------------------------------------------------------------------------------
!> @brief internal microstructure state for all plasticity constitutive models
!----------------------------------------------------------------------------------------------------
submodule(phase) mechanical


  enum, bind(c); enumerator :: &
    PLASTICITY_UNDEFINED_ID, &
    PLASTICITY_NONE_ID, &
    PLASTICITY_ISOTROPIC_ID, &
    PLASTICITY_PHENOPOWERLAW_ID, &
    PLASTICITY_KWNPOWERLAW_ID, &
    PLASTICITY_KINEHARDENING_ID, &
    PLASTICITY_DISLOTWIN_ID, &
    PLASTICITY_DISLOTUNGSTEN_ID, &
    PLASTICITY_NONLOCAL_ID, &
    KINEMATICS_UNDEFINED_ID, &
    KINEMATICS_CLEAVAGE_OPENING_ID, &
    KINEMATICS_SLIPPLANE_OPENING_ID, &
    KINEMATICS_THERMAL_EXPANSION_ID
  end enum

  type(tTensorContainer), dimension(:), allocatable :: &
    ! current value
    phase_mechanical_Fe, &
    phase_mechanical_Fi, &
    phase_mechanical_Fp, &
    phase_mechanical_F, &
    phase_mechanical_Li, &
    phase_mechanical_Lp, &
    phase_mechanical_S, &
    phase_mechanical_P, &
    ! converged value at end of last solver increment
    phase_mechanical_Fi0, &
    phase_mechanical_Fp0, &
    phase_mechanical_F0, &
    phase_mechanical_Li0, &
    phase_mechanical_Lp0, &
    phase_mechanical_S0


  integer(kind(PLASTICITY_undefined_ID)), dimension(:),   allocatable :: &
    phase_plasticity                                                                                !< plasticity of each phase

  integer :: phase_plasticity_maxSizeDotState

  interface

    module subroutine eigendeformation_init(phases)
      class(tNode), pointer :: phases
    end subroutine eigendeformation_init

    module subroutine elastic_init(phases)
      class(tNode), pointer :: phases
    end subroutine elastic_init

    module subroutine plastic_init
    end subroutine plastic_init

    module subroutine phase_hooke_SandItsTangents(S,dS_dFe,dS_dFi,Fe,Fi,ph,en)
      integer, intent(in) :: &
        ph, &
        en
      real(pReal),   intent(in),  dimension(3,3) :: &
        Fe, &                                                                                       !< elastic deformation gradient
        Fi                                                                                          !< intermediate deformation gradient
      real(pReal),   intent(out), dimension(3,3) :: &
        S                                                                                           !< 2nd Piola-Kirchhoff stress tensor in lattice configuration
      real(pReal),   intent(out), dimension(3,3,3,3) :: &
        dS_dFe, &                                                                                   !< derivative of 2nd P-K stress with respect to elastic deformation gradient
        dS_dFi                                                                                      !< derivative of 2nd P-K stress with respect to intermediate deformation gradient
    end subroutine phase_hooke_SandItsTangents

    module subroutine plastic_isotropic_LiAndItsTangent(Li,dLi_dMi,Mi,ph,en)
      real(pReal), dimension(3,3),     intent(out) :: &
        Li                                                                                          !< inleastic velocity gradient
      real(pReal), dimension(3,3,3,3), intent(out)  :: &
        dLi_dMi                                                                                     !< derivative of Li with respect to Mandel stress
      real(pReal), dimension(3,3),     intent(in) :: &
        Mi                                                                                          !< Mandel stress
      integer,                         intent(in) :: &
        ph, &
        en
    end subroutine plastic_isotropic_LiAndItsTangent

    module function plastic_dotState(subdt,co,ip,el,ph,en) result(broken)
      integer, intent(in) :: &
        co, &                                                                                           !< component-ID of integration point
        ip, &                                                                                           !< integration point
        el, &                                                                                           !< element
        ph, &
        en
      real(pReal),  intent(in) :: &
        subdt                                                                                           !< timestep
      logical :: broken
    end function plastic_dotState

    module function plastic_deltaState(ph, en) result(broken)
      integer, intent(in) :: &
        ph, &
        en
      logical :: &
        broken
    end function plastic_deltaState

    module subroutine phase_LiAndItsTangents(Li, dLi_dS, dLi_dFi, &
                                             S, Fi, ph,en)
      integer, intent(in) :: &
        ph,en
      real(pReal),   intent(in),  dimension(3,3) :: &
        S                                                                                               !< 2nd Piola-Kirchhoff stress
      real(pReal),   intent(in),  dimension(3,3) :: &
        Fi                                                                                              !< intermediate deformation gradient
      real(pReal),   intent(out), dimension(3,3) :: &
        Li                                                                                              !< intermediate velocity gradient
      real(pReal),   intent(out), dimension(3,3,3,3) :: &
        dLi_dS, &                                                                                       !< derivative of Li with respect to S
        dLi_dFi

    end subroutine phase_LiAndItsTangents


    module subroutine plastic_LpAndItsTangents(Lp, dLp_dS, dLp_dFi, &
                                               S, Fi, ph,en)
      integer, intent(in) :: &
        ph,en
      real(pReal),   intent(in),  dimension(3,3) :: &
        S, &                                                                                            !< 2nd Piola-Kirchhoff stress
        Fi                                                                                              !< intermediate deformation gradient
      real(pReal),   intent(out), dimension(3,3) :: &
        Lp                                                                                              !< plastic velocity gradient
      real(pReal),   intent(out), dimension(3,3,3,3) :: &
        dLp_dS, &
        dLp_dFi                                                                                         !< derivative of Lp with respect to Fi
    end subroutine plastic_LpAndItsTangents


    module subroutine plastic_isotropic_results(ph,group)
      integer,          intent(in) :: ph
      character(len=*), intent(in) :: group
    end subroutine plastic_isotropic_results

    module subroutine plastic_phenopowerlaw_results(ph,group)
      integer,          intent(in) :: ph
      character(len=*), intent(in) :: group
    end subroutine plastic_phenopowerlaw_results

    module subroutine plastic_kwnpowerlaw_results(ph,group)
      integer,          intent(in) :: ph
      character(len=*), intent(in) :: group
    end subroutine plastic_kwnpowerlaw_results

    module subroutine plastic_kinehardening_results(ph,group)
      integer,          intent(in) :: ph
      character(len=*), intent(in) :: group
    end subroutine plastic_kinehardening_results

    module subroutine plastic_dislotwin_results(ph,group)
      integer,          intent(in) :: ph
      character(len=*), intent(in) :: group
    end subroutine plastic_dislotwin_results

    module subroutine plastic_dislotungsten_results(ph,group)
      integer,          intent(in) :: ph
      character(len=*), intent(in) :: group
    end subroutine plastic_dislotungsten_results

    module subroutine plastic_nonlocal_results(ph,group)
      integer,          intent(in) :: ph
      character(len=*), intent(in) :: group
    end subroutine plastic_nonlocal_results

    module function plastic_dislotwin_homogenizedC(ph,en) result(homogenizedC)
      real(pReal), dimension(6,6) :: homogenizedC
      integer,     intent(in) :: ph,en
    end function plastic_dislotwin_homogenizedC


  end interface
  type :: tOutput                                                                                   !< new requested output (per phase)
    character(len=pStringLen), allocatable, dimension(:) :: &
      label
  end type tOutput
  type(tOutput), allocatable, dimension(:) :: output_constituent

  procedure(integrateStateFPI), pointer :: integrateState

contains


!--------------------------------------------------------------------------------------------------
!> @brief Initialize mechanical field related constitutive models
!> @details Initialize elasticity, plasticity and stiffness degradation models.
!--------------------------------------------------------------------------------------------------
module subroutine mechanical_init(materials,phases)

  class(tNode), pointer :: &
    materials, &
    phases

  integer :: &
    ph, &
    en, &
    Nmembers
  class(tNode), pointer :: &
    num_crystallite, &
    material, &
    constituents, &
    phase, &
    mech

  print'(/,a)', ' <<<+-  phase:mechanical init  -+>>>'

!-------------------------------------------------------------------------------------------------
  allocate(output_constituent(phases%length))

  allocate(phase_mechanical_Fe(phases%length))
  allocate(phase_mechanical_Fi(phases%length))
  allocate(phase_mechanical_Fi0(phases%length))
  allocate(phase_mechanical_Fp(phases%length))
  allocate(phase_mechanical_Fp0(phases%length))
  allocate(phase_mechanical_F(phases%length))
  allocate(phase_mechanical_F0(phases%length))
  allocate(phase_mechanical_Li(phases%length))
  allocate(phase_mechanical_Li0(phases%length))
  allocate(phase_mechanical_Lp0(phases%length))
  allocate(phase_mechanical_Lp(phases%length))
  allocate(phase_mechanical_S(phases%length))
  allocate(phase_mechanical_P(phases%length))
  allocate(phase_mechanical_S0(phases%length))

  do ph = 1, phases%length
    Nmembers = count(material_phaseID == ph)

    allocate(phase_mechanical_Fi(ph)%data(3,3,Nmembers))
    allocate(phase_mechanical_Fe(ph)%data(3,3,Nmembers))
    allocate(phase_mechanical_Fi0(ph)%data(3,3,Nmembers))
    allocate(phase_mechanical_Fp(ph)%data(3,3,Nmembers))
    allocate(phase_mechanical_Fp0(ph)%data(3,3,Nmembers))
    allocate(phase_mechanical_Li(ph)%data(3,3,Nmembers))
    allocate(phase_mechanical_Li0(ph)%data(3,3,Nmembers))
    allocate(phase_mechanical_Lp0(ph)%data(3,3,Nmembers))
    allocate(phase_mechanical_Lp(ph)%data(3,3,Nmembers))
    allocate(phase_mechanical_S(ph)%data(3,3,Nmembers),source=0.0_pReal)
    allocate(phase_mechanical_P(ph)%data(3,3,Nmembers),source=0.0_pReal)
    allocate(phase_mechanical_S0(ph)%data(3,3,Nmembers),source=0.0_pReal)
    allocate(phase_mechanical_F(ph)%data(3,3,Nmembers))
    allocate(phase_mechanical_F0(ph)%data(3,3,Nmembers))

    phase   => phases%get(ph)
    mech    => phase%get('mechanical')
#if defined(__GFORTRAN__)
    output_constituent(ph)%label = output_as1dString(mech)
#else
    output_constituent(ph)%label = mech%get_as1dString('output',defaultVal=emptyStringArray)
#endif
  enddo

  do ph = 1, phases%length
    do en = 1, count(material_phaseID == ph)

      phase_mechanical_Fp0(ph)%data(1:3,1:3,en) = phase_orientation0(ph)%data(en)%asMatrix()        ! Fp reflects initial orientation (see 10.1016/j.actamat.2006.01.005)
      phase_mechanical_Fp0(ph)%data(1:3,1:3,en) = phase_mechanical_Fp0(ph)%data(1:3,1:3,en) &
                                                / math_det33(phase_mechanical_Fp0(ph)%data(1:3,1:3,en))**(1.0_pReal/3.0_pReal)
      phase_mechanical_Fi0(ph)%data(1:3,1:3,en) = math_I3
      phase_mechanical_F0(ph)%data(1:3,1:3,en)  = math_I3

      phase_mechanical_Fe(ph)%data(1:3,1:3,en) = math_inv33(matmul(phase_mechanical_Fi0(ph)%data(1:3,1:3,en), &
                                                                   phase_mechanical_Fp0(ph)%data(1:3,1:3,en)))  ! assuming that euler angles are given in internal strain free configuration
    enddo
    phase_mechanical_Fp(ph)%data = phase_mechanical_Fp0(ph)%data
    phase_mechanical_Fi(ph)%data = phase_mechanical_Fi0(ph)%data
    phase_mechanical_F(ph)%data  = phase_mechanical_F0(ph)%data
  enddo


! initialize elasticity
  call elastic_init(phases)

! initialize plasticity
  allocate(plasticState(phases%length))
  allocate(phase_plasticity(phases%length),source = PLASTICITY_undefined_ID)
  allocate(phase_localPlasticity(phases%length),   source=.true.)

  call plastic_init()

  do ph = 1,phases%length
    plasticState(ph)%state0 = plasticState(ph)%state
  enddo
  phase_plasticity_maxSizeDotState = maxval(plasticState%sizeDotState)


  num_crystallite => config_numerics%get('crystallite',defaultVal=emptyDict)

  select case(num_crystallite%get_asString('integrator',defaultVal='FPI'))

    case('FPI')
      integrateState => integrateStateFPI

    case('Euler')
      integrateState => integrateStateEuler

    case('AdaptiveEuler')
      integrateState => integrateStateAdaptiveEuler

    case('RK4')
      integrateState => integrateStateRK4

    case('RKCK45')
      integrateState => integrateStateRKCK45

    case default
     call IO_error(301,ext_msg='integrator')

  end select


  call eigendeformation_init(phases)


end subroutine mechanical_init


module subroutine mechanical_results(group,ph)

  character(len=*), intent(in) :: group
  integer,          intent(in) :: ph


  call crystallite_results(group,ph)

  select case(phase_plasticity(ph))

    case(PLASTICITY_ISOTROPIC_ID)
      call plastic_isotropic_results(ph,group//'mechanical/')

    case(PLASTICITY_PHENOPOWERLAW_ID)
      call plastic_phenopowerlaw_results(ph,group//'mechanical/')

    case(PLASTICITY_KWNPOWERLAW_ID)
      call plastic_kwnpowerlaw_results(ph,group//'mechanical/')

    case(PLASTICITY_KINEHARDENING_ID)
      call plastic_kinehardening_results(ph,group//'mechanical/')

    case(PLASTICITY_DISLOTWIN_ID)
      call plastic_dislotwin_results(ph,group//'mechanical/')

    case(PLASTICITY_DISLOTUNGSTEN_ID)
      call plastic_dislotungsten_results(ph,group//'mechanical/')

    case(PLASTICITY_NONLOCAL_ID)
      call plastic_nonlocal_results(ph,group//'mechanical/')

  end select


end subroutine mechanical_results


!--------------------------------------------------------------------------------------------------
!> @brief calculation of stress (P) with time integration based on a residuum in Lp and
!> intermediate acceleration of the Newton-Raphson correction
!--------------------------------------------------------------------------------------------------
function integrateStress(F,subFp0,subFi0,Delta_t,co,ip,el) result(broken)

  real(pReal), dimension(3,3), intent(in) :: F,subFp0,subFi0
  real(pReal),                 intent(in) :: Delta_t
  integer, intent(in)::         el, &                                                               ! element index
                                      ip, &                                                         ! integration point index
                                      co                                                            ! grain index

  real(pReal), dimension(3,3)::       Fp_new, &                                                     ! plastic deformation gradient at end of timestep
                                      invFp_new, &                                                  ! inverse of Fp_new
                                      invFp_current, &                                              ! inverse of Fp_current
                                      Lpguess, &                                                    ! current guess for plastic velocity gradient
                                      Lpguess_old, &                                                ! known last good guess for plastic velocity gradient
                                      Lp_constitutive, &                                            ! plastic velocity gradient resulting from constitutive law
                                      residuumLp, &                                                 ! current residuum of plastic velocity gradient
                                      residuumLp_old, &                                             ! last residuum of plastic velocity gradient
                                      deltaLp, &                                                    ! direction of next guess
                                      Fi_new, &                                                     ! gradient of intermediate deformation stages
                                      invFi_new, &
                                      invFi_current, &                                              ! inverse of Fi_current
                                      Liguess, &                                                    ! current guess for intermediate velocity gradient
                                      Liguess_old, &                                                ! known last good guess for intermediate velocity gradient
                                      Li_constitutive, &                                            ! intermediate velocity gradient resulting from constitutive law
                                      residuumLi, &                                                 ! current residuum of intermediate velocity gradient
                                      residuumLi_old, &                                             ! last residuum of intermediate velocity gradient
                                      deltaLi, &                                                    ! direction of next guess
                                      Fe, &                                                         ! elastic deformation gradient
                                      S, &                                                          ! 2nd Piola-Kirchhoff Stress in plastic (lattice) configuration
                                      A, &
                                      B, &
                                      temp_33
  real(pReal), dimension(9) ::        temp_9                                                        ! needed for matrix inversion by LAPACK
  integer,     dimension(9) ::        devNull_9                                                     ! needed for matrix inversion by LAPACK
  real(pReal), dimension(9,9) ::      dRLp_dLp, &                                                   ! partial derivative of residuum (Jacobian for Newton-Raphson scheme)
                                      dRLi_dLi                                                      ! partial derivative of residuumI (Jacobian for Newton-Raphson scheme)
  real(pReal), dimension(3,3,3,3)::   dS_dFe, &                                                     ! partial derivative of 2nd Piola-Kirchhoff stress
                                      dS_dFi, &
                                      dFe_dLp, &                                                    ! partial derivative of elastic deformation gradient
                                      dFe_dLi, &
                                      dFi_dLi, &
                                      dLp_dFi, &
                                      dLi_dFi, &
                                      dLp_dS, &
                                      dLi_dS
  real(pReal)                         steplengthLp, &
                                      steplengthLi, &
                                      atol_Lp, &
                                      atol_Li, &
                                      devNull
  integer                             NiterationStressLp, &                                         ! number of stress integrations
                                      NiterationStressLi, &                                         ! number of inner stress integrations
                                      ierr, &                                                       ! error indicator for LAPACK
                                      o, &
                                      p, &
                                      ph, &
                                      en, &
                                      jacoCounterLp, &
                                      jacoCounterLi                                                 ! counters to check for Jacobian update
  logical :: error,broken


  broken = .true.

  ph = material_phaseID(co,(el-1)*discretization_nIPs + ip)
  en = material_phaseEntry(co,(el-1)*discretization_nIPs + ip)

  call plastic_dependentState(co,ip,el)

  Lpguess = phase_mechanical_Lp(ph)%data(1:3,1:3,en)                                              ! take as first guess
  Liguess = phase_mechanical_Li(ph)%data(1:3,1:3,en)                                              ! take as first guess

  call math_invert33(invFp_current,devNull,error,subFp0)
  if (error) return ! error
  call math_invert33(invFi_current,devNull,error,subFi0)
  if (error) return ! error

  A = matmul(F,invFp_current)                                                                       ! intermediate tensor needed later to calculate dFe_dLp

  jacoCounterLi  = 0
  steplengthLi   = 1.0_pReal
  residuumLi_old = 0.0_pReal
  Liguess_old    = Liguess

  NiterationStressLi = 0
  LiLoop: do
    NiterationStressLi = NiterationStressLi + 1
    if (NiterationStressLi>num%nStress) return ! error

    invFi_new = matmul(invFi_current,math_I3 - Delta_t*Liguess)
    Fi_new    = math_inv33(invFi_new)

    jacoCounterLp  = 0
    steplengthLp   = 1.0_pReal
    residuumLp_old = 0.0_pReal
    Lpguess_old    = Lpguess

    NiterationStressLp = 0
    LpLoop: do
      NiterationStressLp = NiterationStressLp + 1
      if (NiterationStressLp>num%nStress) return ! error

      B  = math_I3 - Delta_t*Lpguess
      Fe = matmul(matmul(A,B), invFi_new)
      call phase_hooke_SandItsTangents(S, dS_dFe, dS_dFi, &
                                        Fe, Fi_new, ph, en)

      call plastic_LpAndItsTangents(Lp_constitutive, dLp_dS, dLp_dFi, &
                                         S, Fi_new, ph,en)

      !* update current residuum and check for convergence of loop
      atol_Lp = max(num%rtol_crystalliteStress * max(norm2(Lpguess),norm2(Lp_constitutive)), &      ! absolute tolerance from largest acceptable relative error
                    num%atol_crystalliteStress)                                                     ! minimum lower cutoff
      residuumLp = Lpguess - Lp_constitutive

      if (any(IEEE_is_NaN(residuumLp))) then
        return ! error
      elseif (norm2(residuumLp) < atol_Lp) then                                                     ! converged if below absolute tolerance
        exit LpLoop
      elseif (NiterationStressLp == 1 .or. norm2(residuumLp) < norm2(residuumLp_old)) then          ! not converged, but improved norm of residuum (always proceed in first iteration)...
        residuumLp_old = residuumLp                                                                 ! ...remember old values and...
        Lpguess_old    = Lpguess
        steplengthLp   = 1.0_pReal                                                                  ! ...proceed with normal step length (calculate new search direction)
      else                                                                                          ! not converged and residuum not improved...
        steplengthLp = num%subStepSizeLp * steplengthLp                                             ! ...try with smaller step length in same direction
        Lpguess      = Lpguess_old &
                     + deltaLp * stepLengthLp
        cycle LpLoop
      endif

      calculateJacobiLp: if (mod(jacoCounterLp, num%iJacoLpresiduum) == 0) then
        jacoCounterLp = jacoCounterLp + 1

        do o=1,3; do p=1,3
          dFe_dLp(o,1:3,p,1:3) = - Delta_t * A(o,p)*transpose(invFi_new)                            ! dFe_dLp(i,j,k,l) = -Delta_t * A(i,k) invFi(l,j)
        enddo; enddo
        dRLp_dLp = math_eye(9) &
                 - math_3333to99(math_mul3333xx3333(math_mul3333xx3333(dLp_dS,dS_dFe),dFe_dLp))
        temp_9 = math_33to9(residuumLp)
        call dgesv(9,1,dRLp_dLp,9,devNull_9,temp_9,9,ierr)                                          ! solve dRLp/dLp * delta Lp = -res for delta Lp
        if (ierr /= 0) return ! error
        deltaLp = - math_9to33(temp_9)
      endif calculateJacobiLp

      Lpguess = Lpguess &
              + deltaLp * steplengthLp
    enddo LpLoop

    call phase_LiAndItsTangents(Li_constitutive, dLi_dS, dLi_dFi, &
                                       S, Fi_new, ph,en)

    !* update current residuum and check for convergence of loop
    atol_Li = max(num%rtol_crystalliteStress * max(norm2(Liguess),norm2(Li_constitutive)), &        ! absolute tolerance from largest acceptable relative error
                  num%atol_crystalliteStress)                                                       ! minimum lower cutoff
    residuumLi = Liguess - Li_constitutive
    if (any(IEEE_is_NaN(residuumLi))) then
      return ! error
    elseif (norm2(residuumLi) < atol_Li) then                                                       ! converged if below absolute tolerance
      exit LiLoop
    elseif (NiterationStressLi == 1 .or. norm2(residuumLi) < norm2(residuumLi_old)) then            ! not converged, but improved norm of residuum (always proceed in first iteration)...
      residuumLi_old = residuumLi                                                                   ! ...remember old values and...
      Liguess_old    = Liguess
      steplengthLi   = 1.0_pReal                                                                    ! ...proceed with normal step length (calculate new search direction)
    else                                                                                            ! not converged and residuum not improved...
      steplengthLi = num%subStepSizeLi * steplengthLi                                               ! ...try with smaller step length in same direction
      Liguess      = Liguess_old &
                   + deltaLi * steplengthLi
      cycle LiLoop
    endif

    calculateJacobiLi: if (mod(jacoCounterLi, num%iJacoLpresiduum) == 0) then
      jacoCounterLi = jacoCounterLi + 1

      temp_33 = matmul(matmul(A,B),invFi_current)
      do o=1,3; do p=1,3
        dFe_dLi(1:3,o,1:3,p) = -Delta_t*math_I3(o,p)*temp_33                                        ! dFe_dLp(i,j,k,l) = -Delta_t * A(i,k) invFi(l,j)
        dFi_dLi(1:3,o,1:3,p) = -Delta_t*math_I3(o,p)*invFi_current
      enddo; enddo
      do o=1,3; do p=1,3
        dFi_dLi(1:3,1:3,o,p) = matmul(matmul(Fi_new,dFi_dLi(1:3,1:3,o,p)),Fi_new)
      enddo; enddo
      dRLi_dLi  = math_eye(9) &
                - math_3333to99(math_mul3333xx3333(dLi_dS,  math_mul3333xx3333(dS_dFe, dFe_dLi) &
                                                          + math_mul3333xx3333(dS_dFi, dFi_dLi)))  &
                - math_3333to99(math_mul3333xx3333(dLi_dFi, dFi_dLi))
      temp_9 = math_33to9(residuumLi)
      call dgesv(9,1,dRLi_dLi,9,devNull_9,temp_9,9,ierr)                                            ! solve dRLi/dLp * delta Li = -res for delta Li
      if (ierr /= 0) return ! error
      deltaLi = - math_9to33(temp_9)
    endif calculateJacobiLi

    Liguess = Liguess &
            + deltaLi * steplengthLi
  enddo LiLoop

  invFp_new = matmul(invFp_current,B)
  call math_invert33(Fp_new,devNull,error,invFp_new)
  if (error) return ! error

  phase_mechanical_P(ph)%data(1:3,1:3,en)  = matmul(matmul(F,invFp_new),matmul(S,transpose(invFp_new)))
  phase_mechanical_S(ph)%data(1:3,1:3,en)  = S
  phase_mechanical_Lp(ph)%data(1:3,1:3,en) = Lpguess
  phase_mechanical_Li(ph)%data(1:3,1:3,en) = Liguess
  phase_mechanical_Fp(ph)%data(1:3,1:3,en) = Fp_new / math_det33(Fp_new)**(1.0_pReal/3.0_pReal)    ! regularize
  phase_mechanical_Fi(ph)%data(1:3,1:3,en) = Fi_new
  phase_mechanical_Fe(ph)%data(1:3,1:3,en) = matmul(matmul(F,invFp_new),invFi_new)
  broken = .false.

end function integrateStress


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with adaptive 1st order explicit Euler method
!> using Fixed Point Iteration to adapt the stepsize
!--------------------------------------------------------------------------------------------------
function integrateStateFPI(F_0,F,subFp0,subFi0,subState0,Delta_t,co,ip,el) result(broken)

  real(pReal), intent(in),dimension(3,3) :: F_0,F,subFp0,subFi0
  real(pReal), intent(in),dimension(:)   :: subState0
  real(pReal), intent(in) :: Delta_t
  integer, intent(in) :: &
    el, &                                                                                            !< element index in element loop
    ip, &                                                                                            !< integration point index in ip loop
    co                                                                                               !< grain index in grain loop
  logical :: &
    broken

  integer :: &
    NiterationState, &                                                                              !< number of iterations in state loop
    ph, &
    en, &
    sizeDotState
  real(pReal) :: &
    zeta
  real(pReal), dimension(phase_plasticity_maxSizeDotState) :: &
    r                                                                                               ! state residuum
  real(pReal), dimension(phase_plasticity_maxSizeDotState,2) :: &
    dotState


  ph = material_phaseID(co,(el-1)*discretization_nIPs + ip)
  en = material_phaseEntry(co,(el-1)*discretization_nIPs + ip)

  broken = plastic_dotState(Delta_t, co,ip,el,ph,en)
  if(broken) return

  sizeDotState = plasticState(ph)%sizeDotState
  plasticState(ph)%state(1:sizeDotState,en) = subState0 &
                                            + plasticState(ph)%dotState (1:sizeDotState,en) * Delta_t

  iteration: do NiterationState = 1, num%nState

    dotState(1:sizeDotState,2) = merge(dotState(1:sizeDotState,1),0.0, nIterationState > 1)
    dotState(1:sizeDotState,1) = plasticState(ph)%dotState(:,en)

    broken = integrateStress(F,subFp0,subFi0,Delta_t,co,ip,el)
    if(broken) exit iteration

    broken = plastic_dotState(Delta_t, co,ip,el,ph,en)
    if(broken) exit iteration

    zeta = damper(plasticState(ph)%dotState(:,en),dotState(1:sizeDotState,1),&
                                                  dotState(1:sizeDotState,2))
    plasticState(ph)%dotState(:,en) = plasticState(ph)%dotState(:,en) * zeta &
                                    + dotState(1:sizeDotState,1) * (1.0_pReal - zeta)
    r(1:sizeDotState) = plasticState(ph)%state(1:sizeDotState,en) &
                      - subState0 &
                      - plasticState(ph)%dotState(1:sizeDotState,en) * Delta_t
    plasticState(ph)%state(1:sizeDotState,en) = plasticState(ph)%state(1:sizeDotState,en) &
                                              - r(1:sizeDotState)
    if (converged(r(1:sizeDotState),plasticState(ph)%state(1:sizeDotState,en),plasticState(ph)%atol(1:sizeDotState))) then
      broken = plastic_deltaState(ph,en)
      exit iteration
    endif

  enddo iteration


  contains

  !--------------------------------------------------------------------------------------------------
  !> @brief calculate the damping for correction of state and dot state
  !--------------------------------------------------------------------------------------------------
  real(pReal) pure function damper(omega_0,omega_1,omega_2)

  real(pReal), dimension(:), intent(in) :: &
    omega_0, omega_1, omega_2

  real(pReal) :: dot_prod12, dot_prod22

  dot_prod12 = dot_product(omega_0-omega_1, omega_1-omega_2)
  dot_prod22 = dot_product(omega_1-omega_2, omega_1-omega_2)

  if (min(dot_product(omega_0,omega_1),dot_prod12) < 0.0_pReal .and. dot_prod22 > 0.0_pReal) then
    damper = 0.75_pReal + 0.25_pReal * tanh(2.0_pReal + 4.0_pReal * dot_prod12 / dot_prod22)
  else
    damper = 1.0_pReal
  endif

  end function damper

end function integrateStateFPI


!--------------------------------------------------------------------------------------------------
!> @brief integrate state with 1st order explicit Euler method
!--------------------------------------------------------------------------------------------------
function integrateStateEuler(F_0,F,subFp0,subFi0,subState0,Delta_t,co,ip,el) result(broken)

  real(pReal), intent(in),dimension(3,3) :: F_0,F,subFp0,subFi0
  real(pReal), intent(in),dimension(:)   :: subState0
  real(pReal), intent(in) :: Delta_t
  integer, intent(in) :: &
    el, &                                                                                            !< element index in element loop
    ip, &                                                                                            !< integration point index in ip loop
    co                                                                                               !< grain index in grain loop
  logical :: &
    broken

  integer :: &
    ph, &
    en, &
    sizeDotState


  ph = material_phaseID(co,(el-1)*discretization_nIPs + ip)
  en = material_phaseEntry(co,(el-1)*discretization_nIPs + ip)

  broken = plastic_dotState(Delta_t, co,ip,el,ph,en)
  if(broken) return

  sizeDotState = plasticState(ph)%sizeDotState
  plasticState(ph)%state(1:sizeDotState,en) = subState0 &
                                            + plasticState(ph)%dotState(1:sizeDotState,en) * Delta_t

  broken = plastic_deltaState(ph,en)
  if(broken) return

  broken = integrateStress(F,subFp0,subFi0,Delta_t,co,ip,el)

end function integrateStateEuler


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with 1st order Euler method with adaptive step size
!--------------------------------------------------------------------------------------------------
function integrateStateAdaptiveEuler(F_0,F,subFp0,subFi0,subState0,Delta_t,co,ip,el) result(broken)

  real(pReal), intent(in),dimension(3,3) :: F_0,F,subFp0,subFi0
  real(pReal), intent(in),dimension(:)   :: subState0
  real(pReal), intent(in) :: Delta_t
  integer, intent(in) :: &
    el, &                                                                                            !< element index in element loop
    ip, &                                                                                            !< integration point index in ip loop
    co                                                                                               !< grain index in grain loop
  logical :: &
    broken

  integer :: &
    ph, &
    en, &
    sizeDotState
  real(pReal), dimension(phase_plasticity_maxSizeDotState) :: residuum_plastic


  ph = material_phaseID(co,(el-1)*discretization_nIPs + ip)
  en = material_phaseEntry(co,(el-1)*discretization_nIPs + ip)

  broken = plastic_dotState(Delta_t, co,ip,el,ph,en)
  if(broken) return

  sizeDotState = plasticState(ph)%sizeDotState

  residuum_plastic(1:sizeDotState) = - plasticState(ph)%dotstate(1:sizeDotState,en) * 0.5_pReal * Delta_t
  plasticState(ph)%state(1:sizeDotState,en) = subState0 &
                                            + plasticState(ph)%dotstate(1:sizeDotState,en) * Delta_t

  broken = plastic_deltaState(ph,en)
  if(broken) return

  broken = integrateStress(F,subFp0,subFi0,Delta_t,co,ip,el)
  if(broken) return

  broken = plastic_dotState(Delta_t, co,ip,el,ph,en)
  if(broken) return

  broken = .not. converged(residuum_plastic(1:sizeDotState) + 0.5_pReal * plasticState(ph)%dotState(:,en) * Delta_t, &
                           plasticState(ph)%state(1:sizeDotState,en), &
                           plasticState(ph)%atol(1:sizeDotState))

end function integrateStateAdaptiveEuler


!---------------------------------------------------------------------------------------------------
!> @brief Integrate state (including stress integration) with the classic Runge Kutta method
!---------------------------------------------------------------------------------------------------
function integrateStateRK4(F_0,F,subFp0,subFi0,subState0,Delta_t,co,ip,el) result(broken)

  real(pReal), intent(in),dimension(3,3) :: F_0,F,subFp0,subFi0
  real(pReal), intent(in),dimension(:)   :: subState0
  real(pReal), intent(in) :: Delta_t
  integer, intent(in) :: co,ip,el
  logical :: broken

  real(pReal), dimension(3,3), parameter :: &
    A = reshape([&
      0.5_pReal, 0.0_pReal, 0.0_pReal, &
      0.0_pReal, 0.5_pReal, 0.0_pReal, &
      0.0_pReal, 0.0_pReal, 1.0_pReal],&
      shape(A))
  real(pReal), dimension(3), parameter :: &
    C = [0.5_pReal, 0.5_pReal, 1.0_pReal]
  real(pReal), dimension(4), parameter :: &
    B = [1.0_pReal/6.0_pReal, 1.0_pReal/3.0_pReal, 1.0_pReal/3.0_pReal, 1.0_pReal/6.0_pReal]


  broken = integrateStateRK(F_0,F,subFp0,subFi0,subState0,Delta_t,co,ip,el,A,B,C)

end function integrateStateRK4


!---------------------------------------------------------------------------------------------------
!> @brief Integrate state (including stress integration) with the Cash-Carp method
!---------------------------------------------------------------------------------------------------
function integrateStateRKCK45(F_0,F,subFp0,subFi0,subState0,Delta_t,co,ip,el) result(broken)

  real(pReal), intent(in),dimension(3,3) :: F_0,F,subFp0,subFi0
  real(pReal), intent(in),dimension(:)   :: subState0
  real(pReal), intent(in) :: Delta_t
  integer, intent(in) :: co,ip,el
  logical :: broken

  real(pReal), dimension(5,5), parameter :: &
    A = reshape([&
      1._pReal/5._pReal,       .0_pReal,             .0_pReal,               .0_pReal,                  .0_pReal, &
      3._pReal/40._pReal,      9._pReal/40._pReal,   .0_pReal,               .0_pReal,                  .0_pReal, &
      3_pReal/10._pReal,       -9._pReal/10._pReal,  6._pReal/5._pReal,      .0_pReal,                  .0_pReal, &
      -11._pReal/54._pReal,    5._pReal/2._pReal,    -70.0_pReal/27.0_pReal, 35.0_pReal/27.0_pReal,     .0_pReal, &
      1631._pReal/55296._pReal,175._pReal/512._pReal,575._pReal/13824._pReal,44275._pReal/110592._pReal,253._pReal/4096._pReal],&
      shape(A))
  real(pReal), dimension(5), parameter :: &
    C = [0.2_pReal, 0.3_pReal, 0.6_pReal, 1.0_pReal, 0.875_pReal]
  real(pReal), dimension(6), parameter :: &
    B = &
      [37.0_pReal/378.0_pReal, .0_pReal, 250.0_pReal/621.0_pReal, &
      125.0_pReal/594.0_pReal, .0_pReal, 512.0_pReal/1771.0_pReal], &
    DB = B - &
      [2825.0_pReal/27648.0_pReal,    .0_pReal,                18575.0_pReal/48384.0_pReal,&
      13525.0_pReal/55296.0_pReal, 277.0_pReal/14336.0_pReal,  1._pReal/4._pReal]


  broken = integrateStateRK(F_0,F,subFp0,subFi0,subState0,Delta_t,co,ip,el,A,B,C,DB)

end function integrateStateRKCK45


!--------------------------------------------------------------------------------------------------
!> @brief Integrate state (including stress integration) with an explicit Runge-Kutta method or an
!! embedded explicit Runge-Kutta method
!--------------------------------------------------------------------------------------------------
function integrateStateRK(F_0,F,subFp0,subFi0,subState0,Delta_t,co,ip,el,A,B,C,DB) result(broken)

  real(pReal), intent(in),dimension(3,3) :: F_0,F,subFp0,subFi0
  real(pReal), intent(in),dimension(:)   :: subState0
  real(pReal), intent(in) :: Delta_t
  real(pReal), dimension(:,:), intent(in) :: A
  real(pReal), dimension(:),   intent(in) :: B, C
  real(pReal), dimension(:),   intent(in), optional :: DB
  integer, intent(in) :: &
    el, &                                                                                            !< element index in element loop
    ip, &                                                                                            !< integration point index in ip loop
    co                                                                                               !< grain index in grain loop
  logical :: broken

  integer :: &
    stage, &                                                                                        ! stage index in integration stage loop
    n, &
    ph, &
    en, &
    sizeDotState
  real(pReal), dimension(phase_plasticity_maxSizeDotState,size(B)) :: plastic_RKdotState


  ph = material_phaseID(co,(el-1)*discretization_nIPs + ip)
  en = material_phaseEntry(co,(el-1)*discretization_nIPs + ip)

  broken = plastic_dotState(Delta_t,co,ip,el,ph,en)
  if(broken) return

  sizeDotState = plasticState(ph)%sizeDotState

  do stage = 1, size(A,1)

    plastic_RKdotState(1:sizeDotState,stage) = plasticState(ph)%dotState(:,en)
    plasticState(ph)%dotState(:,en) = A(1,stage) * plastic_RKdotState(1:sizeDotState,1)

    do n = 2, stage
      plasticState(ph)%dotState(:,en) = plasticState(ph)%dotState(:,en) &
                                      + A(n,stage) * plastic_RKdotState(1:sizeDotState,n)
    enddo

    plasticState(ph)%state(1:sizeDotState,en) = subState0 &
                                              + plasticState(ph)%dotState (1:sizeDotState,en) * Delta_t

    broken = integrateStress(F_0 + (F - F_0) * Delta_t * C(stage),subFp0,subFi0,Delta_t * C(stage),co,ip,el)
    if(broken) exit

    broken = plastic_dotState(Delta_t*C(stage),co,ip,el,ph,en)
    if(broken) exit

  enddo
  if(broken) return


  plastic_RKdotState(1:sizeDotState,size(B)) = plasticState (ph)%dotState(:,en)
  plasticState(ph)%dotState(:,en) = matmul(plastic_RKdotState(1:sizeDotState,1:size(B)),B)
  plasticState(ph)%state(1:sizeDotState,en) = subState0 &
                                            + plasticState(ph)%dotState (1:sizeDotState,en) * Delta_t

  if(present(DB)) &
    broken = .not. converged(matmul(plastic_RKdotState(1:sizeDotState,1:size(DB)),DB) * Delta_t, &
                             plasticState(ph)%state(1:sizeDotState,en), &
                             plasticState(ph)%atol(1:sizeDotState))

  if(broken) return

  broken = plastic_deltaState(ph,en)
  if(broken) return

  broken = integrateStress(F,subFp0,subFi0,Delta_t,co,ip,el)

end function integrateStateRK


!--------------------------------------------------------------------------------------------------
!> @brief writes crystallite results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine crystallite_results(group,ph)

  character(len=*), intent(in) :: group
  integer,          intent(in) :: ph

  integer :: ou
  real(pReal), allocatable, dimension(:,:)   :: selected_rotations


  call results_closeGroup(results_addGroup(group//'/mechanical'))

  do ou = 1, size(output_constituent(ph)%label)

    select case (output_constituent(ph)%label(ou))
      case('F')
        call results_writeDataset(phase_mechanical_F(ph)%data,group//'/mechanical/','F',&
                                 'deformation gradient','1')
      case('F_e')
        call results_writeDataset(phase_mechanical_Fe(ph)%data,group//'/mechanical/','F_e',&
                                 'elastic deformation gradient','1')
      case('F_p')
        call results_writeDataset(phase_mechanical_Fp(ph)%data,group//'/mechanical/','F_p', &
                                 'plastic deformation gradient','1')
      case('F_i')
        call results_writeDataset(phase_mechanical_Fi(ph)%data,group//'/mechanical/','F_i', &
                                 'inelastic deformation gradient','1')
      case('L_p')
        call results_writeDataset(phase_mechanical_Lp(ph)%data,group//'/mechanical/','L_p', &
                                 'plastic velocity gradient','1/s')
      case('L_i')
        call results_writeDataset(phase_mechanical_Li(ph)%data,group//'/mechanical/','L_i', &
                                 'inelastic velocity gradient','1/s')
      case('P')
        call results_writeDataset(phase_mechanical_P(ph)%data,group//'/mechanical/','P', &
                                 'first Piola-Kirchhoff stress','Pa')
      case('S')
        call results_writeDataset(phase_mechanical_S(ph)%data,group//'/mechanical/','S', &
                                 'second Piola-Kirchhoff stress','Pa')
      case('O')
        selected_rotations = select_rotations(phase_orientation(ph)%data)
        call results_writeDataset(selected_rotations,group//'/mechanical',output_constituent(ph)%label(ou),&
                                 'crystal orientation as quaternion','q_0 (q_1 q_2 q_3)')
        call results_addAttribute('lattice',phase_lattice(ph),group//'/mechanical/'//output_constituent(ph)%label(ou))
        if (any(phase_lattice(ph) == ['hP', 'tI'])) &
          call results_addAttribute('c/a',phase_cOverA(ph),group//'/mechanical/'//output_constituent(ph)%label(ou))
    end select
  enddo


  contains

!--------------------------------------------------------------------------------------------------
!> @brief Convert orientation for output: ToDo: implement in HDF5/results
!--------------------------------------------------------------------------------------------------
  function select_rotations(dataset)

    type(rotation), dimension(:), intent(in) :: dataset
    real(pReal), dimension(4,size(dataset,1)) :: select_rotations
    integer :: en

    do en = 1, size(dataset,1)
      select_rotations(:,en) = dataset(en)%asQuaternion()
   enddo

 end function select_rotations

end subroutine crystallite_results


!--------------------------------------------------------------------------------------------------
!> @brief Forward data after successful increment.
! ToDo: Any guessing for the current states possible?
!--------------------------------------------------------------------------------------------------
module subroutine mechanical_forward()

  integer :: ph


  do ph = 1, size(plasticState)
    phase_mechanical_Fi0(ph) = phase_mechanical_Fi(ph)
    phase_mechanical_Fp0(ph) = phase_mechanical_Fp(ph)
    phase_mechanical_F0(ph)  = phase_mechanical_F(ph)
    phase_mechanical_Li0(ph) = phase_mechanical_Li(ph)
    phase_mechanical_Lp0(ph) = phase_mechanical_Lp(ph)
    phase_mechanical_S0(ph)  = phase_mechanical_S(ph)
    plasticState(ph)%state0 = plasticState(ph)%state
  enddo

end subroutine mechanical_forward


!--------------------------------------------------------------------------------------------------
!> @brief calculate stress (P)
!--------------------------------------------------------------------------------------------------
module function crystallite_stress(dt,co,ip,el) result(converged_)

  real(pReal), intent(in) :: dt
  integer, intent(in) :: &
    co, &
    ip, &
    el
  logical :: converged_

  real(pReal) :: &
    formerSubStep
  integer :: &
    ph, en, sizeDotState
  logical :: todo
  real(pReal) :: subFrac,subStep
  real(pReal), dimension(3,3) :: &
    subFp0, &
    subFi0, &
    subLp0, &
    subLi0, &
    subF0, &
    subF
  real(pReal), dimension(:), allocatable :: subState0


  ph = material_phaseID(co,(el-1)*discretization_nIPs + ip)
  en = material_phaseEntry(co,(el-1)*discretization_nIPs + ip)
  sizeDotState = plasticState(ph)%sizeDotState

  subLi0 = phase_mechanical_Li0(ph)%data(1:3,1:3,en)
  subLp0 = phase_mechanical_Lp0(ph)%data(1:3,1:3,en)
  subState0 = plasticState(ph)%State0(:,en)

  if (damageState(ph)%sizeState > 0) &
    damageState(ph)%subState0(:,en) = damageState(ph)%state0(:,en)

  subFp0 = phase_mechanical_Fp0(ph)%data(1:3,1:3,en)
  subFi0 = phase_mechanical_Fi0(ph)%data(1:3,1:3,en)
  subF0  = phase_mechanical_F0(ph)%data(1:3,1:3,en)
  subFrac = 0.0_pReal
  subStep = 1.0_pReal/num%subStepSizeCryst
  todo = .true.
  converged_ = .false.                                                      ! pretend failed step of 1/subStepSizeCryst

  todo = .true.
  cutbackLooping: do while (todo)

    if (converged_) then
      formerSubStep = subStep
      subFrac = subFrac + subStep
      subStep = min(1.0_pReal - subFrac, num%stepIncreaseCryst * subStep)

      todo = subStep > 0.0_pReal                        ! still time left to integrate on?

      if (todo) then
        subF0  = subF
        subLp0 = phase_mechanical_Lp(ph)%data(1:3,1:3,en)
        subLi0 = phase_mechanical_Li(ph)%data(1:3,1:3,en)
        subFp0 = phase_mechanical_Fp(ph)%data(1:3,1:3,en)
        subFi0 = phase_mechanical_Fi(ph)%data(1:3,1:3,en)
        subState0 = plasticState(ph)%state(:,en)
        if (damageState(ph)%sizeState > 0) &
          damageState(ph)%subState0(:,en) = damageState(ph)%state(:,en)

      endif
!--------------------------------------------------------------------------------------------------
!  cut back (reduced time and restore)
    else
      subStep       = num%subStepSizeCryst * subStep
      phase_mechanical_Fp(ph)%data(1:3,1:3,en) = subFp0
      phase_mechanical_Fi(ph)%data(1:3,1:3,en) = subFi0
      phase_mechanical_S(ph)%data(1:3,1:3,en) = phase_mechanical_S0(ph)%data(1:3,1:3,en)          ! why no subS0 ? is S0 of any use?
      if (subStep < 1.0_pReal) then                                                                 ! actual (not initial) cutback
        phase_mechanical_Lp(ph)%data(1:3,1:3,en) = subLp0
        phase_mechanical_Li(ph)%data(1:3,1:3,en) = subLi0
      endif
      plasticState(ph)%state(:,en) = subState0
      if (damageState(ph)%sizeState > 0) &
        damageState(ph)%state(:,en) = damageState(ph)%subState0(:,en)

      todo = subStep > num%subStepMinCryst                          ! still on track or already done (beyond repair)
    endif

!--------------------------------------------------------------------------------------------------
!  prepare for integration
    if (todo) then
      subF = subF0 &
           + subStep * (phase_mechanical_F(ph)%data(1:3,1:3,en) - phase_mechanical_F0(ph)%data(1:3,1:3,en))
      converged_ = .not. integrateState(subF0,subF,subFp0,subFi0,subState0(1:sizeDotState),subStep * dt,co,ip,el)
      converged_ = converged_ .and. .not. integrateDamageState(subStep * dt,co,(el-1)*discretization_nIPs + ip)
    endif

  enddo cutbackLooping

end function crystallite_stress


!--------------------------------------------------------------------------------------------------
!> @brief Restore data after homog cutback.
!--------------------------------------------------------------------------------------------------
module subroutine mechanical_restore(ce,includeL)

  integer, intent(in) :: ce
  logical, intent(in) :: &
    includeL                                                                                        !< protect agains fake cutback

  integer :: &
    co, ph, en


  do co = 1,homogenization_Nconstituents(material_homogenizationID(ce))
    ph = material_phaseID(co,ce)
    en = material_phaseEntry(co,ce)
    if (includeL) then
      phase_mechanical_Lp(ph)%data(1:3,1:3,en) = phase_mechanical_Lp0(ph)%data(1:3,1:3,en)
      phase_mechanical_Li(ph)%data(1:3,1:3,en) = phase_mechanical_Li0(ph)%data(1:3,1:3,en)
    endif                                                                                           ! maybe protecting everything from overwriting makes more sense

    phase_mechanical_Fp(ph)%data(1:3,1:3,en)   = phase_mechanical_Fp0(ph)%data(1:3,1:3,en)
    phase_mechanical_Fi(ph)%data(1:3,1:3,en)   = phase_mechanical_Fi0(ph)%data(1:3,1:3,en)
    phase_mechanical_S(ph)%data(1:3,1:3,en)    = phase_mechanical_S0(ph)%data(1:3,1:3,en)

    plasticState(ph)%state(:,en) = plasticState(ph)%State0(:,en)
  enddo

end subroutine mechanical_restore

!--------------------------------------------------------------------------------------------------
!> @brief Calculate tangent (dPdF).
!--------------------------------------------------------------------------------------------------
module function phase_mechanical_dPdF(dt,co,ce) result(dPdF)

  real(pReal), intent(in) :: dt
  integer, intent(in) :: &
    co, &                                                                                            !< counter in constituent loop
    ce
  real(pReal), dimension(3,3,3,3) :: dPdF

  integer :: &
    o, &
    p, ph, en
  real(pReal), dimension(3,3)     ::   devNull, &
                                       invSubFp0,invSubFi0,invFp,invFi, &
                                       temp_33_1, temp_33_2, temp_33_3
  real(pReal), dimension(3,3,3,3) ::   dSdFe, &
                                       dSdF, &
                                       dSdFi, &
                                       dLidS, &                                                     ! tangent in lattice configuration
                                       dLidFi, &
                                       dLpdS, &
                                       dLpdFi, &
                                       dFidS, &
                                       dFpinvdF, &
                                       rhs_3333, &
                                       lhs_3333, &
                                       temp_3333
  real(pReal), dimension(9,9)::        temp_99
  logical :: error


  ph = material_phaseID(co,ce)
  en = material_phaseEntry(co,ce)

  call phase_hooke_SandItsTangents(devNull,dSdFe,dSdFi, &
                                          phase_mechanical_Fe(ph)%data(1:3,1:3,en), &
                                          phase_mechanical_Fi(ph)%data(1:3,1:3,en),ph,en)
  call phase_LiAndItsTangents(devNull,dLidS,dLidFi, &
                                     phase_mechanical_S(ph)%data(1:3,1:3,en), &
                                     phase_mechanical_Fi(ph)%data(1:3,1:3,en), &
                                     ph,en)

  invFp = math_inv33(phase_mechanical_Fp(ph)%data(1:3,1:3,en))
  invFi = math_inv33(phase_mechanical_Fi(ph)%data(1:3,1:3,en))
  invSubFp0 = math_inv33(phase_mechanical_Fp0(ph)%data(1:3,1:3,en))
  invSubFi0 = math_inv33(phase_mechanical_Fi0(ph)%data(1:3,1:3,en))

  if (sum(abs(dLidS)) < tol_math_check) then
    dFidS = 0.0_pReal
  else
    lhs_3333 = 0.0_pReal; rhs_3333 = 0.0_pReal
    do o=1,3; do p=1,3
      lhs_3333(1:3,1:3,o,p) = lhs_3333(1:3,1:3,o,p) &
                            + matmul(invSubFi0,dLidFi(1:3,1:3,o,p)) * dt
      lhs_3333(1:3,o,1:3,p) = lhs_3333(1:3,o,1:3,p) &
                            + invFi*invFi(p,o)
      rhs_3333(1:3,1:3,o,p) = rhs_3333(1:3,1:3,o,p) &
                            - matmul(invSubFi0,dLidS(1:3,1:3,o,p)) * dt
    enddo; enddo
    call math_invert(temp_99,error,math_3333to99(lhs_3333))
    if (error) then
      call IO_warning(warning_ID=600, &
                      ext_msg='inversion error in analytic tangent calculation')
      dFidS = 0.0_pReal
    else
      dFidS = math_mul3333xx3333(math_99to3333(temp_99),rhs_3333)
    endif
    dLidS = math_mul3333xx3333(dLidFi,dFidS) + dLidS
  endif

  call plastic_LpAndItsTangents(devNull,dLpdS,dLpdFi, &
                                             phase_mechanical_S(ph)%data(1:3,1:3,en), &
                                             phase_mechanical_Fi(ph)%data(1:3,1:3,en),ph,en)
  dLpdS = math_mul3333xx3333(dLpdFi,dFidS) + dLpdS

!--------------------------------------------------------------------------------------------------
! calculate dSdF
  temp_33_1 = transpose(matmul(invFp,invFi))
  temp_33_2 = matmul(phase_mechanical_F(ph)%data(1:3,1:3,en),invSubFp0)
  temp_33_3 = matmul(matmul(phase_mechanical_F(ph)%data(1:3,1:3,en),invFp), invSubFi0)

  do o=1,3; do p=1,3
    rhs_3333(p,o,1:3,1:3)  = matmul(dSdFe(p,o,1:3,1:3),temp_33_1)
    temp_3333(1:3,1:3,p,o) = matmul(matmul(temp_33_2,dLpdS(1:3,1:3,p,o)), invFi) &
                           + matmul(temp_33_3,dLidS(1:3,1:3,p,o))
  enddo; enddo
  lhs_3333 = math_mul3333xx3333(dSdFe,temp_3333) * dt &
           + math_mul3333xx3333(dSdFi,dFidS)

  call math_invert(temp_99,error,math_eye(9)+math_3333to99(lhs_3333))
  if (error) then
    call IO_warning(warning_ID=600, &
                    ext_msg='inversion error in analytic tangent calculation')
    dSdF = rhs_3333
  else
    dSdF = math_mul3333xx3333(math_99to3333(temp_99),rhs_3333)
  endif

!--------------------------------------------------------------------------------------------------
! calculate dFpinvdF
  temp_3333 = math_mul3333xx3333(dLpdS,dSdF)
  do o=1,3; do p=1,3
    dFpinvdF(1:3,1:3,p,o) = - matmul(invSubFp0, matmul(temp_3333(1:3,1:3,p,o),invFi)) * dt
  enddo; enddo

!--------------------------------------------------------------------------------------------------
! assemble dPdF
  temp_33_1 = matmul(phase_mechanical_S(ph)%data(1:3,1:3,en),transpose(invFp))
  temp_33_2 = matmul(phase_mechanical_F(ph)%data(1:3,1:3,en),invFp)
  temp_33_3 = matmul(temp_33_2,phase_mechanical_S(ph)%data(1:3,1:3,en))

  dPdF = 0.0_pReal
  do p=1,3
    dPdF(p,1:3,p,1:3) = transpose(matmul(invFp,temp_33_1))
  enddo
  do o=1,3; do p=1,3
    dPdF(1:3,1:3,p,o) = dPdF(1:3,1:3,p,o) &
                      + matmul(matmul(phase_mechanical_F(ph)%data(1:3,1:3,en),dFpinvdF(1:3,1:3,p,o)),temp_33_1) &
                      + matmul(matmul(temp_33_2,dSdF(1:3,1:3,p,o)),transpose(invFp)) &
                      + matmul(temp_33_3,transpose(dFpinvdF(1:3,1:3,p,o)))
  enddo; enddo

end function phase_mechanical_dPdF


module subroutine mechanical_restartWrite(groupHandle,ph)

  integer(HID_T), intent(in) :: groupHandle
  integer, intent(in) :: ph


  call HDF5_write(plasticState(ph)%state,groupHandle,'omega')
  call HDF5_write(phase_mechanical_Fi(ph)%data,groupHandle,'F_i')
  call HDF5_write(phase_mechanical_Li(ph)%data,groupHandle,'L_i')
  call HDF5_write(phase_mechanical_Lp(ph)%data,groupHandle,'L_p')
  call HDF5_write(phase_mechanical_Fp(ph)%data,groupHandle,'F_p')
  call HDF5_write(phase_mechanical_S(ph)%data,groupHandle,'S')
  call HDF5_write(phase_mechanical_F(ph)%data,groupHandle,'F')

end subroutine mechanical_restartWrite


module subroutine mechanical_restartRead(groupHandle,ph)

  integer(HID_T), intent(in) :: groupHandle
  integer, intent(in) :: ph


  call HDF5_read(plasticState(ph)%state0,groupHandle,'omega')
  call HDF5_read(phase_mechanical_Fi0(ph)%data,groupHandle,'F_i')
  call HDF5_read(phase_mechanical_Li0(ph)%data,groupHandle,'L_i')
  call HDF5_read(phase_mechanical_Lp0(ph)%data,groupHandle,'L_p')
  call HDF5_read(phase_mechanical_Fp0(ph)%data,groupHandle,'F_p')
  call HDF5_read(phase_mechanical_S0(ph)%data,groupHandle,'S')
  call HDF5_read(phase_mechanical_F0(ph)%data,groupHandle,'F')

end subroutine mechanical_restartRead


!----------------------------------------------------------------------------------------------
!< @brief Get first Piola-Kichhoff stress (for use by non-mech physics)
!----------------------------------------------------------------------------------------------
module function mechanical_S(ph,en) result(S)

  integer, intent(in) :: ph,en
  real(pReal), dimension(3,3) :: S


  S = phase_mechanical_S(ph)%data(1:3,1:3,en)

end function mechanical_S


!----------------------------------------------------------------------------------------------
!< @brief Get plastic velocity gradient (for use by non-mech physics)
!----------------------------------------------------------------------------------------------
module function mechanical_L_p(ph,en) result(L_p)

  integer, intent(in) :: ph,en
  real(pReal), dimension(3,3) :: L_p


  L_p = phase_mechanical_Lp(ph)%data(1:3,1:3,en)

end function mechanical_L_p


!----------------------------------------------------------------------------------------------
!< @brief Get elastic deformation gradient (for use by non-mech physics)
!----------------------------------------------------------------------------------------------
module function mechanical_F_e(ph,en) result(F_e)

  integer, intent(in) :: ph,en
  real(pReal), dimension(3,3) :: F_e


  F_e = phase_mechanical_Fe(ph)%data(1:3,1:3,en)

end function mechanical_F_e


!----------------------------------------------------------------------------------------------
!< @brief Get second Piola-Kichhoff stress (for use by homogenization)
!----------------------------------------------------------------------------------------------
module function phase_P(co,ce) result(P)

  integer, intent(in) :: co, ce
  real(pReal), dimension(3,3) :: P


  P = phase_mechanical_P(material_phaseID(co,ce))%data(1:3,1:3,material_phaseEntry(co,ce))

end function phase_P


!----------------------------------------------------------------------------------------------
!< @brief Get deformation gradient (for use by homogenization)
!----------------------------------------------------------------------------------------------
module function phase_F(co,ce) result(F)

  integer, intent(in) :: co, ce
  real(pReal), dimension(3,3) :: F


  F = phase_mechanical_F(material_phaseID(co,ce))%data(1:3,1:3,material_phaseEntry(co,ce))

end function phase_F


!----------------------------------------------------------------------------------------------
!< @brief Set deformation gradient (for use by homogenization)
!----------------------------------------------------------------------------------------------
module subroutine phase_set_F(F,co,ce)

  real(pReal), dimension(3,3), intent(in) :: F
  integer, intent(in) :: co, ce


  phase_mechanical_F(material_phaseID(co,ce))%data(1:3,1:3,material_phaseEntry(co,ce)) = F

end subroutine phase_set_F


end submodule mechanical
