!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, University of Manchester
!> @author Madeleine Bignon, University of Manchester
!> @brief  Powerlaw crystal plasticity coupled KWN formulation
!--------------------------------------------------------------------------------------------------
submodule(phase:plastic) kwnpowerlaw

  real(pReal), parameter :: &
    lnorm  = 1.0e-10_pReal, &
    kBnorm = 1.38064852e-3_pReal, &
    Rnorm  = 8.3145e-10_pReal

  type :: tParameters
    real(pReal) :: &
      dot_gamma_0_sl = 1.0_pReal, &                                                                 !< reference shear strain rate for slip
      dot_gamma_0_tw = 1.0_pReal, &                                                                 !< reference shear strain rate for twin
      n_sl           = 1.0_pReal, &                                                                 !< stress exponent for slip
      n_tw           = 1.0_pReal, &                                                                 !< stress exponent for twin
      f_sat_sl_tw    = 1.0_pReal, &                                                                 !< push-up factor for slip saturation due to twinning
      c_1            = 1.0_pReal, &
      c_2            = 1.0_pReal, &
      c_3            = 1.0_pReal, &
      c_4            = 1.0_pReal, &
      h_0_sl_sl      = 1.0_pReal, &                                                                 !< reference hardening slip - slip
      h_0_tw_sl      = 1.0_pReal, &                                                                 !< reference hardening twin - slip
      h_0_tw_tw      = 1.0_pReal, &                                                                 !< reference hardening twin - twin
      a_sl           = 1.0_pReal
    real(pReal),               allocatable, dimension(:) :: &
      xi_inf_sl, &                                                                                  !< maximum critical shear stress for slip
      h_int, &                                                                                      !< per family hardening activity (optional)
      gamma_char                                                                                    !< characteristic shear for twins
    real(pReal),               allocatable, dimension(:,:) :: &
      h_sl_sl, &                                                                                    !< slip resistance from slip activity
      h_sl_tw, &                                                                                    !< slip resistance from twin activity
      h_tw_sl, &                                                                                    !< twin resistance from slip activity
      h_tw_tw                                                                                       !< twin resistance from twin activity
    real(pReal),               allocatable, dimension(:,:,:) :: &
      P_sl, &
      P_tw, &
      nonSchmid_pos, &
      nonSchmid_neg
    integer :: &
      sum_N_sl, &                                                                                   !< total number of active slip system
      sum_N_tw                                                                                      !< total number of active twin systems
    logical :: &
      nonSchmidActive = .false.
    integer :: &
      kwn_nSteps              ! discretization in r-space
    real(pReal) :: &
      kwn_stepsize, &         ! discretization in r-space
      kwn_step0               ! minimum radius
    real(pReal) :: &
      lattice_param, &        ! lattice parameter in Angstrom
      atomic_volume, &        ! atomic volume in Angstrom^3
      molar_volume, &         ! molar volume in m^3
      misfit_energy, &        ! normalized precipitate misfit energy in J/Angstrom/m^3
      gamma_coherent, &       ! coherent precipitate surface energy in J/m^2
      vacancy_generation, &   ! vacancy generation rate coefficient
      vacancy_sink_spacing, & ! vacancy sink spacing
      vacancy_energy, &       ! normalized vacancy formation energy (Q/kB) in 1/k
      migration_energy, &     ! normalized solute migration energy (Q/kB) in 1/k
      diffusion0, &           ! solute diffusivity in Angstrom^2/s
      vacancy_migration_energy, & ! normalized solute migration energy (Q/kB) in 1/k
      vacancy_diffusion0, &   ! vacancy diffusivity in Angstrom^2/s 
      c0_matrix(2), &            ! initial matrix solute composition in mol fraction
      ceq_matrix(2), &           ! equilibrium matrix composition at flat interface in mol fraction
      ceq_precipitate(2), &         ! equilibrium precipitate composition in mol fraction
      mean_radius_initial,&  ! average radius initial distribution in meters
	  standard_deviation,& ! standard deviation initial distribution (log normal law assumed)
	  volume_fraction_initial, &! initial 
	  stoechiometry(3), &
	  saturation_dislocation_density, &
	  dislocation_arrangement, &
	  jog_formation_energy, &
	  q_dislocation ! normalized pipe diffusion migration energy (Q/kB) in 1/k
    real(pReal) :: &
      shear_modulus, &        ! matrix shear modulus in Pa
      burgers_vec, &          ! Burgers vector in Angstrom                                                                     
      solute_strength , &        ! solute strengthening coefficient
      precipitate_strength ! precipitate strenghthening coefficient
    real(pReal), dimension(:),   allocatable :: &
      bins                    ! size bins, &
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type :: tKwnpowerlawState
    real(pReal), pointer, dimension(:,:) :: &
      xi_slip, &
      xi_twin, &
      gamma_slip, &
      gamma_twin, &
      precipitate_density
    real(pReal), pointer, dimension(  :) :: &
      c_vacancy, &
      time

  end type tKwnpowerlawState

 type :: tKwnpowerlawMicrostructure
   real(pReal),                  dimension(:),   allocatable :: &
     total_precipitate_density, &                                                             
     avg_precipitate_radius, &
     precipitate_volume_frac
  real(pReal),                  dimension(:,:),   allocatable :: &
     c_matrix, &
     interface_concentration
 end type tKwnpowerlawMicrostructure

!--------------------------------------------------------------------------------------------------
! containers for parameters and state
  type(tParameters),                allocatable, dimension(:) :: param
  type(tKwnpowerlawState),          allocatable, dimension(:) :: &
    dotState, &
    state
  type(tKwnpowerlawMicrostructure), allocatable, dimension(:) :: dependentState

contains


!--------------------------------------------------------------------------------------------------
!> @brief Perform module initialization.
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function plastic_kwnpowerlaw_init() result(myPlasticity)

  logical, dimension(:), allocatable :: myPlasticity
  integer :: &
    ph, i, &
    Nmembers, &
    sizeState, sizeDotState, &
    startIndex, endIndex
  integer,     dimension(:), allocatable :: &
    N_sl, N_tw
  real(pReal), dimension(:), allocatable :: &
    xi_0_sl, &                                                                                      !< initial critical shear stress for slip
    xi_0_tw, &                                                                                      !< initial critical shear stress for twin
    a   
	                                                                              !< non-Schmid coefficients

  character(len=pStringLen) :: &
    extmsg = ''
  class(tNode), pointer :: &
    phases, &
    phase, &
    mech, &
    pl

  real(pReal) :: &
    radiusL, radiusR, radiusC, N_0, shape_factor, x
  real(pReal), dimension(:), allocatable :: &
  	normalized_distribution_function
  real(pReal) :: &
    total_precipitate_density, avg_precipitate_radius, precipitate_volume_frac, c_matrix(2)

  
  myPlasticity = plastic_active('kwnpowerlaw')
  if(count(myPlasticity) == 0) return

  print'(/,a)', ' <<<+-  phase:mechanical:plastic:kwnpowerlaw init  -+>>>'
  print'(a,i0)', ' # phases: ',count(myPlasticity); flush(IO_STDOUT)


  phases => config_material%get('phase')
  allocate(param(phases%length))
  allocate(state(phases%length))
  allocate(dotState(phases%length))
  allocate(dependentState(phases%length))

  do ph = 1, phases%length
    if(.not. myPlasticity(ph)) cycle

    associate(prm => param(ph), dot => dotState(ph), stt => state(ph), dst => dependentState(ph))

    phase => phases%get(ph)
    mech  => phase%get('mechanical')
    pl  => mech%get('plastic')

!--------------------------------------------------------------------------------------------------
! slip related parameters
    N_sl         = pl%get_as1dInt('N_sl',defaultVal=emptyIntArray)
    prm%sum_N_sl = sum(abs(N_sl))
    slipActive: if (prm%sum_N_sl > 0) then
      prm%P_sl = lattice_SchmidMatrix_slip(N_sl,phase%get_asString('lattice'),&
                                           phase%get_asFloat('c/a',defaultVal=0.0_pReal))

      if(phase%get_asString('lattice') == 'cI') then
        a = pl%get_as1dFloat('a_nonSchmid',defaultVal=emptyRealArray)
        if(size(a) > 0) prm%nonSchmidActive = .true.
        prm%nonSchmid_pos  = lattice_nonSchmidMatrix(N_sl,a,+1)
        prm%nonSchmid_neg  = lattice_nonSchmidMatrix(N_sl,a,-1)
      else
        prm%nonSchmid_pos  = prm%P_sl
        prm%nonSchmid_neg  = prm%P_sl
      endif
      prm%h_sl_sl   = lattice_interaction_SlipBySlip(N_sl, &
                                                     pl%get_as1dFloat('h_sl_sl'), &
                                                     phase%get_asString('lattice'))

      xi_0_sl             = pl%get_as1dFloat('xi_0_sl',   requiredSize=size(N_sl))
      prm%xi_inf_sl       = pl%get_as1dFloat('xi_inf_sl', requiredSize=size(N_sl))
      prm%h_int           = pl%get_as1dFloat('h_int',     requiredSize=size(N_sl), &
                                            defaultVal=[(0.0_pReal,i=1,size(N_sl))])

      prm%dot_gamma_0_sl  = pl%get_asFloat('dot_gamma_0_sl')
      prm%n_sl            = pl%get_asFloat('n_sl')
      prm%a_sl            = pl%get_asFloat('a_sl')
      prm%h_0_sl_sl       = pl%get_asFloat('h_0_sl_sl')

      ! expand: family => system
      xi_0_sl             = math_expand(xi_0_sl,      N_sl)
      prm%xi_inf_sl       = math_expand(prm%xi_inf_sl,N_sl)
      prm%h_int           = math_expand(prm%h_int,    N_sl)

      ! sanity checks
      if (    prm%dot_gamma_0_sl  <= 0.0_pReal)      extmsg = trim(extmsg)//' dot_gamma_0_sl'
      if (    prm%a_sl            <= 0.0_pReal)      extmsg = trim(extmsg)//' a_sl'
      if (    prm%n_sl            <= 0.0_pReal)      extmsg = trim(extmsg)//' n_sl'
      if (any(xi_0_sl             <= 0.0_pReal))     extmsg = trim(extmsg)//' xi_0_sl'
      if (any(prm%xi_inf_sl       <= 0.0_pReal))     extmsg = trim(extmsg)//' xi_inf_sl'

    else slipActive
      xi_0_sl = emptyRealArray
      allocate(prm%xi_inf_sl,prm%h_int,source=emptyRealArray)
      allocate(prm%h_sl_sl(0,0))
    endif slipActive

!--------------------------------------------------------------------------------------------------
! twin related parameters
    N_tw         = pl%get_as1dInt('N_tw', defaultVal=emptyIntArray)
    prm%sum_N_tw = sum(abs(N_tw))
    twinActive: if (prm%sum_N_tw > 0) then
      prm%P_tw            = lattice_SchmidMatrix_twin(N_tw,phase%get_asString('lattice'),&
                                                      phase%get_asFloat('c/a',defaultVal=0.0_pReal))
      prm%h_tw_tw         = lattice_interaction_TwinByTwin(N_tw,&
                                                           pl%get_as1dFloat('h_tw_tw'), &
                                                           phase%get_asString('lattice'))
      prm%gamma_char      = lattice_characteristicShear_twin(N_tw,phase%get_asString('lattice'),&
                                                             phase%get_asFloat('c/a',defaultVal=0.0_pReal))

      xi_0_tw             = pl%get_as1dFloat('xi_0_tw',requiredSize=size(N_tw))

      prm%c_1             = pl%get_asFloat('c_1',defaultVal=0.0_pReal)
      prm%c_2             = pl%get_asFloat('c_2',defaultVal=1.0_pReal)
      prm%c_3             = pl%get_asFloat('c_3',defaultVal=0.0_pReal)
      prm%c_4             = pl%get_asFloat('c_4',defaultVal=0.0_pReal)
      prm%dot_gamma_0_tw  = pl%get_asFloat('dot_gamma_0_tw')
      prm%n_tw            = pl%get_asFloat('n_tw')
      prm%f_sat_sl_tw     = pl%get_asFloat('f_sat_sl_tw')
      prm%h_0_tw_tw       = pl%get_asFloat('h_0_tw_tw')

      ! expand: family => system
      xi_0_tw       = math_expand(xi_0_tw,N_tw)

      ! sanity checks
      if (prm%dot_gamma_0_tw <= 0.0_pReal)  extmsg = trim(extmsg)//' dot_gamma_0_tw'
      if (prm%n_tw           <= 0.0_pReal)  extmsg = trim(extmsg)//' n_tw'

    else twinActive
      xi_0_tw = emptyRealArray
      allocate(prm%gamma_char,source=emptyRealArray)
      allocate(prm%h_tw_tw(0,0))
    endif twinActive

!--------------------------------------------------------------------------------------------------
! slip-twin related parameters
    slipAndTwinActive: if (prm%sum_N_sl > 0 .and. prm%sum_N_tw > 0) then
      prm%h_0_tw_sl  = pl%get_asFloat('h_0_tw_sl')
      prm%h_sl_tw    = lattice_interaction_SlipByTwin(N_sl,N_tw,&
                                                      pl%get_as1dFloat('h_sl_tw'), &
                                                      phase%get_asString('lattice'))
      prm%h_tw_sl    = lattice_interaction_TwinBySlip(N_tw,N_sl,&
                                                      pl%get_as1dFloat('h_tw_sl'), &
                                                      phase%get_asString('lattice'))
    else slipAndTwinActive
      allocate(prm%h_sl_tw(prm%sum_N_sl,prm%sum_N_tw))                                              ! at least one dimension is 0
      allocate(prm%h_tw_sl(prm%sum_N_tw,prm%sum_N_sl))                                              ! at least one dimension is 0
      prm%h_0_tw_sl = 0.0_pReal
    endif slipAndTwinActive

!--------------------------------------------------------------------------------------------------
! KWN related parameters
    prm%kwn_nSteps = pl%get_asInt('kwn_nsteps',defaultVal=0)
    kwnActive: if (prm%kwn_nSteps > 0) then
      prm%kwn_stepsize = pl%get_asFloat('kwn_stepsize',defaultVal=1.0_pReal)
      prm%kwn_step0 = pl%get_asFloat('kwn_step0',defaultVal=1.0_pReal)
      prm%lattice_param = pl%get_asFloat('lattice_parameter')/lnorm !A
      prm%atomic_volume = pl%get_asFloat('atomic_volume')/lnorm/lnorm/lnorm !A^3/at
      prm%molar_volume = pl%get_asFloat('molar_volume') !m^3/mol
      prm%misfit_energy = pl%get_asFloat('misfit_energy',defaultVal=0.0_pReal)*lnorm !J/m^2/A
      prm%gamma_coherent = pl%get_asFloat('gamma_coherent') !J/m^2
      prm%vacancy_generation = pl%get_asFloat('vacancy_generation')
      prm%vacancy_sink_spacing = pl%get_asFloat('vacancy_sink_spacing') 
      prm%vacancy_energy = pl%get_asFloat('vacancy_energy')/kBnorm/lnorm/lnorm
      prm%migration_energy = pl%get_asFloat('solute_migration_energy')/kBnorm/lnorm/lnorm
      prm%diffusion0 = pl%get_asFloat('solute_diffusion0')/lnorm/lnorm !A^2/s
      prm%vacancy_migration_energy = pl%get_asFloat('vacancy_migration_energy')/kBnorm/lnorm/lnorm !no unit
      prm%vacancy_diffusion0 = pl%get_asFloat('vacancy_diffusion0') 
      prm%c0_matrix = pl%get_as1dFloat('c0_matrix', requiredSize=2)     
      prm%ceq_matrix = pl%get_as1dFloat('ceq_matrix', requiredSize=2)     
      prm%ceq_precipitate = pl%get_as1dFloat('ceq_precipitate', requiredSize=2)     
      prm%stoechiometry = pl%get_as1dFloat('stoechiometry', requiredSize=3)     
      prm%saturation_dislocation_density = pl%get_asFloat('saturation_dislocation_density')
      prm%dislocation_arrangement= pl%get_asFloat('dislocation_arrangement')
      prm%jog_formation_energy=pl%get_asFloat('jog_formation_energy')/kBnorm/lnorm/lnorm!no unit
      prm%shear_modulus = pl%get_asFloat('shear_modulus')*lnorm !J/m^2/A
      prm%burgers_vec = pl%get_asFloat('burgers_vector')/lnorm !A
      prm%solute_strength = pl%get_asFloat('solute_strength',defaultVal=0.0_pReal)
      prm%precipitate_strength = pl%get_asFloat('precipitate_strength_constant',defaultVal=0.0_pReal)
      prm%q_dislocation = pl%get_asFloat('dislocation_migration_energy',defaultVal=1.0e-15_pReal)/kBnorm/lnorm/lnorm !no unit
      

      ! when there is an initial distribution
      prm%mean_radius_initial = pl%get_asFloat('initial_mean_radius') !m
	  prm%standard_deviation = pl%get_asFloat('standard_deviation') !m
	  prm%volume_fraction_initial = pl%get_asFloat('initial_volume_fraction')

      ! sanity checks
      if (prm%kwn_stepsize             <= 0.0_pReal) extmsg = trim(extmsg)//' nkwnstepsize'
      if (prm%lattice_param            <= 0.0_pReal) extmsg = trim(extmsg)//' lattice_param'
      if (prm%atomic_volume            <= 0.0_pReal) extmsg = trim(extmsg)//' atomic_volume'
      if (prm%molar_volume             <= 0.0_pReal) extmsg = trim(extmsg)//' molar_volume'
      if (prm%misfit_energy            <  0.0_pReal) extmsg = trim(extmsg)//' misfit_energy'
      if (prm%gamma_coherent           <= 0.0_pReal) extmsg = trim(extmsg)//' gamma_coherent'
      if (prm%vacancy_generation       <  0.0_pReal) extmsg = trim(extmsg)//' vacancy_generation'
      if (prm%vacancy_sink_spacing     <  0.0_pReal) extmsg = trim(extmsg)//' vacancy_sink_spacing'
      if (prm%vacancy_energy           <  0.0_pReal) extmsg = trim(extmsg)//' vacancy_energy'
      if (prm%migration_energy         <= 0.0_pReal) extmsg = trim(extmsg)//' migration_energy'
      if (prm%diffusion0               <  0.0_pReal) extmsg = trim(extmsg)//' diffusion0'
      if (prm%vacancy_migration_energy <= 0.0_pReal) extmsg = trim(extmsg)//' vacancy_migration_energy'
      if (prm%vacancy_diffusion0       <  0.0_pReal) extmsg = trim(extmsg)//' vacancy_diffusion0'
      if (minval(prm%c0_matrix )       <  0.0_pReal) extmsg = trim(extmsg)//' c0_matrix'
      if (minval(prm%ceq_matrix)       <  0.0_pReal) extmsg = trim(extmsg)//' ceq_matrix'
      if (minval(prm%ceq_precipitate)  <  0.0_pReal) extmsg = trim(extmsg)//' ceq_precipitate'
    
      ! initialize bins
      allocate(prm%bins(0:prm%kwn_nSteps), source=0.0_pReal)
      kwnbins_init: do i = 0, prm%kwn_nSteps
        !prm%bins(i) = 10.0_pReal**(real(i,pReal)*prm%kwn_stepsize + prm%kwn_step0)
        prm%bins(i)=prm%kwn_step0+real(i,pReal)*prm%kwn_stepsize
      enddo kwnbins_init   
    endif kwnActive

!--------------------------------------------------------------------------------------------------
!  output pararameters

#if defined (__GFORTRAN__)
    prm%output = output_as1dString(pl)
#else
    prm%output = pl%get_as1dString('output',defaultVal=emptyStringArray)
#endif

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    Nmembers = count(material_phaseID == ph)
    sizeDotState = size(['xi_sl   ','gamma_sl']) * prm%sum_N_sl &
                 + size(['xi_tw   ','gamma_tw']) * prm%sum_N_tw &
                 + prm%kwn_nSteps + 2
    sizeState = sizeDotState


    call phase_allocateState(plasticState(ph),Nmembers,sizeState,sizeDotState,0)

!--------------------------------------------------------------------------------------------------
! state aliases and initialization
    startIndex = 1
    endIndex   = prm%sum_N_sl
    stt%xi_slip => plasticState(ph)%state   (startIndex:endIndex,:)
    stt%xi_slip =  spread(xi_0_sl, 2, Nmembers)
    dot%xi_slip => plasticState(ph)%dotState(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_xi',defaultVal=1.0_pReal)
    if(any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_xi'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_tw
    stt%xi_twin => plasticState(ph)%state   (startIndex:endIndex,:)
    stt%xi_twin =  spread(xi_0_tw, 2, Nmembers)
    dot%xi_twin => plasticState(ph)%dotState(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_xi',defaultVal=1.0_pReal)
    if(any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_xi'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    stt%gamma_slip => plasticState(ph)%state   (startIndex:endIndex,:)
    dot%gamma_slip => plasticState(ph)%dotState(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)
    if(any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_gamma'
    ! global alias
    plasticState(ph)%slipRate => plasticState(ph)%dotState(startIndex:endIndex,:)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_tw
    stt%gamma_twin => plasticState(ph)%state   (startIndex:endIndex,:)
    dot%gamma_twin => plasticState(ph)%dotState(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)
    if(any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_gamma'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%kwn_nSteps
    stt%precipitate_density => plasticState(ph)%state   (startIndex:endIndex,:)
    dot%precipitate_density => plasticState(ph)%dotState(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_precipitate',defaultVal=1.0_pReal)
    if(any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_precipitate'

    startIndex = endIndex + 1
    endIndex   = endIndex + 1
    stt%c_vacancy => plasticState(ph)%state(startIndex,:)
    stt%c_vacancy = spread(pl%get_asFloat('c0_vacancy',defaultVal=0.0_pReal),1,Nmembers)
    
    dot%c_vacancy => plasticState(ph)%dotState(startIndex,:)
    plasticState(ph)%atol(startIndex) = pl%get_asFloat('atol_solute',defaultVal=1.0_pReal)
    if(plasticState(ph)%atol(startIndex) < 0.0_pReal) extmsg = trim(extmsg)//' atol_solute'

    startIndex = endIndex + 1
    endIndex   = endIndex + 1
    stt%time => plasticState(ph)%state   (startIndex,:)
    dot%time => plasticState(ph)%dotState(startIndex,:)
    plasticState(ph)%atol(startIndex) = 1.0e-6_pReal

    allocate(dst%total_precipitate_density(  Nmembers), source=0.0_pReal)
    allocate(dst%avg_precipitate_radius   (  Nmembers), source=0.0_pReal)
    allocate(dst%precipitate_volume_frac  (  Nmembers), source=0.0_pReal)
    allocate(dst%c_matrix                 (2,Nmembers), source=0.0_pReal)
    allocate(dst%interface_concentration(prm%kwn_nSteps-1,Nmembers), source=prm%ceq_precipitate(1)/2.0)
    
    
    
    !!!!!!!!!!!!!!!!!!!!
    !initialize precipitate density for an already existing distribution
    if (prm%mean_radius_initial>0) then
      allocate(normalized_distribution_function(prm%kwn_nSteps), source=0.0_pReal)
      shape_factor = sqrt(log(1+prm%standard_deviation**2/prm%mean_radius_initial**2))
       
      distribution_function: do i = 1, prm%kwn_nSteps
        x = prm%bins(i)*lnorm/prm%mean_radius_initial
        !log normal distribution
        radiusL = prm%bins(i-1)
        radiusR = prm%bins(i  )
        radiusC = (radiusL+radiusR)/2.0
        normalized_distribution_function(i) = 1.0/sqrt(2.0*PI)/shape_factor/(radiusC) &
                                            * exp(-(log(x)+shape_factor**2/2.0)**2/2.0/shape_factor**2) !at this stage in /A
      enddo distribution_function

      !initialize the number density distribution
      !the normalized distribution function gives the shape of the distribution, it needs to be multiplied by N0 to match the initial precipitate fraction
      N_0 = 0.0_pReal
      do i = 1, prm%kwn_nSteps
        radiusL = prm%bins(i-1)
        radiusR = prm%bins(i  )
        radiusC = (radiusL+radiusR)/2
        N_0 = N_0 + normalized_distribution_function(i)*radiusC**3*(radiusR-RadiusL)*4.0/3.0*PI    
      enddo
      N_0 = prm%volume_fraction_initial/N_0
      normalized_distribution_function = normalized_distribution_function*N_0
      
      
      
      stt%precipitate_density =  spread(normalized_distribution_function, 2, Nmembers)
      precipitate_volume_frac = 0.0_pReal
      total_precipitate_density = 0.0_pReal
      avg_precipitate_radius = 0.0_pReal
      do i = 1, prm%kwn_nSteps
        radiusL = prm%bins(i-1)
        radiusR = prm%bins(i  )
        total_precipitate_density = total_precipitate_density &
                                  + normalized_distribution_function(i) &
                                  * (radiusR - radiusL)
        avg_precipitate_radius = avg_precipitate_radius &
                               + normalized_distribution_function(i) &
                               * (radiusR**2.0_pReal - radiusL**2.0_pReal)/2.0_pReal
        precipitate_volume_frac = precipitate_volume_frac &
                                + 1.0_pReal/6.0_pReal*PI &
                                * (radiusR+ radiusL)**3.0_pReal &
                                * (radiusR - radiusL) &
                                * normalized_distribution_function(i) 
      enddo
      if (total_precipitate_density > 0.0_pReal) &
        avg_precipitate_radius = avg_precipitate_radius &
                               / total_precipitate_density
      dst%total_precipitate_density = total_precipitate_density
      dst%avg_precipitate_radius = avg_precipitate_radius
      dst%precipitate_volume_frac = precipitate_volume_frac
      if (precipitate_volume_frac < 1.0_pReal) &
        c_matrix = (prm%c0_matrix - precipitate_volume_frac*prm%ceq_precipitate) &
                 / (1.0_pReal - precipitate_volume_frac)
      dst%c_matrix =  spread(c_matrix, 2, Nmembers)
    endif
 
    end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(kwnpowerlaw)')

  enddo

end function plastic_kwnpowerlaw_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate plastic velocity gradient and its tangent.
!> @details asummes that deformation by dislocation glide affects twinned and untwinned volume
!  equally (Taylor assumption). Twinning happens only in untwinned volume
!--------------------------------------------------------------------------------------------------
pure module subroutine kwnpowerlaw_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)

  real(pReal), dimension(3,3),     intent(out) :: &
    Lp                                                                                              !< plastic velocity gradient
  real(pReal), dimension(3,3,3,3), intent(out) :: &
    dLp_dMp                                                                                         !< derivative of Lp with respect to the Mandel stress

  real(pReal), dimension(3,3), intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,               intent(in) :: &
    ph, &
    en

  integer :: &
    i,k,l,m,n
  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    gdot_slip_pos,gdot_slip_neg, &
    dgdot_dtauslip_pos,dgdot_dtauslip_neg, &
     tau_slip_pos, tau_slip_neg, tau_kwn
  real(pReal), dimension(param(ph)%sum_N_tw) :: &
    gdot_twin,dgdot_dtautwin

  Lp = 0.0_pReal
  dLp_dMp = 0.0_pReal

  associate(prm => param(ph))

  call kinetics_slip(Mp,ph,en,gdot_slip_pos,gdot_slip_neg,  tau_slip_pos, tau_slip_neg, tau_kwn, dgdot_dtauslip_pos,dgdot_dtauslip_neg)
  slipSystems: do i = 1, prm%sum_N_sl
    Lp = Lp + (gdot_slip_pos(i)+gdot_slip_neg(i))*prm%P_sl(1:3,1:3,i)
    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
      dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                       + dgdot_dtauslip_pos(i) * prm%P_sl(k,l,i) * prm%nonSchmid_pos(m,n,i) &
                       + dgdot_dtauslip_neg(i) * prm%P_sl(k,l,i) * prm%nonSchmid_neg(m,n,i)
  enddo slipSystems

  call kinetics_twin(Mp,ph,en,gdot_twin,dgdot_dtautwin)
  twinSystems: do i = 1, prm%sum_N_tw
    Lp = Lp + gdot_twin(i)*prm%P_tw(1:3,1:3,i)
    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
      dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                       + dgdot_dtautwin(i)*prm%P_tw(k,l,i)*prm%P_tw(m,n,i)
  enddo twinSystems

  end associate

end subroutine kwnpowerlaw_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief Calculate the rate of change of microstructure.
!--------------------------------------------------------------------------------------------------
module subroutine kwnpowerlaw_dotState(Mp,ph,en)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en
  real(pReal) :: &
    T

  real(pReal) :: &
    c_SlipSlip,c_TwinSlip,c_TwinTwin, &
    xi_slip_sat_offset,&
    sumGamma,sumF
  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    left_SlipSlip,right_SlipSlip, &
    gdot_slip_pos,gdot_slip_neg, &
     tau_slip_pos, tau_slip_neg, tau_kwn
  real(pReal) :: &
    deltaGv, &
    interface_energy, &
    radius_crit, &
    nucleation_site_density, &
    zeldovich_factor, &
    beta_star, &
    nucleation_rate,&
    diffusion_coefficient, &
    radiusL, radiusR, radiusC, &
    growth_rate, flux, &
    vacancy_generation, &
    vacancy_annihilation, &
    vacancy_sink_spacing, &
    c_j, &
    c_thermal_vacancy
 
  integer(pInt) :: &
    bin

  associate(prm => param(ph), &
            dot => dotState(ph), &
            stt => state(ph), &
            dst => dependentState(ph))

  T = thermal_T(ph,en)
  !print*, 'T:', T
  sumGamma = sum(stt%gamma_slip(:,en))
  sumF     = sum(stt%gamma_twin(:,en)/prm%gamma_char)

!--------------------------------------------------------------------------------------------------
! system-independent (nonlinear) prefactors to M_Xx (X influenced by x) matrices
  c_SlipSlip = prm%h_0_sl_sl * (1.0_pReal + prm%c_1*sumF** prm%c_2)
  c_TwinSlip = prm%h_0_tw_sl * sumGamma**prm%c_3
  c_TwinTwin = prm%h_0_tw_tw * sumF**prm%c_4

!--------------------------------------------------------------------------------------------------
!  calculate left and right vectors
  left_SlipSlip  = 1.0_pReal + prm%h_int
  xi_slip_sat_offset = prm%f_sat_sl_tw*sqrt(sumF)
  right_SlipSlip = abs(1.0_pReal-stt%xi_slip(:,en) / (prm%xi_inf_sl+xi_slip_sat_offset)) **prm%a_sl &
                 * sign(1.0_pReal,1.0_pReal-stt%xi_slip(:,en) / (prm%xi_inf_sl+xi_slip_sat_offset))

!--------------------------------------------------------------------------------------------------
! shear rates
  call kinetics_slip(Mp,ph,en,gdot_slip_pos,gdot_slip_neg,  tau_slip_pos, tau_slip_neg, tau_kwn)
  dot%gamma_slip(:,en) = abs(gdot_slip_pos+gdot_slip_neg)
  call kinetics_twin(Mp,ph,en,dot%gamma_twin(:,en))

!--------------------------------------------------------------------------------------------------
! hardening
  dot%xi_slip(:,en) = c_SlipSlip * left_SlipSlip * &
                      matmul(prm%h_sl_sl,dot%gamma_slip(:,en)*right_SlipSlip) &
                    + matmul(prm%h_sl_tw,dot%gamma_twin(:,en))

  dot%xi_twin(:,en) = c_TwinSlip * matmul(prm%h_tw_sl,dot%gamma_slip(:,en)) &
                    + c_TwinTwin * matmul(prm%h_tw_tw,dot%gamma_twin(:,en))

!--------------------------------------------------------------------------------------------------
! KWN model (work in progress)
  kwnActive: if (prm%kwn_nSteps > 0) then
    ! vacancy evolution rate
    ! ToDo: Use analytical solution
    interface_energy = prm%gamma_coherent
    
        vacancy_sink_spacing = min( &
                                  0.27*prm%shear_modulus*prm%burgers_vec &
                               /  maxval(stt%xi_slip(:,en)), &
                                  prm%vacancy_sink_spacing &
                               )
    
    
    ! modified version Madeleine, following Joe's paper
    c_thermal_vacancy = 23.0*exp(-prm%vacancy_energy/T) !Murty
    c_j = exp(-prm%jog_formation_energy/T)
    
  !  vacancy_generation = prm%vacancy_generation* sum(dot%gamma_slip(:,en)*(stt%xi_slip(:,en))) &
   ! 				   * prm%atomic_volume*lnorm/prm%vacancy_energy/kBnorm &
	!                   + 0.5*c_j*prm%atomic_volume/4.0/prm%burgers_vec**3*sum(dot%gamma_slip(:,en))
	                   
	!tets
	 vacancy_generation = prm%vacancy_generation* sum(dot%gamma_slip(:,en)*abs(tau_slip_pos(:))) &
    				   * prm%atomic_volume*lnorm/prm%vacancy_energy/kBnorm &
	                   + 0.5*c_j*prm%atomic_volume/4.0/prm%burgers_vec**3*sum(dot%gamma_slip(:,en))
	       

	!print*, 'sum(dot%gamma_slip(:,en)*abs(tau_slip_pos(:)))', sum(dot%gamma_slip(:,en)*abs(tau_slip_pos(:)))
	!print*, 'sum(dot%gamma_slip(:,en)*(stt%xi_slip(:,en)))', sum(dot%gamma_slip(:,en)*(stt%xi_slip(:,en)))
	
	              
	vacancy_annihilation = prm%vacancy_diffusion0*exp(-prm%vacancy_migration_energy/T) &
					     * ( 1.0/prm%vacancy_sink_spacing**2 &
					       ) &  
					     * stt%c_vacancy(en)   
    dot%c_vacancy(en) = vacancy_generation - vacancy_annihilation
    
    !change Madeleine to account for multiple elements
    deltaGv = -Rnorm*T/prm%molar_volume*sum(prm%ceq_precipitate(:)*log(dst%c_matrix(:,en)/prm%ceq_matrix(:))) ! J/m^2/A
    radius_crit = -2.0_pReal*prm%gamma_coherent &
                / (  &
                      prm%misfit_energy + deltaGv &
                   ) !A
                   
    diffusion_coefficient = prm%diffusion0*exp(-(prm%migration_energy)/T)*(1.0+ stt%c_vacancy(en)/c_thermal_vacancy)	&
    						 +2*( maxval(stt%xi_slip(:,en))/0.27/prm%shear_modulus/prm%burgers_vec)**2*prm%atomic_volume/prm%burgers_vec&
  	 						 *prm%diffusion0*exp(-(prm%q_dislocation )/T) !A^2/s
   ! print*, 'T:', T-273, '°C'
   !  print*, 'Diffusion coefficient:', T-273, '°C'
    nucleation_site_density = sum(dst%c_matrix(:,en))/prm%atomic_volume !/A^3
    zeldovich_factor = prm%atomic_volume*sqrt(prm%gamma_coherent/kBnorm/T) &
                     / 2.0_pReal/PI/radius_crit/radius_crit ![]
    beta_star = 4.0_pReal*PI*diffusion_coefficient*sum(dst%c_matrix(:,en)) &
              * radius_crit*radius_crit/(prm%lattice_param**4.0) ![/s]
              
    if (stt%time(en) > 0.0_pReal) then
      nucleation_rate = nucleation_site_density*zeldovich_factor*beta_star &
                      * exp( &
                             - 4.0_pReal*PI*interface_energy*radius_crit*radius_crit/3.0_pReal/kBnorm/T &
                            )
    nucleation_rate=0.0 ! block nucleation rate for now
    else
      nucleation_rate = 0.0_pReal
    endif                        

    ! growth rate in # per Angstrom^3
    dot%precipitate_density(:,en) = 0.0_pReal
    kwnbins_growth: do bin = 1, prm%kwn_nSteps-1
      radiusC = prm%bins(bin  )
      radiusL = prm%bins(bin-1)
      radiusR = prm%bins(bin+1)
      ! ToDo: limit max value of interface_c to ceq
      ! change Madeleine to make the expression of the radius at the interface consistent with that of the critical radius (divide by ceq precipitate)
      ! calculate interface concentration - numerically as it is for a multicomponent system
      ! Calculate the equilibrium at the interface for all bins, only once (when time =0), for it is considered as constant 
      growth_rate = diffusion_coefficient/radiusC &
                  * (dst%c_matrix(1,en)     - dst%interface_concentration(bin,en)) &
                  / (prm%ceq_precipitate(1) - dst%interface_concentration(bin,en)) !A/s
      !print*,'bin:', radiusC
      !print*,  'growth rate: ' ,growth_rate,  'A/s'
      !print*, 'interface concentration', dst%interface_concentration(bin,en)
      !print*, 'diffusion coefficient', diffusion_coefficient
      if (growth_rate > 0.0_pReal) then
        flux = stt%precipitate_density(bin  ,en)*growth_rate
      else
        flux = stt%precipitate_density(bin+1,en)*growth_rate
      endif   
      dot%precipitate_density(bin  ,en) = dot%precipitate_density(bin  ,en) - flux/(radiusC - radiusL)
      dot%precipitate_density(bin+1,en) = dot%precipitate_density(bin+1,en) + flux/(radiusR - radiusC)
      if (radius_crit > radiusL .and. radius_crit < radiusC) &
        dot%precipitate_density(bin,en) = dot%precipitate_density(bin,en) &
                                        + nucleation_rate/(radiusC - radiusL)
    enddo kwnbins_growth

    ! enforce boundary condition, phi(0) = 0
    dot%precipitate_density(1,en) = 0.0_pReal
   
    dot%time(en) = 1.0_pReal
  endif kwnActive

  end associate

end subroutine kwnpowerlaw_dotState


!--------------------------------------------------------------------------------------------------
!> @brief calculates the derived quantities from state
!--------------------------------------------------------------------------------------------------
module subroutine kwnpowerlaw_dependentState(ph,en)

  integer, intent(in) :: &
    ph, &
    en

  real(pReal) :: &
    T
  real(pReal) :: &
    radiusL, radiusR
  integer(pInt) :: &
    bin, iter
  real(pReal) :: &
    interface_c, residual, jacobian, searchdirection, err_current
  			 
  			 
  associate(prm => param(ph), &
            stt => state(ph), &
            dst => dependentState(ph))
  T = thermal_T(ph,en)
  
!--------------------------------------------------------------------------------------------------
! KWN model (work in progress)
  kwnActive: if (prm%kwn_nSteps > 0) then
  
    dst%precipitate_volume_frac(en) = 0.0_pReal
    dst%total_precipitate_density(en) = 0.0_pReal
    dst%avg_precipitate_radius(en) = 0.0_pReal
    kwnbins: do bin = 1, prm%kwn_nSteps
      radiusL = prm%bins(bin-1)
      radiusR = prm%bins(bin  )
      dst%total_precipitate_density(en) = dst%total_precipitate_density(en) &
                                        + stt%precipitate_density(bin,en) &
                                        * (radiusR - radiusL)
      dst%avg_precipitate_radius(en) = dst%avg_precipitate_radius(en) &
                                     + stt%precipitate_density(bin,en) &
                                     * (radiusR**2.0_pReal - radiusL**2.0_pReal)/2.0_pReal
  
      !Madeleine : small change compared to the initial version
      dst%precipitate_volume_frac(en) = dst%precipitate_volume_frac(en) &
                                      + 1.0_pReal/6.0_pReal*PI &
                                      * (radiusR+ radiusL)**3.0_pReal &
                                      * (radiusR - radiusL) &
                                      * stt%precipitate_density(bin,en) 	
    enddo kwnbins
    if (dst%total_precipitate_density(en) > 0.0_pReal) &
      dst%avg_precipitate_radius(en) = dst%avg_precipitate_radius(en) &
                                     / dst%total_precipitate_density(en)
                                     
	
	if (dst%precipitate_volume_frac(en) < 1.0_pReal) &
      dst%c_matrix(:,en) = (prm%c0_matrix(:) - dst%precipitate_volume_frac(en)*prm%ceq_precipitate(:)) &
                         / (1.0_pReal - dst%precipitate_volume_frac(en))
      
    kwnbins_c: do bin = 1, prm%kwn_nSteps-1
      ! interface concentration calculation
      radiusR = prm%bins(bin  )
      interface_c  = dst%interface_concentration(bin,en)
      err_current = huge(1.0_pReal)
 	  iter = 0
      do while (err_current > 1e-12)
        iter = iter + 1
        residual = ((interface_c/prm%ceq_matrix(1))**prm%stoechiometry(1)) &
                 * ( &
                    (prm%c0_matrix(2) + prm%stoechiometry(2)/prm%stoechiometry(1)*(interface_c - prm%c0_matrix(1)))/prm%ceq_matrix(2) &
                   )**prm%stoechiometry(2) &
                 - exp(2.0*prm%molar_volume*prm%gamma_coherent/Rnorm/T/radiusR*(sum(prm%stoechiometry)))
        jacobian = (prm%stoechiometry(1)*(interface_c/prm%ceq_matrix(1))**(prm%stoechiometry(1) - 1))/prm%ceq_matrix(1) &
                 * ( &
                    (prm%c0_matrix(2) + prm%stoechiometry(2)/prm%stoechiometry(1)*(interface_c - prm%c0_matrix(1)))/prm%ceq_matrix(2) &
                   )**prm%stoechiometry(2) &
                 + ((interface_c/prm%ceq_matrix(1))**prm%stoechiometry(1)) &
                 * (prm%stoechiometry(2)*prm%stoechiometry(2)/prm%stoechiometry(1)/prm%ceq_matrix(2)) &
                 * ( &
                    (prm%c0_matrix(2) + prm%stoechiometry(2)/prm%stoechiometry(1)*(interface_c - prm%c0_matrix(1)))/prm%ceq_matrix(2) &
                   )**(prm%stoechiometry(2) - 1)
        searchdirection = residual/jacobian
        interface_c = interface_c - searchdirection
        err_current = abs(searchdirection)
      enddo
      dst%interface_concentration(bin,en) = interface_c
    enddo kwnbins_c
  endif kwnActive

  end associate

end subroutine kwnpowerlaw_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_kwnpowerlaw_results(ph,group)

  integer,          intent(in) :: ph
  character(len=*), intent(in) :: group

  integer :: o

  associate(prm => param(ph), &
            stt => state(ph), &
            dst => dependentState(ph))
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))

      case('xi_sl')
        if(prm%sum_N_sl>0) call results_writeDataset(stt%xi_slip   ,group,trim(prm%output(o)), &
                                                     'resistance against plastic slip','Pa')
      case('gamma_sl')
        if(prm%sum_N_sl>0) call results_writeDataset(stt%gamma_slip,group,trim(prm%output(o)), &
                                                     'plastic shear','1')

	  	
	  
      case('xi_tw')
        if(prm%sum_N_tw>0) call results_writeDataset(stt%xi_twin,   group,trim(prm%output(o)), &
                                                     'resistance against twinning','Pa')
      case('gamma_tw')
        if(prm%sum_N_tw>0) call results_writeDataset(stt%gamma_twin,group,trim(prm%output(o)), &
                                                     'twinning shear','1')
      case ('phi')
        call results_writeDataset(stt%precipitate_density,group,'precipitate_density', &
                                  'precipitate number density','/m^3')
                                 
      case ('phi_total')
        call results_writeDataset(dst%total_precipitate_density,group,'total_precipitate_density', &
                                  'total precipitate number density','/m^3')
                                 ! print*,'total_precipitate_density', dst%total_precipitate_density, '/A^3'
      case ('f')
        call results_writeDataset(dst%precipitate_volume_frac,group,'precipitate_volume_frac', &
                                  'precipitate volume fraction','1')
                                  ! print*,'precipitate volume fraction', dst%precipitate_volume_frac
                                  ! print*, 'precipitate strength', prm%precipitate_strength 
      case ('r_avg')
        call results_writeDataset(dst%avg_precipitate_radius,group,'avg_precipitate_radius', &
                                  'average precipitate radius','m')
                                  !print*,'average precipitate radius',dst%avg_precipitate_radius, 'A'
      case ('solute_c')
        call results_writeDataset(dst%c_matrix,group,'c_matrix', &
                                  'matrix solute composition','1')

 		
      case ('vacancy_c')
        call results_writeDataset(stt%c_vacancy,group,'c_vacancy', &
                                  'matrix vacancy composition','1')
                                  !print*,'Vacancy composition',stt%c_vacancy

    
        

    end select
  enddo outputsLoop
  end associate

end subroutine plastic_kwnpowerlaw_results


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on slip systems and their derivatives with respect to resolved
!         stress.
!> @details Derivatives are calculated only optionally.
! NOTE: Against the common convention, the result (i.e. intent(out)) variables are the last to
! have the optional arguments at the end.
!--------------------------------------------------------------------------------------------------
!pure subroutine kinetics_slip(Mp,ph,en, &
 !                             gdot_slip_pos,gdot_slip_neg,dgdot_dtau_slip_pos,dgdot_dtau_slip_neg)
 pure subroutine kinetics_slip(Mp,ph,en, &
                              gdot_slip_pos,gdot_slip_neg, tau_slip_pos, tau_slip_neg, tau_kwn, dgdot_dtau_slip_pos,dgdot_dtau_slip_neg)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en

  real(pReal),                  intent(out), dimension(param(ph)%sum_N_sl) :: &
    gdot_slip_pos, &
    gdot_slip_neg, &
    tau_slip_pos, &
    tau_slip_neg, &
    tau_kwn

  real(pReal),                  intent(out), optional, dimension(param(ph)%sum_N_sl) :: &
    dgdot_dtau_slip_pos, &
    dgdot_dtau_slip_neg


!  real(pReal), dimension(param(ph)%sum_N_sl) :: &
 !   tau_slip_pos, &
  !  tau_slip_neg, &
   ! tau_kwn
  real(pReal) :: &
    mean_particle_strength, &
    radiusL, radiusR, radius_transition
    
  integer :: i, bin

  associate(prm => param(ph), &
            stt => state(ph), &
            dst => dependentState(ph))

  do i = 1, prm%sum_N_sl
    tau_slip_pos(i) =       math_tensordot(Mp,prm%nonSchmid_pos(1:3,1:3,i))
    tau_slip_neg(i) = merge(math_tensordot(Mp,prm%nonSchmid_neg(1:3,1:3,i)), &
                            0.0_pReal, prm%nonSchmidActive)
  enddo

  tau_kwn = 0.0_pReal
  kwnActive: if (prm%kwn_nSteps > 0) then
    !Madeleine : modified compared to the version of Pratheek  - contains one fitted constant here (fitted with Deschamps 2012 7068) and the transition radius
    if (dst%avg_precipitate_radius(en) > 0.0_pReal) then
      mean_particle_strength = 0.0_pReal
      kwnbins_strength: do bin = 1, prm%kwn_nSteps
        radiusL = prm%bins(bin-1)
        radiusR = prm%bins(bin  )
        radius_transition=33.0 !Fribourg2011 in A
        if (radiusL < radius_transition) then ! particle shearing
          mean_particle_strength = mean_particle_strength  &
                                 + (radiusR**2.0_pReal - radiusL**2.0_pReal)/2.0_pReal &
                                 * stt%precipitate_density(bin,en)
        else                            ! particle looping
          mean_particle_strength = mean_particle_strength &
                                 + radius_transition &
                                 * (radiusR - radiusL) &
                                 * stt%precipitate_density(bin,en)
        endif
      enddo kwnbins_strength
      mean_particle_strength = mean_particle_strength/dst%total_precipitate_density(en) ! A (A^-2/(A-3))
      tau_kwn = prm%precipitate_strength & !fitted constant -test
             *	(prm%shear_modulus/lnorm) &
              * (dst%precipitate_volume_frac(en)*3.0/2.0/PI)**(1.0_pReal/2.0_pReal) &
              / dst%avg_precipitate_radius(en)/sqrt(radius_transition) &
              * mean_particle_strength**(3.0_pReal/2.0_pReal) ! Pa
              
   
         
    endif   
   
            

  endif kwnActive

  where(dNeq0(tau_slip_pos))
    gdot_slip_pos = prm%dot_gamma_0_sl * merge(0.5_pReal,1.0_pReal, prm%nonSchmidActive) &          ! 1/2 if non-Schmid active
                  * sign(abs(tau_slip_pos/(prm%solute_strength*sum(dst%c_matrix(:,en))**(2.0_pReal/3.0_pReal)+(stt%xi_slip(:,en)**2+tau_kwn**2)**(0.5)))**prm%n_sl,  tau_slip_pos)
  else where
    gdot_slip_pos = 0.0_pReal
  end where

  where(dNeq0(tau_slip_neg))
    gdot_slip_neg = prm%dot_gamma_0_sl * 0.5_pReal &                                                ! only used if non-Schmid active, always 1/2
                  * sign(abs(tau_slip_neg/(prm%solute_strength*sum(dst%c_matrix(:,en))**(2.0_pReal/3.0_pReal)+(stt%xi_slip(:,en)**2+tau_kwn**2)**(0.5)))**prm%n_sl,  tau_slip_neg)
  else where
    gdot_slip_neg = 0.0_pReal
  end where

  if (present(dgdot_dtau_slip_pos)) then
    where(dNeq0(gdot_slip_pos))
      dgdot_dtau_slip_pos = gdot_slip_pos*prm%n_sl/tau_slip_pos
    else where
      dgdot_dtau_slip_pos = 0.0_pReal
    end where
  endif
  if (present(dgdot_dtau_slip_neg)) then
    where(dNeq0(gdot_slip_neg))
      dgdot_dtau_slip_neg = gdot_slip_neg*prm%n_sl/tau_slip_neg
    else where
      dgdot_dtau_slip_neg = 0.0_pReal
    end where
  endif
  end associate

end subroutine kinetics_slip


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on twin systems and their derivatives with respect to resolved
!         stress. Twinning is assumed to take place only in untwinned volume.
!> @details Derivatives are calculated only optionally.
! NOTE: Against the common convention, the result (i.e. intent(out)) variables are the last to
! have the optional arguments at the end.
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_twin(Mp,ph,en,&
                              gdot_twin,dgdot_dtau_twin)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en

  real(pReal), dimension(param(ph)%sum_N_tw), intent(out) :: &
    gdot_twin
  real(pReal), dimension(param(ph)%sum_N_tw), intent(out), optional :: &
    dgdot_dtau_twin

  real(pReal), dimension(param(ph)%sum_N_tw) :: &
    tau_twin
  integer :: i

  associate(prm => param(ph), stt => state(ph))

  do i = 1, prm%sum_N_tw
    tau_twin(i)  = math_tensordot(Mp,prm%P_tw(1:3,1:3,i))
  enddo

  where(tau_twin > 0.0_pReal)
    gdot_twin = (1.0_pReal-sum(stt%gamma_twin(:,en)/prm%gamma_char)) &                              ! only twin in untwinned volume fraction
              * prm%dot_gamma_0_tw*(abs(tau_twin)/stt%xi_twin(:,en))**prm%n_tw
  else where
    gdot_twin = 0.0_pReal
  end where

  if (present(dgdot_dtau_twin)) then
    where(dNeq0(gdot_twin))
      dgdot_dtau_twin = gdot_twin*prm%n_tw/tau_twin
    else where
      dgdot_dtau_twin = 0.0_pReal
    end where
  endif

  end associate

end subroutine kinetics_twin

end submodule kwnpowerlaw
