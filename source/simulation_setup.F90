!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module that defines simulation type and procedures, togethere with 
! all subroutines to read simulation related directives
!
! Copyright:  2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author    - i.scivetti Mar 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module simulation_setup

  Use constants,        Only : max_components
  
  Use input_types,      Only : in_integer, &
                               in_integer_array,&
                               in_logic,   &
                               in_string,  &
                               in_param,   &
                               in_param_array, &
                               in_scalar
                               
  Use numprec,          Only : wi, &
                               wp
                               
  Use process_data,     Only : capital_to_lower_case, &
                               remove_symbols,&
                               set_read_status,&
                               check_for_rubbish, &
                               check_end, &
                               get_word_length,&
                               duplication_error
                               
  Use unit_output,      Only : error_stop,&
                               info 

  Implicit None
  Private

  ! Maximum directives for simulations
  Integer(Kind=wi), Parameter, Public  :: max_directives=100 
  ! Maximum number of Boltzmann ions
  Integer(Kind=wi), Parameter, Public  :: max_number_boltzmann_ions=20

  ! Components inherited from model_data
  Type :: component_in_block
    Character(Len=8) :: tag
    Character(Len=2) :: element
    Integer(Kind=wi) :: atomic_number
    Real(Kind=wp)    :: valence_electrons
  End Type component_in_block

  ! NGWF type 
  Type :: type_ngwf
    Character(Len=8) :: tag
    Character(Len=2) :: element
    Integer(Kind=wi) :: ni
    Real(Kind=wp)    :: radius
  End Type type_ngwf

  ! Type for pseudopotentials 
  Type :: type_pseudo
    Character(Len=256) :: file_name
    Character(Len=256) :: potential 
    Character(Len=8)   :: tag
    Character(Len=2)   :: element
    Logical            :: defined
  End Type type_pseudo

  ! Type for extra directives 
  Type, Public :: type_extra 
    Character(Len=256) :: array(max_directives)
    Character(Len=256) :: key(max_directives)
    Character(Len=256) :: set(max_directives)
    Integer(Kind=wi)   :: N0
  End Type type_extra

  ! type for reference data
  Type, Public :: type_ref_data
    Character(Len=256) :: key
    Character(Len=16)  :: keytype
    Character(Len=256) :: msn
    Character(Len=256) :: set_default
    Character(Len=16)  :: units
    Integer(Kind=wi)   :: N0
  End Type type_ref_data

  ! Type for basis_set
  Type :: type_basis
    Character(Len=8)   :: tag
    Character(Len=2)   :: element
    Character(Len=256) :: basis
    Character(Len=256) :: type
    Logical            :: defined
  End Type type_basis

  ! Type for magnetization
  Type :: type_mag
    Character(Len=8)  :: tag
    Character(Len=2)  :: element
    Real(Kind=wp)     :: value
  End Type type_mag

  ! Type for Hubbard corrections 
  Type :: type_hubbard
    Character(Len=8)  :: tag
    Character(Len=2)  :: element
    Integer(Kind=wi)  :: l_orbital  
    Real(Kind=wp)     :: U 
    Real(Kind=wp)     :: J 
  End Type type_hubbard

  ! Type to read quantities from a list of species from two-rows blocks
  Type, Public :: species_list
    Character(Len=8)  :: tag
    Character(Len=2)  :: element
    Real(Kind=wp)     :: value
  End Type species_list

  ! Type for GC-DFT
  Type :: type_gcdft 
    Type(in_logic)  ::  activate     
    Type(in_param)  ::  reference_potential    
    Type(in_param)  ::  electrode_potential
    Type(in_scalar) ::  electron_threshold
    Type(in_param)  ::  target_workfunction
    Type(in_scalar) ::  mixing_coefficient
  End Type
 
  ! Type for boltzmann_ions 
  Type :: boltzmann_ions
    Character(Len=8)  :: tag
    Real(Kind=wp)     :: charge
    Real(Kind=wp)     :: conc
    Real(Kind=wp)     :: necs_shift
  End Type boltzmann_ions 

  ! Type for DFT settings
  Type :: dft_type
    ! Flag to ensure DFT block is not defined more than once 
    Logical  :: generate=.False.
    ! Flag to set spin polarised simulation 
    Type(in_logic)  :: gapw 
    ! Type XC functional
    Type(in_string)  :: xc_level     
    ! Type XC version 
    Type(in_string)  :: xc_version     
    ! Staring XC base approach 
    Character(Len=256) :: xc_base
    ! XC reference
    Character(Len=256) :: xc_ref
    ! vdW
    Type(in_string)    :: vdw   
    ! vdW reference
    Character(Len=256) :: vdw_ref
    ! vdW kernel
    Logical            :: need_vdw_kernel 
    Character(Len=256) :: vdw_kernel_file
    ! Flag to set spin polarised simulation 
    Type(in_logic)  :: spin_polarised
    ! Energy cutoff
    Type(in_param)  :: encut
    ! Precision
    Type(in_string) :: precision   
    ! Smearing
    Type(in_string) :: smear
    ! Width for smearing
    Type(in_param)  :: width_smear
    ! Mixing 
    Type(in_string) :: mixing
    ! SFC steps
    Type(in_integer) :: scf_steps
    ! Energy convergence
    Type(in_param)   :: delta_e
    ! k-point sampling
    Integer(Kind=wi)       :: total_kpoints 
    Type(in_integer_array) :: kpoints
    ! basis set
    Type(in_logic)   :: basis_info 
    Type(type_basis), Allocatable :: basis_set(:)
    ! pseudo-potentials
    Type(in_logic)   :: pp_info 
    Type(type_pseudo), Allocatable :: pseudo_pot(:)
    ! Maximum l_orbital
    Type(in_integer) :: max_l_orbital
    ! Total magnetization 
    Type(in_param)   :: total_magnetization 
    ! magnetization
    Type(in_logic)   :: mag_info 
    Type(type_mag),  Allocatable :: magnetization(:)
    ! hubbard
    Type(in_logic)   :: hubbard_info
    Logical          :: hubbard_all_U_zero
    Type(type_hubbard),  Allocatable :: hubbard(:)
    ! Orbital Transformation (OT), only valid for CP2K
    Type(in_logic)   :: ot
    ! Ensemble DFT (EDFT), only valid for CASTEP and ONETEP
    Type(in_logic)   :: edft 
    ! Bands paralellization, only for VASP
    Type(in_integer) :: npar
    ! kpoints paralellization, only for VASP
    Type(in_integer) :: kpar
    ! Bands
    Type(in_integer) :: bands
    ! NGWF, compulsory only for ONETEP
    Type(in_logic)   :: ngwf_info
    Type(type_ngwf),  Allocatable  :: ngwf(:)
    ! PAW for onetep
    Logical :: onetep_paw
    
    ! GC-DFT
    Type(type_gcdft) :: gc

  End Type

  ! Type for motion settings
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  Type :: motion_type
    ! Flag to ensure motion block is not defined more than once 
    Logical  :: generate=.False.
    ! Relaxation method 
    Type(in_string) :: relax_method
    ! force convergence
    Type(in_param_array) :: delta_f
    ! Time step
    Type(in_param) :: timestep
    ! Number of ionic step, either for relaxation or MD
    Type(in_integer) :: ion_steps
    ! Change simulation cell volume
    Type(in_logic)  :: change_cell_volume
    ! Change simulation cell shape 
    Type(in_logic)  :: change_cell_shape
    ! Ensemble 
    Type(in_string) :: ensemble
    ! Temperature 
    Type(in_param)  :: temperature 
    ! Pressure 
    Type(in_param)  :: pressure 
    ! Thermostat 
    Type(in_string) :: thermostat
    ! Thermostat relaxation time 
    Type(in_param)  :: relax_time_thermostat 
    ! Barostat 
    Type(in_string) :: barostat
    ! Barostat relaxation time 
    Type(in_param)  :: relax_time_barostat 
    ! Masses 
    Type(in_logic)   :: mass_info
    Type(species_list), Allocatable :: mass(:) 
  End Type

  ! Type for solvation 
  Type :: solvation_type
    Character(Len=256) :: bib_epsilon
    Type(in_logic)    :: info
    Type(in_string)   :: cavity_model
    Type(in_string)   :: dielectric_function
    Type(in_scalar)   :: density_threshold
    Type(in_scalar)   :: density_min_threshold
    Type(in_scalar)   :: density_max_threshold
    Type(in_scalar)   :: beta_fg_parameter 
    Type(in_scalar)   :: dielectric_constant
    Type(in_logic)    :: soft_radii_info  
    Type(in_scalar)   :: soft_sphere_scale    
    Type(in_scalar)   :: soft_sphere_delta    
    Type(species_list), Allocatable :: soft_radii(:)
    !Apolar terms
    Type(in_string)   :: apolar_terms
    Type(in_string)   :: sasa_definition
    Type(in_scalar)   :: apolar_scaling 
    Logical           :: both_surfaces
    Type(in_param)    :: smear_ion_width
    Type(in_param)    :: dispersive_pressure      
    Type(in_param_array) :: repulsion_parameter
    Type(in_param_array) :: surface_tension 
  End Type solvation_type

  ! Type for electrolyte
  Type :: electrolyte_type
    Type(in_logic)    :: info 
    ! Counter charge
    Type(in_logic)     :: info_pcc
    Type(in_param)     :: dist_edge
    Type(in_param)     :: gaussian_width
    Character(Len=256) :: plane_orientation
    ! Poisson-Boltzmann 
    Type(in_logic)    :: info_pb
    Type(in_string)   :: solver
    Type(in_string)   :: neutral_scheme 
    Type(in_string)   :: steric_potential 
    Type(in_param)    :: boltzmann_temp
    Type(in_param)    :: steric_isodensity 
    Type(in_param)    :: steric_smearing
    Type(in_param)    :: capping 
    Logical           :: set_necs_shift
    Type(in_logic)    :: boltzmann_ions_info
    Integer           :: number_boltzmann_ions
    Type(in_logic)    :: solvent_radii_info
    Type(species_list), Allocatable :: solvent_radii(:)
    Type(boltzmann_ions) :: boltzmann_ions(max_number_boltzmann_ions)
  End Type electrolyte_type

  ! Type for multi-grid 
  Type :: multigrid
    Type(in_logic)    :: info
  End Type multigrid

  ! Type for the modelling related variables 
  Type, Public :: simul_type
    Private
    ! General
    !!!!!!!!!
    ! Flag to generate simulation files
    Logical, Public  :: generate=.False.
    ! Details for the components
    Type(component_in_block), Public :: component(max_components)
    !number of total tags
    Integer(Kind=wi),   Public  :: total_tags
    ! Code 
    Character(Len=256), Public  :: code_format
    ! Code version 
    Character(Len=256), Public  :: code_version
    ! vector normal to the surface 
    Character(Len=256), Public  :: normal_vector
    ! Simulation cell
    Real(Kind=wp),      Public  :: cell(3,3)
    ! Length of cell vectors
    Real(Kind=wp),      Public  :: cell_length(3)
    ! Large cell
    Logical,            Public  :: large_cell
    
    ! Specific variables 
    !!!!!!!!!!!!!!!!!!!!
    ! Type of the simulation to be performed
    Type(in_string), Public :: simulation     
    ! Level of theory 
    Type(in_string), Public :: theory_level
    ! DFT directives
    Type(dft_type),    Public :: dft 
    ! Ions related variables
    Type(motion_type), Public :: motion
    ! Extra directives
    Type(in_logic),    Public :: extra_info
    Type(type_extra),  Public :: extra_directives
    ! Set directives
    Type(type_extra),  Public :: set_directives

    ! Solvation
    Type(solvation_type), Public :: solvation
    ! Poisson-Boltzmann 
    Type(electrolyte_type), Public :: electrolyte 

  Contains
    Private
    Procedure, Public  :: init_input_dft_variables    =>  allocate_input_dft_variables
    Procedure, Public  :: init_input_motion_variables =>  allocate_input_motion_variables
    Procedure, Public  :: init_input_solvation_variables =>  allocate_input_solvation_variables
    Procedure, Public  :: init_solvent_radii =>  allocate_solvent_radii
    Final              :: cleanup

  End Type simul_type

  Public :: read_simulation_settings
  
Contains

  Subroutine allocate_input_dft_variables(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate essential DFT variables to build input files for simulations 
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(simul_type), Intent(InOut)  :: T

    Integer(Kind=wi)    :: fail(6)
    Character(Len=256)  :: message

    Allocate(T%dft%pseudo_pot(T%total_tags),     Stat=fail(1))
    Allocate(T%dft%magnetization(T%total_tags),  Stat=fail(2))
    Allocate(T%dft%hubbard(T%total_tags),        Stat=fail(3))
    Allocate(T%dft%basis_set(T%total_tags),      Stat=fail(4))
    Allocate(T%dft%kpoints%value(3),             Stat=fail(5))
    Allocate(T%dft%ngwf(T%total_tags),           Stat=fail(6))

    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems of "DFT" variables to build input files&
                               & for simulations'
      Call error_stop(message)
    End If

    !Set to False just in case
    T%dft%basis_info%stat=.False. 
    T%dft%pp_info%stat=.False. 
    T%dft%mag_info%stat=.False.
    T%dft%hubbard_info%stat=.False.
    T%dft%ngwf_info%stat=.False.

  End Subroutine allocate_input_dft_variables

  Subroutine allocate_input_motion_variables(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate essential motion variables to build input files for simulations 
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(simul_type), Intent(InOut)  :: T

    Integer(Kind=wi)    :: fail(3)
    Character(Len=256)  :: message

    Allocate(T%motion%delta_f%value(1),               Stat=fail(1))
    Allocate(T%motion%delta_f%units(2),               Stat=fail(2))
    Allocate(T%motion%mass(T%total_tags),             Stat=fail(3))  

    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems of "motion" variables to build input files&
                               & for simulations'
      Call error_stop(message)
    End If

    !Set to False just in case
    T%motion%mass_info%stat=.False.

  End Subroutine allocate_input_motion_variables

  Subroutine allocate_input_solvation_variables(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate essential solvation input variable 
    !
    ! author    - i.scivetti  March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(simul_type), Intent(InOut)  :: T

    Integer(Kind=wi)    :: fail(5)
    Character(Len=256)  :: message

    Allocate(T%solvation%surface_tension%value(1),     Stat=fail(1))
    Allocate(T%solvation%surface_tension%units(2),     Stat=fail(2))
    Allocate(T%solvation%repulsion_parameter%value(1), Stat=fail(3))
    Allocate(T%solvation%repulsion_parameter%units(2), Stat=fail(4))    
    Allocate(T%solvation%soft_radii(T%total_tags),     Stat=fail(5))  

    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems of "solvation" variables to build input files&
                               & for simulations'
      Call error_stop(message)
    End If

    !Set to False just in case
    T%solvation%soft_radii_info%stat=.False.

  End Subroutine allocate_input_solvation_variables

  Subroutine allocate_solvent_radii(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate essential elctrolyte input variable 
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(simul_type), Intent(InOut)  :: T

    Integer(Kind=wi)    :: fail(1)
    Character(Len=256)  :: message

    Allocate(T%electrolyte%solvent_radii(T%total_tags), Stat=fail(1))  

    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems of "electrolyte" variables to build input files&
                               & for simulations'
      Call error_stop(message)
    End If

  End Subroutine allocate_solvent_radii

  Subroutine cleanup(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deallocate variables
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type) :: T

    If (Allocated(T%motion%delta_f%units)) Then
      Deallocate(T%motion%delta_f%units)
    End If

    If (Allocated(T%motion%delta_f%value)) Then
      Deallocate(T%motion%delta_f%value)
    End If

    If (Allocated(T%motion%mass)) Then
      Deallocate(T%motion%mass)
    End If


    If (Allocated(T%dft%kpoints%value)) Then
      Deallocate(T%dft%kpoints%value)
    End If

    If (Allocated(T%dft%pseudo_pot)) Then
      Deallocate(T%dft%pseudo_pot)
    End If

    If (Allocated(T%dft%hubbard)) Then
      Deallocate(T%dft%hubbard)
    End If

    If (Allocated(T%dft%basis_set)) Then
      Deallocate(T%dft%basis_set)
    End If

    If (Allocated(T%dft%ngwf)) Then
      Deallocate(T%dft%ngwf)
    End If

    If (Allocated(T%solvation%soft_radii)) Then
      Deallocate(T%solvation%soft_radii)
    End If
    
    If (Allocated(T%solvation%repulsion_parameter%value)) Then
      Deallocate(T%solvation%repulsion_parameter%value)
    End If    
    
    If (Allocated(T%solvation%repulsion_parameter%units)) Then
      Deallocate(T%solvation%repulsion_parameter%units)
    End If    
    
    If (Allocated(T%solvation%surface_tension%value)) Then
      Deallocate(T%solvation%surface_tension%value)
    End If     

    If (Allocated(T%solvation%surface_tension%units)) Then
      Deallocate(T%solvation%surface_tension%units)
    End If         

  End Subroutine cleanup

  Subroutine read_simulation_settings(iunit, simulation_data) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read simulation directives, which will be used to generate input 
    ! files, required for atomistic level simulation
    !
    ! author    - i.scivetti Feb 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: simulation_data

    Character(Len=256) :: word
    Integer(Kind=wi)   :: length, io
  
    Character(Len=256)  :: message
    Character(Len=265)  :: set_error

    set_error = '***ERROR in &simulation_settings -'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_simulation_settings" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_simulation_settings') Exit
      Call check_for_rubbish(iunit, '&simulation_settings')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == 'simulation_type') Then 
        Read (iunit, Fmt=*, iostat=io) word, simulation_data%simulation%type
        Call set_read_status(word, io, simulation_data%simulation%fread, simulation_data%simulation%fail, &
                           & simulation_data%simulation%type)

      Else If (word(1:length) == 'theory_level') Then 
        Read (iunit, Fmt=*, iostat=io) word, simulation_data%theory_level%type
        Call set_read_status(word, io, simulation_data%theory_level%fread, simulation_data%theory_level%fail, &
                           & simulation_data%theory_level%type)

      Else If (word(1:length) == '&dft_settings') Then
        Read (iunit, Fmt=*, iostat=io) word 
        If (simulation_data%dft%generate) Then
          Call duplication_error(word)
        End If
        simulation_data%dft%generate=.True.
        Call simulation_data%init_input_dft_variables()
        ! Now, it is ready to read information inside &dft_settings
        Call read_dft_settings(iunit, simulation_data)

      Else If (word(1:length) == '&motion_settings') Then
        Read (iunit, Fmt=*, iostat=io) word 
        If (simulation_data%motion%generate) Then
          Call duplication_error(word)
        End If
        simulation_data%motion%generate=.True.
        ! Now, it is ready to read information inside &dft_settings
        Call simulation_data%init_input_motion_variables()
        Call read_motion_settings(iunit, simulation_data)

      ! Solvation 
      Else If (word(1:length) == '&solvation') Then
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, simulation_data%solvation%info%fread, simulation_data%solvation%info%fail)
        simulation_data%solvation%info%stat = .True.
        Call simulation_data%init_input_solvation_variables()
        Call read_solvation_directives(iunit, simulation_data)

      ! Electrolyte
      Else If (word(1:length) == '&electrolyte') Then
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, simulation_data%electrolyte%info%fread, simulation_data%electrolyte%info%fail)
        simulation_data%electrolyte%info%stat = .True.
        Call read_electrolyte(iunit, simulation_data)
        
      ! Extra directives
      Else If (word(1:length) == '&extra_directives') Then
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, simulation_data%extra_info%fread, simulation_data%extra_info%fail)
        simulation_data%extra_info%stat = .True.
        ! Read extra input 
        Call read_extra_directives(iunit, simulation_data)

      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid simulation settings.',&
                                & ' See manual. Have you properly closed the block with "&end_simulation_settings"?'
        Call error_stop(message)
      End If

    End Do
  End Subroutine read_simulation_settings

  Subroutine read_extra_directives(iunit, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read extra directives defined by the user
    !
    ! author    - i.scivetti Feb 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: simulation_data

    Character(Len=256)  :: word
    Integer(Kind=wi)    :: i, io

    Character(Len=256)  :: message
    Character(Len=265)  :: set_error

    set_error = '***ERROR in &modules (inside &extra_directives):'

    i=1
    Do 
      Read (iunit, Fmt='(a)', iostat=io) word
      If (is_iostat_end(io)) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'End of file? Have you closed the block with "&end_extra_directives"?'
        Call error_stop(message)
      End If

      word=Trim(Adjustl(word))
      Call capital_to_lower_case(word) 

      If (Trim(word) == '&end_extra_directives') Then
        Exit
      Else
        If (word(1:1)=='&') Then
          If (word(1:7)/= '&block ' .And. word(1:10)/='&endblock ') Then
            Write (message,'(2a)') Trim(set_error), ' It appears that &extra_directives has not been closed properly. Please use&
                               & "&end_extra_directives".'
            Call error_stop(message) 
          Else
            Backspace iunit 
            Read (iunit, Fmt='(a)', iostat=io) word 
            simulation_data%extra_directives%array(i)=Trim(Adjustl(word))
            i=i+1
          End If
        Else 
          Backspace iunit 
          Read (iunit, Fmt='(a)', iostat=io) word 
          simulation_data%extra_directives%array(i)=Trim(Adjustl(word))
          i=i+1
        End If  
      End If
    End Do
    simulation_data%extra_directives%N0=i-1

  End Subroutine read_extra_directives

  Subroutine read_electrolyte(iunit, T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read directives from &electrolyte
    !
    ! author    - i.scivetti Feb 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: T

    Character(Len=256) :: word
    Integer(Kind=wi)   :: length, io
 
    Character(Len=256)  :: message
    Character(Len=265)  :: set_error

    set_error = '***ERROR in &electrolyte -'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_electrolyte" to close the block. Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_electrolyte') Then
        Exit
      End If
      Call check_for_rubbish(iunit, '&electrolyte')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == '&poisson_boltzmann') Then 
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, T%electrolyte%info_pb%fread, T%electrolyte%info_pb%fail)
        T%electrolyte%info_pb%stat = .True.
        Call read_electrolyte_PB(iunit, T)

      Else If (word(1:length) == '&planar_counter_charge') Then 
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, T%electrolyte%info_pcc%fread, T%electrolyte%info_pcc%fail)
        T%electrolyte%info_pcc%stat = .True.
        Call read_electrolyte_PCC(iunit, T)
        
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word), '" is not recognised as a valid motion settings.',&
                                & ' See manual. Have you properly closed the block with "&end_electrolyte"?'
        Call error_stop(message)
      End If

     End Do

  End Subroutine read_electrolyte 

  Subroutine read_electrolyte_PCC(iunit, T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read directives from &planar_counter_charge
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: T

    Character(Len=256) :: word
    Integer(Kind=wi)   :: length, io
 
    Character(Len=256)  :: message
    Character(Len=265)  :: set_error

    set_error = '***ERROR in &planar_counter_charge (within &electrolyte) -'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_planar_counter_charge" to close the block. Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_planar_counter_charge') Then
        Exit
      End If
      Call check_for_rubbish(iunit, '&planar_counter_charge')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == 'distance_to_edge') Then
        Read (iunit, Fmt=*, iostat=io) word, T%electrolyte%dist_edge%value, T%electrolyte%dist_edge%units
        Call set_read_status(word, io, T%electrolyte%dist_edge%fread, T%electrolyte%dist_edge%fail)
        Call capital_to_lower_case(T%electrolyte%dist_edge%units)

      Else If (word(1:length) == 'gaussian_width') Then
        Read (iunit, Fmt=*, iostat=io) word, T%electrolyte%gaussian_width%value, T%electrolyte%gaussian_width%units
        Call set_read_status(word, io, T%electrolyte%gaussian_width%fread, T%electrolyte%gaussian_width%fail)
        Call capital_to_lower_case(T%electrolyte%gaussian_width%units)
        
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word), '" is not recognised as a valid motion settings.',&
                                & ' See manual. Have you properly closed the block with "&end_planar_counter_charge"?'
        Call error_stop(message)
      End If

     End Do

  End Subroutine read_electrolyte_PCC   
  
  Subroutine read_electrolyte_PB(iunit, T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read directives from &poisson_boltzmann
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: T

    Character(Len=256) :: word
    Integer(Kind=wi)   :: length, io
 
    Character(Len=256)  :: message
    Character(Len=265)  :: set_error

    set_error = '***ERROR in &poisson_boltzmann (within &electrolyte) -'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_poisson_boltzmann" to close the block. Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_poisson_boltzmann') Then
        Exit
      End If
      Call check_for_rubbish(iunit, '&poisson_boltzmann')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == 'solver') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%electrolyte%solver%type
        Call set_read_status(word, io, T%electrolyte%solver%fread, T%electrolyte%solver%fail, T%electrolyte%solver%type)

      Else If (word(1:length) == 'boltzmann_temperature') Then
        Read (iunit, Fmt=*, iostat=io) word, T%electrolyte%boltzmann_temp%value, T%electrolyte%boltzmann_temp%units
        Call set_read_status(word, io, T%electrolyte%boltzmann_temp%fread, T%electrolyte%boltzmann_temp%fail)
        Call capital_to_lower_case(T%electrolyte%boltzmann_temp%units)

      Else If (word(1:length) == 'steric_isodensity') Then
        Read (iunit, Fmt=*, iostat=io) word, T%electrolyte%steric_isodensity%value
        Call set_read_status(word, io, T%electrolyte%steric_isodensity%fread, T%electrolyte%steric_isodensity%fail)
        Call capital_to_lower_case(T%electrolyte%steric_isodensity%units)

      Else If (word(1:length) == 'steric_smearing') Then
        Read (iunit, Fmt=*, iostat=io) word, T%electrolyte%steric_smearing%value, T%electrolyte%steric_smearing%units
        Call set_read_status(word, io, T%electrolyte%steric_smearing%fread, T%electrolyte%steric_smearing%fail)
        Call capital_to_lower_case(T%electrolyte%steric_smearing%units)
       
      Else If (word(1:length) == 'capping') Then
        Read (iunit, Fmt=*, iostat=io) word, T%electrolyte%capping%value
        Call set_read_status(word, io, T%electrolyte%capping%fread, T%electrolyte%capping%fail)
        Call capital_to_lower_case(T%electrolyte%capping%units)

      Else If (word(1:length) == 'neutral_scheme') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%electrolyte%neutral_scheme%type
        Call set_read_status(word, io, T%electrolyte%neutral_scheme%fread, T%electrolyte%neutral_scheme%fail,&
        T%electrolyte%neutral_scheme%type)

      Else If (word(1:length) == 'steric_potential') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%electrolyte%steric_potential%type
        Call set_read_status(word, io, T%electrolyte%steric_potential%fread, T%electrolyte%steric_potential%fail,&
        T%electrolyte%steric_potential%type)
        
      Else If (word(1:length) == '&boltzmann_ions') Then 
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, T%electrolyte%boltzmann_ions_info%fread, T%electrolyte%boltzmann_ions_info%fail)
        T%electrolyte%boltzmann_ions_info%stat = .True.
        Call read_boltzmann_ions(iunit, T)

      Else If (word(1:length) == '&solvent_radii') Then
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, T%electrolyte%solvent_radii_info%fread, T%electrolyte%solvent_radii_info%fail)
        Call T%init_solvent_radii()
        T%electrolyte%solvent_radii_info%stat = .True.        
        Call read_species_list(iunit, T%total_tags, T%electrolyte%solvent_radii, 'solvent_radii', '&poisson_boltzmann')
        
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word), '" is not recognised as a valid motion settings.',&
                                & ' See manual. Have you properly closed the block with "&end_poisson_boltzmann"?'
        Call error_stop(message)
      End If

     End Do

  End Subroutine read_electrolyte_PB  
  
  Subroutine read_boltzmann_ions(iunit, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the specification for boltzmann ions 
    !
    ! author    - i.scivetti Feb 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: simulation_data

    Integer(Kind=wi) :: io, i, j, k, nions
    Character(Len=256)  :: messages(11), word
    Character(Len=64 )  :: set_error
    Logical  :: endblock, pass
    Logical  :: header(4), error_duplication
    Character(Len=64 )  :: field_ref(3)
    

    set_error = '***ERROR in &boltzmann_ions (inside &poisson_boltzmann):'
    Write (messages(1),'(a)') set_error

    header=.False.
    simulation_data%electrolyte%set_necs_shift=.False.
    error_duplication=.False.
    field_ref(1)='tags'
    field_ref(2)='charge'
    field_ref(3)='conc'
        
    Write (messages(2),'(1x,a)')    'The correct structure must be:'
    Write (messages(3),'(1x,a)')    '&boltzmann_ions'
    Write (messages(4),'(1x,a)')    '  number_boltzmann_ions  Nion'
    Write (messages(5),'(1x,a)')    '  tags         tg_1    tg_2    ....  tg_Nions  (Compulsory)'
    Write (messages(6),'(1x,a)')    '  charge       Q_tg1   Q_tg2   .... Q_tgNions  (Compulsory)'
    Write (messages(7),'(1x,a)')    '  conc         C_tg1   C_tg2   .... C_tgNions  (Compulsory)'
    Write (messages(8),'(1x,a)')    '  necs_shift   x_tg1   x_tg2   .... x_tgNions  (only needed is&
                                    & "neutral_scheme" is set to counterions_fixed)'
    Write (messages(9),'(1x,a)')    '&end_boltzmann_ions'
    Write (messages(10),'(1x,a)')   'where in this case Nions is the value defined by "number_boltzmann_ions"' 
    Write (messages(11),'(1x,a)')   'See manual for details'
    
    pass=.True.
    Do While (pass)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&boltzmann_ions (inside &poisson_boltzmann)')
      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          Call check_for_rubbish(iunit, '&boltzmann_ions (inside &poisson_boltzmann)')
          Call capital_to_lower_case(word)
          If (Trim(word) /= 'number_boltzmann_ions') Then
            Write (messages(2),'(1x,a)') 'The first directive within &boltzmann_ions must be "number_boltzmann_ions"'
            Call info(messages, 2)
            Call error_stop(' ') 
          Else
            Read (iunit, Fmt=*, iostat=io) word, simulation_data%electrolyte%number_boltzmann_ions
            If (io /= 0 .Or. simulation_data%electrolyte%number_boltzmann_ions < 2) Then
              Write (messages(2),'(1x,a)') 'Problems to read directives "number_boltzmann_ions"'
              Call info(messages, 2)
              Call error_stop(' ') 
            End If
            If (simulation_data%electrolyte%number_boltzmann_ions > max_number_boltzmann_ions) Then
              Write (messages(2),'(1x,a,i2)') 'Directives "number_boltzmann_ions" is larger than the max_intra allowed&
                                          & number of Boltzmann ions = ', max_number_boltzmann_ions
              Call info(messages, 2)
              Call error_stop(' ') 
            End If
            pass=.False.
          End If
        Else
           Call info(messages, 11)
           Call error_stop(' ') 
        End If
      End If  
    End Do

    nions=simulation_data%electrolyte%number_boltzmann_ions    

    i=1
    pass=.True.
    Do While (i <= 4 .And. pass)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&boltzmann_ions (inside &poisson_boltzmann)')
      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          Call check_for_rubbish(iunit, '&boltzmann_ions (inside &poisson_boltzmann)')
          Call capital_to_lower_case(word) 

          If (Trim(word)=='tags') Then
            If (.Not. header(1)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%electrolyte%boltzmann_ions(j)%tag, j = 1, nions)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read tags for atoms'
                Call info(messages, 2)
                Call error_stop(' ')
              End If
              Do j=1, nions-1
                 Do k=j+1, nions
                   If (Trim(simulation_data%electrolyte%boltzmann_ions(j)%tag)==&
                       Trim(simulation_data%electrolyte%boltzmann_ions(k)%tag)) Then
                     Write (messages(2),'(3(1x,a))') 'Tag', Trim(simulation_data%electrolyte%boltzmann_ions(j)%tag),&
                                                    & 'is repeated in the list!'
                     Write (messages(3),'((1x,a))') 'All tags for the components of the species must be declared,&
                                                    & each tag only once'
                     Call info(messages, 3)
                     Call error_stop(' ')
                   End If
                 End Do
              End Do  
              header(1)=.True.
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='charge') Then
            If (.Not. header(2)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%electrolyte%boltzmann_ions(j)%charge, j = 1, nions)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the "charge" values for tags. Please check'
                Call info(messages, 2)
                Call error_stop(' ')
              End If  
              header(2)=.True.        
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='conc') Then
            If (.Not. header(3)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%electrolyte%boltzmann_ions(j)%conc, j = 1, nions)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the "conc" (which orbital to apply the correction)&
                                           & for each tags. Please check'
                Call info(messages, 2)
                Call error_stop(' ')
              End If  
              header(3)=.True.        
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='necs_shift') Then
            If (.Not. header(4)) Then 
              simulation_data%electrolyte%set_necs_shift=.True.
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%electrolyte%boltzmann_ions(j)%necs_shift, j = 1, nions)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the "necs_shift" values for tags. Please check'
                Call info(messages, 2)
                Call error_stop(' ')
              End If  
              header(4)=.True.        
            Else
              error_duplication=.True.
            End If
          Else
              Write (messages(2),'(1x,3a)') 'Wrong descriptor "', Trim(word),&
                                         &'". Valid options are: "tags", "charge", "conc" and "necs_shift".&
                                         & Please refer to the manual.'
              Call info(messages, 2)
              Call error_stop(' ')
          End If
        Else

          If (i<4) Then
            Call info(messages, 11)
            Call error_stop(' ')
          Else
            i=i+1
          End If

          If (Trim(word)=='&end_boltzmann_ions') Then
            Backspace iunit
            pass=.False.
          End If
          
        End If
      End If
      If (error_duplication) Then
        Write (messages(2),'(1x,3a)') 'Descriptor "', Trim(word), '" is duplicated within &boltzmann_ions'
        Call info(messages, 2)
        Call error_stop(' ')
      End If
    End Do 

    Do i=1,3
     If(.Not. header(i)) Then
       Write (messages(2),'(3a)') 'Compulsory field "', Trim(field_ref(i)), '" is missing.' 
       Call info(messages,2)
       Call error_stop(' ')
     End If 
    End Do
    
    endblock=.False.
    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&boltzmann_ions (inside &dft_settings)')
      Call capital_to_lower_case(word)
      If (word /= '&end_boltzmann_ions') Then
        If (word(1:1) /= '#') Then
            Write (messages(2),'(3a)') 'Valid descriptors have already been defined. Directive "',&
                                    & Trim(word), '" is not valid. Block must be&
                                    & closed with sentence &end_boltzmann_ions' 
            Call info(messages,2)
            Call error_stop(' ')
        End If
      Else
        endblock=.True.
      End If
    End Do

  End Subroutine read_boltzmann_ions
  

  Subroutine read_solvation_directives(iunit, T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read directives from &solvation
    !
    ! author    - i.scivetti Feb 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: T

    Character(Len=256) :: word
    Integer(Kind=wi)   :: length, io
 
    Integer             :: j
    Character(Len=256)  :: message
    Character(Len=265)  :: set_error

    set_error = '***ERROR in &solvation -'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_solvation" to close the block. Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_solvation') Then
        Exit
      End If
      Call check_for_rubbish(iunit, '&solvation')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == 'cavity_model') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%solvation%cavity_model%type
        Call set_read_status(word, io, T%solvation%cavity_model%fread, T%solvation%cavity_model%fail,&
                             T%solvation%cavity_model%type)

      Else If (word(1:length) == 'dielectric_function') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%solvation%dielectric_function%type
        Call set_read_status(word, io, T%solvation%dielectric_function%fread, T%solvation%dielectric_function%fail,&
                              T%solvation%dielectric_function%type)

       Else If (word(1:length) == 'density_threshold') Then 
         Read (iunit, Fmt=*, iostat=io) word, T%solvation%density_threshold%value
         Call set_read_status(word, io, T%solvation%density_threshold%fread, T%solvation%density_threshold%fail)
                              
       Else If (word(1:length) == 'beta_fg_parameter') Then 
         Read (iunit, Fmt=*, iostat=io) word, T%solvation%beta_fg_parameter%value
         Call set_read_status(word, io, T%solvation%beta_fg_parameter%fread, T%solvation%beta_fg_parameter%fail)
                              
       Else If (word(1:length) == 'dielectric_constant') Then 
         Read (iunit, Fmt=*, iostat=io) word, T%solvation%dielectric_constant%value
         Call set_read_status(word, io, T%solvation%dielectric_constant%fread, T%solvation%dielectric_constant%fail)

       Else If (word(1:length) == 'density_min_threshold') Then 
         Read (iunit, Fmt=*, iostat=io) word, T%solvation%density_min_threshold%value
         Call set_read_status(word, io, T%solvation%density_min_threshold%fread, T%solvation%density_min_threshold%fail)

       Else If (word(1:length) == 'density_max_threshold') Then 
         Read (iunit, Fmt=*, iostat=io) word, T%solvation%density_max_threshold%value
         Call set_read_status(word, io, T%solvation%density_max_threshold%fread, T%solvation%density_max_threshold%fail)

       Else If (word(1:length) == 'solvent_surface_tension') Then
         Read (iunit, Fmt=*, iostat=io) word, T%solvation%surface_tension%value(1),&
                                      & (T%solvation%surface_tension%units(j), j=1,2)
         Call set_read_status(word, io, T%solvation%surface_tension%fread, T%solvation%surface_tension%fail)
         Call capital_to_lower_case(T%solvation%surface_tension%units(1))
         Call capital_to_lower_case(T%solvation%surface_tension%units(2))
       
       Else If (word(1:length) == 'solvent_dispersive_pressure') Then
         Read (iunit, Fmt=*, iostat=io) word, T%solvation%dispersive_pressure%value, T%solvation%dispersive_pressure%units
         Call set_read_status(word, io, T%solvation%dispersive_pressure%fread, T%solvation%dispersive_pressure%fail)
         Call capital_to_lower_case(T%solvation%dispersive_pressure%units)

       ! Pure ONETEP directives  
       Else If (word(1:length) == 'smear_ion_width') Then 
         Read (iunit, Fmt=*, iostat=io) word, T%solvation%smear_ion_width%value, T%solvation%smear_ion_width%units
         Call set_read_status(word, io, T%solvation%smear_ion_width%fread, T%solvation%smear_ion_width%fail)
         Call capital_to_lower_case(T%solvation%smear_ion_width%units)

       Else If (word(1:length) == 'soft_sphere_scale') Then 
         Read (iunit, Fmt=*, iostat=io) word, T%solvation%soft_sphere_scale%value
         Call set_read_status(word, io, T%solvation%soft_sphere_scale%fread, T%solvation%soft_sphere_scale%fail)

       Else If (word(1:length) == 'soft_sphere_delta') Then 
         Read (iunit, Fmt=*, iostat=io) word, T%solvation%soft_sphere_delta%value
         Call set_read_status(word, io, T%solvation%soft_sphere_delta%fread, T%solvation%soft_sphere_delta%fail)
         
       Else If (word(1:length) == '&soft_sphere_radii') Then
         Read (iunit, Fmt=*, iostat=io) word
         Call set_read_status(word, io, T%solvation%soft_radii_info%fread, T%solvation%soft_radii_info%fail)
         T%solvation%soft_radii_info%stat = .True.
         Call read_species_list(iunit, T%total_tags, T%solvation%soft_radii, 'soft_sphere_radii', '&solvation')

       Else If (word(1:length) == 'apolar_terms') Then 
         Read (iunit, Fmt=*, iostat=io) word, T%solvation%apolar_terms%type
         Call set_read_status(word, io, T%solvation%apolar_terms%fread, T%solvation%apolar_terms%fail,&
                              T%solvation%apolar_terms%type)

       Else If (word(1:length) == 'sasa_definition') Then 
         Read (iunit, Fmt=*, iostat=io) word, T%solvation%sasa_definition%type
         Call set_read_status(word, io, T%solvation%sasa_definition%fread, T%solvation%sasa_definition%fail,&
                             T%solvation%sasa_definition%type)
         
       Else If (word(1:length) == 'apolar_scaling') Then 
         Read (iunit, Fmt=*, iostat=io) word, T%solvation%apolar_scaling%value
         Call set_read_status(word, io, T%solvation%apolar_scaling%fread, T%solvation%apolar_scaling%fail)
       
       ! Pure CP2K directives
       Else If (word(1:length) == 'solvent_repulsion_parameter') Then
         Read (iunit, Fmt=*, iostat=io) word, T%solvation%repulsion_parameter%value(1),&
                                      & (T%solvation%repulsion_parameter%units(j), j=1,2)
         Call set_read_status(word, io, T%solvation%repulsion_parameter%fread, T%solvation%repulsion_parameter%fail)
         Call capital_to_lower_case(T%solvation%repulsion_parameter%units(1))
         Call capital_to_lower_case(T%solvation%repulsion_parameter%units(2))
         
        Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word), '" is not recognised as a valid motion settings.',&
                                & ' See manual. Have you properly closed the block with "&end_solvation"?'
        Call error_stop(message)
      End If

     End Do

  End Subroutine read_solvation_directives 
  
  Subroutine read_species_list(iunit, total_tags, quantity, sub_block, inblock)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read values for species in a block. The information read and stored
    ! depends on the input "sub_block". 
    !
    ! This subroutine generalises and replaces read_mass
    !
    ! author    - i.scivetti  March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),   Intent(In   ) :: iunit
    Integer(Kind=wi),   Intent(In   ) :: total_tags
    Type(species_list), Intent(InOut) :: quantity(total_tags)
    Character(Len=*),   Intent(In   ) :: sub_block
    Character(Len=*),   Intent(In   ) :: inblock
    

    Integer(Kind=wi) :: io, i, j, k
    Character(Len=256)  :: messages(8), word
    Character(Len=64 )  :: set_error
    Logical  :: endblock
    Logical  :: header(2), error_duplication
    
    set_error = '***ERROR in sub-block &'//Trim(sub_block)//' (inside '//Trim(inblock)//'):'
    Write (messages(1),'(a)') set_error

    header=.False.
    error_duplication=.False.

    i=1

    Do While (i <= 2)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&'//Trim(sub_block)//' (inside '//Trim(inblock)//')')
      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          Call check_for_rubbish(iunit, '&'//Trim(sub_block)//' (inside '//Trim(inblock)//')')
          Call capital_to_lower_case(word) 
          If (Trim(word)=='tags') Then
            If (.Not. header(1)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (quantity(j)%tag, j = 1, total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read tags for atoms'
                Call info(messages, 2)
                Call error_stop(' ')
              End If
              Do j=1, total_tags-1
                 Do k=j+1, total_tags
                   If (Trim(quantity(j)%tag)==Trim(quantity(k)%tag)) Then
                     Write (messages(2),'(3(1x,a))') 'Tag', Trim(quantity(j)%tag),&
                                                  & 'is repeated in the list!'
                     Write (messages(3),'((1x,a))') 'All tags for the components of the species must be declared,&
                                                   & each tag only once'
                     Call info(messages, 3)
                     Call error_stop(' ')
                   End If
                 End Do
              End Do  
              header(1)=.True.
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='values') Then
            If (.Not. header(2)) Then 
              i=i+1
              Read (iunit,Fmt=*,iostat=io) word, (quantity(j)%value, j=1, total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read values for the defined tags. Please check if tags&
                                            & are consistent with those defined in "&Species_Components".'
                Call info(messages, 2)
                Call error_stop(' ')
              End If  
              header(2)=.True.        
            Else
              error_duplication=.True.
            End If
          Else
            Write (messages(2),'(1x,3a)') 'Wrong descriptor "', Trim(word),'". Valid options are "tags" and "values"&
                                         & Please refer to the manual.'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        Else
          Write (messages(2),'(1x,a)')      'The correct structure must be:'
          Write (messages(3),'(1x,a)')      '&'//Trim(sub_block)
          Write (messages(4),'(1x,a)')      '  tags      tg1        tg2      tg3    ....   tgNsp'
          If (Trim(sub_block)=='masses') Then
            Write (messages(5),'(1x,a)')      '  values  mass_tg1  mass_tg2  mass_tg3 .... mass_tgNsp'
          Else
            Write (messages(5),'(1x,a)')      '  values  R_tg1      R_tg2    R_tg3  ....   R_tgNsp'
          End If  
          Write (messages(6),'(1x,a)')      '&end_'//Trim(sub_block)
          Write (messages(7),'(1x,a,i3,a)') 'where in this case Nsp = ', total_tags,&
                                          & ', which corresponds to the number of tags defined in "&Species_Components".'
          Write (messages(8),'(1x,a)')      'See manual for details'
          Call info(messages, 8)
          Call error_stop(' ')
        End If
      End If
      If (error_duplication) Then
        Write (messages(2),'(1x,3a)') 'Descriptor "', Trim(word), '" is duplicated within &'//Trim(sub_block)
        Call info(messages, 2)
        Call error_stop(' ')
      End If
    End Do 

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&'//Trim(sub_block)//' (inside '//Trim(inblock)//')')
      Call capital_to_lower_case(word)
      If (word /= Trim('&end_'//Trim(sub_block))) Then
        If (word(1:1) /= '#') Then
            Write (messages(2),'(3a)') 'Descriptors "tags" and "values" have already been defined. Directive "',&
                                    & Trim(word), '" is not valid. Block must be&
                                    & closed with sentence &end_'//Trim(sub_block) 
            Call info(messages,2)
            Call error_stop(' ')
        End If
      Else
        endblock=.True.
      End If
    End Do

  End Subroutine read_species_list
  
  
  Subroutine read_motion_settings(iunit, T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read directives from &motion_settings
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: T

    Character(Len=256) :: word
    Integer(Kind=wi)   :: length, io
 
    Integer             :: j
    Character(Len=256)  :: message
    Character(Len=265)  :: set_error

    set_error = '***ERROR in &motion_settings -'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_motion_settings" to close the block. Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_motion_settings') Then
        Exit
      End If
      Call check_for_rubbish(iunit, '&motion_settings')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == 'relax_method') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%motion%relax_method%type
        Call set_read_status(word, io, T%motion%relax_method%fread, T%motion%relax_method%fail, T%motion%relax_method%type)

      Else If (word(1:length) == 'ensemble') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%motion%ensemble%type
        Call set_read_status(word, io, T%motion%ensemble%fread, T%motion%ensemble%fail, T%motion%ensemble%type)

      Else If (word(1:length) == 'change_cell_volume') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%motion%change_cell_volume%stat
        Call set_read_status(word, io, T%motion%change_cell_volume%fread, T%motion%change_cell_volume%fail)

      Else If (word(1:length) == 'change_cell_shape') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%motion%change_cell_shape%stat
        Call set_read_status(word, io, T%motion%change_cell_shape%fread, T%motion%change_cell_shape%fail)

      Else If (word(1:length) == 'force_tolerance') Then
        Read (iunit, Fmt=*, iostat=io) word, T%motion%delta_f%value(1), (T%motion%delta_f%units(j), j=1,2)
        Call set_read_status(word, io, T%motion%delta_f%fread, T%motion%delta_f%fail)
        Call capital_to_lower_case(T%motion%delta_f%units(1))
        Call capital_to_lower_case(T%motion%delta_f%units(2))

      Else If (word(1:length) == 'timestep') Then
        Read (iunit, Fmt=*, iostat=io) word, T%motion%timestep%value, T%motion%timestep%units
        Call set_read_status(word, io, T%motion%timestep%fread, T%motion%timestep%fail)
        Call capital_to_lower_case(T%motion%timestep%units)

      Else If (word(1:length) == 'ion_steps') Then
        Read (iunit, Fmt=*, iostat=io) word, T%motion%ion_steps%value
        Call set_read_status(word, io, T%motion%ion_steps%fread, T%motion%ion_steps%fail)

      Else If (word(1:length) == 'pressure') Then
        Read (iunit, Fmt=*, iostat=io) word, T%motion%pressure%value, T%motion%pressure%units
        Call set_read_status(word, io, T%motion%pressure%fread, T%motion%pressure%fail)
        Call capital_to_lower_case(T%motion%pressure%units)

      Else If (word(1:length) == 'temperature') Then
        Read (iunit, Fmt=*, iostat=io) word, T%motion%temperature%value, T%motion%temperature%units
        Call set_read_status(word, io, T%motion%temperature%fread, T%motion%temperature%fail)
        Call capital_to_lower_case(T%motion%temperature%units)

      Else If (word(1:length) == 'thermostat') Then
        Read (iunit, Fmt=*, iostat=io) word, T%motion%thermostat%type
        Call set_read_status(word, io, T%motion%thermostat%fread, T%motion%thermostat%fail, T%motion%thermostat%type)

      Else If (word(1:length) == 'relax_time_thermostat') Then
        Read (iunit, Fmt=*, iostat=io) word, T%motion%relax_time_thermostat%value, T%motion%relax_time_thermostat%units
        Call set_read_status(word, io, T%motion%relax_time_thermostat%fread, T%motion%relax_time_thermostat%fail)
        Call capital_to_lower_case(T%motion%relax_time_thermostat%units)

      Else If (word(1:length) == 'barostat') Then
        Read (iunit, Fmt=*, iostat=io) word, T%motion%barostat%type
        Call set_read_status(word, io, T%motion%barostat%fread, T%motion%barostat%fail, T%motion%barostat%type)

      Else If (word(1:length) == 'relax_time_barostat') Then
        Read (iunit, Fmt=*, iostat=io) word, T%motion%relax_time_barostat%value, T%motion%relax_time_barostat%units
        Call set_read_status(word, io, T%motion%relax_time_barostat%fread, T%motion%relax_time_barostat%fail)
        Call capital_to_lower_case(T%motion%relax_time_barostat%units)

      ! Masses
      Else If (word(1:length) == '&masses') Then
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, T%motion%mass_info%fread, T%motion%mass_info%fail)
        T%motion%mass_info%stat = .True.
        !Read masses from block &mass
        Call read_species_list(iunit, T%total_tags, T%motion%mass, 'masses', '&motion')

      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word), '" is not recognised as a valid motion settings.',&
                                & ' See manual. Have you properly closed the block with "&end_motion_settings"?'
        Call error_stop(message)
      End If

     End Do

  End Subroutine read_motion_settings


  Subroutine read_dft_settings(iunit, T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read DFT directives from &DFT_settings
    !
    ! author         - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: T 

    Character(Len=256) :: word
    Integer(Kind=wi)   :: length, io
    Integer(Kind=wi)   :: i
  
    Character(Len=256) :: message
    Character(Len=265) :: set_error

    set_error = '***ERROR in &dft_settings -'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&dft_settings')
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_dft_settings" to close the block. Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_dft_settings') Then
        Exit
      End If
      Call check_for_rubbish(iunit, '&dft_settings')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == 'xc_level') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%dft%xc_level%type
        Call set_read_status(word, io, T%dft%xc_level%fread, T%dft%xc_level%fail, T%dft%xc_level%type)

      Else If (word(1:length) == 'xc_version') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%dft%xc_version%type
        Call set_read_status(word, io, T%dft%xc_version%fread, T%dft%xc_version%fail, T%dft%xc_version%type)

      Else If (word(1:length) == 'vdw') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%vdw%type
        Call set_read_status(word, io, T%dft%vdw%fread, T%dft%vdw%fail, T%dft%vdw%type)

      Else If (word(1:length) == 'spin_polarised') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%dft%spin_polarised%stat
        Call set_read_status(word, io, T%dft%spin_polarised%fread, T%dft%spin_polarised%fail)

      Else If (word(1:length) == 'gapw') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%dft%gapw%stat
        Call set_read_status(word, io, T%dft%gapw%fread, T%dft%gapw%fail)

      Else If (word(1:length) == 'energy_cutoff') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%encut%value, T%dft%encut%units
        Call set_read_status(word, io, T%dft%encut%fread, T%dft%encut%fail)
        Call capital_to_lower_case(T%dft%encut%units)
      
      Else If (word(1:length) == 'precision') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%dft%precision%type
        Call set_read_status(word, io, T%dft%precision%fread, T%dft%precision%fail, T%dft%precision%type)

      Else If (word(1:length) == 'smearing') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%dft%smear%type
        Call set_read_status(word, io, T%dft%smear%fread, T%dft%smear%fail, T%dft%smear%type)

      Else If (word(1:length) == 'width_smear') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%width_smear%value, T%dft%width_smear%units
        Call set_read_status(word, io, T%dft%width_smear%fread, T%dft%width_smear%fail)
        Call capital_to_lower_case(T%dft%width_smear%units)

      Else If (word(1:length) == 'mixing_scheme') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%dft%mixing%type
        Call set_read_status(word, io, T%dft%mixing%fread, T%dft%mixing%fail, T%dft%mixing%type)
      
      Else If (word(1:length) == 'scf_energy_tolerance') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%delta_e%value, T%dft%delta_e%units
        Call set_read_status(word, io, T%dft%delta_e%fread, T%dft%delta_e%fail)
        Call capital_to_lower_case(T%dft%delta_e%units)

      Else If (word(1:length) == 'scf_steps') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%scf_steps%value
        Call set_read_status(word, io, T%dft%scf_steps%fread, T%dft%scf_steps%fail)

      Else If (word(1:length) == 'kpoints') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%dft%kpoints%tag, (T%dft%kpoints%value(i), i=1,3)
        Call set_read_status(word, io, T%dft%kpoints%fread, T%dft%kpoints%fail, T%dft%kpoints%tag)

      Else If (word(1:length) == 'max_l_orbital') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%max_l_orbital%value
        Call set_read_status(word, io, T%dft%max_l_orbital%fread, T%dft%max_l_orbital%fail)

      Else If (word(1:length) == 'total_magnetization') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%total_magnetization%value
        Call set_read_status(word, io, T%dft%total_magnetization%fread, T%dft%total_magnetization%fail)

      Else If (word(1:length) == 'ot') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%ot%stat
        Call set_read_status(word, io, T%dft%ot%fread, T%dft%ot%fail)

      Else If (word(1:length) == 'edft') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%edft%stat
        Call set_read_status(word, io, T%dft%edft%fread, T%dft%edft%fail)

      Else If (word(1:length) == 'npar') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%npar%value
        Call set_read_status(word, io, T%dft%npar%fread, T%dft%npar%fail)

      Else If (word(1:length) == 'kpar') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%kpar%value
        Call set_read_status(word, io, T%dft%kpar%fread, T%dft%kpar%fail)

      Else If (word(1:length) == 'bands') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%bands%value
        Call set_read_status(word, io, T%dft%bands%fread, T%dft%bands%fail)

      ! Pseudopotentials
      Else If (word(1:length) == '&pseudo_potentials') Then 
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, T%dft%pp_info%fread, T%dft%pp_info%fail)
        T%dft%pp_info%stat=.True.
        ! Read names of the pseudo potentials 
        Call read_pseudo_poptentials(iunit, T)

      ! Basis set
      Else If (word(1:length) == '&basis_set') Then
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, T%dft%basis_info%fread, T%dft%basis_info%fail)
        T%dft%basis_info%stat=.True.
        ! Read type of basis 
        Call read_basis_set(iunit, T)

      ! Magnetization    
      Else If (word(1:length) == '&magnetization') Then 
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, T%dft%mag_info%fread, T%dft%mag_info%fail)
        T%dft%mag_info%stat = .True.
        ! Read names of the pseudo potentials 
        Call read_dft_magnetization(iunit, T)

      ! Hubbard    
      Else If (word(1:length) == '&hubbard') Then 
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, T%dft%hubbard_info%fread, T%dft%hubbard_info%fail)
        T%dft%hubbard_info%stat = .True.
        ! Read Hubbard corrections
        Call read_dft_hubbard(iunit, T)

      ! NGWF
      Else If (word(1:length) == '&ngwf') Then 
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, T%dft%ngwf_info%fread, T%dft%ngwf_info%fail)
        T%dft%ngwf_info%stat = .True.
        ! Read NGWF setings 
        Call read_dft_ngwf(iunit, T)

      Else If (word(1:length) == '&gcdft') Then
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, T%dft%gc%activate%fread, T%dft%gc%activate%fail)
        T%dft%gc%activate%stat = .True.
        ! Read GCDFT setings 
        Call read_gcdft(iunit, T)
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word), '" is not recognised as a valid DFT settings.',&
                                & ' See manual. Have you properly closed the block with "&end_dft_settings"?'
        Call error_stop(message)
      End If

     End Do

  End Subroutine read_dft_settings

  Subroutine read_gcdft(iunit, T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read sub-block &gcdft with settings
    ! for GC-DFT simulations
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: T

    Character(Len=256) :: word
    Integer(Kind=wi)   :: length, io

    Character(Len=256) :: message
    Character(Len=265) :: set_error

    set_error = '***ERROR in sub-block &gcdft (inside &dft_settings):'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&gcdft (inside &dft_settings)')
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_gcdft" to close the block. Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_gcdft') Then
        Exit
      End If
      Call check_for_rubbish(iunit, '&gcdft')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      ! ONETEP directives
      Else If (word(1:length) == 'reference_potential') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%gc%reference_potential%value, T%dft%gc%reference_potential%units
        Call set_read_status(word, io, T%dft%gc%reference_potential%fread, T%dft%gc%reference_potential%fail)
        Call capital_to_lower_case(T%dft%gc%reference_potential%units)

      Else If (word(1:length) == 'electrode_potential') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%gc%electrode_potential%value, T%dft%gc%electrode_potential%units
        Call set_read_status(word, io, T%dft%gc%electrode_potential%fread, T%dft%gc%electrode_potential%fail)
        Call capital_to_lower_case(T%dft%gc%electrode_potential%units)

      Else If (word(1:length) == 'electron_threshold') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%gc%electron_threshold%value
        Call set_read_status(word, io, T%dft%gc%electron_threshold%fread, T%dft%gc%electron_threshold%fail)

      ! CP2K directives
      Else If (word(1:length) == 'target_workfunction') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%gc%target_workfunction%value, T%dft%gc%target_workfunction%units
        Call set_read_status(word, io, T%dft%gc%target_workfunction%fread, T%dft%gc%target_workfunction%fail)
        Call capital_to_lower_case(T%dft%gc%target_workfunction%units)      

      Else If (word(1:length) == 'mixing_coefficient') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%gc%mixing_coefficient%value
        Call set_read_status(word, io, T%dft%gc%mixing_coefficient%fread, T%dft%gc%mixing_coefficient%fail)
        
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word), '" is not recognised as a valid GC-DFT settings.',&
                                & ' See manual. Have you properly closed the block with "&end_gcdft"?'
        Call error_stop(message)
      End If

     End Do

  End Subroutine 

  Subroutine read_basis_set(iunit, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the type of basis assigned to each atomic site
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: simulation_data

    Character(Len=256) :: word
    Integer(Kind=wi)   :: io
    Integer(Kind=wi)   :: i
    Logical            :: endblock

    Character(Len=256) :: message
    Character(Len=265) :: set_error

    set_error = '***ERROR in &basis_set (inside &dft_settings):'
    simulation_data%dft%basis_set(:)%defined=.False.

    i=1
    Do While (i<= simulation_data%total_tags)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&basis_set (inside &dft_settings)')
      If (word(1:1)/='#') Then
        Call check_for_rubbish(iunit, '&basis_set (inside &dft_settings)')
        Read (iunit, Fmt=*, iostat=io) simulation_data%dft%basis_set(i)%tag, simulation_data%dft%basis_set(i)%type
        Call capital_to_lower_case(simulation_data%dft%basis_set(i)%type)
        simulation_data%dft%basis_set(i)%defined=.True.
        i=i+1  
      End If
    End Do

    If (i-1 < simulation_data%total_tags) Then
      Write (message,'(2a)') Trim(set_error), ' The number of declared atomic tags is less than those defined in&
                          & "&input_composition". Please define all tags with the corresponding&
                          & basis set type'
      Call error_stop(message)
    End If

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&basis_set (inside &dft_settings)')
      Call capital_to_lower_case(word)
      If (word /= '&end_basis_set') Then
        If (word(1:1) /= '#') Then
          Write (message,'(2a)') Trim(set_error), ' It seems the user has provided wrong or additional information.&
                               & Have you defined all the elements? This block must be closed with sentence "&end_basis_set"'
          Call error_stop(message)
        End If
      Else
          endblock=.True.
      End If
    End Do

    Do i=1, simulation_data%total_tags
      If (.Not. simulation_data%dft%basis_set(i)%defined) Then
        Write (message,'(2a)') Trim(set_error), ' The basis set for species "'//Trim(simulation_data%dft%basis_set(i)%tag)//'"& 
                             & has not been defined. Please add the specification and rerun'
        Call error_stop(message)
      End If
    End Do    

  End Subroutine read_basis_set

  Subroutine read_pseudo_poptentials(iunit, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the names of the pseudo potentials files, which
    ! are required to build files for atomistic level simulation
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: simulation_data

    Character(Len=256) :: word
    Integer(Kind=wi)   :: io
    Integer(Kind=wi)   :: i
    Logical            :: endblock

    Character(Len=256) :: message
    Character(Len=265) :: set_error

    set_error = '***ERROR in &pseudo_potentials (inside &dft_settings):'
    simulation_data%dft%pseudo_pot(:)%defined=.False.

    i=1
    Do While (i<= simulation_data%total_tags)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&pseudo_potentials (inside &dft_settings)')
      If (word(1:1)/='#') Then
        Call check_for_rubbish(iunit, '&pseudo_potentials (inside &dft_settings)')
        Read (iunit, Fmt=*, iostat=io) simulation_data%dft%pseudo_pot(i)%tag, simulation_data%dft%pseudo_pot(i)%file_name
        simulation_data%dft%pseudo_pot(i)%defined=.True.
        i=i+1
      End If
    End Do

    If (i-1 < simulation_data%total_tags) Then
      Write (message,'(2a)') Trim(set_error), ' The number of declared atomic tags is less than those defined in&
                          & "&input_composition". Please define all tags with the corresponding file name for the&
                          & pseudopotential'
      Call error_stop(message)
    End If

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&pseudo_potentials (inside &dft_settings)')
      Call capital_to_lower_case(word)
      If (word /= '&end_pseudo_potentials') Then
        If (word(1:1) /= '#') Then
          Write (message,'(2a)') Trim(set_error), ' It seems the user has provided wrong or additional information.&
                               & Have you defined all the elements? This block must be closed with sentence&
                               & "&end_pseudo_potentials"'
          Call error_stop(message)
        End If
      Else
          endblock=.True.
      End If
    End Do

    Do i=1, simulation_data%total_tags
      If (.Not. simulation_data%dft%pseudo_pot(i)%defined) Then
        Write (message,'(2a)') Trim(set_error), ' The pseudopotential for species "'&
                             &//Trim(simulation_data%dft%pseudo_pot(i)%tag)//'"& 
                             & has not been defined. Please add the specification and rerun'
        Call error_stop(message)
      End If
    End Do    

  End Subroutine read_pseudo_poptentials

  Subroutine read_dft_magnetization(iunit, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the initial magnetization for the atomic tags 
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),   Intent(In   ) :: iunit
    Type(simul_type),   Intent(InOut) :: simulation_data

    Integer(Kind=wi) :: io, i, j, k
    Character(Len=256)  :: messages(8), word
    Character(Len=64 )  :: set_error
    Logical  :: endblock
    Logical  :: header(2), error_duplication

    set_error = '***ERROR in &magnetization (inside &dft_settings):'
    Write (messages(1),'(a)') set_error

    header=.False.
    error_duplication=.False.

    i=1

    Do While (i <= 2)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&magnetization (inside &dft_settings)')
      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          Call check_for_rubbish(iunit, '&magnetization (inside &dft_settings)')
          Call capital_to_lower_case(word) 

          If (Trim(word)=='tags') Then
            If (.Not. header(1)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%dft%magnetization(j)%tag, j = 1, simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read tags for atoms'
                Call info(messages, 2)
                Call error_stop(' ')
              End If
              Do j=1, simulation_data%total_tags-1
                 Do k=j+1, simulation_data%total_tags
                   If (Trim(simulation_data%dft%magnetization(j)%tag)==Trim(simulation_data%dft%magnetization(k)%tag)) Then
                     Write (messages(2),'(3(1x,a))') 'Tag', Trim(simulation_data%dft%magnetization(j)%tag),&
                                                  & 'is repeated in the list!'
                     Write (messages(3),'((1x,a))') 'All tags for the components of the species must be declared,&
                                                   & each tag only once'
                     Call info(messages, 3)
                     Call error_stop(' ')
                   End If
                 End Do
              End Do  
              header(1)=.True.
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='values') Then
            If (.Not. header(2)) Then 
              i=i+1
              Read (iunit,Fmt=*,iostat=io) word, (simulation_data%dft%magnetization(j)%value, j=1,simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the magnetization values for the defined tags. Please check if tags&
                                            & are consistent with those defined in "&Species_Components".'
                Call info(messages, 2)
                Call error_stop(' ')
              End If  
              header(2)=.True.        
            Else
              error_duplication=.True.
            End If
          Else
            Write (messages(2),'(1x,3a)') 'Wrong descriptor "', Trim(word),'". Valid options are "tags" and "values"&
                                         & Please refer to the manual.'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        Else
          Write (messages(2),'(1x,a)')      'The correct structure must be:'
          Write (messages(3),'(1x,a)')      '&magnetization'
          Write (messages(4),'(1x,a)')      '  tags      tg1     tg2     tg3    .... tgNsp'
          Write (messages(5),'(1x,a)')      '  values    mu_tg1  mu_tg2  mu_tg3 .... mu_tgNsp'
          Write (messages(6),'(1x,a)')      '&end_magnetization'
          Write (messages(7),'(1x,a,i3,a)') 'where in this case Nsp = ', simulation_data%total_tags,&
                                          & ', which corresponds to the number of tags defined in "&Species_Components".'
          Write (messages(8),'(1x,a)')      'See manual for details'
          Call info(messages, 8)
          Call error_stop(' ')
        End If
      End If
      If (error_duplication) Then
        Write (messages(2),'(1x,3a)') 'Descriptor "', Trim(word), '" is duplicated within &magnetization'
        Call info(messages, 2)
        Call error_stop(' ')
      End If
    End Do 

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call capital_to_lower_case(word)
      If (word /= '&end_magnetization') Then
        If (word(1:1) /= '#') Then
            Write (messages(2),'(3a)') 'Descriptors "tags" and "values" have already been defined. Directive "',&
                                    & Trim(word), '" is not valid. Block must be&
                                    & closed with sentence &end_magnetization.' 
            Call info(messages,2)
            Call error_stop(' ')
        End If
      Else
        endblock=.True.
      End If
    End Do

  End Subroutine read_dft_magnetization

  Subroutine read_dft_hubbard(iunit, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the hubbard specification for each atomic tags 
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: simulation_data

    Integer(Kind=wi) :: io, i, j, k
    Character(Len=256)  :: messages(10), word
    Character(Len=64 )  :: set_error
    Logical  :: endblock
    Logical  :: header(4), error_duplication

    set_error = '***ERROR in &hubbard (inside &dft_settings):'
    Write (messages(1),'(a)') set_error

    header=.False.
    error_duplication=.False.

    i=1

    Do While (i <= 4)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&hubbard (inside &dft_settings)')
      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          Call check_for_rubbish(iunit, '&hubbard (inside &dft_settings)')
          Call capital_to_lower_case(word) 

          If (Trim(word)=='tags') Then
            If (.Not. header(1)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%dft%hubbard(j)%tag, j = 1, simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read tags for atoms'
                Call info(messages, 2)
                Call error_stop(' ')
              End If
              Do j=1, simulation_data%total_tags-1
                 Do k=j+1, simulation_data%total_tags
                   If (Trim(simulation_data%dft%hubbard(j)%tag)==Trim(simulation_data%dft%hubbard(k)%tag)) Then
                     Write (messages(2),'(3(1x,a))') 'Tag', Trim(simulation_data%dft%hubbard(j)%tag),&
                                                    & 'is repeated in the list!'
                     Write (messages(3),'((1x,a))') 'All tags for the components of the species must be declared,&
                                                    & each tag only once'
                     Call info(messages, 3)
                     Call error_stop(' ')
                   End If
                 End Do
              End Do  
              header(1)=.True.
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='l_orbital') Then
            If (.Not. header(2)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%dft%hubbard(j)%l_orbital, j = 1, simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the l_orbital (which orbital to apply the correction)&
                                           & for each tags. Please check'
                Call info(messages, 2)
                Call error_stop(' ')
              End If  
              header(2)=.True.        
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='u') Then
            If (.Not. header(3)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%dft%hubbard(j)%U, j = 1, simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the "U" values for tags. Please check'
                Call info(messages, 2)
                Call error_stop(' ')
              End If  
              header(3)=.True.        
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='j') Then
            If (.Not. header(4)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%dft%hubbard(j)%J, j = 1, simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the "J" values for tags. Please check'
                Call info(messages, 2)
                Call error_stop(' ')
              End If  
              header(4)=.True.        
            Else
              error_duplication=.True.
            End If
          Else
            Write (messages(2),'(1x,3a)') 'Wrong descriptor "', Trim(word),'". Valid options are: "tags", "l_orbital", "U" and "J".&
                                         & Please refer to the manual.'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        Else
          Write (messages(2),'(1x,a)')    'The correct structure must be:'
          Write (messages(3),'(1x,a)')    '&hubbard'
          Write (messages(4),'(1x,a)')    '  tags       tg1     tg2     tg3    .... tgNsp'
          Write (messages(5),'(1x,a)')    '  l_orbital  l_tg1   l_tg2   l_tg3  .... l_tgNsp'
          Write (messages(6),'(1x,a)')    '  U          U_tg1   U_tg2   U_tg3  .... U_tgNsp'
          Write (messages(7),'(1x,a)')    '  J          J_tg1   J_tg2   J_tg3  .... J_tgNsp'
          Write (messages(8),'(1x,a)')    '&end_hubbard'
          Write (messages(9),'(1x,a,i3,a)') 'where in this case Nsp = ', simulation_data%total_tags, &
                                          & ', which corresponds to the number of tags defined in "&Species_Components".'
          Write (messages(10),'(1x,a)')    'See manual for details'
          Call info(messages, 10)
          Call error_stop(' ')
        End If
      End If
      If (error_duplication) Then
        Write (messages(2),'(1x,3a)') 'Descriptor "', Trim(word), '" is duplicated within &hubbard'
        Call info(messages, 2)
        Call error_stop(' ')
      End If
    End Do 

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&hubbard (inside &dft_settings)')
      Call capital_to_lower_case(word)
      If (word /= '&end_hubbard') Then
        If (word(1:1) /= '#') Then
            Write (messages(2),'(3a)') 'Descriptors "tags", "l_orbital", "U" and "J" have already been defined. Directive "',&
                                    & Trim(word), '" is not valid. Block must be&
                                    & closed with sentence &end_hubbard' 
            Call info(messages,2)
            Call error_stop(' ')
        End If
      Else
        endblock=.True.
      End If
    End Do

  End Subroutine read_dft_hubbard


  Subroutine read_dft_ngwf(iunit, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the NGWF specification for each atomic tags (only for ONETEP)
    !
    ! author    - i.scivetti  March 2026 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: simulation_data

    Integer(Kind=wi) :: io, i, j, k
    Character(Len=256)  :: messages(9), word
    Character(Len=64 )  :: set_error
    Logical  :: endblock
    Logical  :: header(4), error_duplication

    set_error = '***ERROR in &ngwf (inside &dft_settings):'
    Write (messages(1),'(a)') set_error

    header=.False.
    error_duplication=.False.

    i=1

    Do While (i <= 3)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&ngwf')
      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          Call check_for_rubbish(iunit, '&ngwf (inside &dft_settings)')
          Call capital_to_lower_case(word) 
          If (Trim(word)=='tags') Then
            If (.Not. header(1)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%dft%ngwf(j)%tag, j = 1, simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read tags for atoms'
                Call info(messages, 2)
                Call error_stop(' ')
              End If
              Do j=1, simulation_data%total_tags-1
                 Do k=j+1, simulation_data%total_tags
                   If (Trim(simulation_data%dft%ngwf(j)%tag)==Trim(simulation_data%dft%ngwf(k)%tag)) Then
                     Write (messages(2),'(3(1x,a))') 'Tag', Trim(simulation_data%dft%ngwf(j)%tag),&
                                                    & 'is repeated in the list!'
                     Write (messages(3),'((1x,a))') 'All tags for the components of the species must be declared,&
                                                    & each tag only once'
                     Call info(messages, 3)
                     Call error_stop(' ')
                   End If
                 End Do
              End Do  
              header(1)=.True.
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='number') Then
            If (.Not. header(2)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%dft%ngwf(j)%ni, j = 1, simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the number of ngwf functions (row labelled as number)&
                                           & for each tags. Please check'
                Call info(messages, 2)
                Call error_stop(' ')
              End If  
              header(2)=.True.        
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='radius') Then
            If (.Not. header(3)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%dft%ngwf(j)%radius, j = 1, simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the valus of "radius" for the tags. Please check'
                Call info(messages, 2)
                Call error_stop(' ')
              End If  
              header(3)=.True.        
            Else
              error_duplication=.True.
            End If
          Else
            Write (messages(2),'(1x,3a)') 'Wrong descriptor "', Trim(word),&
                                        & '". Valid options are: "tags", "number" and "radius". Please refer to the manual.'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        Else
          Write (messages(2),'(1x,a)')    'The correct structure must be:'
          Write (messages(3),'(1x,a)')    '&ngwf'
          Write (messages(4),'(1x,a)')    '  tags         tg1     tg2     tg3    .... tgNsp'
          Write (messages(5),'(1x,a)')    '  number       n_tg1   n_tg2   n_tg3  .... n_tgNsp'
          Write (messages(6),'(1x,a)')    '  radius       r_tg1   r_tg2   r_tg3  .... r_tgNsp'
          Write (messages(7),'(1x,a)')    '&end_ngwf'
          Write (messages(8),'(1x,a,i3,a)') 'where in this case Nsp = ', simulation_data%total_tags, &
                                        & ', which corresponds to the number of tags defined in "&Species_Components".'
          Write (messages(9),'(1x,a)')    'See manual for details.'
          Call info(messages, 9)
          Call error_stop(' ')
        End If
      End If
      If (error_duplication) Then
        Write (messages(2),'(1x,3a)') 'Descriptor "', Trim(word), '" is duplicated within &ngwf'
        Call info(messages, 2)
        Call error_stop(' ')
      End If
    End Do 

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word

      Call capital_to_lower_case(word)
      If (word /= '&end_ngwf') Then
        If (word(1:1) /= '#') Then
            Write (messages(2),'(3a)') 'Descriptors "tags", "number" and "radius" have already been defined. Directive "',&
                                    & Trim(word), '" is not valid. Block must be closed with sentence &end_ngwf' 
            Call info(messages,2)
            Call error_stop(' ')
        End If
      Else
        endblock=.True.
      End If
    End Do

  End Subroutine read_dft_ngwf

  
End Module simulation_setup
