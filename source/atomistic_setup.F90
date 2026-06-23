!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module that allocate arrays and defines atomistic related variables
! 
! Copyright:  2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author       - i.scivetti  Feb 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module atomistic_setup 

  Use constants,          Only : Bohr_to_A,&
                                 twopi, &
                                 chemsymbol, & 
                                 NPTE,&
                                 min_inter_bond, &
                                 min_intra_bond,&
                                 max_intra_bond,&
                                 max_components,&
                                 max_species, &
                                 max_num_species_units, &
                                 max_at_species, &
                                 max_length_name_species
                                 
  Use fileset,            Only : file_type,              &
                                 FILE_SET,               &
                                 FILE_INPUT_ELECTRODE, &
                                 FOLDER_INPUT_GEOM
                                 
  Use input_types,        Only : in_integer, &
                                 in_integer_array, & 
                                 in_logic,   &
                                 in_string,  &
                                 in_param,   &
                                 in_scalar
  Use numprec,            Only : li, &
                                 wi, &
                                 wp
                                 
  Use process_data,       Only : capital_to_lower_case, &
                                 check_for_symbols, &
                                 get_word_length, &
                                 check_end, &
                                 check_for_rubbish, &
                                 duplication_error, &
                                 set_read_status
  
  Use unit_output,        Only : error_stop,&
                                 info 

  Implicit None
  Private

  ! Type to describe the atoms
  Type :: atom_type
     Real(Kind=wp)     :: r(3)
     Character(Len=8)  :: tag
     Character(Len=2)  :: element
     Integer(Kind=wi)  :: atomic_number
     Logical           :: dynamics(3)
     Logical           :: in_species
     Logical           :: vanish
  End Type 
 
  ! Types related to format for VASP files 
  Type :: list_type
    Character(Len=8)  :: tag(max_components)
    Character(Len=2)  :: element(max_components) 
    Integer(Kind=wi)  :: N0(max_components)   
    Integer(Kind=wi)  :: num_elements
    Integer(Kind=wi)  :: net_elements
    Character(Len=32) :: coord_type 
  End Type 
 
  ! Type for the  &input_composition
  Type :: component_in_block
    Character(Len=8)  :: tag(max_components)
    Character(Len=2)  :: element(max_components)
    Integer(Kind=wi)  :: atomic_number(max_components)
    Integer(Kind=wi)  :: N0(max_components)
    Integer(Kind=wi)  :: numtot
  End Type 

  ! Type for the components of a given species 
  Type :: component_type
    Character(Len=8)  :: tag(max_components)
    Character(Len=2)  :: element(max_components)
    Integer(Kind=wi)  :: atomic_number(max_components)
    Integer(Kind=wi)  :: N(max_components)
    Integer(Kind=wi)  :: N0(max_components)
    Logical           :: fread= .False.
  End Type  
  
  ! Type for the definition of geometry of the species   
  Type :: define_species
     Character(Len=8) :: tag(max_at_species)
     Character(Len=2) :: element(max_at_species)
     Integer(Kind=wi) :: atomic_number(max_components)
     Real(Kind=wp)    :: r0(max_at_species,3)
     Logical          :: in_species(max_at_species)
  End Type
 
  ! Type for the species units
  Type :: species_units
     Integer(Kind=wi) :: list(max_at_species)
     Logical          :: vanish
  End Type 

  Type :: species_setup
    !Description
    Character(Len=32) :: tag    
    Character(Len=32) :: vartype0    
    Integer(Kind=wi)  :: N0_target
    Real(Kind=wp)     :: bond_cutoff
    Logical           :: fread= .False.
    Logical           :: fail = .False.
    Integer(Kind=wi)  :: num_components
    Integer(Kind=wi)  :: atoms_per_species 
    Character(Len=32) :: topology
    Type(component_type) :: component
  End Type  
  
  ! Type to describe species within the models
  Type, Public :: species_model
    Character(Len=32)     :: tag
    Character(Len=32)     :: topology
    Integer(Kind=wi)      :: num
    Integer(Kind=wi)      :: num_show
    Integer(Kind=wi)      :: num_extra  
    Integer(Kind=wi)      :: D_num
    Integer(Kind=wi)      :: atoms_per_species 
    Integer(Kind=wi)      :: num_components
    Logical               :: change_content
    Real(Kind=wp)         :: s
    Type(component_type)  :: component
    Type(define_species)  :: definition
    Type(species_units)   :: units(max_num_species_units)
  End Type 

  ! Type for the sample models 
  Type, Public :: sample_type
    ! Path to subfolders
    Character(Len=256) :: path 
    ! normal 
    Real(Kind=wp)      :: normal(3)    
    ! Cell vectors 
    Real(Kind=wp)      :: cell(3,3)
    ! Inverse cell vectors
    Real(Kind=wp)      :: invcell(3,3)
    ! Length of cell vectors
    Real(Kind=wp)      :: cell_length(3)
    ! Scale factor vasp
    Real(Kind=wp)      :: scale_factor_vasp=1.0_wp
    ! Atoms in the model
    Type(atom_type), Allocatable :: atom(:)
    ! Maximum number of atom
    Integer(Kind=wi) :: max_atoms
    ! Total number of atoms  
    Integer(Kind=wi) :: num_atoms
    ! Total number of atoms (needed if inserting species) 
    Integer(Kind=wi) :: num_atoms_extra
    ! list elements
    Type(list_type)  :: list
    ! Minimum number of species units among all the species
    Integer(Kind=wi) :: min_species
    ! Minimum number of atomic  components among all the fixed species
    Integer(Kind=wi) :: min_components
    ! Information for the species  
    Type(species_model), Allocatable :: species(:)
    ! Slab Area
    Real(Kind=wp)             :: slab_area
    ! geometric centre of slab
    Real(Kind=wp)             :: slab_centre(3)
    ! logical to cover both surfaces of slab
    Logical                   :: both_surfaces
    ! Invert atom
    Real(Kind=wp)             :: inverted(3)
    ! Surface shift
    Real(Kind=wp)             :: surface_shift(3)
    ! need_species
    Logical                   :: change_species_number
    ! need_species
    Logical                   :: size_changed    
  End Type  

  ! Type for the modelling related variables 
  Type, Public :: model_type
    Private
    ! Multiple for the amount of output_atom
    Type(in_integer), Public :: multiple_output_atoms 
    ! Multiple for the amount of input_atoms
    Type(in_integer), Public :: multiple_input_atoms 
    ! Constrained dynamics
    Logical, Public          :: selective_dyn=.False.
    ! Flag to rotate the species
    Type(in_logic),  Public  :: rotate_species 
    ! Format of input model 
    Type(in_string), Public  :: input_model_format
    ! Format of output model 
    Type(in_string), Public  :: output_model_format
    ! Vector normal to the slab surface (only for deposited) 
    Type(in_string), Public  :: normal_vector
    ! Flag to insert species
    Logical, Public :: insert_species
    ! Flag to remove species
    Logical, Public :: remove_species
    ! Total retition of the input model 
    Integer(Kind=wi), Public :: input_times
    ! Arrays for sampling the spatial region within the simulation cell 
    Integer(Kind=wi), Public, Dimension(3)  :: scan
    Integer(Kind=wi), Public, Dimension(3)  :: npoints 
    ! Distance cutoff (minium value for inter-species separation)
    Type(in_param), Public   :: distance_cutoff
    ! Repetition of the input model to build the output sample
    Type(in_integer_array), Public  :: repeat_input_model 
    ! Number of different species types
    Integer(Kind=wi), Public :: types_species 
    ! Type for those components defined in input_composition
    Type(component_in_block), Public  ::  component
    ! Block with vectors for the input structure model 
    Type(in_string), Public  :: input_cell 
    ! Block for tags and amount of atoms that compose the input model
    Type(in_string), Public  :: input_composition 
    ! Input model 
    Type(sample_type), Public  :: input 
    ! Generated atomic sample
    Type(sample_type), Public  :: sample
    ! Discretization of the spatial coordinates
    Type(in_param), Public   :: delta_space
    ! Number of species involved 
    Type(in_integer) ,   Public  :: num_species
    ! Type of analysis involved
    Type(in_string) ,   Public  :: analysis 
    ! Directive for the definition of &species
    Type(in_string), Public        :: species
    ! Type for the definition of &species
    Type(species_setup), Allocatable, Public :: species_info(:)
    ! Directive for the definition of &species_components
    Type(in_string), Public        :: species_components
    ! Number of different tags used to name the involved species 
    Integer(Kind=wi),             Public    :: total_tags
    ! instruction of how to build the electrolyte model
    Type(in_string), Public        :: arr_added_species
    ! Decide if the electrode is centered or not
    Type(in_logic), Public    :: centre_electrode
    ! Decide if the cell size will be optimised or not
    Type(in_logic), Public    :: optimise_size
    ! Decide if both surfaces of the electrode slab are considered or not
    Type(in_logic), Public    :: both_surfaces
    ! Level where to start deposition from
    Type(in_param), Public    :: add_species_from
    ! Add extra space
    Type(in_param), Public    :: add_extra_space
    
  Contains
     Private
     Procedure          :: atomic_arrays_input  => allocate_input_atomic_arrays
     Procedure, Public  :: init_species_info    => allocate_species_info
     Procedure, Public  :: init_input_variables => allocate_input_variables
     Procedure, Public  :: species_arrays       => allocate_species_arrays
     Procedure, Public  :: atomic_arrays_model   => allocate_atomic_arrays_model
     Final              :: cleanup
  End Type model_type

  Public :: read_species, read_species_components, read_input_composition, read_input_cell
  Public :: read_input_model
  Public :: check_species, check_components_species, check_atomic_settings
  Public :: about_cell

Contains

  Subroutine allocate_input_variables(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate essential variables to build atomistic models
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(model_type), Intent(InOut)  :: T

    Integer(Kind=wi)     :: fail(1)
    Character(Len=256)   :: message

    Allocate(T%repeat_input_model%value(3), Stat=fail(1))

    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for input variables needed for&
                                & atomistic modelling (subroutine allocate_input_variables)'
      Call error_stop(message)
    End If

  End Subroutine allocate_input_variables
  
  Subroutine allocate_species_info(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate arrays for defined species
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(model_type), Intent(InOut)  :: T

    Integer(Kind=wi) :: fail
    Integer(Kind=wi) :: number_species
    Character(Len=256)   :: message

    number_species=T%num_species%value
    fail =0

    Allocate(T%species_info(number_species), Stat=fail)
    If (fail > 0) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for quantities of the &species'
      Call error_stop(message)
    End If

  End Subroutine allocate_species_info  

  Subroutine allocate_input_atomic_arrays(T,num_atoms)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate atomic arrays for the input model
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(model_type),  Intent(InOut)  :: T
    Integer(Kind=wi),   Intent(In   )  :: num_atoms

    Integer(Kind=wi)    :: fail(1)
    Character(Len=256)  :: message


    T%input%max_atoms=T%multiple_input_atoms%value*num_atoms  
    fail=0

    Allocate(T%input%atom(T%input%max_atoms),   Stat=fail(1))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for the atomic arrays of input model&
                                & (subroutine allocate_input_atomic_arrays). It seems the input model is&
                                & too large. Reduce the value of "multiple_input_atoms" (set to 10000 by default)'
      Call error_stop(message)
    End If

    T%input%atom(:)%in_species=.False.
    T%input%atom(:)%vanish=.False.

  End Subroutine allocate_input_atomic_arrays  
  
  Subroutine allocate_species_arrays(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate species array for characterization of 
    ! the participating species 
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(model_type), Intent(InOut)  :: T

    Integer(Kind=wi)    :: fail(2)
    Integer(Kind=wi)    :: i, j, num_species    
    Character(Len=256)  :: message

    num_species= T%num_species%value
    Allocate(T%input%species(num_species),    Stat=fail(1))
    Allocate(T%sample%species(num_species),   Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for species arrays&
                               & (subroutine allocate_species_arrays)'
      Call error_stop(message)
    End If

    ! Copy variables from the T%species_info
    Do i = 1, num_species
      T%input%species(i)%tag = T%species_info(i)%tag
      T%input%species(i)%topology = T%species_info(i)%topology
      T%input%species(i)%atoms_per_species = T%species_info(i)%atoms_per_species
      T%input%species(i)%num_components = T%species_info(i)%num_components
      T%input%species(i)%component%tag = T%species_info(i)%component%tag
      T%input%species(i)%component%element = T%species_info(i)%component%element
      T%input%species(i)%component%atomic_number = T%species_info(i)%component%atomic_number
      T%input%species(i)%component%N0=T%species_info(i)%component%N0
      If (Trim(T%species_info(i)%vartype0)=='fixed') Then
        T%input%species(i)%change_content=.False.
      Else
        T%input%species(i)%change_content=.True.
      End If 
    End Do
    
    ! Copy from input to sample
    Do i = 1, num_species
      T%sample%species(i)%tag=T%input%species(i)%tag
      T%sample%species(i)%topology=T%input%species(i)%topology
      T%sample%species(i)%atoms_per_species=T%input%species(i)%atoms_per_species
      T%sample%species(i)%component=T%input%species(i)%component
      T%sample%species(i)%num_components=T%input%species(i)%num_components
      T%sample%species(i)%definition=T%input%species(i)%definition
      T%sample%species(i)%change_content=T%input%species(i)%change_content
    End Do

    Do i = 1, num_species
      Do j = 1, max_num_species_units
        T%sample%species(i)%units(j)%vanish=.False.
        T%input%species(i)%units(j)%vanish=.False.
      End Do
    End Do

  End Subroutine allocate_species_arrays

  Subroutine allocate_atomic_arrays_model(T,num_atoms)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate atomic arrays for the sample model
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(model_type),   Intent(InOut)  :: T
    Integer(Kind=wi),    Intent(In   )  :: num_atoms

    Integer(Kind=wi)     :: fail(1)
    Character(Len=256)   :: message

    fail=0
    T%sample%max_atoms=T%multiple_output_atoms%value*T%input_times*num_atoms  

    Allocate(T%sample%atom(T%sample%max_atoms),        Stat=fail(1))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for the atomic arrays&
                               & of output structure model (subroutine allocate_atomic_arrays_model)'
      Call error_stop(message)
    End If

    T%sample%atom(:)%in_species=.False.
    T%sample%atom(:)%vanish=.False.

  End Subroutine allocate_atomic_arrays_model

  Subroutine cleanup(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deallocate variables
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type) :: T

    If (Allocated(T%species_info)) Then
      Deallocate(T%species_info)
    End If     
    
    If (Allocated(T%input%atom)) Then
      Deallocate(T%input%atom)
    End If 
 
    If (Allocated(T%repeat_input_model%value)) Then
      Deallocate(T%repeat_input_model%value)
    End If 

    If (Allocated(T%input%species)) Then
      Deallocate(T%input%species)
    End If 

    If (Allocated(T%sample%species)) Then
      Deallocate(T%sample%species)
    End If 
    
    If (Allocated(T%sample%atom)) Then
      Deallocate(T%sample%atom)
    End If 

  End Subroutine cleanup 
  
  Subroutine read_species(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read settings inside species, needed for building the 
    ! atomistic models
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(model_type),  Intent(InOut) :: model_data 

    Integer(Kind=wi)  ::  io, i

    Character(Len=256)  :: word, word2, messages(3)
    Character(Len=64 )  :: error_species
    Logical  :: error, endblock, fread 

    error= .False.
    error_species = '***ERROR in &species of SETTINGS file'
    Write (messages(1),'(a)') error_species 

    fread= .True.
    Do While (fread)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&species')
      If (word(1:1)/='#') Then
        fread=.False.
        Call check_for_rubbish(iunit, '&species')
      End If
    End Do

    ! Read number of species
    Read (iunit, Fmt=*, iostat=io) model_data%num_species%tag, model_data%num_species%value
    Call set_read_status(word, io, model_data%num_species%fread, model_data%num_species%fail)

    If (Trim(model_data%num_species%tag) /= 'number_species') Then
      Write (messages(2),'(3a)') 'Directive "', Trim(model_data%num_species%tag), &
                         & '" has been found, but directive "number_species" is expected.'
      error=.True.
    End If 

    If (model_data%num_species%fail) Then
      Write (messages(2),'(a)') 'Wrong (or missing) specification for directive "number_species"'
      error=.True.
    Else  
      If (model_data%num_species%value<1) Then
        Write (messages(2),'(a)') 'Number of species MUST be >= 1'
        error=.True.
      ElseIf (model_data%num_species%value>max_species) Then
        Write (messages(2),'(a,i3,a)') 'Are you sure you want to consider a system involving more than ', max_species,&
                                    & ' different types of species?. We have set this value as a sensible limit...'
        error=.True.
      End If
    End If

    If (error) Then
      Call info(messages,2) 
      Call error_stop(' ')
    End If


    ! Initialise arrays
    Call model_data%init_species_info()   

    i=1
    Do While (i <= model_data%num_species%value)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&species')
      If (word(1:1)/='#') Then
        Call check_for_rubbish(iunit, '&species')
  
        Read (iunit, Fmt=*, iostat=io) model_data%species_info(i)%tag, model_data%species_info(i)%N0_target, &
                                       model_data%species_info(i)%vartype0
        If (io/=0) Then
          word2=model_data%species_info(i)%tag
          Call capital_to_lower_case(word2)
          If (Trim(word2) == '&end_species') Then
            Write (messages(2),'(2(a,i2),a)') 'Missing specification for atomic species. Only ',  i-1,&
                                    &' species set out of ', model_data%num_species%value, ' (see number_species)'
            Call info(messages,2) 
            Call error_stop(' ')
          End If 
        End If

        Call capital_to_lower_case(model_data%species_info(i)%vartype0)

        Call set_read_status(word, io, model_data%species_info(i)%fread, model_data%species_info(i)%fail)
        If (model_data%species_info(i)%fail) Then
          Write (messages(2),'(2a)') 'Problems with the specification of species ', Trim(model_data%species_info(i)%tag)
          Call info(messages,2)
          Call error_stop(' ')
        End If 
        i=i+1
      End If
    End Do 

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&species')
      Call capital_to_lower_case(word)
      If (word /= '&end_species') Then
        If (word(1:1) /= '#') Then
          If ((i-1)/=model_data%num_species%value) Then 
            Write (messages(2),'(a)') 'Number of atomic species specified is larger than&
                                & the value given by directive "number_species"'
          Else
            Write (messages(2),'(a)') 'Block must be closed with sentence &end_species. Please check.&
                                     & Is the number of defined species the same as set for directive "number_species"?'
          End If   
          Call info(messages,2) 
          Call error_stop(' ')
        End If
      Else
        endblock=.True.
      End If  
    End Do
    
    i=1
    endblock=.False.
    Do While (i <= model_data%num_species%value)
      If (model_data%species_info(i)%vartype0 == 'fixed') Then
        endblock=.True. 
      End If
      i=i+1
    End Do
  
    If (.Not. endblock) Then
       Write (messages(2),'(a)') 'At least one species must be defined with "fixed" content!'
       Write (messages(3),'(a)') 'SOMETHING IS WRONG...REVIEW THE SETTINGS AND SYNTAX FOR EACH SPECIES'
       Call info(messages,3) 
       Call error_stop(' ')       
    End If

  End Subroutine read_species
  
  
  Subroutine read_species_components(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read atomic composition from &species_components for 
    ! each chemical species defined in &species
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(model_type),  Intent(InOut) :: model_data

    Integer(Kind=wi)   :: io, i, j, k, l, num_comp, length
    Character(Len=256) :: messages(8), word
    Character(Len=64 ) :: error_species_components
    Character(Len=64)  ::  separator
    Logical  :: endblock, species_read, loop, fread

    error_species_components = '***ERROR in &species_components of SETTINGS file'
    Write (messages(1),'(a)') error_species_components
    Write (messages(3),'(1x,a)') ' '
    Write (messages(5),'(1x,a)')  'tag_1  element_1  N0_1  "separator"&
                                 & tag_2  element_2  N0_2 .... "separator"&  
                                 & tag_Nc element_Nc N0_Nc' 
    Write (messages(6),'(1x,a)') ' '
    Write (messages(7),'(1x,a)') 'where N0_i is the amount of atoms for each of the Nc component of the species "i". "separator"&
                                & can be any string to separate the specifications of each components for a better visualization'
    Write (messages(8),'(1x,a)') 'IMPORTANT: All tags must have maximum of 4 characters. Have you missed any of the NC components?&
                                & No comment is allowed in between the two senteces.'

    i=1
    Do While (i <= model_data%num_species%value)

      Read (iunit, Fmt=*, iostat=io) word
      If (io/=0) Then
        Exit 
      End If

      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          species_read=.False.
          j=1
          loop=.True.
          Do While (j <= model_data%num_species%value .And. loop)
            If (Trim(word)==Trim(model_data%species_info(j)%tag)) Then 
              species_read=.True.
              If (.Not. model_data%species_info(j)%component%fread) Then
                model_data%species_info(j)%component%fread = .True.
                loop=.False.
                Call check_for_rubbish(iunit, '&species_components')
                Read (iunit, Fmt=*, iostat=io) word, model_data%species_info(j)%num_components,&
                                                   model_data%species_info(j)%topology

                If (model_data%species_info(j)%num_components==1 .And. Trim(model_data%species_info(j)%topology)=='molecule')Then
                  Write (messages(2),'(1x,3a)') 'Wrong format for the atomic specification of species "', Trim(word), '"' 
                  Write (messages(3),'(1x,a)')  'A molecular species must have a number of atoms larger than 1'
                  Write (messages(4),'(1x,a)')  'If this species has only 1 atom, change its topology to "atom"'
                  Call info(messages, 4)
                  Call error_stop(' ')
                End If
 
 
                Write (messages(2),'(1x,3a)') 'Wrong format for the atomic specification of species "', Trim(word), &
                                             '". Format for each species must be:' 
                Write (messages(4),'(1x,a)')  'Species_tag  number_components(Nc)   topology   &
                                              &bond_cutoff (bond_cutoff only needed if topology is "molecule")' 
                If (io /= 0) Then
                  Call info(messages, 8)
                  Call error_stop(' ')
                End If 

                Call capital_to_lower_case(model_data%species_info(j)%topology) 

                If (model_data%species_info(j)%topology=='molecule') Then
                  Backspace iunit
                  Read (iunit, Fmt=*, iostat=io) word, model_data%species_info(j)%num_components,&
                                                   model_data%species_info(j)%topology, model_data%species_info(j)%bond_cutoff
                  If (io /= 0) Then
                    Write (messages(2),'(1x,3a)') 'Wrong (or missing) setting for bond_cutoff of species "', Trim(word), &
                                             '". Format for each species must be:' 
                        
                    Call info(messages, 8)
                    Call error_stop(' ')
                  End If 
                  ! Check if the bond_cutoff value is sensible
                  If (model_data%species_info(j)%bond_cutoff > max_intra_bond) Then
                    Write (messages(2),'(1x,a,f4.2,3a,f4.2,a)') 'The maximum limit for the bonding distance of ',&
                                                          & model_data%species_info(j)%bond_cutoff,&
                                                          &' Angstrom set for species "', Trim(word), '" MUST NOT be larger than ',&
                                                          & max_intra_bond, ' Angstrom. Please change'
                    Call info(messages,2)
                    Call error_stop(' ')
                  End If
                
                  If (model_data%species_info(j)%bond_cutoff < min_intra_bond) Then
                    Write (messages(2),'(1x,a,f4.2,3a,f4.2,a)') 'The minimum limit for the bonding distance of ',&
                                                          & model_data%species_info(j)%bond_cutoff,&
                                                          &' Angstrom set for species "', Trim(word), '" MUST NOT be lower than ',&
                                                          & min_intra_bond, ' Angstrom. Please change'
                    Call info(messages,2)
                    Call error_stop(' ')
                  End If

                End If

    
                If (model_data%species_info(j)%num_components < 1) Then
                  Write (messages(2),'(1x,3a)') 'The number of atomic components for species ', Trim(word), &
                                               ' must be larger than zero!'
                  Call info(messages, 2)
                  Call error_stop(' ')
                End If 

                num_comp=model_data%species_info(j)%num_components
                fread= .True.
                Do While (fread)
                  Read (iunit, Fmt=*, iostat=io) word
                  If (word(1:1)/='#') Then
                    fread=.False.
                    Call check_for_rubbish(iunit, '&species_components')
                  End If
                End Do

                Read (iunit, Fmt=*, iostat=io) (model_data%species_info(j)%component%tag(k), &
                                                   model_data%species_info(j)%component%element(k), &
                                                   model_data%species_info(j)%component%N0(k), separator, &   
                                                   k = 1, num_comp-1), &
                                                   model_data%species_info(j)%component%tag(num_comp), &
                                                   model_data%species_info(j)%component%element(num_comp), &
                                                   model_data%species_info(j)%component%N0(num_comp)
 
                If (io /= 0) Then
                  Call info(messages, 8)
                  Call error_stop(' ')
                End If 

                If (Trim(model_data%species_info(j)%vartype0)/='fixed') Then
                  If (model_data%species_info(j)%num_components==1 .And. model_data%species_info(j)%component%N0(1)==1 .And. &
                     model_data%species_info(j)%topology/='atom') Then
                    Write (messages(2),'(1x,3a)') 'It appears that the topology of species ', Trim(word), ' should be set&
                                                 & to "atom" in &species_components. Please check'
                    Call info(messages, 2)
                    Call error_stop(' ')
                  End If   

                  If (model_data%species_info(j)%component%N0(1)/=1 .And. model_data%species_info(j)%topology/='molecule') Then
                    Write (messages(2),'(1x,3a)') 'It appears that the topology of species ', Trim(word), ' should be set&
                                                 & to "molecule" in &species_components. Please check'
                    Call info(messages, 2)
                    Call error_stop(' ')
                  End If   
                End If

              Else
                Write (messages(2),'(3(1x,a))') 'Atomic description for species', &
                                               Trim(word), 'has been defined more than once'
                Call info(messages, 2)
                Call error_stop(' ')
              End If
            End If
            j=j+1
          End Do

          If (.Not. species_read) Then
            Write (messages(2),'(3a)') 'Species "', Trim(Adjustl(word)), '" has NOT been defined in&
                             & &species. Please check for typo (capitalization DOES matter&
                             & for the definition of the species)'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
          If (.Not. model_data%species_info(j-1)%component%fread) Then
            If (word(1:1)=="&") Then
              Write (messages(2),'(a)') 'Missing atomic description. Not all species defined&
                                       & in &species have been included' 
              Call info(messages, 2)
              Call error_stop(' ')
            Else 
              Write (messages(2),'(3(1x,a))') 'Atomic description for species', &
                                             Trim(word), 'has been defined more than once'
              Call info(messages, 2)
              Call error_stop(' ')
            End IF
          End If
          i=i+1
        Else
          Write (messages(2),'(1x,a)') 'End of block found! Missing atomistic description for&
                                     & the species defined in &species'
          Call info(messages, 2)
          Call error_stop(' ')
        End If
      End If
    End Do 

    endblock=.False.
    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call capital_to_lower_case(word)
      If (word /= '&end_species_components') Then
        If (word(1:1) /= '#') Then
          If ((i-1)/=model_data%num_species%value) Then
            Write (messages(2),'(a)') 'Atomistic details for the involved species must be closed with&
                                    & sentence &end_species_components. Check 1) syntax 2) if the&
                                    & number of specifications is larger than the value assigned to "number_species"'
          Else
            Write (messages(2),'(a)') 'The number of specifications for the composition of the involved&
                                     & species is larger than the number of species (see&
                                     & value for "number_species" in &species)'
          End If
          Call info(messages,2)
          Call error_stop(' ')
        End If
      Else
          endblock=.True.
      End If
    End Do

    i=1
    endblock=.False.
    Do While (i <= model_data%num_species%value)
      If (model_data%species_info(i)%topology == 'electrode') Then
        endblock=.True. 
      End If
      i=i+1
    End Do
  
    If (.Not. endblock) Then
       Write (messages(2),'(a)') 'There must be at least one species with "electrode" topology'
       Call info(messages,2) 
       Call error_stop(' ')       
    End If    
    
    !Check that tags for compoenents do not exceed 4 characters
    Do i=1, model_data%num_species%value
      Do j=1, model_data%species_info(i)%num_components
        Call get_word_length(model_data%species_info(i)%component%tag(j),length)
        If (length > 4) Then
          Write (messages(2),'(5a)') 'Tag "', Trim(model_data%species_info(i)%component%tag(j)), '" for atomic component&
                                   & of species "', Trim(model_data%species_info(i)%tag), '" must have a maximum of 4 characters.&
                                   & Please use a different tag'
          Call info(messages,2)
          Call error_stop(' ')  
        End If
      End Do
    End Do
  
    ! Within the same species
    Do i=1, model_data%num_species%value
      Do j=1, model_data%species_info(i)%num_components
        Do k=1, model_data%species_info(i)%num_components
          If (j/=k) Then
            If (Trim(model_data%species_info(i)%component%tag(j))==Trim(model_data%species_info(i)%component%tag(k))) Then
              Write (messages(2),'(1x,5a)') 'Component tag "', Trim(model_data%species_info(i)%component%tag(k)), &
                             &'" for species "', Trim(model_data%species_info(i)%tag), '" has been used more than once.&  
                             & Tags for components must be unequivocaly defined.'
              Call info(messages, 2)
              Call error_stop(' ')  
            End If 
          End If
        End Do
      End Do 
    End Do
    ! Check if tags for components are unequivocaly defined
    ! In other species
    Do i=1, model_data%num_species%value
      Do l=1, model_data%num_species%value
        If (i/=l) Then
          Do j=1, model_data%species_info(i)%num_components
            Do k=1, model_data%species_info(l)%num_components
                If (Trim(model_data%species_info(i)%component%tag(j))==Trim(model_data%species_info(l)%component%tag(k))) Then
                  Write (messages(2),'(1x,7a)') 'Component tag "', Trim(model_data%species_info(l)%component%tag(k)), &
                                 &'" has been used for species "', Trim(model_data%species_info(i)%tag), '" and species "',&  
                                 & Trim(model_data%species_info(l)%tag), '". Tags for components must be unequivocaly defined.'
                  Call info(messages, 2)
                  Call error_stop(' ')  
                End If 
            End Do
          End Do 
        End If
      End Do
    End Do

    ! Calculate the total number of different tags
    model_data%total_tags=0
    Do i=1, model_data%num_species%value
      model_data%total_tags=model_data%total_tags+model_data%species_info(i)%num_components
    End Do

    ! Calculate total number of atoms per species
    Do i=1, model_data%num_species%value
      model_data%species_info(i)%atoms_per_species=0
      Do j=1,model_data%species_info(i)%num_components
        model_data%species_info(i)%atoms_per_species = model_data%species_info(i)%atoms_per_species + &
                                                   model_data%species_info(i)%component%N0(j)
      End Do
    End Do

  End Subroutine read_species_components  

  Subroutine read_input_composition(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the amount of atoms of each atomic components of the 
    ! participating species from &input_composition
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(model_type),  Intent(InOut) :: model_data

    Integer(Kind=wi)   :: io, i, j, k
    Character(Len=256) :: messages(8), word
    Character(Len=64 ) :: error_input_composition
    Logical  :: endblock, loop
    Logical  :: header(2), error_duplication

    error_input_composition = '***ERROR in &input_composition of SETTINGS file'
    Write (messages(1),'(a)') error_input_composition

    header=.False.
    error_duplication=.False.

    i=1

    Do While (i <= 2)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&input_composition')
      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          Call capital_to_lower_case(word) 
          If (Trim(word)=='tags') Then
            If (.Not. header(1)) Then 
              i=i+1
              Call check_for_rubbish(iunit, '&input_composition') 
              Read (iunit, Fmt=*, iostat=io) word, (model_data%component%tag(j), j = 1, model_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read tags for atoms'
                Call info(messages, 2)
                Call error_stop(' ')
              End If
              Do j=1, model_data%total_tags-1
                 Do k=j+1, model_data%total_tags
                   If (Trim(model_data%component%tag(j))==Trim(model_data%component%tag(k))) Then
                     Write (messages(2),'(3(1x,a))') 'Tag', Trim(model_data%component%tag(j)), 'is repeated in the list!'
                     Write (messages(3),'((1x,a))') 'All tags for the atomic components of the species must be declared,&
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
          Else If (Trim(word)=='amounts') Then
            If (.Not. header(2)) Then 
              i=i+1
              Call check_for_rubbish(iunit, '&input_composition') 
              Read (iunit, Fmt=*, iostat=io) word, (model_data%component%N0(j), j = 1, model_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the amount of atoms for each atomic tag'
                Call info(messages, 2)
                Call error_stop(' ')
              End If  
              header(2)=.True.        
            Else
              error_duplication=.True.
            End If
          Else
            Write (messages(2),'(1x,3a)') 'Wrong descriptor "', Trim(word),'". Valid options are "tags" and "amounts".&
                                         & Please refer to the manual'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        Else
          Write (messages(2),'(1x,a)')    'The correct structure for the block must be:'
          Write (messages(3),'(1x,a)')    '&input_composition'
          Write (messages(4),'(1x,a)')    '  tags      tg1    tg2    tg3   .... tgNsp'
          Write (messages(5),'(1x,a)')    '  amounts   N_tg1  N_tg2  N_tg3 .... N_tgNsp'
          Write (messages(6),'(1x,a)')    '&end_input_composition'
          Write (messages(7),'(1x,a,i3)') 'where, from &species_components, in this case Nsp = ', model_data%total_tags
          Write (messages(8),'(1x,a)')    'See manual for details'
          Call info(messages, 8)
          Call error_stop(' ')
        End If
      End If
      If (error_duplication) Then
        Write (messages(2),'(1x,3a)') 'Descriptor "', Trim(word), '" is duplicated within the block'
        Call info(messages, 2)
        Call error_stop(' ')
      End If
    End Do 

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&input_composition') 
      Call capital_to_lower_case(word)
      If (word /= '&end_input_composition') Then
        If (word(1:1) /= '#') Then
            Write (messages(2),'(3a)') 'Descriptors "tags" and "amounts" have already been defined. Directive "',&
                                    & Trim(word), '" is not valid. Block must be&
                                    & closed with sentence &end_input_composition.' 
            Call info(messages,2)
            Call error_stop(' ')
        End If
      Else
        endblock=.True.
      End If
    End Do

    ! Check consistency of the tags with those defined in species_components
    ! Assing the chemical elements to model_data%input%element
    Do i=1, model_data%total_tags
      loop=.True.
      j=1
      Do While (j <= model_data%num_species%value .And. loop)
        k=1
        Do While (k <= model_data%species_info(j)%num_components .And. loop)
          If (Trim(model_data%species_info(j)%component%tag(k))==model_data%component%tag(i)) Then
            model_data%component%element(i)=Trim(model_data%species_info(j)%component%element(k))
            loop=.False.
          End If   
          k=k+1
        End Do
        j=j+1        
      End Do 
      If (loop) Then
        Write (messages(2),'(3a)') 'Tag "', Trim(model_data%component%tag(i)), '" has not been defined in&
                                 & &species_components. Please check'
        Call info(messages,2)
        Call error_stop(' ')
      End If
    End Do

 
    ! Check if the number of atoms are correct
    Do i=1, model_data%total_tags
      If (model_data%component%N0(i)< 0) Then
        Write (messages(2),'(3a)') 'Tag "', Trim(model_data%component%tag(i)), '" CANNOT have a negative number of&
                                 & atoms within the input structure! Please correct'
        Call info(messages,2)
        Call error_stop(' ')
      End If
    End Do

   ! Calculate the number of total atoms set in the block
    model_data%component%numtot=0
    Do i=1, model_data%total_tags
      model_data%component%numtot=model_data%component%numtot+model_data%component%N0(i)
    End Do

   ! Assing atomic numbers
   Do i=1, model_data%total_tags
     loop=.True.
     k=1
     Do While (k <= NPTE .And. loop)
       If (Trim(chemsymbol(k))==Trim(model_data%component%element(i))) Then
         loop=.False.
         model_data%component%atomic_number(i)=k
       End If
       k=k+1
     End Do
   End Do

  End Subroutine read_input_composition

  
  Subroutine read_input_cell(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the simulation cell vectors of the input structure
    ! defined in &input_cell
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),   Intent(In   ) :: iunit
    Type(model_type),   Intent(InOut) :: model_data
   
    Integer(Kind=wi)   :: io, i, j
    Character(Len=64 ) :: error_input_cell
    Character(Len=256) :: messages(2), word
    Logical            :: endblock

    error_input_cell = '***ERROR in &input_cell of SETTINGS file'
    Write (messages(1),'(a)') error_input_cell

    i=1
   
    Do While (i <= 3)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&input_cell')
      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          Call check_for_rubbish(iunit, '&input_cell') 
          Read (iunit, Fmt=*, iostat=io) (model_data%input%cell(i,j), j=1,3)
          If (io/=0) Then
            Write (messages(2),'(a,i2)') 'Problems with the definition of cell vector', i
            Call info(messages, 2)
            Call error_stop(' ')
          End If
          i=i+1
        Else
          Write (messages(2),'(1x,a)') 'End of block found! Not all the cell vectors for the&  
                                     & input structure have been defined. Please check.'
          Call info(messages, 2)
          Call error_stop(' ')
        End If
      End If
    End Do
 
    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&input_cell')
      Call capital_to_lower_case(word)
      If (word /= '&end_input_cell') Then
        If (word(1:1) /= '#') Then
          Write (messages(2),'(a)') 'Block for cell vectors must be closed with&
                                  & sentence &end_input_cell.'
          Call info(messages,2)
          Call error_stop(' ')
        End If
      Else
          endblock=.True.
      End If
    End Do

  End Subroutine read_input_cell
  

  Subroutine check_species(files, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the format and directives of &species
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(InOut) :: model_data

    Integer(Kind=wi)  ::  i, length
    Logical  ::  error
    Character(Len=256)  :: messages(3)
    Character(Len=64 )  :: error_species
    Character(Len=64 )  :: error_set

    error_set = '***ERROR in file '//Trim(files(FILE_SET)%filename)//' -'
    error_species = '***ERROR in &species of file '//Trim(files(FILE_SET)%filename)

    error=.False.

    ! Check paramentes in species
    If (.Not. model_data%species%fread) Then
      Write (messages(1),'(4a)') Trim(error_set), ' &species, needed for "', Trim(model_data%analysis%type),&
                               &'" analysis, not found'
      error=.True.
      Call info(messages,1)
      Call error_stop(' ')
    End If
    
    Write (messages(1),'(a)') error_species

    Do i=1, model_data%num_species%value
      Call check_for_symbols(model_data%species_info(i)%tag, '{[()]}', model_data%species_info(i)%fail)
      If (model_data%species_info(i)%fail) Then
        Write (messages(2),'(3a)') 'Species name ', Trim(model_data%species_info(i)%tag), ' contains brackets that MUST be&
                                 & avoided. Please remove all bracket or rename the species species.'

        error=.True.
        Call info(messages,2)
      End If

      Call get_word_length(model_data%species_info(i)%tag,length)
      If (length > max_length_name_species) Then
        Write (messages(2),'(3a,i2,a)') 'Species name ', Trim(model_data%species_info(i)%tag), &
                                       & ' exceeds the maximum number of ', max_length_name_species,&
                                       & ' set by default. Please rename the species using a shorter string.'

        error=.True.
        Call info(messages,2)
      End If

      If (model_data%species_info(i)%N0_target < 1) Then
        Write (messages(2),'(3a)') 'The number of species ', Trim(model_data%species_info(i)%tag), &
                                & ' must be equal or larger than 1.'
        error=.True.
        Call info(messages,2)
      End If

      If (Trim(model_data%species_info(i)%vartype0) /= 'fixed' .And. Trim(model_data%species_info(i)%vartype0) /= 'variable') Then
        Write (messages(2),'(3a)') 'The content type for species ', Trim(model_data%species_info(i)%tag),&
                                 & ' must be either "fixed" or "variable".'
        error=.True.
        Call info(messages,2)
      End If
    End Do


    If (error) Then
      Call info('REVIEW THE STRUCURE OF THE BLOCK!', 1)
      Call error_stop('')
    End If


  End Subroutine check_species
  
  Subroutine check_components_species(files, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the format and directives for the atoms defined for 
    ! each component in &species_components
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(InOut) :: model_data

    Integer(Kind=wi)  ::  i, j, k, m
    Character(Len=256)  :: messages(4)
    Character(Len=64 )  :: error_species_components
    Logical :: loop

    error_species_components = '***ERROR in &species_components of file '//Trim(files(FILE_SET)%filename)

    Write (messages(1),'(a)') error_species_components

    ! Check correctness in the definition for the atomic components inside the same species
    Do i=1, model_data%num_species%value
      Do j= 1, (model_data%species_info(i)%num_components-1)
       Do k=j+1, model_data%species_info(i)%num_components
          If (Trim(model_data%species_info(i)%component%tag(j)) == &
            Trim(model_data%species_info(i)%component%tag(k)) ) Then
            Write (messages(2),'(2(a,i2),3a)') 'Tags for components ', j, ' and ', k, ' of species "', &
                                             & Trim(model_data%species_info(i)%tag), '" are identical.'
            Write (messages(3),'(a)')          'Please use different local tags for components of the same species'
            Call info(messages, 3)
            Call error_stop(' ')
          End If 
          If (Trim(model_data%species_info(i)%component%element(j)) == &
            Trim(model_data%species_info(i)%component%element(k)) ) Then
            Write (messages(2),'(2(a,i2),3a)') 'Chemical element for components ', j, ' and ', k, ' of species "', &
                                             & Trim(model_data%species_info(i)%tag), '" are identical.'
            Write (messages(3),'(a)')          'Please group atoms of the same type in only one component'
            Call info(messages, 3)
            Call error_stop(' ')
          End If 
        End Do
      End Do

      If (Trim(model_data%species_info(i)%topology) /= 'electrode' .And. & 
        Trim(model_data%species_info(i)%topology) /= 'molecule' .And. &  
        Trim(model_data%species_info(i)%topology) /= 'atom' ) Then
        Write (messages(2),'(4a)') 'Incorrect definition for the topology of species "', Trim(model_data%species_info(i)%tag), &
                                 & '": ', Trim(model_data%species_info(i)%topology)
                      
        Write (messages(3),'(a)')  'Valid options are: electrode, molecule or atom'
        Call info(messages, 3)
        Call error_stop(' ')
      End If

      If (Trim(model_data%species_info(i)%topology) == 'electrode') Then
        If (Trim(model_data%species_info(i)%vartype0)/='fixed') Then
          Write (messages(2),'(5a)') 'Inconsistent definition for species "', Trim(model_data%species_info(i)%tag), &
                                 & '".  If the topology is chosen to be "', Trim(model_data%species_info(i)%topology),&
                                 & '", it must correspond to a "fixed" variable during the reaction'
                      
          Call info(messages, 2)
          Call error_stop(' ')
        End If
      End If

 
      !Check if global tags are consistent with the elements of the periodic table  
      Do j= 1, model_data%species_info(i)%num_components
        loop=.True.
        k=1
        Do While (k <= NPTE .And. loop)
          If (Trim(chemsymbol(k))==Trim(model_data%species_info(i)%component%element(j))) Then
           model_data%species_info(i)%component%atomic_number(j)=k
           loop=.False.
          End If
          k=k+1
        End Do
        If (loop) Then 
          Write (messages(2),'(3a,i2,3a)') 'Chemical element "' , Trim(model_data%species_info(i)%component%element(j)), &
                                         & '" defined for component ',  j, ' of species "', Trim(model_data%species_info(i)%tag), &
                                         & '" does not correspond to any element of the Periodic Table. Please use a valid chemical&
                                         & element'
          Call info(messages, 2)
          Call error_stop(' ')
        End If
      End Do
    End Do

    ! Check correctness in the definition for the atomic components inside the same species
    Do i=1, model_data%num_species%value-1
      Do j = i+1, model_data%num_species%value 
        Do k=1, model_data%species_info(i)%num_components 
          Do m=1, model_data%species_info(j)%num_components
            If (Trim(model_data%species_info(i)%component%tag(k)) == &
              Trim(model_data%species_info(j)%component%tag(m)) ) Then
              Write (messages(2),'(3a,i2,a,i2,5a)')  'Tag ', Trim(model_data%species_info(i)%component%tag(k)), &
                                          & ' is the same for components ', k, ' and ' , m , ' of species ',            &
                                          & Trim(model_data%species_info(i)%tag), ' and ', Trim(model_data%species_info(j)%tag),&
                                          & ', respectively.' 
              Write (messages(3),'(a)')   'Please ammend! Tags for the components of species must be unique and not used&
                                          & in multiple species.'
              Call info(messages, 3)
              Call error_stop(' ')      
            End If
          End Do
        End Do        
      End Do
    End Do
 
  End Subroutine check_components_species

  Subroutine check_atomic_settings(files, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the format and directives for building atomistic models
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(InOut) :: model_data

    Character(Len=256)  :: messages(8), message
    Character(Len=64 )  :: error_set
    Integer(Kind=wi)    :: j
    Logical             :: fprint

    error_set = '***ERROR in file '//Trim(files(FILE_SET)%filename)//' -'

    ! block species components
    If (model_data%species_components%fread) Then
      Call check_components_species(files, model_data) 
    Else 
      Write (message,'(2(1x,a))') Trim(error_set), 'Missing &species_components, which is required to specify the&
                                                      & atomistic details for the components of each chemical species.' 
      Call error_stop(message)
    End If

    ! check &input_composition
    If (.Not. model_data%input_composition%fread) Then
      Write (messages(1),'(2a)')  Trim(error_set),&
                              & 'Building atomistic models requires the specification of &input_composition'
      Call info(messages,1)
      Call error_stop(' ')
    End If
    
    ! Check "input_model_format"
    If (model_data%input_model_format%fread) Then
      If (model_data%input_model_format%fail) Then
        Write (messages(1),'(2a)')  Trim(error_set), ' Wrong specification for directive "input_model_format"'
        Call info(messages,1)
        Call error_stop(' ') 
      Else
        If (Trim(model_data%input_model_format%type) /= 'vasp'    .And. &
          Trim(model_data%input_model_format%type) /= 'onetep'   .And. &
          Trim(model_data%input_model_format%type) /= 'castep'   .And. &
          Trim(model_data%input_model_format%type) /= 'xyz') Then
          Write (messages(1),'(2a)') Trim(error_set), ' Specification for directive "input_model_format" should either be:'
          Write (messages(2),'(1x,a)') '- vasp'
          Write (messages(3),'(1x,a)') '- onetep'
          Write (messages(4),'(1x,a)') '- castep'
          Write (messages(5),'(1x,a)') '- xyz'
          Call info(messages, 5)
          Call error_stop(' ') 
        End If
      End If
    Else
      Write (messages(1),'(2a)')  Trim(error_set), 'Atomistic analysis requires specification of directive&
                                                      & "input_model_format".'
      Call info(messages,1)
      Call error_stop(' ') 
    End If

    ! Check "output_model_format"
    If (model_data%output_model_format%fread) Then
      If (model_data%output_model_format%fail) Then
        Write (messages(1),'(2a)')  Trim(error_set), ' Wrong specification for directive "output_model_format"'
        Call info(messages,1)
        Call error_stop(' ') 
      Else
        If (Trim(model_data%output_model_format%type) /= 'vasp' .And. &
          Trim(model_data%output_model_format%type) /= 'onetep'  .And. &
          Trim(model_data%output_model_format%type) /= 'castep'   .And. &
          Trim(model_data%output_model_format%type) /= 'cp2k'   .And. &
          Trim(model_data%output_model_format%type) /= 'xyz') Then
          Write (messages(1),'(2a)') Trim(error_set), ' Specification for directive "output_model_format" should either be:'
          Write (messages(2),'(1x,a)') '- vasp'
          Write (messages(3),'(1x,a)') '- onetep'
          Write (messages(4),'(1x,a)') '- castep'
          Write (messages(5),'(1x,a)') '- cp2k'
          Write (messages(6),'(1x,a)') '- xyz'
          Call info(messages, 6)
          Call error_stop(' ') 
        End If
      End If
    Else
      Write (messages(1),'(2a)')  Trim(error_set), ' Requested analysis needs the specification of&
                                                      & directive "output_model_format"'
      Call info(messages,1)
      Call error_stop(' ') 
    End If

    ! Need to read &input_cell only if input_model_format is "xyz"
    If (Trim(model_data%input_model_format%type) == 'xyz' ) Then
      If (.Not. model_data%input_cell%fread) Then
        Write (messages(1),'(4a)')  Trim(error_set), ' Format ', Trim(model_data%input_model_format%type),&
                                                      & ' for the input model&
                                                      & requires the specification of the cell vectors via &input_cell'
        Call info(messages,1)
        Call error_stop(' ') 
      End If
    Else
      If (model_data%input_cell%fread) Then 
        Write (messages(1),'(4a)') ' Specification for cell vectors in &input_cell will be neglected. Cell vectors must&
                                 & be specified in file ', Trim(files(FILE_SET)%filename), ' according to the format',&
                                 & Trim(model_data%input_model_format%type)  
        Call info(messages,1)
      End If
    End If

    ! Check if directive "normal_vector" has been defined correctly
    If (model_data%normal_vector%fread) Then
      If (model_data%normal_vector%fail) Then
        Write (messages(1),'(2a)')  Trim(error_set), ' Wrong specification for directive "normal_vector"'
        Call info(messages,1)
        Call error_stop(' ')
      Else
        If (Trim(model_data%normal_vector%type) /= 'c1' .And. &
            Trim(model_data%normal_vector%type) /= 'c2' .And. &
            Trim(model_data%normal_vector%type) /= 'c3') Then
           Write (messages(1),'(2a)') Trim(error_set), ' Specification for directive "normal_vector" should either be &
                                   &"c1", "c2" or "c3", which refer to the three cell vectors'
          Call info(messages,1)
          Call error_stop(' ')
        End If
      End If
    Else
      Write (messages(1),'(2a)')  Trim(error_set), ' The user must specify the "normal_vector" directive, which corresponds to&
                                & the axis normal to the electrode surface. Options: c1, c2, or c3'
      Call info(messages,1)
      Call error_stop(' ')
    End If

    ! repeat_input_model
    If (model_data%repeat_input_model%fread) Then
      If (model_data%repeat_input_model%fail) Then
        Write (message,'(2(1x,a))') Trim(error_set),&
                                & 'Wrong (or missing) settings for "repeat_input_model" directive (see manual).'
        Call error_stop(message)
      Else
        Do j = 1, 3
          If (model_data%repeat_input_model%value(j) < 1) Then
            Write (message,'(1x,2a, i2,a)') Trim(error_set), &
                                    &' Input value', j, ' for directive "repeat_input_model" must be positive'
            Call error_stop(message)
          End If
        End Do 
      End If
    Else  
      model_data%repeat_input_model%value=1
    End If

    If (model_data%repeat_input_model%fread) Then
      ! Check consistency between normal vector and repeat_input_model
      If (Trim(model_data%normal_vector%type) == 'c1' .And. model_data%repeat_input_model%value(1) /= 1) Then 
        Write (message,'(2a)')  Trim(error_set), ' Since the normal_vector directive is set to "c1", the FIRST component&
                                 & "repeat_input_model" must be equal to 1!'    
        Call error_stop(message)                         
      Else If (Trim(model_data%normal_vector%type) == 'c2' .And. model_data%repeat_input_model%value(2) /= 1) Then 
         Write (message,'(2a)')  Trim(error_set), ' Since the "normal_vector" directive is set to "c2", the SECOND component&
                                 & of "repeat_input_model" must be equal to 1!'    
         Call error_stop(message)
      Else If (Trim(model_data%normal_vector%type) == 'c3' .And. model_data%repeat_input_model%value(3) /= 1) Then 
         Write (message,'(2a)')  Trim(error_set), ' Since the "normal_vector" directive is set to "c3", the THIRD component&
                                 & of "repeat_input_model" must be equal to 1!'    
         Call error_stop(message)         
      End If
    End If
    
    
    ! multiple_input_atoms 
    If (model_data%multiple_input_atoms%fread) Then
      If (model_data%multiple_input_atoms%fail) Then
        Write (message,'(2(1x,a))') Trim(error_set),&
                                & 'Wrong (or missing) settings for "multiple_input_atoms" directive.&
                                & Is it not too large?'
        Call error_stop(message)
      Else
        If (model_data%multiple_input_atoms%value <= 0) Then
          Write (message,'(1x,2a)') Trim(error_set), ' Directive "multiple_input_atoms" must be larger than 0.'
          Call error_stop(message)
        End If
      End If
    Else
      model_data%multiple_input_atoms%value=10000
    End If

    ! multiple_output_atoms 
    If (model_data%multiple_output_atoms%fread) Then
      If (model_data%multiple_output_atoms%fail) Then
        Write (message,'(2(1x,a))') Trim(error_set),&
                                & 'Wrong (or missing) settings for "multiple_output_atoms" directive.&
                                & Is it not too large?'
        Call error_stop(message)
      Else
        If (model_data%multiple_output_atoms%value < 2) Then
          Write (message,'(1x,2a)') Trim(error_set), ' Directive "multiple_output_atoms" must be larger than 2'
          Call error_stop(message)
        End If
      End If
    Else
      model_data%multiple_output_atoms%value=2
    End If

    ! Set if centering the structure or not
    If (model_data%centre_electrode%fread) Then
      If (model_data%centre_electrode%fail) Then
        Write (message,'(2(1x,a))') Trim(error_set), 'Wrong settings for "centre_electrode" directive.'
        Call error_stop(message)
      End If
    Else
      model_data%centre_electrode%stat=.True.
    End If

    ! Both surfaces 
    If (model_data%both_surfaces%fread) Then
      If (model_data%both_surfaces%fail) Then
        Write (message,'(2(1x,a))') Trim(error_set), 'Wrong settings for "both_surfaces" directive.'
        Call error_stop(message)
      End If
    Else
      model_data%both_surfaces%stat=.True.
    End If
    
    ! Set if size will be optimised or not
    If (model_data%optimise_size%fread) Then
      If (model_data%optimise_size%fail) Then
        Write (message,'(2(1x,a))') Trim(error_set), 'Wrong settings for "optimise_size" directive.'
        Call error_stop(message)
      End If
      If (model_data%optimise_size%stat .And. (.Not. model_data%both_surfaces%stat)) Then
        Write (message,'(2(1x,a))') Trim(error_set), 'To optimise the cell size the user must set&
                                   & "both_surfaces" to .True.'
        Call error_stop(message)      
      End If
    Else
      model_data%optimise_size%stat=.True.
    End If    
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    If (model_data%analysis%type == 'build_model') Then
      ! check directive arrangement_added_species
      If (model_data%arr_added_species%fread) Then
        If (model_data%arr_added_species%fail) Then
          Write (messages(1),'(2a)')  Trim(error_set), ' Wrong (or missing) specification for directive "arrangement_added_species"'
          Call error_stop(messages(1))
        Else
          If (Trim(model_data%arr_added_species%type) /= 'deposited' .And. &
              Trim(model_data%arr_added_species%type) /= 'random') Then
            Write (messages(1),'(2a)') Trim(error_set), ' Wrong specification for directive "arrangement_added_species". It should&
                                     & be either "deposited" or "random"'
            Call error_stop(messages(1))
          End If
        End If
      Else
        Write (messages(1),'(a,1x,a)')  Trim(error_set), 'Specification of directive "arrangement_added_species" is required.'
        Call error_stop(messages(1))
      End If    
      
      ! rotate_species
      If (model_data%rotate_species%fread) Then
        If (model_data%rotate_species%fail) Then
          Write (message,'(2(1x,a))') Trim(error_set), 'Wrong specification for directive "rotate_species". It must be &
                                    & either .True. or .False. (.True. by default)'
          Call error_stop(message)
        End If
      Else
        model_data%rotate_species%stat=.True.
      End If
      
      ! delta_space
      If (model_data%delta_space%fread) Then
        If (model_data%delta_space%fail) Then
          Write (message,'(2(1x,a))') Trim(error_set), 'Wrong settings for "delta_space" directive.'
          Call error_stop(message)
        Else
          If (model_data%delta_space%value < epsilon(model_data%delta_space%value)) Then
            If (fprint) Call info(' ',1)
            Write (message,'(2(1x,a))') Trim(error_set), &
                                      &'Input value for delta_space MUST be larger than zero'
            Call error_stop(message)
          End If
          If (Trim(model_data%delta_space%units) == 'angstrom' ) Then
            model_data%delta_space%convert=1.0_wp
          Else If (Trim(model_data%delta_space%units) == 'bohr' ) Then
            model_data%delta_space%convert= Bohr_to_A
          Else
            Write (message,'(2a)')  Trim(error_set), 'Invalid units for directive "delta_space". Options: Angstrom or Bohr'
            Call info(message, 1)
            Call error_stop(' ') 
          End If
          model_data%delta_space%units='Angstrom' 
          model_data%delta_space%value=model_data%delta_space%value*&
                                                     model_data%delta_space%convert
        End If
      Else
        model_data%delta_space%units='Angstrom'
        model_data%delta_space%value= 0.1_wp 
      End If
      
      ! distance cutoff
      If (model_data%distance_cutoff%fread) Then
        If (model_data%distance_cutoff%fail) Then
          Write (message,'(2(1x,a))') Trim(error_set), 'Wrong settings for "distance_cutoff" directive.'
          Call error_stop(message)
        Else
          If (model_data%distance_cutoff%value < epsilon(model_data%distance_cutoff%value)) Then
            If (fprint) Call info(' ',1)
            Write (message,'(2(1x,a))') Trim(error_set), &
                                      &'Input value for "distance_cutoff" MUST be larger than zero'
            Call error_stop(message)
          End If
          If (Trim(model_data%distance_cutoff%units) == 'angstrom' ) Then
            model_data%distance_cutoff%convert=1.0_wp
          Else If (Trim(model_data%distance_cutoff%units) == 'bohr' ) Then
            model_data%distance_cutoff%convert= Bohr_to_A
          Else
            Write (message,'(2a)')  Trim(error_set), 'Invalid units for directive "distance_cutoff".&
                                  & Options: Angstrom or Bohr'
            Call info(message, 1)
            Call error_stop(' ') 
          End If
          model_data%distance_cutoff%units='Angstrom' 
          model_data%distance_cutoff%value=model_data%distance_cutoff%value*&
                                                     model_data%distance_cutoff%convert
          If (model_data%distance_cutoff%value < min_inter_bond) Then
            Write (message,'(1x,a,2(a,f4.2),a)') Trim(error_set), &
                                            & ' The minimum limit for the inter-species distance of ',&
                                            & model_data%distance_cutoff%value,&
                                            &' Angstrom (set via "distance_cutoff" directive) is unphysical.&
                                            & The value MUST NOT be lower than ', min_inter_bond, ' Angstrom. Please change'
            Call info(message,1)
            Call error_stop(' ')
          End If
        End If
      Else
        Write (messages(1),'(2a)')  Trim(error_set), ' The user must specify the "distance_cutoff" directive, which&
                                  & corresponds to the minimum distance between atoms of different species'
        Call info(messages,1)
        Call error_stop(' ')
      End If
      
      ! level where to start adding the electrolyte from
      If (model_data%add_species_from%fread) Then
        If (model_data%add_species_from%fail) Then
          Write (message,'(2(1x,a))') Trim(error_set), 'Wrong settings for "add_species_from" directive.'
          Call error_stop(message)
        Else
          If (Trim(model_data%add_species_from%units) == 'angstrom' ) Then
            model_data%add_species_from%convert=1.0_wp
          Else If (Trim(model_data%add_species_from%units) == 'bohr' ) Then
            model_data%add_species_from%convert= Bohr_to_A
          Else
            Write (message,'(2a)')  Trim(error_set), 'Invalid units for directive "add_species_from".&
                                  & Options: Angstrom or Bohr'
            Call info(message, 1)
            Call error_stop(' ') 
          End If
          model_data%add_species_from%units='Angstrom' 
          model_data%add_species_from%value=&
                    & model_data%add_species_from%value*model_data%add_species_from%convert
        End If
      End If

      ! add extra space...
      If (model_data%add_extra_space%fread) Then
        If (model_data%add_extra_space%fail) Then
          Write (message,'(2(1x,a))') Trim(error_set), 'Wrong settings for "add_extra_space" directive.'
          Call error_stop(message)
        Else
          If (Trim(model_data%add_extra_space%units) == 'angstrom' ) Then
            model_data%add_extra_space%convert=1.0_wp
          Else If (Trim(model_data%add_extra_space%units) == 'bohr' ) Then
            model_data%add_extra_space%convert= Bohr_to_A
          Else
            Write (message,'(2a)')  Trim(error_set), 'Invalid units for directive "add_extra_space".&
                                  & Options: Angstrom or Bohr'
            Call info(message, 1)
            Call error_stop(' ') 
          End If
          model_data%add_extra_space%units='Angstrom' 
          model_data%add_extra_space%value=&
                    & model_data%add_extra_space%value*model_data%add_extra_space%convert
        End If
        If (.Not. model_data%both_surfaces%stat) Then
          Write (message,'(2a)')  Trim(error_set), 'Definition of "add_extra_space" is only valid if&
                                & "both_surfaces" is set to True. '
          Call info(message, 1)
          Call error_stop(' ') 
        End If 
        If (.Not. model_data%optimise_size%stat) Then
          Write (message,'(2a)')  Trim(error_set), 'Definition of "add_extra_space" is only valid if&
                                & "optimise_size" is set to True. '
          Call info(message, 1)
          Call error_stop(' ') 
        End If 
      Else
        If (model_data%both_surfaces%stat) Then
          model_data%add_extra_space%units='Angstrom' 
          model_data%add_extra_space%value=6.0_wp
        End If
      End If
    End If
   
  End Subroutine check_atomic_settings

  Subroutine read_input_model(files, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read input model from file INPUT_GEOM/INPUT_ELECTRODE
    ! Format Options: 
    ! - VASP
    ! - xyz 
    ! - ONETEP
    ! - CASTEP  
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(model_type),  Intent(InOut) :: model_data

    Logical            :: safe
    Integer(Kind=wi)   :: ifolder
    Character(Len=256) :: message
    Character(Len=256) :: messages(3)
    Character(Len=32 ) :: input_file, set_error, path_to_file

    input_file=Trim(files(FILE_INPUT_ELECTRODE)%filename)

    ! Check if folder INPUT_GEOM exists
    Call execute_command_line('[ -d '//Trim(FOLDER_INPUT_GEOM)//' ]', exitstat=ifolder)
    If (ifolder/=0) Then
      Call info(' ', 1)
      Write (messages(1), '(1x,3a)') '***ERROR: folder ', Trim(FOLDER_INPUT_GEOM), ' cannot be found.'
      Write (messages(2), '(1x,3a)') 'This folder must contain file ', Trim(input_file),' and the xyz&
                                   & files only for the participant MOLECULAR species (if any).'
      Write (messages(3), '(1x,a)') 'The requested analysis cannot be conducted. Please create the folder&
                                   & and add the required information.'
      Call info(messages, 3)
      Call error_stop(' ')
    End If
 
    path_to_file = Trim(FOLDER_INPUT_GEOM)//'/'//Trim(files(FILE_INPUT_ELECTRODE)%filename)
    set_error = '***ERROR -'

    Inquire(File=path_to_file, Exist=safe)

    If (.not.safe) Then
      Call info(' ', 1)
      Write (message,'(4(1x,a))') Trim(set_error), 'File', Trim(path_to_file), ' not found'
      Call error_stop(message)
    End If

    ! Open the INPUT_ELECTRODE file
    Open(Newunit=files(FILE_INPUT_ELECTRODE)%unit_no, File=Trim(path_to_file),Status='old')

    ! Select the calling of the subroutine according to the specification of input_model_format directive
    If (Trim(model_data%input_model_format%type) == 'vasp') Then
        Call read_input_vasp_format(files, model_data)
    Else If (Trim(model_data%input_model_format%type) == 'xyz') Then
        Call read_input_xyz_format(files, model_data)
    Else If (Trim(model_data%input_model_format%type) == 'castep' .Or. & 
             Trim(model_data%input_model_format%type) == 'onetep') Then
        Call read_input_geom_format(files, model_data)  
    End If

    Call about_cell(model_data%input%cell,model_data%input%invcell,model_data%input%cell_length)

    ! Close file
    Close(files(FILE_INPUT_ELECTRODE)%unit_no) 
 
  End Subroutine read_input_model  

  Subroutine read_input_geom_format(files, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read file INPUT_GEOM/INPUT_ELECTRODE according to either 
    ! CASTEP or ONETEP format 
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(InOut) :: model_data

    Integer(Kind=wi) :: io, num_geo, num, iunit
    Integer(Kind=wi) :: internal
    Integer(Kind=wi) :: i, j, k, m
    Character(Len=256)    :: messages(2)
    
    Character(Len=256)    :: set_error, input_file, path_to_file
    Character(Len=256)    :: exec_numgeo, word

    input_file = Trim(files(FILE_INPUT_ELECTRODE)%filename)
    path_to_file=Trim(FOLDER_INPUT_GEOM)//'/'//Trim(input_file)
    iunit=files(FILE_INPUT_ELECTRODE)%unit_no

    set_error = '***ERROR in file '//Trim(path_to_file)//&
               &' (inconsistency with CASTEP/ONETEP format)'

    ! Check the correctness of the file

    exec_numgeo='grep "<-- E" '//Trim(path_to_file)//' > number_geom.dat'
    Call execute_command_line(exec_numgeo) 
    Call execute_command_line('wc -l number_geom.dat > nlines.dat')
    Open(Newunit=internal, File='nlines.dat' ,Status='old')
    Read (internal, Fmt=*, iostat=io) num_geo
    If (num_geo == 0) Then
      Call error_stop(set_error)
    Else If (num_geo > 1) Then
      Write (messages(1),'(2a,i3,a)') Trim(set_error), ': number of geometrical configurations indentified&
                                   & is equal to', num_geo, '. Please provide only one configuration.'
      Call error_stop(messages(1)) 
    End If
    Close(internal)
    Call execute_command_line('rm nlines.dat number_geom.dat')

    exec_numgeo='grep "<-- R" '//Trim(path_to_file)//' > number_atoms.dat'
    Call execute_command_line(exec_numgeo) 
    Call execute_command_line('wc -l number_atoms.dat > nat.dat')
    Open(Newunit=internal, File='nat.dat' ,Status='old')
    Read (internal, Fmt=*, iostat=io) model_data%input%num_atoms
    If (model_data%input%num_atoms == 0) Then
      Call error_stop(set_error)
    End If
    Close(internal)
    Call execute_command_line('rm nat.dat number_atoms.dat')

    If (model_data%input%num_atoms /= model_data%component%numtot) Then
     Write (messages(1),'(3a)') 'ERROR***: Inconsistency between the number of atoms in file ',&
                             & Trim(path_to_file),&
                             &' and the amount of atoms defined in &input_composition. Please check'   
     Call error_stop(messages(1))
    End If

    ! Read the file
    Read (iunit, Fmt=*, iostat=io) word  
    Read (iunit, Fmt=*, iostat=io) word  
    Do i= 1, 3  
      Read (iunit, Fmt=*, iostat=io) (model_data%input%cell(i,j), j=1,3) 
      If (io/=0) Then
        Write (messages(1),'(2a,i1,a)') Trim(set_error), ': Problems with the specification of cell vector ', i,&
                                     & '. See those lines that end with "<-- h"'  
        Call error_stop(messages(1))
      End If 
      model_data%input%cell(i,:)=Bohr_to_A*model_data%input%cell(i,:) 
    End Do

    ! Allocate atomic arrays for the input model   
    Call model_data%atomic_arrays_input(model_data%input%num_atoms)  

    ! Read atomic coordinates
    j=1; k=0; i=0
    Do While (i < model_data%input%num_atoms)
      If (model_data%component%N0(j)/=0) Then
        i=i+1
        Read (iunit, Fmt=*, iostat=io) model_data%input%atom(i)%element, num, (model_data%input%atom(i)%r(m), m=1,3)
        If (io/=0) Then
          Write (messages(1),'(2a,i5)') Trim(set_error), ' Wrong specification for atom', i
          Call info(messages, 1)
          Call error_stop(' ')
        Else
          If (Trim(model_data%input%atom(i)%element)/=Trim(model_data%component%element(j))) Then
            Write (messages(1),'(a,i5,5a)') '***ERROR: Chemical element of atom ', i, &
                                          & ' is set to "', Trim(model_data%input%atom(i)%element), '", but& 
                                          & according to the definition of &input_composition it must be "',&
                                          & Trim(model_data%component%element(j)), '".'
            Write (messages(2),'(3a)')      'Check the consistency between the data in ', &
                                          &  Trim(path_to_file), ' and &input_composition'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        End If
        ! Transform to Angstom
        model_data%input%atom(i)%r=Bohr_to_A*model_data%input%atom(i)%r 
        k=k+1
        model_data%input%atom(i)%tag=model_data%component%tag(j)
        If (k==model_data%component%N0(j)) Then
          k=0
          j=j+1
        End If
      Else
        k=0
        j=j+1
      End If
    End Do

  End Subroutine read_input_geom_format 


  Subroutine read_input_xyz_format(files, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read file INPUT_ELECTRODE according to xyz format 
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(InOut) :: model_data

    Integer(Kind=wi) :: io, iunit
    Integer(Kind=wi) :: i, j, k, m
    Character(Len=256)    :: title
    Character(Len=256)    :: messages(5)
    
    Character(Len=256)    :: set_error, error_format(10), input_file, path_to_file

    input_file = Trim(files(FILE_INPUT_ELECTRODE)%filename) 
    path_to_file = Trim(FOLDER_INPUT_GEOM)//'/'//Trim(input_file)    
    iunit=files(FILE_INPUT_ELECTRODE)%unit_no

    set_error = '***ERROR in file '//Trim(path_to_file)//' (inconsistency with xyz format):'
    error_format(1) ='The structure of the '//Trim(input_file)//' file in xyz format must be:'
    error_format(2) = ' '
    error_format(3) = 'Number of atoms (Nat)'
    error_format(4) = 'Comment with the description of the system (compulsory)'
    error_format(5) = 'Element_1         X_1      Y_1      Z_1'
    error_format(6) = 'Element_2         X_2      Y_2      Z_2'
    error_format(7) = '...........       .....    .....    .....'
    error_format(8) = 'Element_Nat       X_Nat    Y_Nat    Z_Nat'
    error_format(9) = ' ' 
    error_format(10) = 'Please check consistency between the structure of the file and&
                       & directive "output_model_format".'
    ! Start reading the file
    !!!!!!!!!!!!!!!!!!!!!!!!

    ! Read number of atoms 
    Read (iunit, Fmt=*, iostat=io) model_data%input%num_atoms
    If (io/=0) Then
      Write (messages(1),'(2a)') Trim(set_error), ' Invalid specification for the number of atoms, which must be&
                               & specified in the first line.'
      Call info(messages,1)
      Call info(error_format,10)
      Call error_stop(' ')
    End If    

    If (model_data%input%num_atoms /= model_data%component%numtot) Then
     Write (messages(1),'(3a)') 'ERROR***: Inconsistency between the number of atoms in file ', Trim(path_to_file),&
                             &' and the amount of atoms defined in &input_composition. Please check.'   
     Call error_stop(messages(1))
    End If

    ! Read title
    Read (iunit, Fmt=*, iostat=io) title

    ! Allocate atomic arrays for the input model   
    Call model_data%atomic_arrays_input(model_data%input%num_atoms)  

    ! Read atomic coordinates
    j=1; k=0; i=0
    Do While (i < model_data%input%num_atoms)
      If (model_data%component%N0(j)/=0) Then
        i=i+1
        Read (iunit, Fmt=*, iostat=io) model_data%input%atom(i)%element, (model_data%input%atom(i)%r(m), m=1,3)
        If (io/=0) Then
          If (i== model_data%input%num_atoms) Then
            Write (messages(1),'(2a)') Trim(set_error), ' Missing input coordinates somewhere in the list.'
          Else
            Write (messages(1),'(2a,i5)') Trim(set_error), ' Wrong specification for atom', i
          End If
          Call info(messages, 1)
          Call info(error_format,8)
          Call error_stop(' ')
        Else
          If (Trim(model_data%input%atom(i)%element)/=Trim(model_data%component%element(j))) Then
            Write (messages(1),'(a,i5,5a)') '***ERROR: Chemical element of atom ', i, &
                                          & ' is set to "', Trim(model_data%input%atom(i)%element), '", but& 
                                          & according to the definition of &input_composition it must be "',&
                                          & Trim(model_data%component%element(j)), '".'
            Write (messages(2),'(3a)')  'Check the consistency between the data in ', Trim(path_to_file),&
                                          & ' and &input_composition.'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        End If
        k=k+1
        model_data%input%atom(i)%tag=model_data%component%tag(j)
        If (k==model_data%component%N0(j)) Then
          k=0
          j=j+1
        End If
      Else
        k=0
        j=j+1
      End If
    End Do

  End Subroutine read_input_xyz_format


  Subroutine read_input_vasp_format(files, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read file INPUT_ELECTRODE according to VASP/POSCAR format 
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(InOut) :: model_data

    Integer(Kind=wi) :: io, iunit
    Integer(Kind=wi) :: i, j, k, m
    Character(Len=256)  :: title, word
    Character(Len=256)  :: messages(5)
    Logical             :: loop, loop2, error_vasp, error_elements
    Character(Len=  2)  :: buffer
    
    Character(Len=256)  :: input_file, path_to_file, set_error, webpage

    Real(Kind=wp)       :: v_cart(3)

    Character(Len= 2), Allocatable :: element_file(:)
    Integer(Kind=wi),  Allocatable :: amount_file(:)
    Integer(Kind=wi)  :: fail(2)

    error_vasp = .False.
    input_file = Trim(files(FILE_INPUT_ELECTRODE)%filename)
    path_to_file= Trim(FOLDER_INPUT_GEOM)//'/'//Trim(input_file)   
    iunit=files(FILE_INPUT_ELECTRODE)%unit_no

    set_error = '***ERROR in file '//Trim(path_to_file)//' (inconsistency with POSCAR format):'
    webpage='For correct format and details see the following link: https://www.vasp.at/wiki/index.php/POSCAR'

    ! Determine vasp elements and amounts of atoms from the data in &input_composition and &species_components
    model_data%input%list%num_elements=1 
    j=1
    model_data%input%list%N0(j)=model_data%component%N0(j)
    model_data%input%list%element(j)=model_data%component%element(j)   
    buffer=model_data%component%element(j)
    Do i = 2,  model_data%total_tags
      If (Trim(model_data%component%element(i)) /= Trim(buffer)) Then
        buffer=model_data%component%element(i)
        If (model_data%component%N0(i)/=0) Then
          model_data%input%list%num_elements=model_data%input%list%num_elements+1
          j=j+1
          model_data%input%list%N0(j)=model_data%component%N0(i)
          model_data%input%list%element(j)=model_data%component%element(i)
        End If
      Else
        model_data%input%list%N0(j)=model_data%input%list%N0(j)+model_data%component%N0(i)
      End If
    End Do 

    model_data%input%num_atoms=0
    Do i=1, model_data%input%list%num_elements
       model_data%input%num_atoms=model_data%input%num_atoms + model_data%input%list%N0(i)
    End Do

    ! Allocate atomic arrays for the input model   
    Call model_data%atomic_arrays_input(model_data%input%num_atoms)  

    ! Allocate vasp related variables for reading the file (deallocated below)
    Allocate(element_file(model_data%input%list%num_elements), Stat=fail(1))
    Allocate(amount_file(model_data%input%list%num_elements),  Stat=fail(2))

    If (Any(fail > 0)) Then
      Write (messages(1),'(1x,1a)') '***ERROR: Allocation problems for arrays in subroutine "read_input_vasp_format"'
      Call error_stop(messages(1))
    End If

    ! Start reading the file
    !!!!!!!!!!!!!!!!!!!!!!!!
    ! Read elements or title
    Read (iunit, Fmt=*, iostat=io) (element_file(i) , i=1, model_data%input%list%num_elements)
    loop2=.True.
    i=1
    Do While (i<= model_data%input%list%num_elements .And. loop2)
      loop=.True.
      k=1
      Do While (k <= NPTE .And. loop)
        If (Trim(chemsymbol(k))==Trim(element_file(i))) Then
          loop=.False.
        End If
        k=k+1
      End Do
      If (loop) Then
        loop2=.False.
      End If
      i=i+1
    End Do

    If (loop2) Then
      error_elements=.False.
    Else
      error_elements=.True.
      Rewind iunit 
      Read (iunit, Fmt=*, iostat=io) title
    End If   
 
    ! Read scale factor
    Read (iunit, Fmt=*, iostat=io) model_data%input%scale_factor_vasp
    If (io/=0) Then
      Write (messages(1),'(2a)') Trim(set_error), ' Scale factor for cell vectors is invalid.'
      Write (messages(2),'(a)') webpage
      Call info(messages,2)
      Call error_stop(' ')
    End If 

    Do i= 1, 3
      Read (iunit, Fmt=*, iostat=io) (model_data%input%cell(i,j), j=1,3)
      If (io/=0) Then
        Write (messages(1),'(2a,i1,a)') Trim(set_error), ' Definition for cell vector ', i, ' is incorrect.'
        Write (messages(2),'(a)') webpage
        Call info(messages,2)
        Call error_stop(' ')
      End If 
    End Do

    If (error_elements) Then
      Read (iunit, Fmt=*, iostat=io) (element_file(i) , i=1, model_data%input%list%num_elements)
      If (io/=0) Then
        Write (messages(1),'(2a)') Trim(set_error), ' Invalid specification for the list of atomic species.'
        Write (messages(2),'(a)') webpage
        Call info(messages,2)
        Call error_stop(' ') 
      End If
    End If

    Do i=1, model_data%input%list%num_elements
      loop=.True.
      k=1
      Do While (k <= NPTE .And. loop)
        If (Trim(chemsymbol(k))==Trim(element_file(i))) Then
          loop=.False.
        End If
        k=k+1
      End Do
      If (loop) Then
        Write (messages(1),'(4a,i2,a)') Trim(set_error), ' Chemical element "' , Trim(element_file(i)), &
                                       & '" defined for component ',  i, ' of the list does not correspond to an element&
                                       & of the Periodic Table. Please use a valid element.'
        Write (messages(2),'(a)')   webpage
        Write (messages(3),'(3a)') 'IMPORTANT: The user should also check for inconsistencies between the list/number of atoms&
                                & in file ', Trim(input_file),' and the settings in &input_composition.&
                                & Have you missed info?'
        Call info(messages,3)
        Call error_stop(' ')
      End If

      If (Trim(element_file(i)) /= Trim(model_data%input%list%element(i))) Then
        error_vasp=.True.
      End If 
    End Do


    Read (iunit, Fmt=*, iostat=io) (amount_file(i) , i=1, model_data%input%list%num_elements)
    If (io/=0) Then
      Read (iunit, Fmt=*, iostat=io) (amount_file(i) , i=1, model_data%input%list%num_elements)
      If (io/=0) Then
        Write (messages(1),'(2a)') Trim(set_error), ' Invalid specification for the list with the number of atoms&
                                               & for each atomic species'
        Write (messages(2),'(a)') webpage
        Call info(messages,2)
        Call error_stop(' ')
      End If
    End If

    Write (messages(1),'(3a)') '***ERROR: Inconsistent between the list of atomic species in ',  Trim(path_to_file),&
                             & ' (VASP format) and the data specified in &input_composition.'
    Write (messages(2),'(a)') 'The list of elements and amount of atoms per element resulting from the definition&
                             & in &input_composition should be:'
          
    Write (messages(3),'(*(4x,a2))') (model_data%input%list%element(j), j=1, model_data%input%list%num_elements)
    Write (messages(4),'(*(1x,i5))') (model_data%input%list%N0(j), j=1, model_data%input%list%num_elements)
    Write (messages(5),'(3a)') 'Which DO NOT MATCH the order/amount of atoms in file ', Trim(path_to_file), & 
                            & '. Please review the settings of &input_composition (or the geometry of the input structure)'

    Do i=1, model_data%input%list%num_elements
      If (amount_file(i) /= model_data%input%list%N0(i)) Then
        error_vasp=.True.
      End If
    End Do

    If (error_vasp) Then
      Call info(messages, 5)
      Call error_stop(' ')
    End If

    Read (iunit, Fmt=*, iostat=io) word
    Call capital_to_lower_case(word)
    If (word(1:4)=='sele') Then
      model_data%selective_dyn=.True. 
    Else       
      Backspace iunit
    End If

    Read (iunit, Fmt=*, iostat=io) model_data%input%list%coord_type
    If (io/=0) Then
      Write (messages(1),'(2a)') Trim(set_error), ' Invalid specification for the type of coordinates. Valid options are either& 
                                             & "Cartesian" or "Direct"'
      Write (messages(2),'(a)') webpage
      Call info(messages,2)
      Call error_stop(' ')
    End If
    Call capital_to_lower_case(model_data%input%list%coord_type)   
  
    If (Trim(model_data%input%list%coord_type) /= 'cartesian' .And. &
      Trim(model_data%input%list%coord_type) /= 'direct') Then
      Write (messages(1),'(2a)') Trim(set_error), ' Wrong option for the type of coordinates. Valid options are either& 
                                             & "Cartesian" or "Direct"'
      Write (messages(2),'(a)') webpage
      Call info(messages,2)
      Call error_stop(' ')
    End If

    ! Read atomic coordinates
    j=1; k=0; i=0
    Do While (i < model_data%input%num_atoms)
      If (model_data%component%N0(j)/=0) Then
        i=i+1
        If (model_data%selective_dyn) Then
          Read (iunit, Fmt=*, iostat=io) (model_data%input%atom(i)%r(m), m=1,3), &
                                              (model_data%input%atom(i)%dynamics(m), m=1,3)
        Else 
          Read (iunit, Fmt=*, iostat=io) (model_data%input%atom(i)%r(m), m=1,3)
        End If
        If (io/=0) Then
          If (i== model_data%input%num_atoms) Then
            Write (messages(1),'(2a)') Trim(set_error), ' Missing input coordinates somewhere in the list'
          Else
            Write (messages(1),'(2a,i5)') Trim(set_error), ' Wrong specification for the input coordinates of atom ', i
          End If
          Call info(messages, 1)
          Call error_stop(' ')
        End If
        k=k+1
        model_data%input%atom(i)%tag=model_data%component%tag(j)
        model_data%input%atom(i)%element=model_data%component%element(j)
        If (k==model_data%component%N0(j)) Then
          k=0
          j=j+1 
        End If
      Else
        k=0
        j=j+1
      End If
    End Do

    ! Deallocate arrays
    Deallocate(amount_file)
    Deallocate(element_file)

    ! Transform using the scaling factor
    If (Trim(model_data%input%list%coord_type) == 'cartesian') Then
      Do i = 1 , model_data%input%num_atoms
        model_data%input%atom(i)%r=model_data%input%scale_factor_vasp*model_data%input%atom(i)%r
      End Do
    Else If (Trim(model_data%input%list%coord_type) == 'direct') Then
      Do i = 1, model_data%input%num_atoms 
        v_cart=MatMul(model_data%input%atom(i)%r, model_data%input%cell)
        model_data%input%atom(i)%r=model_data%input%scale_factor_vasp*v_cart
      End Do
    End If
    
    model_data%input%cell= model_data%input%scale_factor_vasp * model_data%input%cell

  End Subroutine read_input_vasp_format  

  Subroutine about_cell(A,invA,length)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to invert a general 3x3 matrix and compute the length of the cell
    ! vectors.
    ! Matrix A is the input matrix.
    ! Matrix invA is the output matrix (inverse of A) 
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp), Intent(In   )  :: A(3,3)
    Real(Kind=wp), Intent(  Out)  :: invA(3,3)
    Real(Kind=wp), Intent(  Out)  :: length(3)

    Real(Kind=wp) :: Det
    Real(Kind=wp) :: Cofactor(3,3)
   
    Integer(Kind=wi) :: i, j

    length = 0.0_wp     

    Do i = 1, 3 
      Do j= 1, 3
        length(i) = length(i)+ A(i,j)**2   
      End Do
      length(i)=sqrt(length(i))
    End Do 

    Det =   A(1,1)*A(2,2)*A(3,3)  &
          - A(1,1)*A(2,3)*A(3,2)  &
          - A(1,2)*A(2,1)*A(3,3)  &
          + A(1,2)*A(2,3)*A(3,1)  &
          + A(1,3)*A(2,1)*A(3,2)  &
          - A(1,3)*A(2,2)*A(3,1)

    If (Abs(Det) <= epsilon(det)) Then
      Call error_stop('***ERROR: The determinant of the simulation cell is zero. Please check the definition of the cell vectors')
    End If

    Cofactor(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
    Cofactor(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    Cofactor(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    Cofactor(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
    Cofactor(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
    Cofactor(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
    Cofactor(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
    Cofactor(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
    Cofactor(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

    invA = Transpose(Cofactor) / Det

  End Subroutine about_cell  
  
End Module atomistic_setup
