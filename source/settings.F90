!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Module that:
! - reads the SETTINGS file 
! - defines and checks the correctness of directives and blocks
!
! Copyright:  2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author     - i.scivetti March 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module settings 

  Use atomistic_setup,   Only : model_type, &
                                read_species, &
                                read_species_components, &
                                read_species_components, & 
                                read_input_composition, &
                                read_input_cell, &
                                check_species, &
                                check_components_species,&
                                check_atomic_settings

  Use constants,         Only : Bohr_to_A,  &
                                chemsymbol, & 
                                NPTE
                                
  Use fileset,           Only : file_type, &
                                FILE_SET, &  
                                FILE_OUT, &
                                FILE_RECORD_MODELS, &
                                refresh_out
                                
  Use hpc,               Only : hpc_type, &
                                read_hpc_settings, &    
                                check_hpc_settings
                                
  Use numprec,           Only : wi, &
                                wp
                                
  Use process_data,      Only : capital_to_lower_case, &
                                check_for_rubbish, &
                                get_word_length, &
                                duplication_error,     &
                                set_read_status
                                
  Use simulation_setup,  Only : simul_type, &
                                read_simulation_settings
                                
  Use unit_output,       Only : error_stop,&
                                info

  Use simulation_files_builder, Only : check_simulation_settings                              
                                 

  Implicit None
  Private

  Public :: read_settings, check_all_settings

Contains

  Subroutine read_settings(files, model_data, simulation_data, hpc_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the SETTINGS file
    ! For each directive, subroutine "set_read_fail" assigns fread=.True. and assigns 
    ! fail=.True. (fail=.False.) if the format/syntax for the directive is correct (incorrect)
    ! If the directive is repeated the execution is aborted via the duplication_error subroutine 
    ! Subroutines "set_read_fail" and "duplication_error" are defined in module process_data
    ! 
    ! author        - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),       Intent(InOut) :: files(:)
    Type(model_type),      Intent(InOut) :: model_data 
    Type(simul_type),      Intent(InOut) :: simulation_data 
    Type(hpc_type),        Intent(InOut) :: hpc_data 
 
    Logical            :: safe
    Character(Len=256) :: word
    Integer(Kind=wi)   :: i, j
    Integer(Kind=wi)   :: length, io, iunit
  
    Character(Len=256)  :: message

    Character(Len=32 )  :: set_file
    Character(Len=32 )  :: set_error

    set_file = Trim(files(FILE_SET)%filename)
    set_error = '***ERROR -'

  ! Initialise relevant arrays
    Call model_data%init_input_variables()    
    
    ! Open the file to read instructions
    Inquire(File=files(FILE_SET)%filename, Exist=safe)
    
    If (.not.safe) Then
      Call info(' ', 1)
      Write (message,'(4(1x,a))') Trim(set_error), 'File', Trim(set_file), ' not found'
      Call error_stop(message)
    Else
      Open(Newunit=files(FILE_SET)%unit_no, File=Trim(set_file), Status='old')
      iunit=files(FILE_SET)%unit_no 
    End If

     Read (iunit, Fmt=*, iostat=io) word
     ! If nothing is found, complain and abort
     If (is_iostat_end(io)) Then
       Write (message,'(3(1x,a))') Trim(set_error), Trim(set_file), 'file seems to be empty?. Please check'
       Call error_stop(message)
     End If
     ! Check header has "#" as the first character 
     If (word(1:1)/='#') Then
       Write (message,'(4(1x,a))') Trim(set_error), 'Heading comment in file', Trim(set_file), & 
                                  'is required and MUST be preceded with the symbol "#"'
       Call error_stop(message)
     End If

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Exit
      end If
      Call check_for_rubbish(iunit, Trim(set_file)) 
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
        ! Do nothing if line is a comment of we have an empty line
        Read (iunit, Fmt=*, iostat=io) word
        
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Specifications of settings to build atomic models 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
      Else If (word(1:length) == 'analysis') Then
        Read (iunit, Fmt=*, iostat=io) word, model_data%analysis%type
        Call set_read_status(word, io, model_data%analysis%fread, model_data%analysis%fail, model_data%analysis%type)

      Else If (word(1:length) == 'arrangement_added_species') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%arr_added_species%type
        Call set_read_status(word, io, model_data%arr_added_species%fread, model_data%arr_added_species%fail,&
                           & model_data%arr_added_species%type)
        
      Else If (word(1:length) == '&species') Then
        Read (iunit, Fmt=*, iostat=io) model_data%species%type
        Call set_read_status(word, io, model_data%species%fread, model_data%species%fail)
        ! Read information inside the block
        Call read_species(iunit, model_data)  

      Else If (word(1:length) == '&species_components') Then
        If (.Not. model_data%analysis%fread) Then
          Write (message,'(3(1x,a))') Trim(set_error), 'Directive "Analysis" must be specified before', &
                                     '&species_components'
          Call error_stop(message) 
        End If
        If (.Not.model_data%species%fread) Then
          Write (message,'(3(1x,a))') Trim(set_error), '&species_components must be specified after', &
                                     '&species. Have you specified &species?'
          Call error_stop(message) 
        End If

        Read (iunit, Fmt=*, iostat=io) model_data%species_components%type
        Call set_read_status(word, io, model_data%species_components%fread, model_data%species_components%fail)
        ! Read information inside the block
        Call read_species_components(iunit, model_data)  

      Else If (word(1:length) == '&input_composition') Then
        If (.Not. model_data%species_components%fread) Then
          Write (message,'(3(1x,a))') Trim(set_error), '&input_composition must be specified after', &
                                     '&species_components. Have you specified "&species_components"?'
          Call error_stop(message)
        End If
        Read (iunit, Fmt=*, iostat=io) model_data%input_composition%type
        Call set_read_status(word, io, model_data%input_composition%fread, model_data%input_composition%fail)
        ! Read information inside the block
        Call read_input_composition(iunit, model_data)

      Else If (word(1:length) == '&input_cell') Then
        Read (iunit, Fmt=*, iostat=io) model_data%input_cell%type
        Call set_read_status(word, io, model_data%input_cell%fread, model_data%input_cell%fail)
        ! Read information inside the block
        Call read_input_cell(iunit, model_data)  
 
      Else If (word(1:length) == 'input_model_format') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%input_model_format%type
        Call set_read_status(word, io, model_data%input_model_format%fread, model_data%input_model_format%fail, &
                           & model_data%input_model_format%type)

      Else If (word(1:length) == 'output_model_format') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%output_model_format%type
        Call set_read_status(word, io, model_data%output_model_format%fread, model_data%output_model_format%fail, &
                           & model_data%output_model_format%type)

      Else If (word(1:length) == 'delta_space') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%delta_space%value, model_data%delta_space%units 
        Call set_read_status(word, io, model_data%delta_space%fread, model_data%delta_space%fail)
        Call capital_to_lower_case(model_data%delta_space%units)

      Else If (word(1:length) == 'distance_cutoff') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%distance_cutoff%value, model_data%distance_cutoff%units 
        Call set_read_status(word, io, model_data%distance_cutoff%fread, model_data%distance_cutoff%fail)
        Call capital_to_lower_case(model_data%distance_cutoff%units)

      Else If (word(1:length) == 'repeat_input_model') Then 
        Read (iunit, Fmt=*, iostat=io) word, (model_data%repeat_input_model%value(j), j=1,3)
        Call set_read_status(word, io, model_data%repeat_input_model%fread, model_data%repeat_input_model%fail)

      Else If (word(1:length) == 'rotate_species') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%rotate_species%stat
        Call set_read_status(word, io, model_data%rotate_species%fread, model_data%rotate_species%fail)

      Else If (word(1:length) == 'multiple_input_atoms') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%multiple_input_atoms%value
        Call set_read_status(word, io, model_data%multiple_input_atoms%fread, model_data%multiple_input_atoms%fail)

      Else If (word(1:length) == 'multiple_output_atoms') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%multiple_output_atoms%value
        Call set_read_status(word, io, model_data%multiple_output_atoms%fread, model_data%multiple_output_atoms%fail)

      Else If (word(1:length) == 'normal_vector') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%normal_vector%type
        Call set_read_status(word, io, model_data%normal_vector%fread, model_data%normal_vector%fail,&
                           & model_data%normal_vector%type)

      Else If (word(1:length) == 'both_surfaces') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%both_surfaces%stat
        Call set_read_status(word, io, model_data%both_surfaces%fread, model_data%both_surfaces%fail)

        Else If (word(1:length) == 'add_species_from') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%add_species_from%value, model_data%add_species_from%units 
        Call set_read_status(word, io, model_data%add_species_from%fread, model_data%add_species_from%fail)
        Call capital_to_lower_case(model_data%add_species_from%units)

      Else If (word(1:length) == 'centre_electrode') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%centre_electrode%stat
        Call set_read_status(word, io, model_data%centre_electrode%fread, model_data%centre_electrode%fail)

      Else If (word(1:length) == 'add_extra_space') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%add_extra_space%value, model_data%add_extra_space%units 
        Call set_read_status(word, io, model_data%add_extra_space%fread, model_data%add_extra_space%fail)
        Call capital_to_lower_case(model_data%add_extra_space%units)
        
      Else If (word(1:length) == 'optimise_size') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%optimise_size%stat
        Call set_read_status(word, io, model_data%optimise_size%fread, model_data%optimise_size%fail)
        
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Specifications of settings to build simulation files 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Else If (word(1:length) == '&simulation_settings') Then
        If (.Not.model_data%input_composition%fread) Then
          Write (message,'(2(1x,a))') Trim(set_error), '"&simulation_settings" must be specified after&
                                    & "&input_composition"'
          Call error_stop(message)
        End If
 
        Read (iunit, Fmt=*, iostat=io) word 
        If (simulation_data%generate) Then
          Call duplication_error(word)
        End If
        simulation_data%generate=.True.
        simulation_data%total_tags=model_data%total_tags
      
        ! Asssign variables before reading simulation settings
        Do i=1, simulation_data%total_tags
          simulation_data%component(i)%tag=model_data%component%tag(i)
          simulation_data%component(i)%element=model_data%component%element(i)
          simulation_data%component(i)%atomic_number=model_data%component%atomic_number(i)
        End Do 
        ! Now, it is ready to read information inside &simulation
        Call read_simulation_settings(iunit, simulation_data)
 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Specifications of HPC settings 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Else If (word(1:length) == '&hpc_settings') Then
        Read (iunit, Fmt=*, iostat=io) word 
        If (hpc_data%generate) Then
          Call duplication_error(word)
        End If
        hpc_data%generate=.True.
        Call hpc_data%init_input_variables()

        ! Now, it is ready to read information inside &hpc_settings
        Call read_hpc_settings(iunit, hpc_data)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! Directive not recognised. Inform and kill 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      Else
        If (word(1:1)=='&') Then
          Write (message,'(1x,6a)') Trim(set_error), ' unknown directive found in ', Trim(set_file), &
                                  &' file: ', Trim(word),'. Do you use "&" to define a block? If so,&
                                  & make sure the block is located in the right place, as it might be&
                                  & a sub-block actually. Please also check for the right syntax.'
        Else
          Write (message,'(5(1x,a))') Trim(set_error), 'unknown directive found in', Trim(set_file), &
                                  &'file: ', Trim(word)
        End If 
        Call error_stop(message)
      End If

    End Do

    ! Close file
    Close(files(FILE_SET)%unit_no)

  End Subroutine read_settings
  
  
  Subroutine check_all_settings(files, model_data, simulation_data, hpc_data) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the correctness of all directive defined in SETTINGS files
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),      Intent(InOut) :: files(:)
    Type(model_type),     Intent(InOut) :: model_data 
    Type(simul_type),     Intent(InOut) :: simulation_data
    Type(hpc_type),       Intent(InOut) :: hpc_data  
    
    Character(Len=256) :: messages(7)
    Character(Len=64 ) :: error_set 

    error_set = '***ERROR in file '//Trim(files(FILE_SET)%filename)//' -'   
  
    !!!!!!!!!!!!!!!!!!!!!!! 
    ! Check the type of analysis selected
    If (model_data%analysis%type /= 'build_model'       .And.  &
       model_data%analysis%type /= 'only_simulation_directives' ) Then 
      Write (messages(1),'(2(1x,a))')  Trim(error_set), 'No (or wrong) specification for directive "analysis"' 
      Call info(messages, 1)
      Write (messages(1),'(1x,a)')  'Possible settings for option "analysis" are:'
      Write (messages(2),'(1x,a)')  '- Build_model'
      Write (messages(3),'(1x,a)')  '- Only_simulation_directives'
      Call info(messages,3)
      Call error_stop(' ') 
    Else
      Call info(' ',1)
      Write (messages(1),'(1x,2a)') 'Requested analysis: ', Trim(model_data%analysis%type)
      Call info(messages,1)
              Call info(' ----------------------------------------------------------------------------------------', 1)
      If (model_data%analysis%type == 'build_model' ) Then
        Write (messages(1),'(1x,a)')  'This option generates an atomistic model, whose composition will target the'
        Write (messages(2),'(1x,a)')  'values defined in &species. Input files for simulations will also be generated'
        Write (messages(3),'(1x,a)')  'only if the &simulation_settings is defined.'
        Call info(messages,3)
      Else If (model_data%analysis%type == 'only_simulation_directives' ) Then
        Write (messages(1),'(1x,a)')  'This option generates (or regenerates):'
        Write (messages(2),'(1x,a)')  ' - input files for simulation (if "&simulation_settings" is present)'
        Write (messages(3),'(1x,a)')  ' - script files for HPC submission (if "&hpc_settings" is present)'
        Write (messages(4),'(1x,2a)') 'without the need to generate the atomistic models again. This analysis&
                                    & strictly needs the file ', Trim(files(FILE_RECORD_MODELS)%filename)
        Write (messages(5),'(1x,a)')  'which is created after the generation of atomistic models and contains&
                                    & all the relevant settings.'
        Write (messages(6),'(1x,a)')  'The output with information of the atomistic models has been saved to&
                                    & the OUT_BACKUP file.'
        Write (messages(7),'(1x,3a)') 'Thus, the new ', Trim(files(FILE_OUT)%filename), ' file must be&
                                    & interpreted as a complement of the OUT_BACKUP file.' 
        Call info(messages,7)
      End If
              Call info(' ----------------------------------------------------------------------------------------', 1)
    End If    

    Call check_species(files, model_data)
    Call check_atomic_settings(files, model_data)
    
    If (simulation_data%generate) Then
      simulation_data%code_format=model_data%output_model_format%type
      simulation_data%normal_vector=model_data%normal_vector%type
      simulation_data%solvation%both_surfaces=model_data%both_surfaces%stat 
      Call check_simulation_settings(files, model_data, simulation_data)
    End If
    If (hpc_data%generate) Then
      Call check_hpc_settings(files, model_data%output_model_format%type, hpc_data)
    End If

    If (Trim(model_data%analysis%type) == 'only_simulation_directives') Then
      If ((.Not. simulation_data%generate) .And. (.Not. hpc_data%generate)) Then  
        Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                   & 'The user has requested the generation of files for simulations (dft_settings)&
                                   & and/or the generation of HPC files (HPC_settings).'
        Write (messages(2),'(1x,a)') 'However, none of these blocks is defined. Requested analysis is not possible.&
                                   & Please correct.'
        Call info(messages, 2)
        Call error_stop(' ')
      End If
    End If
    
    ! Refresh out
    Call refresh_out(files) 

  End Subroutine check_all_settings

End module settings

