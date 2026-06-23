!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module that generates atomistic models of electrochemical interfaces from the input
! from the input reference structure and the information provided in the SETTINGS
! file. This module also generates:
!
! - backup information of the atomistic models
! - input files for the simulation of the atomistic models
! - scripts for job submission to HPC 
!
! In addition, the generated atomistic model can be used to rebuild 
! input files for simulations and/or HPC scripts if the user decides any 
! parameter must be adjusted (without the need to re-generate atomistic models)    
!
! Copyright:  2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author       - i.scivetti  March 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module atomistic_builder 

  Use atomistic_setup,    Only : model_type,  &
                                 sample_type, &
                                 read_input_model

  Use atomistic_tools,    Only : check_consistency_input_model, &
                                 compute_number_input_species, &
                                 identify_species_input, &
                                 input_model_species_vs_target, &
                                 check_orthorhombic_cell,&
                                 check_add_species_from, &
                                 define_model_cell, &
                                 compute_area_slab, &
                                 match_target_number_species, &
                                 init_random_seed, &
                                 insert_species, &
                                 remove_species, &
                                 define_repeated_model, &
                                 create_list_net_elements, &
                                 set_electrode_boundary, &
                                 centre_electrode, &
                                 check_pcc_vs_generated_model, &
                                 optimise_perp_cell_size, &
                                 surface_shift,&
                                 normal_along_vector

   Use code_castep,       Only : print_castep_settings
   Use code_cp2k,         Only : print_cp2k_settings
   Use code_onetep,       Only : print_onetep_settings
   Use code_vasp,         Only : print_vasp_settings
                                
                              
   Use fileset,           Only : file_type,              &
                                 FILE_RECORD_MODELS,&
                                 FILE_MODEL_SUMMARY, &
                                 FILE_SET_SIMULATION, &
                                 FILE_OUTPUT_STRUCTURE,  &
                                 FILE_KPOINTS, &                                  
                                 FILE_HPC_SETTINGS, &
                                 FOLDER_DFT, &
                                 FOLDER_SIMULATION, &
                                 refresh_out   

   Use hpc,               Only : hpc_type, &
                                 summary_hpc_settings 
                                 
  Use numprec,            Only : wi, &
                                 wp
                          
  Use simulation_setup,   Only : simul_type
                          
  Use unit_output,        Only : error_stop,&
                                 info 

  Use simulation_files_builder, Only : summary_simulation_settings, &
                                       warning_simulation_settings
                                 
  Implicit none                               
                                 
  Public ::  build_atomistic_model

Contains
 
  Subroutine build_atomistic_model(files, model_data, simulation_data, hpc_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to build atomistic models and input files for simulation
    !
    ! author    - i.scivetti Jan 2026
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(InOut) :: model_data
    Type(simul_type),    Intent(InOut) :: simulation_data
    Type(hpc_type),      Intent(In   ) :: hpc_data

    Character(Len=256) :: messages(5)
    Character(Len=256) :: message
    Real(Kind=wp), Dimension(3)  :: v1, v2
    Integer(Kind=wi)   :: i,  record_unit
    Logical :: fortho
    Logical :: activate_random, activate_loop, alloc_sample_arrays
   
    activate_random=.True.
    activate_loop=.True.
    alloc_sample_arrays=.True.

    model_data%input_times=1
    Do i = 1,3 
      model_data%input_times=model_data%input_times*model_data%repeat_input_model%value(i)
    End Do

    Call info(' ', 1)
    Write (messages(1),'(1x,a)') '====================='
    Write (messages(2),'(1x,a)') 'Build atomistic model'
    Write (messages(3),'(1x,a)') '====================='
    Call info(messages, 3)

    ! Refresh out
    Call refresh_out(files)
      
    Call print_atomistic_settings(model_data)
    
    Call read_input_model(files, model_data)

    If (simulation_data%solvation%info%stat) Then
      Call check_orthorhombic_cell(model_data%input%cell, fortho) 
      If (.Not.fortho) Then
        Write (message,'(1x,1a)') '***ERROR: computation with implicit solvent is only possible for&
                                 & orthorhombic cells.'
        Call error_stop(message)
      End If        
    End If

    If (.Not. model_data%add_species_from%fread) Then
      Call set_electrode_boundary(model_data)
    End If
    
    ! Centre structure 
    If (model_data%centre_electrode%stat) Then
      Call centre_electrode(model_data, 'shift_initial')
    End If
    
    ! Check where to start adding species
    Call check_add_species_from(model_data, simulation_data)

    Call check_consistency_input_model(model_data)
    Call model_data%species_arrays()  
    Call compute_number_input_species(model_data)

    model_data%types_species=model_data%num_species%value
     
    Call identify_species_input(model_data)

    If (Trim(model_data%analysis%type)=='build_model') Then
      Call input_model_species_vs_target(files, model_data)
    End If

    If (model_data%both_surfaces%stat .And. model_data%input%change_species_number) Then
      Call surface_shift(model_data)
    End If    
    
    Call info(' ', 1)
    Call info('Atomistic details', 1)
    Call info('=================', 1)
    
    ! Compute the surface area in case of deposited 
    If (Trim(model_data%normal_vector%type)=='c3') Then
      v1(:)=model_data%input%cell(1,:)
      v2(:)=model_data%input%cell(2,:)
    ElseIf (Trim(model_data%normal_vector%type)=='c2') Then
      v1(:)=model_data%input%cell(1,:)
      v2(:)=model_data%input%cell(3,:)
    ElseIf (Trim(model_data%normal_vector%type)=='c1') Then
      v1(:)=model_data%input%cell(2,:)
      v2(:)=model_data%input%cell(3,:)
    End If
    Call compute_area_slab(v1,v2,model_data%input%slab_area)

    ! Opening file for record of relevant modelling settings
    Open(Newunit=files(FILE_RECORD_MODELS)%unit_no, File=files(FILE_RECORD_MODELS)%filename,Status='Replace')
    record_unit=files(FILE_RECORD_MODELS)%unit_no
    Write (record_unit, '(2(3x,a))') Trim(model_data%analysis%type), Trim(model_data%output_model_format%type) 

    model_data%input%num_atoms_extra=model_data%input%num_atoms

    ! Settings for extra atoms
    Do i=1, model_data%types_species  
      model_data%input%species(i)%num_extra=model_data%input%species(i)%num
    End Do
    
    Call match_target_number_species(model_data, 'input')
    ! Set initialization for random number: This is done only once when activate_random=.True.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    If ((model_data%remove_species .Or. model_data%insert_species) .And. activate_random) Then
      Call init_random_seed()         
      activate_random=.False.
    End If
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    ! Remove species in the input model to match target composition
    If (model_data%remove_species) Then
      Call remove_species(model_data%input, model_data%types_species, model_data%arr_added_species%type, 'input')
    End If
    ! Insert species to match target composition. This is done after having removed species to make some space
    If (model_data%insert_species) Then
      Call insert_species(model_data%arr_added_species%type, model_data)
    End If
    
    If (model_data%both_surfaces%stat) Then
      If (model_data%optimise_size%stat) then
        Call optimise_perp_cell_size(model_data, simulation_data)
      End If
    End If
    
    ! Check compliance with planar counter charge electrolyte
    If (simulation_data%electrolyte%info_pcc%stat) Then
      Call check_pcc_vs_generated_model(model_data, simulation_data) 
    End If

    ! Define simulation cell
    Call define_model_cell(model_data, simulation_data)

    ! Calculate the amount of atoms of the modelled sample
    model_data%sample%num_atoms=model_data%input_times*model_data%input%num_atoms_extra

    ! Allocate atomic arrays for the output model, only once
    If (alloc_sample_arrays) Then
      Call model_data%atomic_arrays_model(model_data%sample%num_atoms)  
      alloc_sample_arrays=.False.
    End If
    
    ! define model
    Call define_repeated_model(model_data)
    
    ! Print composition of generated models
    Call generated_model_summary(files, model_data)

    ! Create list
    Call create_list_net_elements(model_data%sample, model_data%types_species,model_data%both_surfaces%stat)

    ! Print output file
    Call print_output_files(files, model_data, simulation_data, hpc_data)
    
    ! Record relevant settings to file RECORD_MODELS
    Call record_models(record_unit, model_data)

   ! Refresh out
    Call refresh_out(files)
    
    ! If the user has requested the generation of simulation files, the following will provide summary of the 
    ! settings as well as aspects to be taken into consideration
    ! Print summary of HPC settings requested by the user
    If (simulation_data%generate) Then
      Call info(' ', 1)
      Call info(' =======================================================', 1)
      Call info(' Summary of the generated input settings for simulations', 1)
      Call info(' =======================================================', 1)
      Call summary_simulation_settings(simulation_data)
    Else
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'INFO: No input files for atomistic level simulations have been generated.'
      Write (messages(2),'(1x,a)') 'To generate input files, the user must define&
                                 & "&simulation_settings" in the SETTINGS file.'
      Call info(messages,2)
    End If
    
    If (hpc_data%generate) Then
      Write (messages(1),'(1x,a)')  'In addition, the user has requested to build HPC script files.'
      Write (messages(2),'(1x,3a)') 'Each sub-folder contains the generated "', Trim(hpc_data%script_name),&
                                   &'" file for job submission.' 
      Call info(messages,2)
      Call summary_hpc_settings(hpc_data)
    Else
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'INFO: No scripts for HPC job submission has been generated.'
      Write (messages(2),'(1x,a)') 'To generate HPC scripts, the user must define "&hpc_settings"&
                                 & in the SETTINGS file.'
      Call info(messages,2)
    End If

    If ((.Not.simulation_data%generate) .Or. (.Not. hpc_data%generate)) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)')  'INFO: Input files for atomistic level simulations and/or HPC scripts can be&
                                  & generated/updated without' 
      Write (messages(2),'(1x,a)')  'the need to re-build atomistic models using the option "only_simulation_directives"&
                                  & for directive "analysis".'
      Call info(messages,2)
    End If

    If (simulation_data%generate) Then
      ! Print warnings
      Call warning_simulation_settings(simulation_data)
    End If
    
    ! Close RECORD_MODELS
    Close(record_unit)

  End Subroutine build_atomistic_model

  
  Subroutine print_atomistic_settings(model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print relevant settings used to building atomistic models
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),    Intent(InOut) :: model_data

    Character(Len=256)  :: message

    If (model_data%analysis%type == 'build_model') Then
      Write (message,'(1x,2a)') 'Selected option to arrange the electrolyte over the electrode: ',&
                               & Trim(model_data%arr_added_species%type)
      Call info(message,1)  
          
      ! Print delta_space 
      If (model_data%delta_space%fread) Then
        Write (message,'(1x,a,f5.3,1x,a)') 'Spatial discretization of the simulation cell: ',&
                                             & model_data%delta_space%value, '[Angstrom]'
      Else
        Write (message,'(1x,a,f5.3,1x,a)') 'Spatial discretization of the simulation cell (default): ', &
                                            & model_data%delta_space%value, '[Angstrom]'
      End If
      Call info(message,1) 
      
      ! Print distance_cutoff
      If (model_data%distance_cutoff%fread) Then
        Write (message,'(1x,a,f4.2,1x,a)') 'Minimum separation distance between species: ',&
                                             & model_data%distance_cutoff%value, '[Angstrom]'
      Else
        Write (message,'(1x,a,f4.2,1x,a)') 'Minimum separation distance between species (default): ', &
                                            & model_data%distance_cutoff%value, '[Angstrom]'
      End If
      Call info(message,1) 
      
      ! Print in case rotations are prevented 
      If (.Not. model_data%rotate_species%stat) Then
        Write (message, '(1x,a)') 'In case there are molecular species to being incorparated to the model, they will&
                                 & not be randomly rotated but keep their orientation, as defined in the corresponding xyz files'
        Call info(message,1) 
      End If
 
    End If
    
    ! Print if the structure is centered 
    If (model_data%centre_electrode%stat) Then
      Write (message, '(1x,a)') 'The whole model is centered within the simulation cell along the direction&
                                & perpendicular to the surface'
      Call info(message,1)
    End If

    ! Print format for input/output
    Write (message, '(2(1x,a))') 'Format for input model: ', Trim(model_data%input_model_format%type)
    Call info(message,1) 
     
    Write (message, '(2(1x,a))') 'Format for output model:', Trim(model_data%output_model_format%type)
    Call info(message,1) 

  End Subroutine print_atomistic_settings   

  Subroutine generated_model_summary(files, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print compositional details of the generated model 
    !
    ! author       - i.scivetti Feb 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),  Intent(InOut) :: files(:)
    Type(model_type), Intent(InOut) :: model_data

    Integer(Kind=wi)   :: i, num, summary
    Character(Len=256) :: message, line
    Character(Len= 64) :: fmt0, fmt1, fmt2

    fmt0='(15x,a))'
    fmt1='(6x,a,2x,a)'
    fmt2='(1x,a12,4x,i7)'
   
    ! Open MODEL_SUMMARY  
    Open(Newunit=files(FILE_MODEL_SUMMARY)%unit_no, File=files(FILE_MODEL_SUMMARY)%filename, Status='Replace')
    summary=files(FILE_MODEL_SUMMARY)%unit_no
 
    Write (message,'(a)') 'Total number of species'
    Call info(message,1);  Write (summary,'(a)') Trim(message)
    Write (line,'(a)') '--------------------------'
    Call info(line,1);     Write (summary,'(a)') Trim(line)
    Write (message, fmt0)            '|  Total  |'
    Call info(message,1);  Write (summary,'(a)') Trim(message)
    Write (message, fmt1) 'Species', '|  number |'
    Call info(message,1);  Write (summary,'(a)') Trim(message)
    Call info(line,1);     Write (summary,'(a)') Trim(line) 

    Do i=1, model_data%types_species 
      num=model_data%sample%species(i)%D_num+ model_data%sample%species(i)%num_show
      If (model_data%sample%species(i)%change_content .And. model_data%both_surfaces%stat) Then
        num=2*num
      End If
      Write (message, fmt2)  Trim(model_data%sample%species(i)%tag), num
      Call info(message,1); Write (summary,'(a)') Trim(message)
    End Do
    Call info(line,1); Write (summary,'(a)') Trim(line)

    If (model_data%both_surfaces%stat) Then
      If (model_data%input%change_species_number) Then
        Write (message, '(a)')    ' ***IMPORTANT: The electrode model has species at both surfaces.'
        Call info(message, 1); Write (summary,'(a)') Trim(message)
      End If  
      
      If (model_data%optimise_size%stat) Then
        If (model_data%input%size_changed) Then
          Write (message, '(a)') ' ***IMPORTANT: The perpendicular size of the cell has been optimised.'
          Call info(message, 1); Write (summary,'(a)') Trim(message)
        End If
      Else
        Write (message, '(a)') ' ***IMPORTANT: The user might need to enlarge the cell size perpendicular to the surface'
        Call info(message, 1); Write (summary,'(a)') Trim(message)      
      End If
      
      If (model_data%selective_dyn .And. (Trim(model_data%output_model_format%type)=='vasp')) Then
        Write (message, '(a)')    ' ***IMPORTANT: The user should review the correctness of fixing atomic&
                                 & coordinates for the simulation of this model.'
        Call info(message, 1); Write (summary,'(a)') Trim(message)
      End If
      Call info(' ', 1)
    End If
   
    Close(summary)

  End Subroutine generated_model_summary
  
  
  Subroutine print_output_files(files, model_data, simulation_data, hpc_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print output files with:
    !  - the generated sample models
    !  - all related files with settings for the simulation
    !
    ! author         - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(model_type),  Intent(InOut) :: model_data
    Type(simul_type),  Intent(InOut) :: simulation_data
    Type(hpc_type),    Intent(In   ) :: hpc_data

    Character(Len=256) :: exec_mkdir, exec_cp, exec_mv, exec_cat
    Character(Len=256) :: fileformat
    Character(Len=256) :: filename, path
    
    Real(Kind=wp)      :: l(3)
    
    Integer(Kind=wi)   :: nfixed , i, j, isp, indx
    
    fileformat=model_data%output_model_format%type
    filename='SAMPLE.'//Trim(fileformat)

    If (Trim(model_data%normal_vector%type)=='c3') Then
      indx=3
    ElseIf (Trim(model_data%normal_vector%type)=='c2') Then
      indx=2
    ElseIf (Trim(model_data%normal_vector%type)=='c1') Then
      indx=1
    End If

    l(:)=model_data%sample%cell(indx,:)
    Call normal_along_vector(l,model_data%sample%normal)    

    ! compute the geometrical centre of the slab, only if the "both_surfaces" directive is activated
    If (model_data%both_surfaces%stat) Then
      model_data%sample%slab_centre=0.0_wp
      nfixed=0
      Do isp = 1, model_data%types_species
        Do j=1, model_data%sample%species(isp)%num_components
          Do i=1, model_data%sample%num_atoms
            If ((.Not. model_data%sample%atom(i)%vanish) .And.&
               (Trim(model_data%sample%atom(i)%tag)==Trim(model_data%sample%species(isp)%component%tag(j))) ) Then
              If (model_data%sample%species(isp)%topology=='electrode') Then
                model_data%sample%slab_centre=model_data%sample%slab_centre+model_data%sample%atom(i)%r
                nfixed=nfixed+1
              End If
            End If
          End Do
        End Do
      End Do
      model_data%sample%both_surfaces=.True.
      model_data%sample%slab_centre=model_data%sample%slab_centre/Real(nfixed,Kind=wp)
    End If
    
    If (Trim(fileformat)=='vasp') Then
      Call print_vasp_output(files, model_data%sample, model_data%types_species, model_data%selective_dyn,  simulation_data)
    Else If (Trim(fileformat)=='cp2k') Then
      Call print_cp2k_output(files, model_data%sample, model_data%types_species, simulation_data)
    Else If (Trim(fileformat)=='castep') Then
      Call print_castep_output(files, model_data%sample, model_data%types_species, simulation_data)
    Else If (Trim(fileformat)=='onetep') Then
      Call print_onetep_output(files, model_data%sample, model_data%types_species, simulation_data)
    Else If (Trim(fileformat)=='xyz') Then 
      Call print_xyz_output(files, model_data%sample, model_data%types_species)
    End If 

    Call execute_command_line('[ ! -d '//Trim(FOLDER_SIMULATION)//' ] && '//'mkdir '//Trim(FOLDER_SIMULATION))

   ! If (Trim(model_data%analysis%type)=='build_model') Then
      Write (path,'(a)') Trim(FOLDER_SIMULATION)
      exec_mkdir='mkdir '//Trim(path)
      Call execute_command_line('[ ! -d '//Trim(path)//' ] && '//Trim(exec_mkdir))
    !End If

    Call info(' Atomistic model is printed to file '//Trim(path)//'/'//Trim(filename), 1)

    ! Moving simulation setting files to the corresponding directory
    If (simulation_data%generate) Then
      Call info(' Files for the simulation of the generated model are printed to folder '//Trim(path), 1)
      If (Trim(fileformat)=='vasp') Then
        ! Generate POSCAR
        exec_cp='cp '//Trim(files(FILE_OUTPUT_STRUCTURE)%filename)//' '//Trim(path)//'/POSCAR'
        Call execute_command_line(exec_cp)
      Else If (model_data%output_model_format%type=='castep') Then
        exec_cat='cat '//Trim(files(FILE_SET_SIMULATION)%filename)//' '//&
              &Trim(files(FILE_OUTPUT_STRUCTURE)%filename)//' > model.cell'
        Call execute_command_line(exec_cat)
        exec_mv ='mv  '//Trim(files(FILE_SET_SIMULATION)%filename)//' CI-castep.cell' 
        Call execute_command_line(exec_mv)
      Else If (model_data%output_model_format%type=='onetep') Then
        exec_cat='cat '//Trim(files(FILE_SET_SIMULATION)%filename)//' '//&
              &Trim(files(FILE_OUTPUT_STRUCTURE)%filename)//' > model.dat'
        Call execute_command_line(exec_cat)
        exec_mv ='mv  '//Trim(files(FILE_SET_SIMULATION)%filename)//' CI-onetep.dat'  
        Call execute_command_line(exec_mv)
      End If

      Call print_simulation_files(path, files, model_data, simulation_data, hpc_data)
    End If

    ! Moving structure
    exec_mv='mv '//Trim(files(FILE_OUTPUT_STRUCTURE)%filename)//' '//Trim(path)//'/'//Trim(filename)
    Call execute_command_line(exec_mv)
      
    ! Move MODEL_SUMMARY files
    exec_mv='mv '//Trim(files(FILE_MODEL_SUMMARY)%filename)//' '//Trim(path)//'/'//Trim(files(FILE_MODEL_SUMMARY)%filename)
    Call execute_command_line(exec_mv)

    model_data%sample%path=path

  End Subroutine print_output_files 
 
  Subroutine print_simulation_files(path, files, model_data, simulation_data, hpc_data) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print output files with settings for the simulation of the models
    ! This abstraction is convenient to generated files with simulation settings
    ! without the need of re-generating atomistic models (the most expensive bit)
    !
    ! author         - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=256), Intent(In   ) :: path
    Type(file_type),    Intent(InOut) :: files(:)
    Type(model_type),   Intent(InOut) :: model_data
    Type(simul_type),   Intent(In   ) :: simulation_data
    Type(hpc_type),     Intent(In   ) :: hpc_data

    Character(Len=256) :: exec_cp, exec_mv, fileformat
    Character(Len=256) :: path_pp

    fileformat=model_data%output_model_format%type

    !!!! VASP
    If (Trim(fileformat)=='vasp') Then
      ! Generate INCAR
      exec_mv='mv '//Trim(files(FILE_SET_SIMULATION)%filename)//' '//Trim(path)//'/INCAR'   
      Call execute_command_line(exec_mv)
      ! Generate KPOINTS 
      exec_mv='mv '//Trim(files(FILE_KPOINTS)%filename)//' '//Trim(path)//'/KPOINTS'   
      Call execute_command_line(exec_mv)
      ! Generate POTCAR
      If (simulation_data%dft%pp_info%stat)Then
        exec_mv='mv POTCAR'//' '//Trim(path)//'/POTCAR'   
        Call execute_command_line(exec_mv)
      End If
    !!!! CP2K
    Else If (Trim(fileformat)=='cp2k') Then
      ! Generated input.cp2k 
      exec_mv='mv '//Trim(files(FILE_SET_SIMULATION)%filename)//' '//Trim(path)//'/input.cp2k'
      Call execute_command_line(exec_mv)
      ! copy the potential
      If (simulation_data%dft%pp_info%stat)Then
        path_pp=Trim(FOLDER_DFT)//'/PPs/'//Trim(simulation_data%dft%pseudo_pot(1)%file_name)  
        exec_cp='cp '//Trim(path_pp)//' '//Trim(path)
        Call execute_command_line(exec_cp)
      End If
      ! copy basis set 
      If (simulation_data%dft%basis_info%stat)Then
        path_pp=Trim(FOLDER_DFT)//'/BASIS_SET'
        exec_cp='cp '//Trim(path_pp)//' '//Trim(path)
        Call execute_command_line(exec_cp)
      End If
    !!!! CASTEP
    Else If (Trim(fileformat)=='castep') Then
      exec_mv='mv '//'model.param'//' '//Trim(path) 
      Call execute_command_line(exec_mv)
      exec_mv='mv '//'model.cell'//' '//Trim(path) 
      Call execute_command_line(exec_mv)
      exec_mv='mv '//'CI-castep.cell'//' '//Trim(path) 
      Call execute_command_line(exec_mv)
      If (simulation_data%dft%pp_info%stat)Then
        exec_cp='cp DFT/PPs/* '//Trim(path)
        Call execute_command_line(exec_cp)
      End If
    !!!! ONETEP
    Else If (Trim(fileformat)=='onetep') Then
      exec_mv='mv '//'model.dat'//' '//Trim(path) 
      Call execute_command_line(exec_mv)
      exec_mv='mv '//'CI-onetep.dat'//' '//Trim(path) 
      Call execute_command_line(exec_mv)
      If (simulation_data%dft%pp_info%stat) Then
        exec_cp='cp DFT/PPs/* '//Trim(path)
        Call execute_command_line(exec_cp)
      End If
    End If

    ! If vdW correction uses kernel, copy the file
    If (simulation_data%dft%need_vdw_kernel) Then
      exec_cp='cp '// Trim(FOLDER_DFT)//'/'//Trim(simulation_data%dft%vdw_kernel_file)//' '//Trim(path)//'/'   
      Call execute_command_line(exec_cp)
    End If

    ! If hpc settings are defined, copy to the corresponding HPC file
    If (hpc_data%generate) Then   
      exec_cp='cp '//Trim(files(FILE_HPC_SETTINGS)%filename)//' '//Trim(path)//'/'//Trim(hpc_data%script_name)
      Call execute_command_line(exec_cp)
    End If

  End Subroutine print_simulation_files 

  Subroutine place_at_opposite_surface(T, i)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Place an atom at the opposite electrode surface 
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    Type(sample_type), Intent(InOut) :: T
    Integer(Kind=wi),  Intent(In   ) :: i

    Real(Kind=wp)    :: l(3), dist

    l=T%atom(i)%r-T%slab_centre
    dist=Dot_product(l,T%normal)
    T%inverted=T%atom(i)%r-2.0_wp*dist*T%normal+T%surface_shift
    
  End Subroutine place_at_opposite_surface
  
 
  Subroutine print_vasp_output(files, T, types_species, dynamics, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print files in vasp format 
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(sample_type), Intent(InOut) :: T
    Integer(Kind=wi),  Intent(In   ) :: types_species     
    Logical,           Intent(In   ) :: dynamics
    Type(simul_type),  Intent(InOut) :: simulation_data

    Integer(Kind=wi) :: isp, i, j
    Integer(Kind=wi) :: iunit

    ! Open OUTPUT_STRUCTURE file
    Open(Newunit=files(FILE_OUTPUT_STRUCTURE)%unit_no, File=files(FILE_OUTPUT_STRUCTURE)%filename,Status='Replace')
    iunit=files(FILE_OUTPUT_STRUCTURE)%unit_no

    Write (iunit,'(a)') 'Model for sample' 
    Write (iunit,'(f19.16)') T%scale_factor_vasp
    Do i= 1, 3
      Write (iunit,'(3f20.12)') (T%cell(i,j), j=1,3)
    End Do
    Write (iunit,'(*(6x,a2))') (Trim(T%list%element(j)), j=1, T%list%net_elements)   
    Write (iunit,'(*(i8))') (T%list%N0(j), j=1, T%list%net_elements)
    If (dynamics) Then
      Write (iunit, '(a)')  'Selective dynamics'
    End If   
    Write (iunit, '(a)')    'Cartesian'
    Do isp = 1, types_species
      Do j=1, T%species(isp)%num_components
        Do i=1, T%num_atoms
          If ((.Not. T%atom(i)%vanish) .And. (Trim(T%atom(i)%tag)==Trim(T%species(isp)%component%tag(j)))) Then
            If (dynamics) Then
              Write (iunit, '(3f20.12,3l3)') T%atom(i)%r, T%atom(i)%dynamics
              If (T%species(isp)%change_content .And. T%both_surfaces) Then
                Call place_at_opposite_surface(T, i)
                Write (iunit, '(3f20.12,3l3)') T%inverted, T%atom(i)%dynamics
              End If
            Else
              Write (iunit, '(3f20.12,3l3)') T%atom(i)%r
              If (T%species(isp)%change_content .And. T%both_surfaces) Then
                Call place_at_opposite_surface(T, i)
                Write (iunit, '(3f20.12)') T%inverted
              End If
            End If
          End If
        End Do
      End Do
    End Do

    ! Close file
    Close(iunit)

    If (simulation_data%generate) Then
      Call print_vasp_settings(files, T%list%net_elements, T%list%element, T%list%tag, T%list%N0, simulation_data) 
    End If

  End Subroutine print_vasp_output

  Subroutine print_cp2k_output(files, T, types_species, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print files in CP2K format 
    !
    ! author        - i.scivetti March  2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(sample_type), Intent(InOut) :: T
    Integer(Kind=wi),  Intent(In   ) :: types_species     
    Type(simul_type),  Intent(In   ) :: simulation_data

    Integer(Kind=wi) :: isp, i, j, ntot
    Integer(Kind=wi) :: iunit

    ntot=0
    Do isp = 1, types_species
      Do j=1, T%species(isp)%num_components
        Do i=1, T%num_atoms
          If ((.Not. T%atom(i)%vanish) .And. &
            & Trim(T%atom(i)%tag)==Trim(T%species(isp)%component%tag(j))) Then
            ntot=ntot+1
            If(T%species(isp)%change_content .And. T%both_surfaces) Then
              ntot=ntot+1
            End If
          End If
        End Do
      End Do
    End Do

    ! Open OUTPUT_STRUCTURE file
    Open(Newunit=files(FILE_OUTPUT_STRUCTURE)%unit_no, File=files(FILE_OUTPUT_STRUCTURE)%filename,Status='Replace')
    iunit=files(FILE_OUTPUT_STRUCTURE)%unit_no

    ! print atomic position
    Write (iunit,'(i10)') ntot
    Write (iunit,*) '  ' 

    Do isp = 1, types_species
      Do j=1, T%species(isp)%num_components
        Do i=1, T%num_atoms
          If ((.Not. T%atom(i)%vanish) .And. &
            & Trim(T%atom(i)%tag)==Trim(T%species(isp)%component%tag(j))) Then
            Write (iunit,'(a,2x,3f20.12)') Trim(T%atom(i)%tag), T%atom(i)%r
              If (T%species(isp)%change_content .And. T%both_surfaces) Then
                Call place_at_opposite_surface(T, i)
                Write (iunit,'(a,2x,3f20.12)') Trim(T%atom(i)%tag), T%inverted
              End If
          End If
        End Do
      End Do
    End Do

    Close(iunit)

    If (simulation_data%generate) Then
      Call print_cp2k_settings(files, T%list%net_elements, T%list%element, T%list%tag, T%list%N0, simulation_data)
    End If

  End Subroutine print_cp2k_output 

  Subroutine print_castep_output(files, T, types_species, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print files in CASTEP format 
    !
    ! author        - i.scivetti June 2026 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(sample_type), Intent(InOut) :: T
    Integer(Kind=wi),  Intent(In   ) :: types_species     
    Type(simul_type),  Intent(InOut) :: simulation_data

    Integer(Kind=wi) :: isp, i, j, k
    Integer(Kind=wi) :: iunit
    Logical :: loop 
    Character(Len=8) :: atom_tag  

    ! Open OUTPUT_STRUCTURE file
    Open(Newunit=files(FILE_OUTPUT_STRUCTURE)%unit_no, File=files(FILE_OUTPUT_STRUCTURE)%filename,Status='Replace')
    iunit=files(FILE_OUTPUT_STRUCTURE)%unit_no

    Write (iunit,'(a)') '%BLOCK LATTICE_CART'
    Do i= 1, 3
      Write (iunit,'(3f20.12)') (T%cell(i,j), j=1,3)
    End Do  
    Write (iunit,'(a)') '%ENDBLOCK LATTICE_CART' 
    Write (iunit, *) '  '

    ! print atomic position
    Write (iunit,'(a)') '%BLOCK POSITIONS_ABS'

    Do isp = 1, types_species
      Do j=1, T%species(isp)%num_components
        Do i=1, T%num_atoms
          If ((.Not. T%atom(i)%vanish) .And. Trim(T%atom(i)%tag)==Trim(T%species(isp)%component%tag(j))) Then
            If (simulation_data%dft%mag_info%fread) Then
                k=1
                loop=.True.
                Do While (k <= simulation_data%total_tags .And. loop)
                  If (Trim(T%atom(i)%tag)==Trim(simulation_data%dft%magnetization(k)%tag)) Then
                    If (simulation_data%dft%pp_info%fread) Then
                      atom_tag=Trim(T%atom(i)%tag)
                    Else
                      atom_tag=Trim(T%atom(i)%element)
                    End If
                    Write (iunit,'(a,2x,3f20.12,4x,a,f5.2)') Trim(atom_tag), T%atom(i)%r, 'spin=', &
                                                        & simulation_data%dft%magnetization(k)%value 
                    loop=.False.
                    If (T%species(isp)%change_content .And. T%both_surfaces) Then
                      Call place_at_opposite_surface(T, i)
                      Write (iunit,'(a,2x,3f20.12,4x,a,f5.2)') Trim(atom_tag), T%inverted, 'spin=', &
                                                        & simulation_data%dft%magnetization(k)%value 
                    End If
                  End If
                  k=k+1
                End Do
            Else
              If (simulation_data%dft%pp_info%fread) Then
                atom_tag=Trim(T%atom(i)%tag)
              Else
                atom_tag=Trim(T%atom(i)%element)
              End If
              Write (iunit,'(a,2x,3f20.12)') Trim(atom_tag), T%atom(i)%r
              If (T%species(isp)%change_content .And. T%both_surfaces) Then
                Call place_at_opposite_surface(T, i)
                Write (iunit,'(a,2x,3f20.12)') Trim(atom_tag), T%inverted
              End If
            End If
          End If
        End Do
      End Do
    End Do
    Write (iunit,'(a)') '%ENDBLOCK POSITIONS_ABS'


    Close(iunit)

    If (simulation_data%generate) Then
      Call print_castep_settings(files, T%list%net_elements, T%list%element, T%list%tag, T%list%N0, simulation_data)
    End If

  End Subroutine print_castep_output 

  Subroutine print_onetep_output(files, T, types_species, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print files in CASTEP format
    !
    ! author        - i.scivetti Jan  2026 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(sample_type), Intent(InOut) :: T
    Integer(Kind=wi),  Intent(In   ) :: types_species     
    Type(simul_type),  Intent(InOut) :: simulation_data

    Integer(Kind=wi) :: isp, i, j
    Integer(Kind=wi) :: iunit

    ! Open OUTPUT_STRUCTURE file
    Open(Newunit=files(FILE_OUTPUT_STRUCTURE)%unit_no, File=files(FILE_OUTPUT_STRUCTURE)%filename,Status='Replace')
    iunit=files(FILE_OUTPUT_STRUCTURE)%unit_no

    ! Print simulation cell
    Write (iunit,'(a)') '%block lattice_cart'
    Write (iunit,'(a)') 'ang'
    Do i= 1, 3
      Write (iunit,'(3f12.6)') (T%cell(i,j), j=1,3)
    End Do  
    Write (iunit,'(a)') '%endblock lattice_cart' 
    Write (iunit, *) '  '

    ! print atomic position
    Write (iunit,'(a)') '%block positions_abs'
    Write (iunit,'(a)') 'ang'
    Do isp = 1, types_species
      Do j=1, T%species(isp)%num_components
        Do i=1, T%num_atoms
          If ((.Not. T%atom(i)%vanish) .And. &
            & Trim(T%atom(i)%tag)==Trim(T%species(isp)%component%tag(j))) Then
              Write (iunit,'(a,2x,3f12.6)') Trim(T%atom(i)%tag), T%atom(i)%r
              If (T%species(isp)%change_content .And. T%both_surfaces) Then
                Call place_at_opposite_surface(T, i)
                Write (iunit,'(a,2x,3f12.6)') Trim(T%atom(i)%tag), T%inverted
              End If
          End If
        End Do
      End Do
    End Do
    Write (iunit,'(a)') '%endblock positions_abs' 

    Close(iunit)

    If (simulation_data%generate) Then
      Call print_onetep_settings(files, T%list%net_elements, T%list%tag, T%list%N0, simulation_data)
    End If

  End Subroutine print_onetep_output 
 
  Subroutine print_xyz_output(files, T, types_species)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Prints model in xyz format. 
    !
    ! author        - i.scivetti Jan  2026 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(sample_type), Intent(InOut) :: T
    Integer(Kind=wi),  Intent(In   ) :: types_species     

    Integer(Kind=wi) :: isp, i, j, ntot
    Integer(Kind=wi) :: iunit

    ntot=0
    Do isp = 1, types_species
      Do j=1, T%species(isp)%num_components
        Do i=1, T%num_atoms
          If ((.Not. T%atom(i)%vanish) .And. &
            & Trim(T%atom(i)%tag)==Trim(T%species(isp)%component%tag(j))) Then
            ntot=ntot+1
            If(T%species(isp)%change_content .And. T%both_surfaces) Then
              ntot=ntot+1
            End If
          End If
        End Do
      End Do
    End Do

    ! Open OUTPUT_STRUCTURE file
    Open(Newunit=files(FILE_OUTPUT_STRUCTURE)%unit_no, File=files(FILE_OUTPUT_STRUCTURE)%filename,Status='Replace')
    iunit=files(FILE_OUTPUT_STRUCTURE)%unit_no

    ! print atomic position
    Write (iunit,'(i10)') ntot
    Write (iunit,'(a, 9f10.4, a)') 'Lattice = "' , ((T%cell(i,j), j=1,3), i=1,3), '"' 
    Do isp = 1, types_species
      Do j=1, T%species(isp)%num_components
        Do i=1, T%num_atoms
          If ((.Not. T%atom(i)%vanish) .And. &
            & Trim(T%atom(i)%tag)==Trim(T%species(isp)%component%tag(j))) Then
              Write (iunit,'(a,2x,3f12.6)') Trim(T%atom(i)%element), T%atom(i)%r
              If (T%species(isp)%change_content .And. T%both_surfaces) Then
                Call place_at_opposite_surface(T, i)
                Write (iunit,'(a,2x,3f12.6)') Trim(T%atom(i)%element), T%inverted
              End If
          End If
        End Do
      End Do
    End Do

    Close(iunit)

  End Subroutine print_xyz_output
 
  Subroutine record_models(record_unit, model_data) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    ! Subroutine to keep relevant settings from the generation                                  
    ! of atomistic models 
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(InOut) :: record_unit
    Type(model_type),  Intent(InOut) :: model_data

    Integer(Kind=wi) :: i

    Write (record_unit,*) Trim(model_data%sample%path)
    Write (record_unit,*) model_data%sample%list%net_elements
    Write (record_unit,'(*(6x,a4))') (Trim(model_data%sample%list%tag(i)), i=1, model_data%sample%list%net_elements)
    Write (record_unit,'(*(6x,a4))') (model_data%sample%list%element(i), i=1, model_data%sample%list%net_elements)
    Write (record_unit,'(*(2x,i8))')    (model_data%sample%list%N0(i),  i=1, model_data%sample%list%net_elements)

  End Subroutine record_models
  
End Module atomistic_builder  
