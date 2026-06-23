!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!                    Welcome to ALC_MOSES
! A code to automatically generate atomistic models of electrochemical 
! interfaces together with input files for simulations using the 
! Grand-Canonical DFT scheme.
!
! Copyright:  2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC, UKRI)  
!               
! Author:     Ivan Scivetti (i.scivetti)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Program alc_moses

  Use atomistic_setup,   Only: model_type

  Use atomistic_builder, Only: build_atomistic_model 

  Use fileset,           Only: file_type, &
                               NUM_FILES, &
                               print_header_out, &
                               set_system_files, &
                               FILE_HPC_SETTINGS, &
                               FILE_OUT, &
                               FILE_RECORD_MODELS,&
                               FOLDER_RESTART,& 
                               wrapping_up

  Use hpc,               Only: hpc_type, &
                               build_hpc_script 
                               
  Use input_types,       Only: in_string                              
                                
  Use numprec,           Only: wi,& 
                               wp
  
  Use settings,          Only: read_settings, &
                               check_all_settings

  Use unit_output,       Only: info

  Use simulation_setup,  Only: simul_type 
  
  Use simulation_files_builder, Only : generate_simulation_directives_only

Implicit None

! all simulation variables
  Type(file_type)      :: files(NUM_FILES)
  Type(model_type)     :: model_data
  Type(simul_type)     :: simulation_data
  Type(hpc_type)       :: hpc_data

  !Time related variables
  Character(Len=256) :: message
  Integer(kind=wi)   :: start,finish,rate

  Call system_clock(count_rate=rate)
  ! Record initial time
  Call system_clock(start)
  ! Initialise settings for input/output files
  Call set_system_files(files)
  ! Print header of OUT
  Call print_header_out(files) 
  ! Read settings from SET
  Call read_settings(files, model_data, simulation_data, hpc_data)
  ! Check the specification of settings
  Call check_all_settings(files, model_data, simulation_data, hpc_data)
  
  ! HPC related directives
  If (hpc_data%generate) Then
    Call build_hpc_script(files, model_data%output_model_format%type, hpc_data)
  End If

  If (Trim(model_data%analysis%type) /= 'only_simulation_directives') Then
       Call build_atomistic_model(files, model_data, simulation_data, hpc_data) 
  Else
    ! Building simulations files and/or hpc settings only
     Call generate_simulation_directives_only(model_data%output_model_format%type, files, &
                                      & model_data, simulation_data, hpc_data)
  End If

  ! Record final time
  Call system_clock(finish)

  Call info(' ', 1)
  Call info(' ==========================================', 1)
  Write (message, '(1x,a,f9.3,a)') 'Total execution time = ',  Real(finish-start,Kind=wp)/rate,  ' seconds.' 
  Call info(message, 1)
  Call info(' ==========================================', 1)

  ! Print appendix to OUT file
  Call wrapping_up(files)

  ! Create restart for subsequent runs, only if atomic models have been generated
  If (Trim(model_data%analysis%type) /= 'only_simulation_directives') Then
    Call execute_command_line('[ ! -d '//Trim(FOLDER_RESTART)//' ] && '//'mkdir '//Trim(FOLDER_RESTART))
    Call execute_command_line('cp '//Trim(files(FILE_OUT)%filename)//' '//Trim(FOLDER_RESTART)//'/OUT_BACKUP') 
    Call execute_command_line('mv '//Trim(files(FILE_RECORD_MODELS)%filename)//' '//Trim(FOLDER_RESTART)) 
  End If
  
  ! Remove temporary HPC_SETTINGS file
  If (hpc_data%generate) Then
    Call execute_command_line('rm '//Trim(files(FILE_HPC_SETTINGS)%filename))
  End If

End Program alc_moses
