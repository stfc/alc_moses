!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Module for input/output files and related subroutines
!
! Copyright: 2026 Ada Lovelace Centre (ALC)
!            Scientific Computing Department (SCD)
!            The Science and Technology Facilities Council (STFC)
!
! Author -   i.scivetti  Jan 2026
!!!!!!!!!!!!!!!!!!!!11!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module fileset

  Use constants,    Only: code_name,&
                          code_VERSION, &
                          date_RELEASE
  Use numprec,      Only: wi
  use unit_output,  Only: info, &
                          set_output_unit 

  Implicit None
  Private

  ! File data
  Type, Public :: file_type
    Private
    ! Filename
    Character(Len=256), Public :: filename
    ! Fortran unit number, set with newunit=T%unit_no
    Integer(Kind=wi), Public   :: unit_no = -2
  Contains
    Procedure, Public :: init => file_type_init
    Procedure, Public :: Close => close_file
  End Type file_type

  ! SETTINGS file
  Integer(Kind=wi), Parameter, Public :: FILE_SET = 1
  ! OUT_EQC file
  Integer(Kind=wi), Parameter, Public :: FILE_OUT = 2 
  ! INPUT STRUCTURE
  Integer(Kind=wi), Parameter, Public :: FILE_INPUT_ELECTRODE = 3
  ! OUTPUT VASP STRUCTURE 
  Integer(Kind=wi), Parameter, Public :: FILE_OUTPUT_STRUCTURE = 4
  ! SET_SIMULATION file
  Integer(Kind=wi), Parameter, Public :: FILE_SET_SIMULATION  = 5 
  ! KPOINTS file
  Integer(Kind=wi), Parameter, Public :: FILE_KPOINTS  = 6 
  ! HPC_SETTINGS file
  Integer(Kind=wi), Parameter, Public :: FILE_HPC_SETTINGS  = 7 
  ! RECORD_MODELS file 
  Integer(Kind=wi), Parameter, Public :: FILE_RECORD_MODELS = 8
  ! MODEL_SUMMARY file 
  Integer(Kind=wi), Parameter, Public :: FILE_MODEL_SUMMARY = 9
  
  ! Size of filename array
  Integer(Kind=wi), Parameter, Public :: NUM_FILES = 9

  ! Folder data
  Character(Len=256), Public :: FOLDER_SIMULATION = "SIMULATION_FILES"
  Character(Len=256), Public :: FOLDER_DFT        = "DFT"
  Character(Len=256), Public :: FOLDER_RESTART    = "RESTART"
  Character(Len=256), Public :: FOLDER_INPUT_GEOM = "INPUT_GEOM"

  Public :: set_system_files, print_header_out, wrapping_up, refresh_out

Contains

  Subroutine refresh_out(files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to refresh the output
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES) 

    Call files(FILE_OUT)%close ()
    Open (Newunit=files(FILE_OUT)%unit_no, File=files(FILE_OUT)%filename, Position='Append')
    Call set_output_unit(files(FILE_OUT)%unit_no)

  End Subroutine refresh_out

  Subroutine file_type_init(T, filename)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to initialise files
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    Class(file_type)                :: T
    Character(Len=*), Intent(In   ) :: filename

    T%filename = Trim(filename)
  End Subroutine file_type_init


  Subroutine set_names_files(files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set default names for files
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES)

    Character(Len=256), Dimension(NUM_FILES)   :: set_names
    Integer(Kind=wi)                           :: file_no

    ! Default file names array
    ! Populate default names array
    set_names(FILE_SET)         = "SETTINGS"
    set_names(FILE_OUT)         = "OUTPUT"
    set_names(FILE_INPUT_ELECTRODE)  = "INPUT_ELECTRODE"
    set_names(FILE_OUTPUT_STRUCTURE) = "OUTPUT_STRUCTURE"
    set_names(FILE_SET_SIMULATION)   = "SET_SIMULATION"
    set_names(FILE_HPC_SETTINGS)     = "HPC_SETTINGS"
    set_names(FILE_RECORD_MODELS)    = "RECORD_MODELS"
    set_names(FILE_MODEL_SUMMARY)    = "MODEL_SUMMARY"
    set_names(FILE_KPOINTS)          = "SET_KPOINTS"

    ! Set default filenames
    Do file_no = 1, NUM_FILES
      Call files(file_no)%init(set_names(file_no))
    End Do
    
  End Subroutine set_names_files


  Subroutine close_file(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to close files
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    Class(file_type) :: T

    Logical :: is_open

    Inquire (T%unit_no, opened=is_open)
    If (is_open) Then
      Close (T%unit_no)
      T%unit_no = -2
    End If

  End Subroutine close_file

  Subroutine set_system_files(files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to open OUTPUT file 
    ! 
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES)

    Call set_names_files(files)   
    Open (Newunit=files(FILE_OUT)%unit_no, File=files(FILE_OUT)%filename, Status='replace')
    Call set_output_unit(files(FILE_OUT)%unit_no)

  End Subroutine set_system_files   

  Subroutine print_header_out(files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print the header to OUTPUT file 
    !  
    ! author        - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES)

    Character(Len=*), Parameter :: fmt1 = '(a)'
    Character(Len=*), Parameter :: fmt2 = '(3a)'
    Character(Len=*), Parameter :: fmt3 = '(4a)'
    Character(Len=128)          :: header(14)

    Write (header(1), fmt1)   Repeat("#", 72)
    Write (header(2), fmt2)  "#                      WELCOME TO ", Trim(code_name),  Repeat(" ", 28)//"#"
    Write (header(3), fmt1)  "#  A software to build atomistic models of electrochemical interfaces  #"
    Write (header(4), fmt1)  "#  together with input files for Grand-Canonical DFT simulations       #"
    Write (header(5), fmt3)  "#  version:  ", Trim(code_VERSION), Repeat(' ',55),                   "#"
    Write (header(6), fmt3)  "#  release:  ", Trim(date_RELEASE), Repeat(' ',48),                   "#"
    Write (header(7), fmt1)  "#                                                                      #"
    Write (header(8), fmt1)  "#  Copyright:  2026  Ada Lovelace Centre (ALC)                         #"
    Write (header(9), fmt1)  "#              Scientific Computing Department (SCD)                   #"
    Write (header(10), fmt1) "#              Science and Technology Facilities Council (STFC)        #"
    Write (header(11), fmt1) "#                                                                      #"
    Write (header(12), fmt1) "#  Author:     Ivan Scivetti (SCD-STFC)                                #"
    Write (header(13), fmt1) "#                                                                      #"
    Write (header(14), fmt1)  Repeat("#", 72)
    Call info(header, 14)
    
    ! Refresh OUT
    Call refresh_out(files)

  End Subroutine print_header_out

  Subroutine wrapping_up(files)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to print final remarks to OUT file 
  ! and close the file 
  !  
  ! author    - i.scivetti Jan 2026
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES)
 
    Character(Len=*), Parameter :: fmt1 = '(1x,a)'
    Character(Len=*), Parameter :: fmt2 = '(1x,3a)'
    Character(Len=128)          :: appex(7)
     
    Write (appex(1), fmt1)   Repeat(" ", 1)
    Write (appex(2), fmt1)   Repeat("#", 35)
    Write (appex(3), fmt1)  "#                                 #" 
    Write (appex(4), fmt1)  "#  Job has finished successfully  #"
    Write (appex(5), fmt2)  "#  Thanks for using ", Trim(code_name), " !!! #"
    Write (appex(6), fmt1)  "#                                 #" 
    Write (appex(7), fmt1)   Repeat("#", 35)
    Call info(appex, 7)

    Close(files(FILE_OUT)%unit_no)    

  End Subroutine wrapping_up  

End Module fileset
