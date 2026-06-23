!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module with tools to build atomistic models
!
! Copyright: 2026 Ada Lovelace Centre (ALC)
!            Scientific Computing Department (SCD)
!            The Science and Technology Facilities Council (STFC)
!
! Author     - i.scivetti  March 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module atomistic_tools 

  Use atomistic_setup,    Only : model_type, &
                                 about_cell, &
                                 sample_type,&
                                 species_model
                                 

  Use constants,          Only : length_tol, &
                                 min_intra_bond,&
                                 max_nn, &
                                 large_cell_limit,&
                                 max_num_species_units, &
                                 max_at_species, &
                                 times_number_attempts, &
                                 twopi
                                 
   Use fileset,            Only : file_type,              &
                                  FILE_INPUT_ELECTRODE, & 
                                  FOLDER_INPUT_GEOM
                                  
  Use input_types,         Only : in_param
  Use numprec,             Only : wi, &
                                  wp, &
                                  li
                                 
  Use simulation_setup,    Only : simul_type
                                  
  Use unit_output,         Only : error_stop,&
                                  info 

  implicit none                               
                                 
  Public :: set_electrode_boundary, centre_electrode
  Public :: check_consistency_input_model, normal_along_vector
  Public :: compute_number_input_species, identify_species_input, input_model_species_vs_target
  Public :: check_orthorhombic_cell, check_add_species_from, define_model_cell
  Public :: compute_area_slab, match_target_number_species, init_random_seed
  Public :: remove_species, insert_species, define_repeated_model, surface_shift
  Public :: create_list_net_elements, check_pcc_vs_generated_model, optimise_perp_cell_size

Contains

  Subroutine init_random_seed()
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to generate the seed for random number generation 
    !
    ! author    - i.scivetti Feb 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi) :: i, n, clock
    Integer(Kind=wi), DIMENSION(:), Allocatable :: seed
                                      
    Call random_seed(size = n)
    
    Allocate(seed(n))
    Call system_clock(Count=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    Call random_seed(put = seed)
    Deallocate(seed)

  End Subroutine init_random_seed

  Subroutine match_target_number_species(T, option)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Find the difference in the number of species between the target values and
    ! those for the input model
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),  Intent(InOut) ::  T
    Character(Len=*),  Intent(In   ) ::  option
  
    Integer(Kind=wi) :: i

    T%remove_species=.False.
    T%insert_species=.False.

    ! Find the number of species to add/remove from the input model to match the units of the target model 
    Do i = 1, T%types_species
      If (Trim(option)=='input') Then
        T%input%species(i)%D_num=T%species_info(i)%N0_target-T%input%species(i)%num
        If (T%input%species(i)%D_num<0) Then
           T%remove_species=.True.
        End If
        If (T%input%species(i)%D_num>0) Then
           T%insert_species=.True.
        End If
      End If
    End Do   

  End Subroutine match_target_number_species

  Subroutine input_model_species_vs_target(files, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print the comparison between the target composition of the 
    ! model to be generated and the composition of the input model (table format) 
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(InOut) :: model_data

    Integer(Kind=wi)   :: i, delta
    Character(Len=256) :: message
    Character(Len= 64) :: fmt
    
    fmt='(1x,a12,x,a1,3(x,i7))'
    
    model_data%input%change_species_number=.False.

    Call info(' ', 1)
    Write (message,'(1x,a)') 'Number of species identified from the input electrode model (file "'&
                            &//Trim(FOLDER_INPUT_GEOM)//'/'//Trim(files(FILE_INPUT_ELECTRODE)%filename)//&
                            &'") against the target values defined in &species'
    Call info(message, 1)
  
    Call info(' --------------------------------------', 1)
    Call info('      Species |   Input  Target   Delta', 1)
    Call info(' --------------------------------------', 1)
    Do i=1, model_data%num_species%value
      delta=model_data%species_info(i)%N0_target-model_data%input%species(i)%num
      Write (message, fmt)  Trim(model_data%species_info(i)%tag), '|', model_data%input%species(i)%num, &
                            & model_data%species_info(i)%N0_target, delta
      Call info(message,1)
      If (delta /= 0) Then
        model_data%input%change_species_number=.True.
      End If
    End Do
    Call info(' --------------------------------------', 1)

    
    Do i = 1, model_data%num_species%value
      If (Trim(model_data%input%species(i)%topology)=='electrode') Then
        If (model_data%input%species(i)%num /= model_data%species_info(i)%N0_target) Then
          Write (message,'(3a)') '***ERROR: the number of fixed, electrode species "',&
                                 & Trim(model_data%input%species(i)%tag),&   
                                 & '" in the input model is different than the target value&
                                 & assigned in &species. Review the settings and the model.'
          Call info(message,1)
          Call error_stop(' ')    
        End If
      End If
    End Do  
    
    If (.Not. model_data%input%change_species_number) Then
      Write (message, '(a)')    ' ***IMPORTANT: No need to insert/remove species to/from the input model'
      Call info(message, 1)
    End If  
    
  End Subroutine input_model_species_vs_target

  Subroutine compute_number_input_species(model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to identify the number of the species (defined in &species)
    ! within the input model 
    !
    ! author    - i.scivetti Feb 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),    Intent(InOut) :: model_data

    Integer(Kind=wi) :: i, j, k
    Character(Len=256) :: messages(2)
    Logical            :: loop
    Real(Kind=wp)      :: nsp, nsp_round
    Integer(Kind=wi)   :: nsp_int, imin

    ! Determine the number of components for each species of the input model
    Do i=1, model_data%total_tags
      loop=.True.
      j=1
      Do While (j <= model_data%num_species%value .And. loop)
        k=1
        Do While (k <= model_data%species_info(j)%num_components .And. loop)
          If (Trim(model_data%species_info(j)%component%tag(k))==model_data%component%tag(i)) Then
            model_data%input%species(j)%component%N(k)=model_data%component%N0(i)
            loop=.False.
          End If
          k=k+1
        End Do
        j=j+1
      End Do
    End Do
    
    ! For each species:
    ! 1) check that the number of atoms provided makes an interget number of units
    ! 2) find the number of units. 
    Do j=1, model_data%num_species%value
      model_data%input%species(j)%num=0
      Do k= 1, model_data%species_info(j)%num_components
        If ((Trim(model_data%input%species(j)%topology) == 'electrode') .And. &
            model_data%input%species(j)%component%N(k) == 0) Then
          Write (messages(1),'(5a)') '***ERROR: Component "',  Trim(model_data%input%species(j)%component%tag(k)),&
                                   & '" of "fixed" species "', Trim(model_data%input%species(j)%tag), &
                                   & '" cannot have zero atoms in the input structure. Please review settings.'          
          Call info(messages,1)
          Call error_stop(' ')    
        End If
        model_data%input%species(j)%num=model_data%input%species(j)%num+model_data%input%species(j)%component%N(k)
      End Do
      nsp=Real(model_data%input%species(j)%num,Kind=wp)/model_data%species_info(j)%atoms_per_species
      nsp_int=Int(nsp)
      nsp_round=Real(nsp_int,Kind=wp)
      If (Abs(nsp-nsp_round) > epsilon(nsp)) Then
        Write (messages(1),'(3a,f12.5)') '***ERROR: Amount of atoms provided for species "',&
                                 & Trim(model_data%input%species(j)%tag),&
                                 & '" does not make an integer number of units for the input model.&
                                 & Computed number of species is ', nsp    
        Write (messages(2),'(a)')  'Please check consistent between the definition of &species_component and&
                                  & &input_composition, as well as the atoms provided in input structure'  
        Call info(messages,2)
        Call error_stop(' ')  
      End If
        model_data%input%species(j)%num=model_data%input%species(j)%num/model_data%species_info(j)%atoms_per_species
    End Do

    ! For all the species of topology "electrode" find that particular species with minimum number of units 
    model_data%input%min_species=Huge(1)
    Do i = 1, model_data%num_species%value
      If (Trim(model_data%input%species(i)%topology)=='electrode') Then
        If (model_data%input%species(i)%num==0) Then
          Write (messages(1),'(3a)') '***ERROR: the input number for "fixed" species "', Trim(model_data%input%species(i)%tag),&   
                                 & '" is zero. This is a wrong setting for a "fixed" species. Please change.'
          Call info(messages,1)
          Call error_stop(' ')    
        End If
        If (model_data%input%species(i)%num < model_data%input%min_species) Then
          imin=i
          model_data%input%min_species=model_data%input%species(i)%num
        End If
      End If
    End Do     

  End Subroutine compute_number_input_species


  Subroutine identify_species_input(model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to group the species by creating a list with the involved atoms
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),    Intent(InOut) :: model_data

    Integer(Kind=wi) :: i, j, k, m
    Logical               :: loop

    Do i = 1, model_data%types_species
      If (model_data%input%species(i)%change_content) Then
        Do j =1, model_data%input%species(i)%num
          loop=.True.
          k=1
          Do While (k <= model_data%input%num_atoms .And. loop)
            m=1
            Do While (m <= model_data%input%species(i)%num_components .And. loop)
              If (Trim(model_data%input%atom(k)%tag)==Trim(model_data%input%species(i)%component%tag(m)) .And. &
                   (.Not. model_data%input%atom(k)%in_species)) Then
                loop=.False.
                model_data%input%species(i)%units(j)%list(1)=k
                If (model_data%input%species(i)%atoms_per_species > 1) Then
                  Call find_neighbours(model_data%species_info(i)%bond_cutoff, model_data, k, j, i)
                Else
                  model_data%input%atom(k)%in_species=.True. 
                End If
              End If
              m=m+1
            End Do
            k=k+1
          End Do
        End Do
      End If
    End Do

  End Subroutine identify_species_input 


  Subroutine find_neighbours(bond_cutoff, model_data, kin, jin, iin)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to find the all the neighbours-bonded atoms that form a species 
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp),     Intent(In   ) :: bond_cutoff
    Type(model_type),  Intent(InOut) :: model_data
    Integer(Kind=wi),  Intent(In   ) :: kin, jin, iin

    Integer(Kind=wi) :: m, k, ka, il
    Integer(Kind=wi) :: knn, nacc, num_nn, num_tot, nlist, knni
    Logical             :: loop, loop_sp, is_nn, is_list
    Real(Kind=wp)       :: a(3), b(3)
    Real(Kind=wp)       :: dist
    Character(Len=256)  :: messages(6)

    Integer(Kind=wi), Dimension(max_nn) :: list_nn, list_nn_next

    num_nn=1
    list_nn(1)=kin

    loop=.True. 
    num_tot=0
    list_nn_next=0 
    nlist=1

    Do While (loop)
      nacc=0
      Do knn=1, num_nn
        ka=list_nn(knn)
        Do k=1, model_data%input%num_atoms
          loop_sp=.True.
          m=1
          Do While (m <= model_data%input%species(iin)%num_components .And. loop_sp)
            If (Trim(model_data%input%atom(k)%tag)==Trim(model_data%input%species(iin)%component%tag(m)) .And. k/=ka .And.&
               (.Not. model_data%input%atom(k)%in_species) .And. (.Not. model_data%input%atom(ka)%in_species) ) Then
               ! Compute distances
               a=model_data%input%atom(k)%r(:)
               b=model_data%input%atom(ka)%r(:)
               Call compute_distance_PBC(a, b, model_data%input%cell, model_data%input%invcell, dist)
               ! Check if distances are within a given cutoff
               If (dist<bond_cutoff) Then 
                 If (dist<min_intra_bond) Then
                   Write (messages(1),'(a)') '***ERROR in the input structure: '
                   Write (messages(2),'(3(a,i4),3a)') 'Intermolecular distance between atoms ', k, ' and ', ka, &
                                                    ' in unit ', jin, ' of species "',&
                                                    & Trim(model_data%input%species(iin)%tag),'"'
                   Write (messages(3),'(a,f4.2,a)') 'is shorter than the input minimum distance criteria of ', &
                                                   &  min_intra_bond,&
                                                   & ' Angstrom for bonding. Please review the input geometry.'
                   Call info(messages,3)
                   Call error_stop(' ')
                 Else
                   loop_sp=.False.
                   is_nn=.True.
                   il=1
                   Do While (il <= max_nn .And. is_nn) 
                     If (k==list_nn_next(il)) Then
                       is_nn=.False.
                     End If
                     il=il+1
                   End Do 
                   If (is_nn) Then    
                     nacc=nacc+1
                     If (nacc==max_nn+1) Then
                        Write (messages(1),'(a,i4,3a)') '***ERROR: Trouble to identify unit ', jin, ' of species "',&
                                                 & Trim(model_data%input%species(iin)%tag),&
                                                 & '" from the input coordinates of the input structure.'    
                        Write (messages(2),'(a)') ' It is likely the species represents a fixed crystal structure with bonds&
                                                 & periodically repeated. If this is the case, the user must define the&
                                                 & species as "fixed" in &species'
                        Call info(messages, 2)
                        Call error_stop(' ')
             
                     End If
                     list_nn_next(nacc)=k
                     il=1
                     is_list=.True.
                     Do While (il <= model_data%input%species(iin)%atoms_per_species .And. is_list) 
                       If (k==model_data%input%species(iin)%units(jin)%list(il)) Then
                         is_list=.False.
                       End If 
                       il=il+1
                     End Do
                     If (is_list) Then
                       nlist=nlist+1
                       model_data%input%species(iin)%units(jin)%list(nlist)=k
                     End If
                   End If
                 End If  
               End If
            End If
            m=m+1
          End Do
        End Do
        If (.Not. model_data%input%atom(ka)%in_species) Then
          model_data%input%atom(ka)%in_species=.True.
          num_tot=num_tot+1
          If (num_tot+nacc > model_data%input%species(iin)%atoms_per_species) Then
            Do knni= 1, nacc
              model_data%input%atom(list_nn_next(knni))%in_species=.True.
            End Do
            loop=.False.
            num_tot=model_data%input%species(iin)%atoms_per_species
          End If
        End If
      End Do
      num_nn=nacc  
      list_nn=list_nn_next
      list_nn_next=0 
      If (nacc==0) loop=.False.
    End Do

    If (model_data%input%species(iin)%atoms_per_species/=num_tot) Then
      Write (messages(1),'(a,i4,3a)') '***ERROR: Trouble to identify unit ', jin, ' of species "',&
                               & Trim(model_data%input%species(iin)%tag),&
                               &'" from the input coordinates of the input structure.'
      Write (messages(2),'(a)')  'Possible reasons:' 
      Write (messages(3),'(a)')  ' 1) wrong input geometry, where one(or more) species are dissociated (e.g. H3O instead of H2O).&
                           & If this is the case, the user must provide an input model consistent with the definition of species.'

      Write (messages(4),'(a,f4.2,a)')  ' 2) the value for the bond cutoff (', bond_cutoff, &
                                    & ' Angstrom) as specified in &species_components&
                                    & is not adequate. The user must adjust this value'

      Write (messages(5),'(a)') ' 3) the species represents a fixed crystal structure with bonds periodically repeated.&
                               & If this is the case, the user must define the species as "fixed" in &species'
      Write (messages(6),'(a)') ' 4) the cell vectors are not strictly consistent with the geometry of the input structure.'
      Call info(messages, 6)
      Call error_stop(' ')
    End If

  End Subroutine find_neighbours
  
  Subroutine normal_along_vector(a,normal)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the normal corresponding to a lattice vector 
    !
    ! author    - i.scivetti Feb 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In   )   :: a(3)
    Real(Kind=wp), Intent(  Out)   :: normal(3)
 
    Real(Kind=wp) :: norm
    Integer(Kind=wi) :: i

    norm=0.0_wp

    Do i=1, 3
      norm=norm+a(i)**2
    End Do

    norm=sqrt(norm)
    normal=a/norm

  End Subroutine normal_along_vector 
  
 
  Subroutine check_consistency_input_model(model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the consistency between the simulation cell and the
    ! atomic coordinates of the input model 
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type), Intent(InOut) :: model_data

    Integer(Kind=wi) :: i, j
    Real(Kind=wp)    :: dist, dist0
    Logical          :: error 

    Real(Kind=wp)       :: a(3), b(3), r(3)
    Character(Len=256)  :: message, messages(2)
    Real(Kind=wp)       :: r0(model_data%input%num_atoms,3)
    Logical             :: changed_geo(model_data%input%num_atoms)

    Do i=1, model_data%input%num_atoms
      r0(i,:)=model_data%input%atom(i)%r(:)
    End Do

    error=.False.

    Write (messages(1),'(a)') '***ERROR: Inconsistency found between the atomic coordinates&
                            & and the simulation cell vectors. Please verify input coordinates and cell'

    ! Move all inside the cell
    Do i = 1, model_data%input%num_atoms
      r(:)=model_data%input%atom(i)%r(:)
      Call check_PBC(r, model_data%input%cell, model_data%input%invcell, 1.0_wp, changed_geo(i))
      If (changed_geo(i)) Then
        model_data%input%atom(i)%r(:)=r(:)
      End If
    End Do

    ! Try to move atoms again....if any atom is now moved, there is an inconsistency
    ! between the coordinates and the cell 
    Do i =  1, model_data%input%num_atoms
      If (changed_geo(i)) Then
        r(:)=model_data%input%atom(i)%r(:)
        Call check_PBC(r, model_data%input%cell, model_data%input%invcell, 1.0_wp, changed_geo(i))
        If (changed_geo(i)) Then
          ! Once again. If any atom is now moved, there is an inconsistency
          Call check_PBC(r, model_data%input%cell, model_data%input%invcell, 1.0_wp, changed_geo(i))
        End If  
        If (changed_geo(i)) Then
          Write (message,'(a,i4)') '***PROBLEMS with atom',  i
          Call info(message, 1)
          error=.True.
        End If
      End If
    End Do

    If (error) Then
       If (.Not. model_data%centre_electrode%stat) Then                 
         Write (messages(2),'(a)')    '*** The user should try by setting the "centre_electrode"&
                                     & directive to .True.'
         Call info(messages,2)
       Else
         Call info(messages,1)
       End If
       Call error_stop(' ')
    End If

    Do i = 1, model_data%input%num_atoms-1
      Do j = i+1, model_data%input%num_atoms
        a(:)=model_data%input%atom(i)%r(:)
        b(:)=model_data%input%atom(j)%r(:)
        Call compute_distance_PBC(a, b, model_data%input%cell, model_data%input%invcell, dist)
        a(:)=r0(i,:)
        b(:)=r0(j,:)
        Call compute_distance_PBC(a, b, model_data%input%cell, model_data%input%invcell, dist0)
        If (Abs(dist-dist0)>length_tol) Then
          Write (message,'(a,2(i7,a))') '***PROBLEMS: Distance between atom ', i, ' and ', j, &
                                   &' does not comply with the crystal symmetry imposed by the cell vectors.'
          Call info(message,1)
          Call info(messages,1)
          Call error_stop(' ')
        Else  
          If (dist < min_intra_bond) Then
            Write (message,'(a,2(i7,a))') '***PROBLEMS: Distance between atom ', i, ' and ', j, &
                                     &' is too short. It is likely that the cell dimensions are&
                                     & shorter than the system size!'
            Call info(message,1)
            Call info(messages,1)
            Call error_stop(' ')
          End If
        End If
      End Do
    End Do
 
  End Subroutine check_consistency_input_model
  
  Subroutine check_PBC(v_cart, basis, inv_basis, ratio, changed_geo)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the components of given vector "v_cart" 
    ! (given in cartesian coordinates of the Euclidean space) in terms of a general
    ! 3D vector basis, specified by the 3x3 matrix "basis".
    ! The inverse of the basis matrix is named "inv_basis".
    ! The factor "ratio" is used to evaluate how large these components are, 
    ! and modify the "vect" accordingly to accound for PBC.
    !
    ! The input value of "ratio" depends on the quantity to evalue. Thus,
    ! - ratio=0.5 is used for nearest-neighbour distances 
    ! - ratio=1.0 is used to evaluate if atomic positions lie within the volume 
    !   defined by the basis.
    !
    ! Logical variable "changed_geo" is used to check is the vector has been modified 
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp), Intent(InOut) :: v_cart(3)
    Real(Kind=wp), Intent(In   ) :: basis(3,3)
    Real(Kind=wp), Intent(In   ) :: inv_basis(3,3)
    Real(Kind=wp), Intent(In   ) :: ratio 
    Logical,       Intent(  Out) :: changed_geo                     

    Real(Kind=wp) :: v_direct(3), limit1, limit2
    Integer(Kind=wi) :: ir, i
    Logical          :: flag

    changed_geo=.False.
 
    If (Abs(ratio-0.5_wp) < epsilon(ratio)) Then
      limit1=ratio
      limit2=-ratio
    Else If (Abs(ratio-1.0_wp) < epsilon(ratio)) Then
      limit1=ratio+length_tol
      limit2=-length_tol
    End If

    ! Express vector difference in terms of the cell vectors
    v_direct= MatMul(v_cart, inv_basis)

    ! PCB effect
    i=1
    flag=.True.
    Do While (i< 4 .And. flag)
      Do ir = 1, 3
        If (v_direct(ir) > limit1) Then
           v_cart(:)= v_cart(:) - basis(ir,:)
           changed_geo=.True.
        Else If (v_direct(ir) < limit2) Then
           v_cart(:)= v_cart(:) + basis(ir,:)
           changed_geo=.True.
        End If
        If (changed_geo) Then
          v_direct= MatMul(v_cart, inv_basis)
        End If
      End Do
      If (.Not. changed_geo) Then
        flag=.False.      
      End If        
      i=i+1
    End Do

  End Subroutine check_PBC   
    
  Subroutine compute_distance_PBC(a, b, cell, invcell, dist)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the distance between atoms with PCB 
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In   ) :: a(3), b(3)
    Real(Kind=wp), Intent(In   ) :: cell(3,3)
    Real(Kind=wp), Intent(In   ) :: invcell(3,3)
    Real(Kind=wp), Intent(  Out) :: dist

    Real(Kind=wp) :: Dr_cart(3)
    Logical :: modified
    Integer :: ir

    ! Vector difference
    Do ir=1,3
      Dr_cart(ir)=a(ir)-b(ir) 
    End Do

    ! Find the vector difference for the nearest neighbours (NN)
    Call check_PBC(Dr_cart, cell, invcell, 0.5_wp, modified)
    ! Calculate norm
    dist=norm2(Dr_cart)

  End Subroutine compute_distance_PBC

  Subroutine check_orthorhombic_cell(A, flag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check if the simulation cell 'A' is orthorhombic or not. 
    ! vectors. If A is not orthorhombic, flag will be set to False
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In   )  :: A(3,3)
    Logical,       Intent(  Out)  :: flag
   
    Integer(Kind=wi) :: i, j

    i=1
    flag=.True.
    Do i = 1, 2
     Do j = i+1, 3
       If (Abs(Dot_product(A(i,:), A(j,:)))>epsilon(1.0_wp)) Then
         flag=.False.      
       End If        
     End Do
    End Do 

  End Subroutine check_orthorhombic_cell

  Subroutine define_model_cell(model_data, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to define the cell of the atomistic model 
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),  Intent(InOut) :: simulation_data
    Type(model_type),  Intent(InOut) :: model_data
    
    Character(Len=256) :: message
    Integer(Kind=wi)   :: i

    ! Set size for the output sample 
    Do i = 1, 3 
       model_data%sample%cell(i,:)=model_data%input%cell(i,:)*model_data%repeat_input_model%value(i) 
    End Do

    !Compute cell properties and inverse
    Call about_cell(model_data%sample%cell,model_data%sample%invcell,model_data%sample%cell_length)    

    If (model_data%analysis%type /= 'only_simulation_directives' ) Then 
      ! Print a message with the size of the model
      If (model_data%repeat_input_model%fread) Then
        If (model_data%repeat_input_model%value(1)==1 .And. &
            model_data%repeat_input_model%value(2)==1 .And. &
            model_data%repeat_input_model%value(3)==1) Then
            Write (message,'(a,3(i3,a))') 'The surface size for the generated model is the same as for the input structure'
        Else
          If (Trim(model_data%normal_vector%type) == 'c1') Then
            Write (message,'(a,2(i3,a))') 'The surface size for the generated model is (', &
                                        & model_data%repeat_input_model%value(2), 'x',&
                                        & model_data%repeat_input_model%value(3), &
                                        & ') times the surface size of the input model.'
          Else If (Trim(model_data%normal_vector%type) == 'c2') Then
            Write (message,'(a,2(i3,a))') 'The surface size for the generated model is (', &
                                        & model_data%repeat_input_model%value(1), 'x',&
                                        & model_data%repeat_input_model%value(3), &
                                        & ') times the surface size of the input model.'
          Else If (Trim(model_data%normal_vector%type) == 'c3') Then
            Write (message,'(a,2(i3,a))') 'The surface size for the generated model is (', &
                                        & model_data%repeat_input_model%value(1), 'x',&
                                        & model_data%repeat_input_model%value(2), &
                                        & ') times the surface size of the input model.'
          End If
                  End If 
      Else
        Write (message,'(a)') 'By default, the surface size of the generated model is the same than the surface&
                                 & size of the input model. The user can increase the surface size using the&
                                 & "repeat_input_model" directive.'
      End If 
      Call info(message,1)
    End If
    
    ! In case the user wants to generate input files for simulations
    If (simulation_data%generate) Then
     !Call large_cell(model_data, simulation_data) 
      ! copy cell info to simulation_data%cell 
      simulation_data%cell=model_data%sample%cell
      simulation_data%cell_length=model_data%sample%cell_length
      ! If the size of the cell is large...
      simulation_data%large_cell=.False.
      Do i=1,3
        If (model_data%sample%cell_length(i)>large_cell_limit) Then
          simulation_data%large_cell=.True.
        End If
      End Do
    End If

  End Subroutine define_model_cell  

  
  Subroutine compute_area_slab(a, b, area)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the surface area from the 
    ! two vectors that define the 2D periodicity 
    !
    ! author    - i.scivetti Feb 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In   ) :: a(3)
    Real(Kind=wp), Intent(In   ) :: b(3)
    Real(Kind=wp), Intent(InOut) :: area 

    Real(Kind=wp) :: cross(3)

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)

    ! Calculate area
    area= norm2(cross)

  End Subroutine compute_area_slab  

  Subroutine remove_species(T, types_species, arr_added_species, stage)  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to remove atoms from the model  
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(sample_type), Intent(InOut) :: T 
    Integer(Kind=wi),  Intent(in   ) :: types_species
    Character(Len=*),  Intent(In   ) :: arr_added_species
    Character(Len=*),  Intent(In   ) :: stage 

    Integer(Kind=wi) :: i, j, k, m, l, na

    Integer(Kind=wi) :: list_rm(max_num_species_units)
    Real(Kind=wp) :: rn
    Logical       :: hit(max_at_species)

    list_rm=0
    Do i = 1, types_species
      If (T%species(i)%change_content) Then
        If (T%species(i)%D_num<0) Then
          If (Trim(arr_added_species)=='random'    .Or. &
             (Trim(arr_added_species)=='deposited' .And. Trim(stage)=='sample')) Then
            hit=.False.
            k=1 
            Do While (k <= (-T%species(i)%D_num)) 
              Call random_number(rn)
              na=floor(T%species(i)%num*rn)+1
              If (.Not. T%species(i)%units(na)%vanish) Then
                list_rm(k)=na 
                T%species(i)%units(na)%vanish=.True.
                k=k+1
              End If
            End Do
          ElseIf (Trim(arr_added_species)=='deposited' .And. Trim(stage)=='input') Then
            k=1
            Do j=T%species(i)%num, T%species(i)%num+T%species(i)%D_num+1, -1
              list_rm(k)=j
              T%species(i)%units(j)%vanish=.True.
              k=k+1
            End Do
          End If
          ! remove 
          Do k =1, -T%species(i)%D_num
             j=list_rm(k)
             Do m= 1, T%species(i)%atoms_per_species
               l=T%species(i)%units(j)%list(m)
               T%atom(l)%vanish=.True.
             End Do  
          End Do

        End If
      End If
    End Do

  End Subroutine remove_species  
  
  Subroutine insert_species(arr_added_species, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to insert species into the model  
    !
    ! author    - i.scivetti Feb 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=32),   Intent(In   ) :: arr_added_species
    Type(model_type),    Intent(InOut) :: model_data

    Logical :: flag
    Integer(Kind=wi) :: i, j, index, loop(3)
    Real(Kind=wp)    :: min_diff, diff, length, shift_level, scaling

    Do i=1, model_data%types_species  
      model_data%input%species(i)%num_extra=model_data%input%species(i)%num
    End Do

    ! We first try to insert molecules, then atoms
    Do i=1, model_data%types_species 
      If (model_data%input%species(i)%D_num>0) Then
        If (model_data%input%species(i)%topology=='molecule') Then
          Call read_molecular_units(model_data%input%species(i), model_data%species_info(i)%bond_cutoff)  
        ElseIf (model_data%input%species(i)%topology=='atom') Then
          model_data%input%species(i)%definition%r0(1,:)=0.0_wp 
          model_data%input%species(i)%definition%tag(1)=model_data%input%species(i)%component%tag(1)
          model_data%input%species(i)%definition%element(1)=model_data%input%species(i)%component%element(1)
          model_data%input%species(i)%definition%atomic_number(1)=model_data%input%species(i)%component%atomic_number(1)
        End If
      End If
    End Do

    ! Find the total number of grid points based on the delta_space
    Do i=1,3
      model_data%npoints(i)=Int(model_data%input%cell_length(i)/model_data%delta_space%value)
      model_data%scan(i)=1 
    End Do
     
    If (Trim(model_data%normal_vector%type)=='c3') Then
      index=3
      length=model_data%input%cell_length(3)
    ElseIf (Trim(model_data%normal_vector%type)=='c2') Then
      index=2
      length=model_data%input%cell_length(2)
    ElseIf (Trim(model_data%normal_vector%type)=='c1') Then
      index=1
      length=model_data%input%cell_length(1)
    End If
      
    min_diff=Huge(1.0_wp)
    Do i=1, model_data%npoints(index)
      diff=Abs((i-1)*length/model_data%npoints(index)-Abs(model_data%add_species_from%value))
      If (diff < min_diff) Then
        min_diff=diff
        model_data%scan(index)=i
      End If
    End Do
      
    If (Trim(arr_added_species)=='random') Then
      model_data%npoints(index)=Int((model_data%input%cell_length(index)-&
                             & model_data%add_species_from%value)/model_data%delta_space%value)
      shift_level=model_data%add_species_from%value/model_data%input%cell_length(index)
      scaling=Abs(model_data%input%cell_length(index)-model_data%add_species_from%value)/model_data%input%cell_length(index)
    Else If (Trim(arr_added_species)=='deposited') Then
      If (Trim(model_data%normal_vector%type)=='c3') Then
        loop(1)=1
        loop(2)=2
        loop(3)=3
      Else If (Trim(model_data%normal_vector%type)=='c2') Then
        loop(1)=1
        loop(2)=3
        loop(3)=2      
      Else If (Trim(model_data%normal_vector%type)=='c1') Then
        loop(1)=2
        loop(2)=3
        loop(3)=1      
      End If    
    End If

    j=1
    flag=.True.
    Do While (flag)
      flag=.False.
      Do i=1, model_data%types_species 
        If (model_data%input%species(i)%change_content) Then
           If (j<= model_data%input%species(i)%D_num) Then  
              flag=.True.
               If (Trim(arr_added_species)=='random') Then 
                 Call fitting_species_random(model_data%input, model_data%npoints,&
                                   & i, j, index, shift_level, scaling, model_data%distance_cutoff, model_data%rotate_species%stat)
               ElseIf (Trim(arr_added_species)=='deposited') Then
                 Call fitting_species_deposited(model_data%input, loop, model_data%scan, model_data%npoints, &
                                   & i, j, model_data%distance_cutoff, model_data%rotate_species%stat)
               End If
           End If
        End If
      End Do
      j=j+1
    End Do

  End Subroutine insert_species
  
  Subroutine fitting_species_deposited(T, loop, scan, npoints, isp, junit, distance_cutoff, rotate)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to fit species over the input electrode surface
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(sample_type),  Intent(InOut) :: T
    Integer(Kind=wi),   Intent(In   ) :: loop(3)
    Integer(Kind=wi),   Intent(InOut) :: scan(3)
    Integer(Kind=wi),   Intent(InOut) :: npoints(3)
    Integer(Kind=wi),   Intent(In   ) :: isp, junit
    Type(in_param),     Intent(In   ) :: distance_cutoff
    Logical,            Intent(In   ) :: rotate 

    Character(Len=256)   :: messages(6)
    Integer(Kind=wi) :: ic(3)
    Integer(Kind=wi) :: i, j, l
    Real(Kind=wp) :: centre(3), a(3), b(3)
    Real(Kind=wp) :: r0_rot(max_at_species,3)
    Real(Kind=wp) :: dist
    Logical       :: nofit, too_short

    nofit=.True.
    ic=scan

    Do While (ic(loop(3)) < npoints(loop(3)) .And. nofit)
      Do While (ic(loop(2)) < npoints(loop(2)) .And. nofit)
        Do While (ic(loop(1)) < npoints(loop(1)) .And. nofit)
          Do l=1,3
            centre(l)=(ic(l)-1.0_wp)/(npoints(l))
          End Do
          centre=MatMul(centre,T%cell)
          If (T%species(isp)%topology=='molecule') Then
            If (rotate) Then  
              Call rotate_species(T%species(isp)%definition%r0, T%species(isp)%atoms_per_species, r0_rot)
            Else
              r0_rot=T%species(isp)%definition%r0
            End If
          Else
           r0_rot(1,:)=T%species(isp)%definition%r0(1,:)
          End If

          Do l=1,3
           Do i=1, T%species(isp)%atoms_per_species
             r0_rot(i,l)=r0_rot(i,l)+centre(l)
           End Do
          End Do

          i=1
          too_short=.False.
          Do While (i <= T%species(isp)%atoms_per_species .And. .Not.(too_short))
            j=1
            Do While (j <= T%num_atoms_extra .And. .Not.(too_short))
              If (.Not.(T%atom(j)%vanish)) Then
              a(:)=r0_rot(i,:)
              b(:)=T%atom(j)%r(:)
              Call compute_distance_PBC(a, b, T%cell, T%invcell, dist)
              If (dist<distance_cutoff%value) Then
                too_short=.True.
              End If
              End If
              j=j+1
            End Do
            i=i+1
          End Do

          If (.Not. too_short) Then
            scan=ic
            ! Assign arrays
            nofit=.False.
            Do i=1,T%species(isp)%atoms_per_species
              T%atom(T%num_atoms_extra+i)%r(:)=r0_rot(i,:)
              T%atom(T%num_atoms_extra+i)%element=T%species(isp)%definition%element(i)
              T%atom(T%num_atoms_extra+i)%atomic_number=T%species(isp)%definition%atomic_number(i)
              T%atom(T%num_atoms_extra+i)%tag=T%species(isp)%definition%tag(i)
              T%atom(T%num_atoms_extra+i)%vanish=.False.
              T%atom(T%num_atoms_extra+i)%in_species=.True.
              T%atom(T%num_atoms_extra+i)%dynamics=.True.        
            End Do
            Do i=1, T%species(isp)%atoms_per_species
              T%species(isp)%units(T%species(isp)%num_extra+1)%list(i)=T%num_atoms_extra+i
            End Do
            T%num_atoms_extra=T%num_atoms_extra+T%species(isp)%atoms_per_species
            T%species(isp)%num_extra=T%species(isp)%num_extra+1
          End If
          ic(loop(1))=ic(loop(1))+1
        End Do
        If (nofit) then
          ic(loop(1))=1
          ic(loop(2))=ic(loop(2))+1
        End If
      End Do
        If (nofit) then
          ic(loop(2))=1
          ic(loop(3))=ic(loop(3))+1
        End If
    End Do

    scan=ic
    
    If (nofit) Then
      Write (messages(1),'(2(a,i3),3a)') '***ERROR: fail to deposit unit ', junit, ' (out of ',&
                               T%species(isp)%D_num, ') for species "', Trim(T%species(isp)%tag), '"'
      Write (messages(2),'(a)') 'This is probably due to:'
      Write (messages(3),'(a)') '1) an incorrect geometry for the input model. Either the surface area is not large enough to&
                              & accommodate the required species or there is not enough vacuum region. In case of the latter,&
                              & check the value of "add_species_from" if defined.'
      If (distance_cutoff%fread) Then
        Write (messages(4),'(a,f6.2,a)') '2) a rather large value for the cutoff distance between species, set to ',&
                                       & distance_cutoff%value, ' Angstrom.&
                                       & Please reduce the value of directive "distance_cutoff".'
      Else
        Write (messages(4),'(a,f4.2,a)') '2) the default cutoff distance between species, set to ',&
                                       & distance_cutoff%value, ' is not enough to fit the species. If the is the case, &
                                       & please increase the vaccum region.'
      End If
      If (.Not. rotate) Then
        Write (messages(5),'(a)') '3) for large molecular species, there could be an inconsistency between the orientation/size&
                                & of the species to deposit and the surface of the modelled substrate. Please check.'
      Else
        Write (messages(5),'(a)') '3) for large molecular species, it is convenient to set "rotate_species" to .False. and adapt&
                                & the orientation of the molecule to the substrate.'
      End If
      Write (messages(6),'(a)') '4) the option set for "normal_vector" is not consistent with the input model.'
      
      Call info(messages,6)
      Call execute_command_line('rm  RECORD_MODELS')
      Call error_stop(' ')
    End If

  End Subroutine fitting_species_deposited

  Subroutine fitting_species_random(T, npoints, isp, junit, index, shift_level, scaling, distance_cutoff, rotate)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to randomly fit electrolyte species over the electrode 
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(sample_type),  Intent(InOut) :: T
    Integer(Kind=wi),   Intent(InOut) :: npoints(3)
    Integer(Kind=wi),   Intent(In   ) :: isp
    Integer(Kind=wi),   Intent(In   ) :: junit
    Integer(Kind=wi),   Intent(In   ) :: index
    Real(Kind=wp),      Intent(In   ) :: shift_level
    Real(Kind=wp),      Intent(In   ) :: scaling
    Type(in_param),     Intent(In   ) :: distance_cutoff
    Logical,            Intent(In   ) :: rotate 
    
    Integer(Kind=wi), Dimension(3) :: ic
    Integer(Kind=wi) :: i, j, l
    Real(Kind=wp) :: centre(3), centre2(3), a(3), b(3)
    Real(Kind=wp) :: r0_rot(max_at_species,3)
    Real(Kind=wp) :: dist, rn
    Integer(Kind=li) :: nlimit, ncount 
    Logical       :: nofit, too_short
    Character(Len=256)   :: messages(6)

    nofit=.True.
    nlimit=1
    Do i=1,3 
      nlimit=nlimit*npoints(i)
    End Do
    nlimit=nlimit*times_number_attempts
    ncount=0

    Do While (nofit .And. ncount<=nlimit)
      ncount=ncount+1
      Do l=1, 3
        Call random_number(rn)
        ic(l)=floor(npoints(l)*rn+1)
      End Do
      Do l=1,3
        If (l==index) Then
          centre2(l)=(ic(l)-1.0_wp)/(npoints(l))*scaling+shift_level
        Else
          centre2(l)=(ic(l)-1.0_wp)/(npoints(l))
        End If
      End Do  
      centre=MatMul(centre2, T%cell)
      If (T%species(isp)%topology=='molecule') Then
        If (rotate) Then
          Call rotate_species(T%species(isp)%definition%r0, T%species(isp)%atoms_per_species, r0_rot)
        Else
          r0_rot=T%species(isp)%definition%r0
        End If
      Else
       r0_rot(1,:)=T%species(isp)%definition%r0(1,:)
      End If
   
      Do i=1, T%species(isp)%atoms_per_species
        Do l=1,3
          r0_rot(i,l)=r0_rot(i,l)+centre(l)
        End Do
      End Do
  
      i=1
      too_short=.False. 
      Do While (i <= T%species(isp)%atoms_per_species .And. .Not.(too_short))
        j=1
        Do While (j <= T%num_atoms_extra .And. .Not.(too_short))
          If (.Not.(T%atom(j)%vanish)) Then
          a(:)=r0_rot(i,:)
          b(:)=T%atom(j)%r(:)
          Call compute_distance_PBC(a, b, T%cell, T%invcell, dist)
          If (dist<distance_cutoff%value) Then
            too_short=.True.
          End If
          End If
          j=j+1
        End Do
        i=i+1
      End Do
      
      If (.Not. too_short) Then
        ! Assign arrays  
        nofit=.False.
        Do i=1,T%species(isp)%atoms_per_species
          If (T%num_atoms_extra+i > T%max_atoms) Then 
            Write(messages(1),'(a)')    '***ERROR: it seems the user wants to build a model with a large number of atoms.' 
            Write(messages(2),'(a,i5)') '   The maximum number of atoms to build a model using the input structure is set to ', &
                                     & T%max_atoms
            Write(messages(3),'(a)')    '   However, the required number of atoms as specified in the settings is larger&
                                     & than this limit.'
            Write(messages(4),'(a)')    '   The user must define a larger value for "multiple_input_atoms" (10000 by default)' 
            Call info(messages, 4) 
            Call error_stop(' ') 
          End If 
          T%atom(T%num_atoms_extra+i)%r(:)=r0_rot(i,:)
          T%atom(T%num_atoms_extra+i)%element=T%species(isp)%definition%element(i)
          T%atom(T%num_atoms_extra+i)%atomic_number=T%species(isp)%definition%atomic_number(i)
          T%atom(T%num_atoms_extra+i)%tag=T%species(isp)%definition%tag(i)
          T%atom(T%num_atoms_extra+i)%vanish=.False.
          T%atom(T%num_atoms_extra+i)%in_species=.True.
          T%atom(T%num_atoms_extra+i)%dynamics=.True.
        End Do    
        Do i=1, T%species(isp)%atoms_per_species
          T%species(isp)%units(T%species(isp)%num_extra+1)%list(i)=T%num_atoms_extra+i
        End Do
        T%num_atoms_extra=T%num_atoms_extra+T%species(isp)%atoms_per_species
        T%species(isp)%num_extra=T%species(isp)%num_extra+1
      End If
      ic(1)=ic(1)+1    
    End Do

    If (ncount>nlimit) Then
      Write (messages(1),'(2(a,i6),3a)') '***ERROR: fail to insert unit ', junit, ' (out of ',&
                               T%species(isp)%D_num, ') for species "', Trim(T%species(isp)%tag), '"'
      Write (messages(2),'(a)') 'The maximum limit of attemps to insert the species has been reached. Possible reasons:' 
      Write (messages(3),'(a)') '1) an incorrect geometry for the input model. Either the surface area is not large enough to&
                              & accommodate the required species or there is not enough vacuum region. In case of the latter,&
                              & check the value of "add_species_from" if defined.'
      If (distance_cutoff%fread) Then
        Write (messages(4),'(a,f6.2,a)') '2) a rather large value for the cutoff distance between species, set to ',&
                                       & distance_cutoff%value, ' Angstrom.&
                                       & Please reduce the value of directive "distance_cutoff".'
      Else
        Write (messages(4),'(a,f4.2,a)') '2) the default cutoff distance between species, set to ',&
                                       & distance_cutoff%value, ' is not enough to fit the species. If the is the case, &
                                       & please increase the vaccum region.'
      End If
      If (.Not. rotate) Then
        Write (messages(5),'(a)') '3) for large molecular species, there could be an inconsistency between the orientation/size&
                                & of the species to deposit and the surface of the modelled substrate. Please check.'
      Else
        Write (messages(5),'(a)') '3) for large molecular species, it is convenient to set "rotate_species" to .False. and adapt&
                                & the orientation of the molecule to the substrate.'
      End If
      Write (messages(6),'(a)') '4) the option set for "normal_vector" is not consistent with the input model.'
      
      Call info(messages,6)      
      Call execute_command_line('rm  RECORD_MODELS')
      Call error_stop(' ')
    End If

  End Subroutine fitting_species_random  

  Subroutine read_molecular_units(T, bond_cutoff)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read species units from file 
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(species_model),     Intent(InOut) :: T 
    Real(Kind=wp),           Intent(In   ) :: bond_cutoff 
 
    Logical  :: safe, error, loop
    Integer(Kind=wi)  :: i, j, m
    Integer(Kind=wi)  :: iunit, io, num_atoms
    Character(Len=256)   :: set_error, title
    Character(Len=256)   :: messages(4), error_format(10)
    Character(Len= 64)   :: filename 
   
    Integer(Kind=wi)  :: ncount(T%num_components) 

    Integer(Kind=wi) :: num_nn, num_tot
    Integer(Kind=wi), Dimension(max_nn) :: list_nn, list_nn_next
    Integer(Kind=wi) :: ka, k, nacc, knni, il, knn, nlist
    Real(Kind=wp), Dimension(3) :: a, b, dr, geom_centre
    Real(Kind=wp) :: dist
    Logical       :: loop_sp, is_nn 

    error=.False.
    filename=Trim(FOLDER_INPUT_GEOM)//'/'//Trim(T%tag)//'.xyz' 
    set_error='***ERROR in file '//Trim(filename)//':'

    error_format(1) ='The setting for each molecular species unit must be:'
    error_format(2) = ' '
    error_format(3) = 'Number of atoms per species (Nsp)'
    error_format(4) = 'Tag/description of the unit (compulsory)'
    error_format(5) = 'Element_1         X_1      Y_1      Z_1       Tag_1'
    error_format(6) = 'Element_2         X_2      Y_2      Z_2       Tag_2'
    error_format(7) = '...........       .....    .....    .....     ...... '
    error_format(8) = 'Element_Nsp       X_Nsp    Y_Nsp    Z_Nsp     Tag_Nsp'    
    error_format(9) = ' '
    error_format(10) = 'where the Tag correspond to the tag used for the component of each species,&
                     & as defined in &Species_Components'

    ! Check if the SETTINGS file exists and open it
    Inquire(File=filename, Exist=safe)
    If (.not.safe) Then
      Write (messages(1),'(4(1x,a))') '***ERROR - File', Trim(filename), 'not found for the&
                                 & specification of species unit', Trim(T%tag)
      Call info(messages,1)
      Call error_stop(' ')
    Else
      Open(Newunit=iunit, File=filename,Status='old')
    End If


    Read (iunit, Fmt=*, iostat=io) num_atoms
    If (io/=0) Then
      Write (messages(1),'(4a)') Trim(set_error), ' First line must be the number of atoms (integer) that constitute species "',&
                                 Trim(T%tag), '". Please check'
      Call info(messages,1)
      error=.True.
    Else
      If (num_atoms/=T%atoms_per_species) Then
        Write (messages(1),'(4a,i3)') Trim(set_error), ' First line must be the number of atoms that constitute species "',&
                               & Trim(T%tag), '", and equal to ', T%atoms_per_species
        Call info(messages,1)
        error=.True.
      End If
    End If

    If (error) Then
      Call info(error_format,10)
      Call error_stop(' ')
    End If

    ! Read title
    Read (iunit, Fmt=*, iostat=io) title


    ! Read specification for the species unit 
    i=0
    Do While (i < T%atoms_per_species)
        i=i+1
        Read (iunit, Fmt=*, iostat=io) T%definition%element(i), (T%definition%r0(i,m), m=1,3), T%definition%tag(i)
        If (io/=0) Then
          If (i== T%atoms_per_species) Then
            Write (messages(1),'(2a,i5)') Trim(set_error), ' Missing specification for atom ', i
          Else
            Write (messages(1),'(2a,i5)') Trim(set_error), ' Wrong specification for atom', i
          End If
          Call info(messages, 1)
          Call info(error_format,10)
          Call error_stop(' ')
        End If
    End Do

    Close(iunit)

    ! Count the number of atoms por component
    ncount=0
    Do i=1, T%atoms_per_species
      loop=.True.
      j=1
      Do While (j <= T%num_components .And. loop)
        If (Trim(T%definition%tag(i))==Trim(T%component%tag(j))) Then
          T%definition%atomic_number(i)=T%component%atomic_number(j)
          ncount(j)=ncount(j)+1
          loop=.False.
          If (Trim(T%definition%element(i))/=Trim(T%component%element(j))) Then
            Write (messages(1),'(4a,i3,3a,i3,3a)') Trim(set_error), ' Setting "', Trim(T%definition%element(i)), &
                                        & '" for the species element of atom ', i, ' does not agree with&
                                        & setting "', Trim(T%component%element(j)), '" specified for component ', &
                                        & j, ' of species "', Trim(T%tag), '" in  &Species_Components'
            Call info(messages,1)
            Call error_stop(' ')  
          End If
        End If
        j=j+1
      End Do
      If (loop) Then
        Write (messages(1),'(4a,i3,3a)') Trim(set_error), ' Setting "', Trim(T%definition%tag(i)), &
                                    & '" for the tag of atom ', i, ' has not been assigned any component tag&
                                    & of species "', Trim(T%tag), '" in  &Species_Components'
        Call info(messages,1)
        Call error_stop(' ')  
      End If
    End Do
    
    Do j=1, T%num_components
        If (ncount(j)/=T%component%N0(j)) Then
           Write (messages(1),'(4a,i3,a,i3,a)') Trim(set_error), ' The number of components with tag "', Trim(T%component%tag(j)),&
                                         & '" is equal to ', ncount(j), ', which is different from the amout of ', &
                                         & T%component%N0(j), ' specified in &Species_Components. Please check the&
                                         & definition for the molecular geometry.'    
          Call info(messages,1)
          Call error_stop(' ')  
        End If
    End Do

    ! Check it is a sensible geometry
    T%definition%in_species=.False.
    num_nn=1
    list_nn(1)=1
    ka=1
    loop=.True. 
    num_tot=0
    list_nn_next=0 
    nlist=1


    Do While (loop)
      nacc=0
      Do knn=1, num_nn
        ka=list_nn(knn)
        Do k=1, T%atoms_per_species
          loop_sp=.True.
          If (k/=ka .And. (.Not. T%definition%in_species(k)) .And. (.Not. T%definition%in_species(ka))) Then
             ! Compute distances
             a=T%definition%r0(k,:)
             b=T%definition%r0(ka,:)
             dr=a-b
             dist=Norm2(dr)
             ! Check if distances are within a given cutoff
             If (dist<bond_cutoff) Then 
               If (dist<min_intra_bond) Then
                 Write (messages(1),'(2a)') Trim(set_error), ' wrong input geometry'
                 Write (messages(2),'(2(a,i4))') 'Intermolecular distance between atoms ', k, ' and ', ka
                 Write (messages(3),'(a,f4.2,a)') 'is shorter than the input minimum distance criteria of ', min_intra_bond,&
                                                 & ' Angstrom for bonding. Please review the input geometry.'  
                 Call info(messages,3)
                 Call error_stop(' ')
               Else
                 loop_sp=.False.
                 is_nn=.True.
                 il=1
                 Do While (il <= max_nn .And. is_nn) 
                   If (k==list_nn_next(il)) Then
                     is_nn=.False.
                   End If
                   il=il+1
                 End Do 
                 If (is_nn) Then    
                   nacc=nacc+1
                   list_nn_next(nacc)=k
                 End If
               End If  
             End If
          End If
        End Do
        If (.Not. T%definition%in_species(ka)) Then
          T%definition%in_species(ka)=.True.
          num_tot=num_tot+1
          If (num_tot+nacc == T%atoms_per_species) Then
            Do knni= 1, nacc
              T%definition%in_species(list_nn_next(knni))=.True.
            End Do
            loop=.False.
            num_tot=T%atoms_per_species
          End If
        End If
      End Do
      num_nn=nacc  
      list_nn=list_nn_next
      list_nn_next=0 
      If (nacc==0) loop=.False.
    End Do

    If (T%atoms_per_species/=num_tot) Then
      Write (messages(1),'(2a)') Trim(set_error), ' Problems with the input geometry'
      Write (messages(2),'(a)')  'Possible reasons:' 
      Write (messages(3),'(a)')          ' 1) wrong input geometry, where species are dissociated (e.g. H3O instead of H2O).&
                                        & If this is the case, the user must provide an input model consistent with the&
                                        & definition of species.' 
      Write (messages(4),'(a,f4.2,3a)')  ' 2) the value for the bond cutoff (', bond_cutoff, &
                                        & ' Angstrom) as specified for species ',  Trim(T%tag) ,' in  &species_components&
                                        & is not sufficient. The user must adjust this value' 
      Call info(messages, 4)
      Call error_stop(' ')
    End If

    ! Compute the geometric centre of the species and refer atomic coordinates to such point in space 
    geom_centre=0.0_wp
    Do i =1, T%atoms_per_species
      Do k=1,3
        geom_centre(k)=geom_centre(k)+T%definition%r0(i,k)
      End Do
    End Do

    geom_centre=geom_centre/T%atoms_per_species

    Do i =1, T%atoms_per_species
      Do k=1,3
        T%definition%r0(i,k)=T%definition%r0(i,k)-geom_centre(k)
      End Do
    End Do  
    
  End Subroutine read_molecular_units  

  Subroutine rotate_species(r0, nat, r0_rot)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to rotate molecular species 
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp),    Intent(In   ) :: r0(max_at_species,3)
    Integer(Kind=wi), Intent(In   ) :: nat
    Real(Kind=wp),    Intent(  Out) :: r0_rot(max_at_species,3)

    Integer(Kind=wi) :: i
    Real(Kind=wp) :: rn, alpha, beta, gamma
    Real(Kind=wp) :: rot_matrix(3,3), a(3), b(3)

    !Define the Euler angles
    Call random_number(rn)
    alpha=twopi*rn
    Call random_number(rn)
    beta =twopi*rn
    Call random_number(rn)
    gamma=twopi*rn
    ! Build rotation matrix
    rot_matrix(1,1)=cos(alpha)*cos(beta)
    rot_matrix(2,1)=sin(alpha)*cos(beta)
    rot_matrix(3,1)=-sin(beta)
    rot_matrix(1,2)=cos(alpha)*sin(beta)*sin(gamma)-sin(alpha)*cos(gamma)
    rot_matrix(2,2)=sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma)
    rot_matrix(3,2)=cos(beta)*sin(gamma)
    rot_matrix(1,3)=cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma)
    rot_matrix(2,3)=sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma)
    rot_matrix(3,3)=cos(beta)*cos(gamma)

    Do i= 1, nat      
      a(:)=r0(i,:)
      b=MatMul(a, rot_matrix)
      r0_rot(i,:)=b(:)
    End Do

  End Subroutine rotate_species   
  
  Subroutine define_repeated_model(model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Uses the modified input model to build the output sample
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),    Intent(InOut) :: model_data

    Integer(Kind=wi) :: i, j, k
    Integer(Kind=wi) :: ic(3), ic1, ic2, ic3
    Integer(Kind=wi) :: na, nat_extra, nsp_extra
    Integer(Kind=wi) :: ind_at, ind_sp, isp 
    Integer(Kind=wi) :: ris(3)
    Real(Kind=wp)    :: shift(3)
    Character(Len=256) :: messages(3)

    ris=model_data%repeat_input_model%value
    nat_extra=model_data%input%num_atoms_extra

    model_data%sample%min_species=model_data%input%min_species*model_data%input_times
    model_data%sample%min_components=model_data%input%min_components

    ! Copy species
    Do i = 1, model_data%types_species    
      model_data%sample%species(i)%num=model_data%input%species(i)%num_extra*model_data%input_times
    End Do

    ! Copy the first num_atoms_extra atoms
    If (model_data%input%num_atoms_extra>model_data%sample%max_atoms) Then
       Write (messages(1),'(1x,1a)') '***ERROR: the initial allocated dimensions for the maximum number of atoms&
                                    & to generate models has been exceeded. '      
       Write (messages(2),'(1x,1a)') '   The user must define a larger value for the "multiple_output_atoms"&
                                    & directive, which is set to 2 by default.'
       Write (messages(3),'(1x,1a)') '   Try first by increasing this directive to 5 first. Increase it more if needed'
       Call info(messages, 3)
       Call error_stop(' ')
    End If            
    Do j =1 , model_data%input%num_atoms_extra
      model_data%sample%atom(j)=model_data%input%atom(j)      
    End Do
  
    Do ic3 = 1, ris(3)
      ic(3)=ic3
      Do ic2 = 1, ris(2)
        ic(2)=ic2
        Do ic1 = 1, ris(1)
          ic(1)=ic1
          ind_at=nat_extra*(ic1-1)+nat_extra*ris(1)*(ic2-1)+nat_extra*ris(1)*ris(2)*(ic3-1)
          ! Set atomic properties
          shift=0.0_wp
          Do j=1,3
            Do k=1,3
              shift(j)=shift(j)+model_data%input%cell(k,j)*(ic(k)-1)
            End Do
          End Do
          Do i=1, nat_extra
            na=i+ind_at
            model_data%sample%atom(na)%r=model_data%input%atom(i)%r+shift
            model_data%sample%atom(na)%element=model_data%input%atom(i)%element
            model_data%sample%atom(na)%atomic_number=model_data%input%atom(i)%atomic_number
            model_data%sample%atom(na)%tag=model_data%input%atom(i)%tag
            model_data%sample%atom(na)%in_species=model_data%input%atom(i)%in_species
            model_data%sample%atom(na)%vanish=model_data%input%atom(i)%vanish
            model_data%sample%atom(na)%dynamics=model_data%input%atom(i)%dynamics
          End Do
          !Set species
          Do isp = 1, model_data%types_species
            nsp_extra=model_data%input%species(isp)%num_extra
            ind_sp=nsp_extra*(ic1-1)+nsp_extra*ris(1)*(ic2-1)+nsp_extra*ris(1)*ris(2)*(ic3-1)
            If (model_data%input%species(isp)%change_content) Then
              Do j=1, model_data%input%species(isp)%num_extra
                If ((j+ind_sp) > max_num_species_units) Then
                  Call error_stop('***ERROR: requested size for the model structure exceeeds the stipulated limits.&
                                 & Please reduce the size using directive "repeat_input_model"') 
                Else
                   model_data%sample%species(isp)%units(j+ind_sp)%vanish = &
                        model_data%input%species(isp)%units(j)%vanish
                End If
                Do k= 1, model_data%input%species(isp)%atoms_per_species
                  model_data%sample%species(isp)%units(j+ind_sp)%list(k) = &
                                                    model_data%input%species(isp)%units(j)%list(k)+ind_at 
                End Do
              End Do
            End If
          End Do 
  
        End Do  
      End Do  
    End Do 
 
    Do isp = 1, model_data%types_species
      model_data%sample%species(isp)%num_show=0
      Do j = 1, model_data%sample%species(isp)%num
        If (.Not. model_data%sample%species(isp)%units(j)%vanish) Then
          model_data%sample%species(isp)%num_show=model_data%sample%species(isp)%num_show+1
        End If
      End Do
    End Do

  End Subroutine define_repeated_model

  Subroutine create_list_net_elements(T, types_species, both_surfaces)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Create list for the atoms 
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(sample_type), Intent(InOut) :: T
    Integer(Kind=wi),  Intent(In   ) :: types_species
    Logical,           Intent(In   ) :: both_surfaces    

    Integer(Kind=wi) :: icount, isp, j, i, inat

    ! Define the amount of net elements
    icount=0
    T%list%net_elements=0
    Do isp = 1, types_species
      Do j=1, T%species(isp)%num_components
        inat=0
        Do i=1, T%num_atoms
          If ((.Not. T%atom(i)%vanish) .And. (Trim(T%atom(i)%tag)==Trim(T%species(isp)%component%tag(j))) ) Then
            inat=inat+1
          End If
        End Do
        If (inat/=0) Then
          icount=icount+1
          T%list%N0(icount)=inat
          If (T%species(isp)%change_content .And. both_surfaces) Then
            T%list%N0(icount)=2*T%list%N0(icount)
          End If
          T%list%tag(icount)=Trim(T%species(isp)%component%tag(j))
          T%list%element(icount)=Trim(T%species(isp)%component%element(j))   
          T%list%net_elements=T%list%net_elements+1
        End If
      End Do
    End Do
    
  End Subroutine create_list_net_elements

  Subroutine set_electrode_boundary(model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to define the coordinate value from which the species will 
    ! be included
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),    Intent(InOut) :: model_data

    Integer(Kind=wi) :: i, j, k, indx
    Real(Kind=wp)    :: max_dist
    Real(Kind=wp)    :: dist
    Real(Kind=wp)    :: r(3), normal(3), atom(3)
    Logical          :: flag

    If (Trim(model_data%normal_vector%type)=='c3') Then
      indx=3
    ElseIf (Trim(model_data%normal_vector%type)=='c2') Then
      indx=2
    ElseIf (Trim(model_data%normal_vector%type)=='c1') Then
      indx=1
    End If

    r(:)=model_data%input%cell(indx,:)
    Call normal_along_vector(r,normal)
    ! Find the atom with the largest component along the normal vector to the surface 
    max_dist=-Huge(1.0_wp)
    Do i = 1, model_data%input%num_atoms
      flag=.True.
      j=1 
      Do While (j <= model_data%num_species%value .And. flag)
         If (Trim(model_data%species_info(j)%topology)=='electrode') Then
           k=1
           Do While (k <= model_data%species_info(j)%num_components .And. flag)
             If (Trim(model_data%input%atom(i)%tag) == Trim(model_data%species_info(j)%component%tag(k))) Then
               r(:)=model_data%input%atom(i)%r(:)
               dist=Dot_product(r,normal)
               If (dist > max_dist) Then
                 max_dist = dist
                 atom=r
               End If
               flag=.False.
             End If
             k=k+1
           End Do
         End If
        j=j+1
      End Do    
    End Do
    
   model_data%add_species_from%value=atom(indx)    
    
  End Subroutine set_electrode_boundary
  
  Subroutine surface_shift(model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to find the surface shift from the topmost and lowest atoms of
    ! the electrode model
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),    Intent(InOut) :: model_data

    Logical          :: flag
    Integer(Kind=wi) :: i, j, k, indx
    Real(Kind=wp)    :: min_dist, max_dist
    Real(Kind=wp)    :: dist
    Real(Kind=wp)    :: r(3), normal(3)
    Real(Kind=wp)    :: topmost(3), lowest(3)
    Character(Len=8)   :: tag_top, tag_low 
    Character(Len=256) :: messages(6)

    If (Trim(model_data%normal_vector%type)=='c3') Then
      indx=3
    ElseIf (Trim(model_data%normal_vector%type)=='c2') Then
      indx=2
    ElseIf (Trim(model_data%normal_vector%type)=='c1') Then
      indx=1
    End If

    r(:)=model_data%input%cell(indx,:)
    Call normal_along_vector(r,normal)
    
    ! Find the topmost with the largest component along the normal vector to the surface 
    max_dist=-Huge(1.0_wp)
    Do i = 1, model_data%input%num_atoms
      flag=.True.
      j=1 
      Do While (j <= model_data%num_species%value .And. flag)
         If (Trim(model_data%species_info(j)%topology)=='electrode') Then
           k=1
           Do While (k <= model_data%species_info(j)%num_components .And. flag)
             If (Trim(model_data%input%atom(i)%tag) == Trim(model_data%species_info(j)%component%tag(k))) Then
               r(:)=model_data%input%atom(i)%r(:)
               dist=Dot_product(r,normal)
               If (dist > max_dist) Then
                 max_dist = dist
                 topmost=r
                 tag_top=model_data%input%atom(i)%tag
               End If
               flag=.False.
             End If
             k=k+1
           End Do
         End If
        j=j+1
      End Do    
    End Do
    
    ! Find the topmost with the largest component along the normal vector to the surface 
    min_dist=Huge(1.0_wp)
    Do i = 1, model_data%input%num_atoms
      flag=.True.
      j=1 
      Do While (j <= model_data%num_species%value .And. flag)
         If (Trim(model_data%species_info(j)%topology)=='electrode') Then
           k=1
           Do While (k <= model_data%species_info(j)%num_components .And. flag)
             If (Trim(model_data%input%atom(i)%tag) == Trim(model_data%species_info(j)%component%tag(k))) Then
               r(:)=model_data%input%atom(i)%r(:)
               dist=Dot_product(r,normal)
               If (dist < min_dist) Then
                 min_dist = dist
                 lowest=r
                 tag_low=model_data%input%atom(i)%tag
               End If
               flag=.False.
             End If
             k=k+1
           End Do
         End If
        j=j+1
      End Do    
    End Do    
    
    model_data%sample%surface_shift=topmost-lowest
    model_data%sample%surface_shift(indx)=0.0_wp
    
    If (Trim(tag_top) /= Trim(tag_low)) Then
         Write (messages(1),'(1x,a)') ' '
         Write (messages(2),'(1x,a)') '*********************************************************************************'
         Write (messages(3),'(1x,a)') '*** WARNING: the user should corroborate if both sides of the electrode      ***'
         Write (messages(4),'(1x,a)') '*** model are equivalent. Adding species at both sides might not be correct, ***'
         Write (messages(5),'(1x,a)') '*** in which case the "both_surfaces" directive should be set to .False.     ***'
         Write (messages(6),'(1x,a)') '********************************************************************************'
         Call info(messages, 6)      
    End If
    
  End Subroutine surface_shift
  
  Subroutine check_add_species_from(model_data, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check the value set for the deposition level againts the cell dimensions
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    Type(model_type),    Intent(In   ) :: model_data
    Type(simul_type),    Intent(In   ) :: simulation_data

    Integer(Kind=wi) :: indx
    Real(Kind=wp)    :: ltop
    Real(Kind=wp)    :: r(3), normal(3), level(3)
    Character(Len=256) :: message
    Character(Len=256) :: messages(2)
    
    If (Trim(model_data%normal_vector%type)=='c3') Then
      indx=3
    ElseIf (Trim(model_data%normal_vector%type)=='c2') Then
      indx=2
    ElseIf (Trim(model_data%normal_vector%type)=='c1') Then
      indx=1
    End If

    r(:)=model_data%input%cell(indx,:)
    Call normal_along_vector(r,normal)

    level=0.0_wp
    level(indx)=model_data%add_species_from%value
    r(:)=model_data%input%cell(indx,:)
    If (Abs(model_data%add_species_from%value) >= norm2(r) .Or. Dot_product(level,normal) < 0.0_wp ) Then
       If (model_data%add_species_from%fread) Then
         Write (message, '(1x,a)') '***ERROR: innapropriate value of "add_species_from" for the input system.&
                                & Please check the input model structure and simulation box.'
       Else
         If (model_data%centre_electrode%stat) Then
           Write (message, '(1x,a)') '***ERROR: Wrong input structure and/or simulation box. Please check'
         Else
           Write (message, '(1x,a)') '***ERROR: Wrong input structure and/or simulation box. The user&
                                    & should set the "centre_electrode" directive to True'
         End If
       End If
       Call error_stop(message)      
    End If 
  
    If (simulation_data%dft%gc%activate%stat) Then
      If (simulation_data%electrolyte%info_pcc%fread) Then
        Write (messages(1), '(1x,a)') '***ERROR: inconsistency in the setting of the planar counter charge for GC-DFT.'
        ltop=model_data%input%cell_length(indx)-simulation_data%electrolyte%dist_edge%value
        If (model_data%add_species_from%value > ltop) Then 
          If (model_data%add_species_from%fread) Then
            Write (messages(2), '(1x,a)') 'Either the input value of "add_species_from" is inappropriate&
                                   & or the cell size perpendicular to the surface is not large enough.'
          
          Else
            Write (messages(2), '(1x,a)') 'Either the plane position is too close to the slab (large value for "distance_to_edge")&
                                   & or the cell size perpendicular to the surface is not large enough.'
          End If
          Call info(messages, 2)
          Call error_stop(' ')
        End If
        ltop=ltop-simulation_data%electrolyte%gaussian_width%value
        If (model_data%add_species_from%value > ltop) Then 
          If (model_data%add_species_from%fread) Then
            Write (messages(2), '(1x,a)') 'Either the input value of "add_species_from" is inappropriate&
                                   & or the cell size perpendicular to the surface is not large enough.'
          
          Else
            Write (messages(2), '(1x,a)') 'Either the width of the gaussian charge distributiion is wrong&
                                          & (unphysical "gaussian_width") or the cell size perpendicular&
                                          & to the surface is not large enough.'
          End If
          Call info(messages, 2)
          Call error_stop(' ')
        End If
        
      End If
    End If
  
  End Subroutine check_add_species_from

  
  Subroutine optimise_perp_cell_size(model_data, simulation_data) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Optimise perpendicular size of the simulation cell 
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),    Intent(InOut) :: model_data
    Type(simul_type),    Intent(In   ) :: simulation_data 

    Integer(Kind=wi) :: i, indx
    Real(Kind=wp)    :: min_dist, max_dist, lpcc
    Real(Kind=wp)    :: dist, fact
    Real(Kind=wp)    :: r(3), normal(3)

    min_dist=Huge(1.0_wp)
    normal=0.0_wp

    If (Trim(model_data%normal_vector%type)=='c3') Then
      indx=3
    ElseIf (Trim(model_data%normal_vector%type)=='c2') Then
      indx=2
    ElseIf (Trim(model_data%normal_vector%type)=='c1') Then
      indx=1
    End If

    r(:)=model_data%input%cell(indx,:)
    Call normal_along_vector(r,normal)
    ! Find the atom with the largest component along the normal vector to the surface 
    max_dist=-Huge(1.0_wp)
    Do i = 1, model_data%input%num_atoms_extra
      r(:)=model_data%input%atom(i)%r(:)
      dist=Dot_product(r,normal)
      If (dist > max_dist) Then
        max_dist = dist
      End If
    End Do

    If (simulation_data%electrolyte%info_pcc%stat) Then
      lpcc=simulation_data%electrolyte%dist_edge%value+simulation_data%electrolyte%gaussian_width%value/2.0_wp
      fact=(max_dist+model_data%add_extra_space%value+lpcc)/model_data%input%cell_length(indx) 
    Else
      fact=(max_dist+model_data%add_extra_space%value)/model_data%input%cell_length(indx)
    End If
    
    If (Abs(fact-1.0_wp)>length_tol) then
      model_data%input%size_changed=.True.
    Else
     model_data%input%size_changed=.False.
    End If

    model_data%input%cell(indx,:)=(2.0_wp*fact-1.0_wp)*model_data%input%cell(indx,:)
    Call about_cell(model_data%input%cell,model_data%input%invcell,model_data%input%cell_length)
    
    Call centre_electrode(model_data, 'shift_optimised')
  
  End Subroutine optimise_perp_cell_size
 
  Subroutine check_pcc_vs_generated_model(model_data, simulation_data) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Cheks 
    !
    ! author    - i.scivetti Jan 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),    Intent(In   ) :: model_data
    Type(simul_type),    Intent(In   ) :: simulation_data 

    Integer(Kind=wi) :: i, indx
    Real(Kind=wp)    :: min_dist, max_dist
    Real(Kind=wp)    :: dist, ltop
    Real(Kind=wp)    :: r(3), normal(3)
    Character(Len=256) :: messages(4)

    min_dist=Huge(1.0_wp)
    normal=0.0_wp

    If (Trim(model_data%normal_vector%type)=='c3') Then
      indx=3
    ElseIf (Trim(model_data%normal_vector%type)=='c2') Then
      indx=2
    ElseIf (Trim(model_data%normal_vector%type)=='c1') Then
      indx=1
    End If

    r(:)=model_data%input%cell(indx,:)
    Call normal_along_vector(r,normal)
    ! Find the atom with the largest component along the normal vector to the surface 
    max_dist=-Huge(1.0_wp)
    Do i = 1, model_data%input%num_atoms_extra
      r(:)=model_data%input%atom(i)%r(:)
      dist=Dot_product(r,normal)
      If (dist > max_dist) Then
        max_dist = dist
      End If
    End Do
    
    Write (messages(1), '(1x,a)') '***ERROR: problems with setting the planar counter charge for GC-DFT.'
    ltop=model_data%input%cell_length(indx)-simulation_data%electrolyte%dist_edge%value
    ltop=ltop-simulation_data%electrolyte%gaussian_width%value/2.0_wp
    If (max_dist > ltop) Then 
      Write (messages(2), '(1x,a)') 'The planar charge distribution for the electrolyte overlaps with atomic species of&
                                    & in the generated model.'
      Write (messages(3), '(1x,a)') 'The user should set sufficient room to prevent this to happen. Review the&
                                    & "&planar_counter_charge"'
      Write (messages(4), '(1x,a)') 'directives, set a larger input cell size (perp. to the surface) and/or add less&
                                    & species in the model.'
      Call info(messages, 4)
      If (.Not. model_data%optimise_size%stat) Then
        Call info(' ', 1)
        Write (messages(1), '(1x,a)') '***ADVISE*** If  the above does not work, set the "optimise_size"&
                                     & directive to .True.'
        Call info(messages, 1)
      End If
      Call error_stop(' ')
    End If
    
   End  Subroutine check_pcc_vs_generated_model

  Subroutine centre_electrode(model_data, option)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to shift the set of atomic coordinates to the centre of the cell
    !
    ! author    - i.scivetti March 2026
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),   Intent(InOut) :: model_data
    Character(Len=*),   Intent(In   ) :: option 

    Integer(Kind=wi) :: i, j, k, indx, nat_elect, net_nat
    Real(Kind=wp)    :: shift(3), suma(3), centre_cell(3)
    Logical          :: flag

    
    If (Trim(model_data%normal_vector%type)=='c3') Then
      indx=3
    ElseIf (Trim(model_data%normal_vector%type)=='c2') Then
      indx=2
    ElseIf (Trim(model_data%normal_vector%type)=='c1') Then
      indx=1
    End If

    centre_cell=0.0_wp
    Do i = 1, 3
      centre_cell(:)=centre_cell(:)+0.5*model_data%input%cell(i,:)    
    End Do
    
    nat_elect=0
    suma=0.0_wp
    
    If (Trim(option)=='shift_initial') Then
      net_nat=model_data%input%num_atoms
    Else If (Trim(option)=='shift_optimised') Then
      net_nat=model_data%input%num_atoms_extra
    End if

    Do i = 1, net_nat
      flag=.True.
      j=1 
      Do While (j <= model_data%num_species%value .And. flag)
         If (Trim(model_data%species_info(j)%topology)=='electrode') Then
           k=1
           Do While (k <= model_data%species_info(j)%num_components .And. flag)
             If (Trim(model_data%input%atom(i)%tag) == Trim(model_data%species_info(j)%component%tag(k))) Then
               suma(:)=suma(:)+model_data%input%atom(i)%r(:)
               nat_elect=nat_elect+1
               flag=.False.
             End If
             k=k+1
           End Do
         End If
        j=j+1
      End Do
    End Do
 
    shift=centre_cell-suma/nat_elect

    ! Centre model
    Do i = 1, net_nat
      model_data%input%atom(i)%r(:)=model_data%input%atom(i)%r(:)+shift(:)
    End Do
    
    ! Shift add_species_from
    If (Trim(option)=='shift_initial') Then    
      model_data%add_species_from%value=model_data%add_species_from%value+shift(indx)
    End If
    
  End Subroutine centre_electrode

   
   
End Module atomistic_tools
