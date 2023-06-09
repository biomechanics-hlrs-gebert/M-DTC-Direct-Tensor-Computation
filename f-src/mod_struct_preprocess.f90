!------------------------------------------------------------------------------
!>  Program for generating a Hexaedra mesh from a 3D scalar field
!------------------------------------------------------------------------------
!>  \section written Written by:
!>  Ralf Schneider
!>
!>  \section modified Last modified:
!>  by: Johannes Gebert \n
!>  on : 20.03.2022
!------------------------------------------------------------------------------
Module gen_geometry

Use vtkio             
Use gen_quadmesh
Use vtk_meta_data
Use ser_binary
Use write_deck

Implicit None

Contains

Subroutine generate_geometry(root, ddc_nn, job_dir, typeraw, glob_success)

Type(tBranch), Intent(InOut) :: root
Character(LEN=*), Intent(in) :: job_dir
integer(ik), Intent(in) :: ddc_nn
Logical, Intent(Out) :: glob_success
    
! Chain Variables
Character(Len=*), Parameter :: link_name = 'gen_quadmesh'

! Puredat Variables
Type(tBranch) :: phi_desc

! Decomp Variables 
Type(tBranch), Pointer :: loc_ddc, ddc, bounds_b

! Branch pointers
Type(tBranch), Pointer :: meta_para, domain_branch

    ! Mesh Variables 
    Real(rk)    , Dimension(:,:) , Allocatable :: nodes, displ
    Integer(ik) , Dimension(:)   , Allocatable :: elem_col,node_col, cref
    Integer(ik) , Dimension(:,:) , Allocatable :: elems
    Character(mcl) :: elt_micro, desc, filename, typeraw
    Character(scl) :: restart
    Integer(ik)   :: elo_macro,alloc_stat
    
    ! Parted Mesh -
    Type(tBranch), Pointer :: PMesh
    INTEGER(ik):: parts, multiplier, nts=0

    Integer(1), Dimension(:,:,:), Allocatable :: Phi_ik1
    Integer(2), Dimension(:,:,:), Allocatable :: Phi_ik2
    Integer(4), Dimension(:,:,:), Allocatable :: Phi_ik4

    Character(*), Parameter :: inpsep = "('#',79('='))"

    Character, Dimension(:), Allocatable :: char_arr
    Real(rk) , Dimension(:), Allocatable :: delta
    Integer(ik), Dimension(:), Allocatable :: bpoints,x_D,nn_D, vdim
    Integer(ik), Dimension(3) :: xa_n, xe_n, xa_n_ext, xe_n_ext
    Integer(ik) :: nn_1,nn_2,nn_3, ii,jj, pos_f, first_bytes, llimit, no_nodes=0, no_elems=0

    Logical :: success, fex
    
    Character(len=9) :: nn_char

    write(nn_char,'(I0)')ddc_nn
    glob_success = .TRUE.

    !----------------------------------------------------------------------------
    ! Get global DDC parameters from root
    !----------------------------------------------------------------------------
    Call Search_branch("Global domain decomposition", root, ddc, success)

    call pd_get(ddc,"nn_D",nn_D)
    call pd_get(ddc,"bpoints", bpoints)
    call pd_get(ddc,"x_D",x_D)

    !----------------------------------------------------------------------------
    ! Calculate special parameters of domain decomposition
    !----------------------------------------------------------------------------
    nn_3 = INT( ddc_nn / ( nn_D(1)*nn_D(2) ))
    nn_2 = INT(( ddc_nn - nn_D(1)*nn_D(2)*nn_3 ) / ( nn_D(1) ))
    nn_1 = ( ddc_nn - nn_D(1)*nn_D(2)*nn_3 - nn_D(1)*nn_2 )

    xa_n = [x_D(1) * nn_1 + bpoints(1) + 1, &
            x_D(2) * nn_2 + bpoints(2) + 1, &
            x_D(3) * nn_3 + bpoints(3) + 1  ]

    xe_n = [x_D(1) * (nn_1 + 1) + bpoints(1), &
            x_D(2) * (nn_2 + 1) + bpoints(2), & 
            x_D(3) * (nn_3 + 1) + bpoints(3)  ] 

    xa_n_ext = xa_n - bpoints
    xe_n_ext = xe_n + bpoints

    !----------------------------------------------------------------------------
    ! Add domain branch to root
    !----------------------------------------------------------------------------
    desc=''
    Write(desc,'(A,I0)')"Domain ", ddc_nn
    call add_branch_to_branch(root, domain_branch)
    
    call raise_branch(trim(desc), 0_pd_ik, 0_pd_ik, domain_branch)

    !----------------------------------------------------------------------------
    ! Add branch for local ddc params to domain branch
    !----------------------------------------------------------------------------
    desc=''
    Write(desc,'(A,I0)') "Local domain Decomposition of domain no ", ddc_nn

    call add_branch_to_branch(domain_branch, loc_ddc)
    call raise_branch(trim(desc), 0_pd_ik, 0_pd_ik, loc_ddc)

    call add_leaf_to_branch(loc_ddc, "nn"      , 1_pd_ik, [ddc_nn])
    call add_leaf_to_branch(loc_ddc, "nn_i"    , 3_pd_ik, [nn_1, nn_2, nn_3])
    call add_leaf_to_branch(loc_ddc, "xa_n"    , 3_pd_ik, xa_n    )
    call add_leaf_to_branch(loc_ddc, "xe_n"    , 3_pd_ik, xe_n    )
    call add_leaf_to_branch(loc_ddc, "xa_n_ext", 3_pd_ik, xa_n_ext)
    call add_leaf_to_branch(loc_ddc, "xe_n_ext", 3_pd_ik, xe_n_ext)
    
    !----------------------------------------------------------------------------
    Select Case (timer_level)
      Case (3)
         call start_timer("  +-- Initialisation of Phi "//trim(nn_char))
      Case (2)
         call start_timer("  +-- Initialisation of Phi "//trim(nn_char))
      Case default
         continue
    End Select

    !----------------------------------------------------------------------------
    ! Load branch
    !----------------------------------------------------------------------------
    Call Search_branch("Input parameters", root, meta_para, success)
    Call pd_get(meta_para,"muCT puredat pro_path",char_arr)

    pro_path = char_to_str(char_arr)
    deallocate(char_arr)
    
    Call pd_get(meta_para,"muCT puredat pro_name",char_arr)

    pro_name = char_to_str(char_arr)
    deallocate(char_arr)

    call pd_get(meta_para, "Grid spacings", delta)  
    Call pd_get(meta_para,"Number of voxels per direction",vdim)

    if ( out_amount /= "PRODUCTION" ) then
       write(un_lf,FMT_MSG_AxF0) "Grid spacings", delta
       write(un_lf,FMT_MSG_AxI0) "Number of voxels per direction", vdim
    End if
    
    !----------------------------------------------------------------------------
    ! Open raw data streams
    !----------------------------------------------------------------------------
    phi_desc = read_tree()

    !----------------------------------------------------------------------------
    ! Read PHI (scalar binary values) from file
    !----------------------------------------------------------------------------
    SELECT CASE(TRIM(ADJUSTL(typeraw)))
    CASE('ik1')
        multiplier = 1_ik
        allocate(Phi_ik1(xa_n(1):xe_n(1), xa_n(2):xe_n(2), xa_n(3):xe_n(3)), stat=alloc_stat)
    CASE('ik2')
        multiplier = 2_ik
        allocate(Phi_ik2(xa_n(1):xe_n(1), xa_n(2):xe_n(2), xa_n(3):xe_n(3)), stat=alloc_stat)
    CASE('ik4')
        multiplier = 4_ik
        allocate(Phi_ik4(xa_n(1):xe_n(1), xa_n(2):xe_n(2), xa_n(3):xe_n(3)), stat=alloc_stat)
    END SELECT

    call alloc_err("phi", alloc_stat)

    DO ii=1, phi_desc%no_leaves
        IF (TRIM(phi_desc%leaves(ii)%desc) == "Scalar data") THEN
            first_bytes = INT(phi_desc%leaves(ii)%lbound, 8)
            EXIT
        END IF 
    END DO 
    
    !----------------------------------------------------------------------------
    ! Retrieve basic geometric information
    !----------------------------------------------------------------------------
    call pd_get(meta_para, "Grid spacings", delta)  
    Call pd_get(meta_para, "Number of voxels per direction",vdim)

    if ( out_amount /= "PRODUCTION" ) then
       write(un_lf,FMT_MSG_AxF0) "Grid spacings", delta
       write(un_lf,FMT_MSG_AxI0) "Number of voxels per direction", vdim
    End if

    !----------------------------------------------------------------------------
    ! Set datatype
    ! intent(inout) from open_stream_files did not copy the units(!)
    !----------------------------------------------------------------------------
    SELECT CASE(TRIM(ADJUSTL(typeraw)))
        CASE('ik1'); nts=1
        CASE('ik2'); nts=2
        CASE('ik4'); nts=3
        CASE DEFAULT; CALL print_err_stop(std_out, "Data type failed.", 1)
    END SELECT
    !
    phi_desc%streams%units(nts) = pd_give_new_unit()
    Open(unit=phi_desc%streams%units(nts), &
        file=phi_desc%streams%stream_files(nts), status="old", &
        action="read", access='stream', form='unformatted', &
        position="REWIND")
    phi_desc%streams%ifopen(nts) = .TRUE.

    Do jj = xa_n(3), xe_n(3)
       Do ii = xa_n(2), xe_n(2)

            pos_f = INT((INT(vdim(1)*vdim(2)*(jj         - 1_8), 8) + &
                         INT(vdim(1)*        (ii         - 1_8), 8) + &
                         INT(                (xa_n(1)    - 1_8), 8) + &
                         first_bytes - 1_8) * multiplier + 1_8 , 8)

            SELECT CASE(nts)
                CASE(1); Read(phi_desc%streams%units(1), pos=pos_f) Phi_ik1(:,ii,jj)
                CASE(2); Read(phi_desc%streams%units(2), pos=pos_f) Phi_ik2(:,ii,jj)
                CASE(3); Read(phi_desc%streams%units(3), pos=pos_f) Phi_ik4(:,ii,jj)
            END SELECT
        End Do
    End Do

    ! call close_stream_files(phi_desc)
    SELECT CASE(nts)
        CASE(1); close(phi_desc%streams%units(1))
        CASE(2); close(phi_desc%streams%units(2))
        CASE(3); close(phi_desc%streams%units(3))
    END SELECT


    Select Case (timer_level)
    Case (3)
       call end_timer("  +-- Initialisation of Phi "//trim(nn_char))
    Case (2)
       call end_timer("  +-- Initialisation of Phi "//trim(nn_char))
    Case default
       continue
    End Select

    !============================================================================
    ! Generate Quadmesh from Phi
    !============================================================================
    Call pd_get(meta_para, 'Lower limit of iso value', llimit)
    Call pd_get(meta_para, 'Element type  on micro scale', char_arr)
    elt_micro = char_to_str(char_arr)
    deallocate(char_arr)

    !------------------------------------------------------------------------------
    ! Get global DDC parameters from root
    !------------------------------------------------------------------------------
    Call Search_branch("Global domain decomposition", root, ddc, success)

    SELECT CASE(TRIM(ADJUSTL(typeraw)))
        CASE('ik1')
            call gen_quadmesh_from_phi(delta, ddc, loc_ddc, llimit, elt_micro, &
             nodes, elems, node_col, elem_col, no_nodes, no_elems, typeraw, phi_ik1=Phi_ik1)

        CASE('ik2')
            call gen_quadmesh_from_phi(delta, ddc, loc_ddc, llimit, elt_micro, &
             nodes, elems, node_col, elem_col, no_nodes, no_elems, typeraw, phi_ik2=Phi_ik2)

        CASE('ik4')
            call gen_quadmesh_from_phi(delta, ddc, loc_ddc, llimit, elt_micro, &
             nodes, elems, node_col, elem_col, no_nodes, no_elems, typeraw, phi_ik4=Phi_ik4)

    END SELECT

    If (out_amount == "DEBUG") THEN
       Write(un_lf, fmt_dbg_sep)
       Write(un_lf, fmt_MSG_xAI0)"Root pointer after exec_single_domain"
       Call log_tree(root,un_lf,.True.)
       Write(un_lf, fmt_dbg_sep)
    END If

    !------------------------------------------------------------------------------
    ! Store PHI as vtk structured points
    !------------------------------------------------------------------------------
    if (out_amount == "DEBUG") then

       filename=''
       write(filename,'(A,I0,A)')trim(job_dir)//trim(project_name)//"_",ddc_nn,"_phi.vtk"
       !------------------------------------------------------------------------------
       ! Check existance of vtk file 
       !------------------------------------------------------------------------------
       Call pd_get(meta_para, 'Restart', char_arr)
       restart = char_to_str(char_arr)
       deallocate(char_arr)


       INQUIRE(FILE=TRIM(filename), EXIST=fex)

        IF((fex) .AND. ((restart=="Y") .OR. (restart=="YES"))) THEN
            CALL execute_command_line("rm -r "//TRIM(filename))
        ELSE IF ((fex).AND. ((restart=="N") .OR. (restart=="NO"))) THEN
            CALL print_err_stop(std_out, TRIM(filename)//" already exists and &
                &restart -> No", 1)
        END IF 

       !------------------------------------------------------------------------------
       ! Write structured points 
       !------------------------------------------------------------------------------
        SELECT CASE(TRIM(ADJUSTL(typeraw)))
            CASE('ik1')
                CALL write_ser_vtk(filename, typeraw, delta, x_D, &
                    (Real(xa_n,rk)-[0.5_rk,0.5_rk,0.5_rk])*delta, Phi_ik1)
            CASE('ik2')
                CALL write_ser_vtk(filename, typeraw, delta, x_D, &
                    (Real(xa_n,rk)-[0.5_rk,0.5_rk,0.5_rk])*delta, Phi_ik2)
            CASE('ik4')
                CALL write_ser_vtk(filename, typeraw, delta, x_D, &
                    (Real(xa_n,rk)-[0.5_rk,0.5_rk,0.5_rk])*delta, Phi_ik4)
        END SELECT


    !     call write_vtk_structured_points(filename = trim(filename),  &
    !         extend = x_D, spacing = delta, &
    !         origin = (Real(xa_n,rk)-[0.5_rk,0.5_rk,0.5_rk])*delta)
    !    !------------------------------------------------------------------------------
    !    ! Write data to vtk file 
    !    !------------------------------------------------------------------------------
    !    call write_vtk_data_int4_scalar_1D(&
    !         matrix = reshape(phi,[x_D(1)*x_D(2)*x_D(3)]), &
    !         filename = trim(filename),  &
    !         desc = "PHI", head = .TRUE.)
    End if

    SELECT CASE(nts)
        CASE(1); Deallocate(Phi_ik1)
        CASE(2); Deallocate(Phi_ik2)
        CASE(3); Deallocate(Phi_ik4)
    END SELECT


    !------------------------------------------------------------------------------
    ! Break if no structure is generated
    !------------------------------------------------------------------------------
    If ((no_nodes == 0) .or. (no_elems == 0)) then
       glob_success = .FALSE.
       Goto 1000
    End If

    !------------------------------------------------------------------------------
    ! Store Quadmesh as vtk unstructured grid and ABAQUS input deck
    !------------------------------------------------------------------------------
    if (out_amount == "DEBUG") then

       filename=''
       write(filename,'(A,I0,A)')trim(job_dir)//trim(project_name)//"_",ddc_nn,"_usg.vtk"
       
       INQUIRE(FILE=TRIM(filename), EXIST=fex)

        IF((fex) .AND. ((restart=="Y") .OR. (restart=="YES"))) THEN
            CALL execute_command_line("rm -r "//TRIM(filename))
            
        ELSE IF ((fex).AND. ((restart=="N") .OR. (restart=="NO"))) THEN
            CALL print_err_stop(std_out, TRIM(filename)//" already exists and &
                &restart -> No", 1)
        END IF 


       if (elt_micro == "HEX08") then
          
          Call write_vtk_unstructured_grid(trim(filename), &
               nodes, elems(1:8, 1:no_elems))
          
       else if (elt_micro == "HEX20") then
       
          Call write_vtk_unstructured_grid(trim(filename), &
               nodes, &
               elems([1,3,5,7, 13,15,17,19, 2,4,6,8, 14,16,18,20, 9,10,11,12],1:no_elems))
       end if

    End if
    
    call pd_get(meta_para,"No of mesh parts per subdomain", parts)
    Call pd_get(meta_para,'Element order on macro scale', elo_macro)
    
    !------------------------------------------------------------------------------
    ! Perform domain decomposition by metis
    !------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------
    ! Add branch for mesh to domain branch
    !------------------------------------------------------------------------------
    desc=''
    Write(desc,'(A,I0)')'Mesh info of '//trim(project_name)//'_',ddc_nn
   
    call add_branch_to_branch(domain_branch, PMesh)
    call raise_branch(trim(desc), parts, 0_pd_ik, PMesh)

    Call part_mesh(nodes, elems, no_nodes, no_elems, parts, PMesh, job_dir, ddc_nn)
   
    call add_leaf_to_branch(PMesh, "No of nodes in mesh",  1_pd_ik, [no_nodes])
    call add_leaf_to_branch(PMesh, "No of elements in mesh",  1_pd_ik, [no_elems])
    
    !------------------------------------------------------------------------------
    ! Retrive valid branch pointers to local and global ddc
    !------------------------------------------------------------------------------
    Call Search_branch("Global domain decomposition", root, ddc, success)
    desc=''
    Write(desc,'(A,I0)')"Local domain Decomposition of domain no ",ddc_nn
    Call Search_branch(trim(desc), root, loc_ddc, success)
    
    Call generate_boundaries(PMesh, ddc, loc_ddc, elo_macro)

    !------------------------------------------------------------------------------
    ! Store boundary displacements of LC1 as point data in parted vtk meshes
    !------------------------------------------------------------------------------
    if (out_amount == "DEBUG") then
       
       Do ii = 1, parts

          Call Search_branch("Boundaries_"//trim(nn_char)//"_1",&
             PMesh%branches(ii),bounds_b, success)

          allocate(cref( &
               minval(PMesh%branches(ii)%leaves(3)%p_int8(1:)):&
               maxval(PMesh%branches(ii)%leaves(3)%p_int8(1:))))
          cref = 0_ik
          
          do jj = 1, PMesh%branches(ii)%leaves(3)%dat_no
             cref(PMesh%branches(ii)%leaves(3)%p_int8(jj))=jj
          End do

          allocate(displ(3,PMesh%branches(ii)%leaves(3)%dat_no))
          displ=0._rk
          
          Do jj = 1, bounds_b%leaves(1)%dat_no
             displ(:,cref(bounds_b%leaves(1)%p_int8(jj))) = &
                  bounds_b%leaves(2)%p_real8((jj-1)*3+1:(jj-1)*3+3)
          End Do
          
        !   filename = ""
        !   Write(filename,'(A,A,I0,A,I0,A)')trim(job_dir),"Part-",ii,"_",ddc_nn,".vtk"

        !   INQUIRE(FILE=TRIM(filename), EXIST=fex)

        !   Call write_vtk_data_real8_vector_1D ( &
        !        displ, Trim(filename), "BoundDispl", .FALSE., "POINT_DATA")
          
          deallocate(cref,displ)
          
       End Do
       
    End if
    
    call add_leaf_to_branch(PMesh, "Coordinates", no_nodes*3, reshape(nodes,[3*no_nodes]))
    
    !======================================================
    !== Break if less than 6 DOFs are constrained =========
    If (sum(PMesh%leaves(3)%p_int8) < 6) then
       glob_success = .FALSE.
    End If
    
    1000 continue 

  End Subroutine generate_geometry

  !**************************************************************************
  !>  Program for boundary application                                       **
  !                                                                        **
  ! ---------- **
  !>  \section written Written by:
  !>  Ralf Schneider
  !>
  !>  \section modified Last modified:
  !>  by: Ralf Schneider \n
  !>  on : 16.09.2015
  ! ------------ *
  Subroutine generate_boundaries(PMesh, ddc, loc_ddc, elo_macro)

    ! Parameters 
    Type(tBranch), Intent(InOut) :: PMesh
    Type(tBranch), Intent(In) :: ddc, loc_ddc
    Integer, Intent(In) :: elo_macro

    Type(tBranch), Pointer :: part_b, bounds_b
    
    Real(rk), Dimension(:,:), Allocatable :: bnode_vals
    Real(rk), Dimension(:), Allocatable :: dim_c, delta
    Real(rk), Dimension(3) :: min_c, max_c, coor

    Integer(ik) :: parts, ddc_nn, no_nodes_macro, ii, jj, no_bnodes, b_items, nnodes
    Integer(ik), Dimension(:), Allocatable :: &
        no_nodes_all, no_elems_all, no_cdofs_all, xa_n, xe_n, bnode_ids

    !------------------------------------------------------------------------------
    ! Get Parameters of domain decomposition
    !------------------------------------------------------------------------------
    call pd_get(loc_ddc, "nn",      ddc_nn)
    call pd_get(ddc,     "x_D_phy", dim_c )
    call pd_get(loc_ddc, "xa_n",    xa_n  )
    call pd_get(loc_ddc, "xe_n",    xe_n  )
    call pd_get(ddc,     "delta",   delta )

    min_c = Real(xa_n - 1,rk) * delta
    max_c = Real(xe_n    ,rk) * delta

    if ( out_amount /= "PRODUCTION" ) then
       Write(un_lf,FMT_MSG_AxF0)'Minimum Coordinates',min_c
       Write(un_lf,FMT_MSG_AxF0)'Maximum Coordinates',max_c
       Write(un_lf,*)
       
       Write(un_lf,FMT_MSG_AxF0)'Cube dimensions',dim_c
       Write(un_lf,*)
    End if
    
    parts = PMesh%no_branches

    IF (elo_macro == 1) THEN
        no_nodes_macro = 8
    ELSE IF (elo_macro == 2) THEN
        no_nodes_macro = 20
    END IF 

    Allocate(no_nodes_all(parts))
    Allocate(no_elems_all(parts))
    Allocate(no_cdofs_all(parts))
    
    !------------------------------------------------------------------------------
    ! For all Parts
    !------------------------------------------------------------------------------
    Do jj = 1, parts

        part_b => PMesh%branches(jj)
        nnodes =  part_b%leaves(1)%dat_no

        no_nodes_all(jj) = nnodes
        no_elems_all(jj) = part_b%leaves(4)%dat_no
       
        !------------------------------------------------------------------------------
        ! Determine number of boundary nodes
        !------------------------------------------------------------------------------
        no_bnodes = 0

        if ( out_amount /= "PRODUCTION" ) then
            Write(un_lf,FMT_MSG_xAI0)'Number of nodes in Part ',jj," during preprocessing: ",nnodes
            Write(un_lf,FMT_MSG_xAI0)'Size of p_real8 in Part ',jj,": ",size(part_b%leaves(2)%p_real8)
        End if

        !------------------------------------------------------------------------------
        ! For all nodes in part
        !------------------------------------------------------------------------------
        Do ii = 1, nnodes

            coor = part_b%leaves(2)%p_real8((ii-1_ik)*3_ik+1_ik:ii*3_ik)
            
            ! Nodes in facet 1
            if ( (coor(3) - min_c(3)) <= (dim_c(3) * delta_b) ) then
                no_bnodes = no_bnodes + 1
                
            ! Nodes in facet 6
            Else if ( (max_c(3) - coor(3)) <= (dim_c(3) * delta_b) ) then
                no_bnodes = no_bnodes + 1
                
            ! Nodes in facet 2
            Else if ( (coor(2) - min_c(2)) <= (dim_c(2) * delta_b) ) then
                no_bnodes = no_bnodes + 1
                
            ! Nodes in facet 4
            Else if ( (max_c(2) - coor(2)) <= (dim_c(2) * delta_b) ) then
                no_bnodes = no_bnodes + 1
                
            ! Nodes in facet 5
            Else if ( (coor(1) - min_c(1)) <= (dim_c(1) * delta_b) ) then
                no_bnodes = no_bnodes + 1
                
            ! Nodes in facet 3
            Else if ( (max_c(1) - coor(1)) <= (dim_c(1) * delta_b) ) then
                no_bnodes = no_bnodes + 1
                
            End if

        End Do

        if ( out_amount /= "PRODUCTION" ) then
            Write(un_lf,FMT_MSG_xAI0)'Number of constrained nodes in Part ',jj," : ",no_bnodes

            IF (no_bnodes == 0) THEN
                WRITE(un_lf,FMT_MSG_AI0AxF0)'get_part_coords: min_c    of Part ',jj," : ",min_c
                WRITE(un_lf,FMT_MSG_AI0AxF0)'get_part_coords: max_c    of Part ',jj," : ",max_c
                WRITE(un_lf,FMT_MSG_AI0AxF0)'get_part_coords: min coor of Part ',jj," : ",minval(part_b%leaves(2)%p_real8)
                WRITE(un_lf,FMT_MSG_AI0AxF0)'get_part_coords: max coor of Part ',jj," : ",maxval(part_b%leaves(2)%p_real8)
            END IF 
        End if
        
        Allocate(bnode_ids(no_bnodes))
        Allocate(bnode_vals(no_nodes_macro*3,no_bnodes*3))

        !------------------------------------------------------------------------------
        ! Boundary application     
        !------------------------------------------------------------------------------
        b_items = 0

        Do ii = 1, nnodes

            coor = part_b%leaves(2)%p_real8((ii-1_ik)*3_ik+1_ik:ii*3_ik)
            
            ! Nodes in facet 1
            if ( (coor(3) - min_c(3)) <= (dim_c(3) * delta_b) ) then
                b_items = b_items + 1
                bnode_ids(b_items) = part_b%leaves(3)%p_int8(ii)

            ! Nodes in facet 6
            Else if ( (max_c(3) - coor(3)) <= (dim_c(3) * delta_b) ) then
                b_items = b_items + 1
                bnode_ids(b_items) = part_b%leaves(3)%p_int8(ii)

            ! Nodes in facet 2
            Else if ( (coor(2) - min_c(2)) <= (dim_c(2) * delta_b) ) then
                b_items = b_items + 1
                bnode_ids(b_items) = part_b%leaves(3)%p_int8(ii)

            ! Nodes in facet 4
            Else if ( (max_c(2) - coor(2)) <= (dim_c(2) * delta_b) ) then
                b_items = b_items + 1
                bnode_ids(b_items) = part_b%leaves(3)%p_int8(ii)
        
            ! Nodes in facet 5
            Else if ( (coor(1) - min_c(1)) <= (dim_c(1) * delta_b) ) then
                b_items = b_items + 1
                bnode_ids(b_items) = part_b%leaves(3)%p_int8(ii)
        
            ! Nodes in facet 3
            Else if ( (max_c(1) - coor(1)) <= (dim_c(1) * delta_b) ) then
                b_items = b_items + 1
                bnode_ids(b_items) = part_b%leaves(3)%p_int8(ii)

            End if
        
        End Do

        b_items = 0

        Do ii = 1, nnodes

            coor = part_b%leaves(2)%p_real8((ii-1_ik)*3_ik+1_ik:ii*3_ik)
        
            ! Nodes in facet 1
            if ( (coor(3) - min_c(3)) <= (dim_c(3) * delta_b) ) then
                Call determine_prescribed_displ(coor(:), min_c, max_c, &
                    no_nodes_macro, bnode_vals(:,b_items+1:b_items+3))
                b_items = b_items+3

            ! Nodes in facet 6
            Else if ( (max_c(3) - coor(3)) <= (dim_c(3) * delta_b) ) then
                Call determine_prescribed_displ(coor(:), min_c, max_c, &
                    no_nodes_macro, bnode_vals(:,b_items+1:b_items+3))
                b_items = b_items+3
                
            ! Nodes in facet 2
            Else if ( (coor(2) - min_c(2)) <= (dim_c(2) * delta_b) ) then
                Call determine_prescribed_displ(coor(:), min_c, max_c, &
                    no_nodes_macro, bnode_vals(:,b_items+1:b_items+3))
                b_items = b_items+3
                
            ! Nodes in facet 4
            Else if ( (max_c(2) - coor(2)) <= (dim_c(2) * delta_b) ) then
                Call determine_prescribed_displ(coor(:), min_c, max_c, &
                    no_nodes_macro, bnode_vals(:,b_items+1:b_items+3))
                b_items = b_items+3
                
            ! Nodes in facet 5
            Else if ( (coor(1) - min_c(1)) <= (dim_c(1) * delta_b) ) then
                Call determine_prescribed_displ(coor(:), min_c, max_c, &
                    no_nodes_macro, bnode_vals(:,b_items+1:b_items+3))
                b_items = b_items+3
                
            ! Nodes in facet 3
            Else if ( (max_c(1) - coor(1)) <= (dim_c(1) * delta_b) ) then
                Call determine_prescribed_displ(coor(:), min_c, max_c, &
                    no_nodes_macro, bnode_vals(:,b_items+1:b_items+3))
                b_items = b_items+3
                
            End if
            
        End Do

        Write(un_lf,FMT_MSG_xAI0)'Number of constrained DOF',b_items

        no_cdofs_all(jj) = b_items

        call add_branch_to_branch(part_b,bounds_b)
        call raise_branch("Boundaries", no_nodes_macro*3, 0, bounds_b)

        Do ii = 1, no_nodes_macro*3
            
            Write(bounds_b%branches(ii)%desc,'(A,I0,A,I0)') "Boundaries"//'_',ddc_nn,'_',ii

            ! no_bnodes == 0?! > dat_no = 0 for some specific nodes
            call add_leaf_to_branch(bounds_b%branches(ii), &
                                    "Boundary_Ids", no_bnodes, bnode_ids)
            call add_leaf_to_branch(bounds_b%branches(ii), &
                                    "Boundary_Values", no_bnodes*3,  bnode_vals(ii,:))
        End Do

        DeAllocate(bnode_ids)
        DeAllocate(bnode_vals)

        End Do

        call add_leaf_to_branch(PMesh, "No of nodes in parts", parts, no_nodes_all)
        call add_leaf_to_branch(PMesh, "No of elems in parts", parts, no_elems_all)
        call add_leaf_to_branch(PMesh, "No of cdofs in parts", parts, no_cdofs_all)
        
        ! DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        if (out_amount == "DEBUG") then
        Write(un_lf,fmt_dbg_sep)
        Write(un_lf,'(A)')"PMesh after Boundary application"
        call log_tree(PMesh,un_lf)
        Write(un_lf,fmt_dbg_sep)
        End if
        ! DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
    End Subroutine generate_boundaries
    
End Module gen_geometry

