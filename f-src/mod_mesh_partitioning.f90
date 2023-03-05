module mesh_partitioning

  Use decomp
  Use chain_routines
  Use puredat
  Use vtkio
  USE metis

  use, intrinsic :: iso_c_binding

  Implicit None 
  
  Type T_PMesh

     Integer(kind=ik) :: nnodes, nelems, nouter_nds

     INTEGER(KIND=IK), DIMENSION(:) , Allocatable :: NN, ncolor, cneigh, bnodes
     REAL(KIND=RK)   , DIMENSION(:,:), Allocatable :: COOR
     INTEGER(KIND=IK), DIMENSION(:,:), Allocatable :: EIND,neigh

  End type T_PMesh

contains

  Subroutine part_mesh(nodes, eind, HU_magnitudes, nnodes, ne, parts, PMesh, job_dir, ddc_nn)

    !**************************************************************************
    ! Declarations ************************************************************
   
    ! Parted Mesh *************************************************************
    Type(tBranch)   , intent(Inout) :: PMesh
    Character(LEN=*), Intent(in) :: job_dir
    Integer(kind=ik), intent(in) :: ddc_nn

    ! Metis variables *********************************************************
    Integer(kind=C_INT64_T), Intent(In) :: ne
    Integer(kind=C_INT64_T), intent(In) :: nnodes

    Integer(kind=C_INT64_T)                              :: nn
    Integer(kind=C_INT64_T), Dimension(:),   Allocatable :: eptr
    Integer(kind=C_INT64_T), Dimension(:,:), Allocatable :: eind
    Integer(kind=C_INT64_T), Dimension(:)  , Allocatable :: HU_magnitudes
    Integer(kind=C_INT64_T), Dimension(:),   Allocatable :: vwgt, vsize
    Integer(kind=C_INT64_T)                              :: ncommon
    Integer(kind=C_INT64_T), Intent(in)                  :: parts
    Real   (kind=c_double) , Dimension(:)  , Allocatable :: tpwgts
    Integer(kind=C_INT64_T), Dimension(40)               :: options
    Integer(kind=C_INT64_T)                              :: objval
    Integer(kind=C_INT64_T), Dimension(:)  , Allocatable :: depart
    Integer(kind=C_INT64_T), Dimension(:)  , Allocatable :: dnpart

    Real(kind=C_double), Dimension(:,:), intent(in) :: nodes

    !**************************************************************************
    Integer(kind=ik), Dimension(:), Allocatable :: nnodes_pp, nelems_pp, nouter_nds_pp
    INTEGER(kind=ik)                            :: ii,jj,kk,idum, tmp_i8
    INTEGER(kind=ik), Dimension(:), Allocatable :: cref
    Character(Len=mcl)                          :: desc, filename
    INTEGER(kind=ik), Dimension(1,2)            :: bounds
    Character(len=mcl)                          :: vtk_file
    Integer(kind=ik), Dimension(:,:), Allocatable :: elems
    Integer(kind=ik)                              :: nn_el
    
    !== Code ====================================================================

    nn_el = size(eind(:,1))

    nn = ne * nn_el

    If (out_amount /= "PRODUCTION" ) then
       Write(un_lf, fmt_msg_xAI0) "Number of Nodes          :",nnodes
       Write(un_lf, fmt_msg_xAI0) "Number of Elements       :",ne
       Write(un_lf, fmt_msg_xAI0) "Length of Element descr. :",nn
       Write(un_lf, fmt_msg_xAI0) "Number of Nodes per Elem.:",nn_el
    End If
    
    Do ii = 1, Parts
          
       Write(desc,'(A,I0)')"Part_",ii
       call raise_branch(trim(desc), 1_pd_ik, 6_pd_ik, PMesh%branches(ii))
       Call raise_leaves(6,&
            ["Node Numbers       ", "Coordinates        ", &
             "Global Node Numbers", "Element Numbers    ", &
             "Topology           ", "HU Magnitudes      "],&
            [4_1,     5_1,     4_1,     4_1,     4_1,     4_1    ],&
            [0_pd_ik, 0_pd_ik, 0_pd_ik, 0_pd_ik, 0_pd_ik, 0_pd_ik],&
            [0_pd_ik, 0_pd_ik, 0_pd_ik, 0_pd_ik, 0_pd_ik, 0_pd_ik],&
            [0_pd_ik, 0_pd_ik, 0_pd_ik, 0_pd_ik, 0_pd_ik, 0_pd_ik],&
            PMesh%branches(ii))
       
       call raise_branch("Connections", 0_pd_ik, 2_pd_ik, PMesh%branches(ii)%branches(1))
       
    End Do

    Allocate(nnodes_pp(parts), nelems_pp(parts), nouter_nds_pp(parts))
    nnodes_pp     = 0
    nelems_pp     = 0
    nouter_nds_pp = 0
    
    !************************************************************************
    ! If we should have more than one part **********************************
    !************************************************************************
    If (parts > 1) then
       
       call start_timer("  +-- Prep EPTR")
    
       nn = ne * nn_el
       
       eind = eind-1
       

        !-------------------------------------------------------------------------------
        ! From the METIS manual section 5.6
        ! Creating the array eptr is comparatively simple, because all elements are of 
        ! the same type and therefore of the same interval regarding eptr/eind
        !-------------------------------------------------------------------------------
        ! 5.6 Mesh data structure
        ! All of the mesh partitioning and mesh conversion routines in METIS take as input the element node array of a mesh.
        ! This element node array is stored using a pair of arrays called eptr and eind, which are similar to the xadj and
        ! adjncy arrays used for storing the adjacency structure of a graph. The size of the eptr array is n+ 1, where n is the
        ! number of elements in the mesh. The size of the eind array is of size equal to the sum of the number of nodes in all
        ! the elements of the mesh. The list of nodes belonging to the ith element of the mesh are stored in consecutive locations
        ! of eind starting at position eptr[i] up to (but not including) position eptr[i+1]. This format makes it easy
        ! to specify meshes of any type of elements, including meshes with mixed element types that have different number of
        ! nodes per element. As it was the case with the format of the mesh file described in Section 4.1.2, the ordering of the
        ! nodes in each element is not important.
        !-------------------------------------------------------------------------------
        Allocate(eptr(ne+1))
        Do ii = 1, ne+1
            eptr(ii) = (ii-1) * nn_el
        End Do
       
       call end_timer("  +-- Prep EPTR")
    
       ! EXEC METIS *********************************************************
    
       ! Allocate other METIS parameters ************************
       ! These should be fitted in metis_interface.c !! *********
       Allocate(vwgt(ne), vsize(ne), tpwgts(parts), & 
                depart(ne), dnpart(nn) )

       vwgt    = 0_C_INT64_T
       vsize   = 0_C_INT64_T
       
       ncommon = 4
       
       tpwgts  = 0._c_double
       options = 0_C_INT64_T
       
       objval  =  0_C_INT64_T
       
       depart  =  0_C_INT64_T
       dnpart  =  0_C_INT64_T
       
       call start_timer("  +-- Metis")

       call F_Metis_PartMeshDual( ne, nn, eptr, reshape(eind,[size(eind)]),&
            vwgt, vsize,  ncommon, parts, tpwgts, options, objval, &
            depart, dnpart)

       If (out_amount /= "PRODUCTION" ) then
          Write(un_lf,fmt_msg_xAI0)"objval of dual  partitioning =",objval
       End If
       
       call end_timer("  +-- Metis")
    
        !-------------------------------------------------------------------------------
        !** Shift C-Style indicees *********************************************
        !-------------------------------------------------------------------------------
        ! This might be replaced with correctly setting up METIS OPTION NUMBERING
        ! Not done yet to avoid possible roots of erros.
        !-------------------------------------------------------------------------------
        call start_timer("  +-- Shift C-Style indicees")
        depart = depart + 1
        eind   = eind   + 1

       If (out_amount /= "PRODUCTION" ) then
          write(un_lf,fmt_msg_xAI0)"MinVal in depart:",minval(depart)
          write(un_lf,fmt_msg_xAI0)"MaxVal in depart:",maxval(depart)
       End If

       call end_timer("  +-- Shift C-Style indicees")
  
       ! Count elements per domain ******************************************
       call start_timer("  +-- Count elements per domain")
       
       do ii = 1, ne
          nelems_pp(depart(ii)) = nelems_pp(depart(ii)) + 1
       End do      
       
       call end_timer("  +-- Count elements per domain")
       
       !============================================================================
       ! Store partitioning as data on vtk unstructured grid
       !============================================================================
        !-------------------------------------------------------------------------------
        ! Colors of the elements are the number or the part of the mesh
        !-------------------------------------------------------------------------------
        if (out_amount == "DEBUG") then

            filename=''
            write(filename,'(A,I0,A)')trim(job_dir)//trim(project_name)//"_",ddc_nn,"_usg.vtk"
            Call write_vtk_data_int4_scalar_1D(&
                matrix = int(depart,4), &
                filename = trim(filename),  &
                desc = "PartNo", head = .TRUE.,location="CELL_DATA")
        End if
       
       ! Split Topology *****************************************************       
       call start_timer("  +-- Split Topology")

       Do ii = 1, Parts

            !------------------------------------------------------------------------------
            ! nn_el     -> Number of nodes per element
            ! nelems_pp -> Number of Elements per part
            !------------------------------------------------------------------------------
            ! PMesh%branches(ii)%leaves(5) = "Topology"
            !------------------------------------------------------------------------------
            PMesh%branches(ii)%leaves(5)%dat_no = nn_el*nelems_pp(ii)
            PMesh%branches(ii)%leaves(5)%pstat  = 1
            Allocate(PMesh%branches(ii)%leaves(5)%p_int8(nn_el*nelems_pp(ii)))

            !------------------------------------------------------------------------------
            ! PMesh%branches(ii)%leaves(6) = Hounsfield units of the voxels (HU Magnitude)
            ! As many entries as elements in part
            !------------------------------------------------------------------------------
            PMesh%branches(ii)%leaves(6)%dat_no = nelems_pp(ii)
            PMesh%branches(ii)%leaves(6)%pstat  = 1
            Allocate(PMesh%branches(ii)%leaves(6)%p_int8(nelems_pp(ii)))
            

            If (out_amount /= "PRODUCTION" ) then
                write(un_lf,fmt_msg_xAI0)"No Elems in part",ii,"=",nelems_pp(ii)
            End If

       End Do

       !------------------------------------------------------------------------------
       ! Count Elements per part
       ! Assemble node list
       !------------------------------------------------------------------------------
       nelems_pp = 0      
       Do ii = 1, ne
            !------------------------------------------------------------------------------
            ! depart(ii) -> depart is an array with length "elements of domain"
            ! array entries are the number of the part, the elements belongs to.
            ! While iterating over ne, the elements and nodes are sorted to their 
            ! corresponding parts. 
            !
            ! PMesh%branches(depart(ii))%leaves(6)%p_int8(nelems_pp(depart(ii)))
            ! For example, nelems_pp(part xy) = HU_magnitudes(corresponding element)
            !------------------------------------------------------------------------------
            nelems_pp(depart(ii)) = nelems_pp(depart(ii)) + 1
            
            !------------------------------------------------------------------------------
            ! Assign nodes to parts
            !------------------------------------------------------------------------------
            PMesh%branches(depart(ii))%leaves(5)%p_int8( &
                (nelems_pp(depart(ii))-1)*nn_el+1:nelems_pp(depart(ii))*nn_el &
                ) = eind(:,ii)


            !------------------------------------------------------------------------------
            ! Assign corresponding HUs to parts
            !------------------------------------------------------------------------------
            PMesh%branches(depart(ii))%leaves(6)%p_int8(nelems_pp(depart(ii))) = HU_magnitudes(ii)

       End do

       ! Renumber nodes in each part starting at 1 **************************
       Do ii = 1, Parts
          
          Allocate(cref( &
               minval(PMesh%branches(ii)%leaves(5)%p_int8) : &
               maxval(PMesh%branches(ii)%leaves(5)%p_int8))  )
          cref = 0
          idum = 0

          Do jj = 1, nelems_pp(ii)
             
             Do kk = 1, nn_el
                tmp_i8 = PMesh%branches(ii)%leaves(5)%p_int8((jj-1)*nn_el+kk)
                
                If ( cref(tmp_i8) == 0 ) Then
                   idum = idum + 1
                   cref(tmp_i8) = idum
                   tmp_i8       = idum
                Else
                   tmp_i8       = cref(tmp_i8)
                End If
                
                !PMesh%branches(ii)%leaves(5)%p_int8((jj-1)*20+kk) = tmp_i8
             End Do
             
          End Do

          ! Set Node pointer sizes and allocate memory
          ! PMesh%branches(ii)%leaves(1) : Node Numbers
          PMesh%branches(ii)%leaves(1)%dat_no = idum
          PMesh%branches(ii)%leaves(1)%pstat  = 1
          Allocate(PMesh%branches(ii)%leaves(1)%p_int8(idum))

          ! PMesh%branches(ii)%leaves(2) : Coordinates
          PMesh%branches(ii)%leaves(2)%dat_no = idum*3
          PMesh%branches(ii)%leaves(2)%pstat  = 1
          Allocate(PMesh%branches(ii)%leaves(2)%p_real8(idum*3))

          ! PMesh%branches(ii)%leaves(3) : Global Node numbers
          !> \todo fix inconsistent allocation of p_int8 starting from 0
          !> with dat_no not being idum+1
          PMesh%branches(ii)%leaves(3)%dat_no = idum
          PMesh%branches(ii)%leaves(3)%pstat  = 1
          Allocate(PMesh%branches(ii)%leaves(3)%p_int8(0:idum))
          
          write(un_lf,fmt_msg_xai0)"NODES in part",ii,"--",&
               PMesh%branches(ii)%leaves(1)%dat_no

          nnodes_pp(ii) = idum
          
          bounds(:,1) = lbound(cref)
          bounds(:,2) = ubound(cref)

          Do jj = bounds(1,1), bounds(1,2)
             PMesh%branches(ii)%leaves(3)%p_int8(cref(jj)) = jj
          End Do
          PMesh%branches(ii)%leaves(3)%p_int8(0) = 0
          
          Do jj = 1, PMesh%branches(ii)%leaves(3)%dat_no
             PMesh%branches(ii)%leaves(1)%p_int8(jj) = jj
             
             PMesh%branches(ii)%leaves(2)%p_real8((jj-1)*3+1:jj*3) = &
                  nodes(:,PMesh%branches(ii)%leaves(3)%p_int8(jj))
          end Do
          
          !*********************************************************
          ! DEBUG OUTPUT
          If (out_amount == "DEBUG") THEN 
             
             ! Output of partitioning in vtk format
             Call start_timer("  +-- VTK Output")
             
             ! Dynamic Allocation and reallocation takes place
             elems = reshape(cref(PMesh%branches(ii)%leaves(5)%p_int8),[nn_el,nelems_pp(ii)])
             
             Write(vtk_file,'(A,A,I0,A,I0,A)')trim(job_dir),"Part-",ii,"_",ddc_nn,".vtk"

             if (nn_el == 8) then
                Call write_vtk_unstructured_grid(Trim(vtk_file), &
                     reshape(PMesh%branches(ii)%leaves(2)%p_real8,[3,nnodes_pp(ii)]), &
                     elems(1:8,:))
             else if (nn_el == 20) then
                Call write_vtk_unstructured_grid(Trim(vtk_file), &
                     reshape(PMesh%branches(ii)%leaves(2)%p_real8,[3,nnodes_pp(ii)]), &
                     elems([1,3,5,7, 13,15,17,19, 2,4,6,8, 14,16,18,20, 9,10,11,12],:))
             end if
                
             Call write_vtk_data_int4_scalar_1D ( &
                  Int(PMesh%branches(ii)%leaves(3)%p_int8(1:nnodes_pp(ii)),4), &
                  Trim(vtk_file), "Glob_NIds", .True., "POINT_DATA")

             Call end_timer("  +-- VTK Output")

          End If
          ! DEBUG OUTPUT
          !*********************************************************

          deallocate(cref)
          
       End Do
       call end_timer("  +-- Split Topology")

    !************************************************************************
    ! If we should have exactly one part ************************************
    !************************************************************************
    else

       ! Node Numbers
       PMesh%branches(1)%leaves(1)%dat_no   = nnodes
       PMesh%branches(1)%leaves(1)%pstat = 1
       Allocate(PMesh%branches(1)%leaves(1)%p_int8(nnodes))
       
       ! Coordinates 
       PMesh%branches(1)%leaves(2)%dat_no   = nnodes*3_ik
       PMesh%branches(1)%leaves(2)%pstat = 1
       Allocate(PMesh%branches(1)%leaves(2)%p_real8(nnodes*3_ik))

       ! Global Node Number 
       PMesh%branches(1)%leaves(3)%dat_no   = nnodes
       PMesh%branches(1)%leaves(3)%pstat = 1
       Allocate(PMesh%branches(1)%leaves(3)%p_int8(nnodes))

       ! Element Numbers 
       PMesh%branches(1)%leaves(4)%dat_no   = ne
       PMesh%branches(1)%leaves(4)%pstat = 1
       Allocate(PMesh%branches(1)%leaves(4)%p_int8(ne))

       ! Topology 
       PMesh%branches(1)%leaves(5)%dat_no   = ne*nn_el
       PMesh%branches(1)%leaves(5)%pstat = 1
       Allocate(PMesh%branches(1)%leaves(5)%p_int8(ne*nn_el))

       !------------------------------------------------------------------------------
       ! PMesh%branches(ii)%leaves(6) = Hounsfield units of the voxels (HU Magnitude)
       !------------------------------------------------------------------------------
       PMesh%branches(1)%leaves(6)%dat_no = ne
       PMesh%branches(1)%leaves(6)%pstat = 1
       Allocate(PMesh%branches(1)%leaves(6)%p_int8(ne))

       ! Fill in nodes numbers and global node numbers
       Do ii = 1, nnodes
          PMesh%branches(1)%leaves(1)%p_int8(ii) = ii
          PMesh%branches(1)%leaves(3)%p_int8(ii) = ii
       End Do

       ! Fill in element numbers 
       Do ii = 1, ne
          PMesh%branches(1)%leaves(4)%p_int8(ii) = ii
       End Do

       ! Fill in coordinates 
       PMesh%branches(1)%leaves(2)%p_real8 = reshape(nodes,[nnodes*3_ik])

       ! Fill in topology 
       PMesh%branches(1)%leaves(5)%p_int8 = reshape(eind, [ne*nn_el])

       !------------------------------------------------------------------------------
       ! Fill in Hounsfield units of the voxels (HU Magnitude)
       !------------------------------------------------------------------------------
       PMesh%branches(1)%leaves(6)%p_int8 = reshape(HU_magnitudes, [ne])

      
       nnodes_pp     = nnodes
       nelems_pp     = ne
       nouter_nds_pp = 0
       
    end If
    
  end Subroutine part_mesh

end module mesh_partitioning
