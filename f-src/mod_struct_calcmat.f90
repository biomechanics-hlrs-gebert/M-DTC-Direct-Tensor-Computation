!==============================================================================
!> This Program calculates the effective numerical stiffness of
!> microstructured Volume elements
!>
Module calcmat

  use tensors
  Use decomp
  use mat_matrices
  use chain_routines
  use linFE
  use mpi
  
  implicit none

contains
 
  subroutine calc_effective_material_parameters(root, nn, ddc_nn, job_dir, fh_mpi) 

    !> \todo Make scaling and backrotation process correct!
    
    !==========================================================================
    !** Declarations **********************************************************
    Type(tBranch), Intent(InOut)          :: root
    Character(LEN=*), Intent(in)          :: job_dir
    integer(Kind=ik), Intent(in)          :: ddc_nn, nn
    Integer(kind=mpi_ik), Dimension(no_streams), Intent(in) :: fh_mpi
    
    Integer(kind=mpi_ik)                                    :: ierr
    Integer(kind=mpi_ik), Dimension(MPI_STATUS_SIZE)        :: status_mpi
    
    Character(len=*), Parameter                :: link_name="struct_calcmat_fmps"

    Integer(Kind=8)                                    :: no_dat, no_nodes,no_cnodes
    Integer(Kind=8), Dimension(:), Allocatable         :: no_cnodes_pp
    Integer(Kind=8)                                    :: macro_order
    CHARACTER(LEN=8)                                   :: el_type

    integer(Kind=8), Dimension(:), allocatable         :: cref_cnodes, tmp_ar_int8
    Real(kind=8)   , Dimension(:)  , allocatable       :: tmp_nn
    Real(kind=8)   , Dimension(:,:), allocatable       :: nodes, vv, ff, stiffness

    Real(kind=8)   , Dimension(:,:,:), allocatable     :: uu, rforces

    Real(kind=8)   , Dimension(:,:,:)    , allocatable :: calc_rforces
    Real(kind=8)   , Dimension(6)                :: ro_stress
    Real(kind=8)   , Dimension(6,6)              :: CC
    Real(kind=8)   , Dimension(6,6)              :: cc_mean, tmp_r6x6, EE, fv,meps
    Real(kind=8)   , Dimension( 8)               :: tmp_r8 
    Real(kind=8)   , Dimension(12)               :: tmp_r12

    Real(Kind=rk)  , Dimension(6,24)             :: int_strain, int_stress
    Real(kind=8)                                 :: E_div, nu_div, G_div
    Real(kind=8)                                 :: E_Modul, nu, rve_strain
    Real(kind=8)                                 :: v_elem, v_cube
    integer(Kind=ik)                             :: ii, jj, kk, ll
    integer(Kind=ik)                             :: no_elem_nodes,micro_elem_nodes
    integer(Kind=ik)                             :: no_lc
    Real(Kind=rk),    Dimension(3)               :: min_c, max_c
    Type(tBranch), Pointer                       :: ddc, loc_ddc, params, db
    Character(Len=mcl)                           :: desc

    Real(kind=rk), Dimension(1)                  :: tmp_real_fd1

    Type(tBranch),pointer                        :: mesh, res_tree
    Type(tLeaf),  pointer                        :: node_leaf_pointer
    Real(kind=rk)                                   :: alpha, phi, eta
    Integer                                         :: ii_phi,ii_eta,kk_phi,kk_eta
    Real(kind=rk), Dimension (3)                    :: n
    Real(kind=rk), Dimension (3,3)                  :: aa
    Real(kind=rk), Dimension (6,6)                  :: ee_orig, BB

    Integer      , Dimension (:,:,:,:), Allocatable :: ang
    Real(kind=rk), Dimension (:,:,:)  , Allocatable :: crit_1,crit_2
    Real(kind=rk), Dimension (0 :16)                :: crit_min
    Integer      , Dimension (3)                    :: s_loop,e_loop, mlc

    Real(kind=rk)                                   :: div_10_exp_jj
    Real(kind=rk)                                   :: cos_alpha, sin_alpha           
    Real(kind=rk)                                   :: One_Minus_cos_alpha 
    Real(kind=rk)                                   :: n12, n13, n23

    Real(kind=rk), Dimension(:), Allocatable :: delta, x_D_phy
    Character    , Allocatable, Dimension(:) :: EL_TYPE_ARR
    Integer                                  :: num_leaves, alloc_stat
    Type(tLeaf), Allocatable, Dimension(:)   :: leaf_list
    Real(kind=8)   , Dimension(:,:,:), allocatable      :: edat
    Integer(Kind=ik), Dimension(:), Allocatable         :: xa_n, xe_n
    Logical                                             :: success

    Character(len=9)                                    :: nn_char
    !  Type(tPs_ortho) :: oc
    !==========================================================================
    write(nn_char,'(I0)')ddc_nn

    Allocate(ang(3,0:180,0:180,0:90))
    Allocate(crit_1(0:180,0:180,0:90), crit_2(0:180,0:180,0:90))

    Select Case (timer_level)
    Case (3)
       call start_timer("  +-- Loading input data "//trim(nn_char))
    Case (2)
       call start_timer("  +-- Loading input data "//trim(nn_char))
    Case default
       continue
    End Select

    !** Load input params *****************************************************
    Call Search_branch("Input parameters", root, params, success)
    Call pd_get(params, "Young_s modulus", E_Modul)
    Call pd_get(params, "Poisson_s ratio", nu)
    Call pd_get(params, "Average strain on RVE", rve_strain)
    Call pd_get(params, "Element order on macro scale", macro_order)

    !**************************************************************************
    !** Load domain decomposition parameters **********************************

    !** Get global DDC parameters from root *********
    Call Search_branch("Global domain decomposition", root, ddc, success)

    call pd_get(ddc,"delta",delta)
    call pd_get(ddc,"x_D_phy",x_D_phy)

    !** Get local DDC parameters from root **********
    Write(desc,'(A,I0)')"Local domain Decomposition of domain no ",ddc_nn

    Call Search_branch(trim(desc), root, loc_ddc, success)
    call pd_get(loc_ddc,"xa_n",xa_n)
    call pd_get(loc_ddc,"xe_n",xe_n)

    min_c = Real(xa_n - 1,rk) * delta
    max_c = Real(xe_n    ,rk) * delta

    !** Set node number of macro element **************************************
    If (macro_order == 1) then
       no_elem_nodes = 8
       no_lc         = 24
    Else
       Write(*,*)"Element order other than 1 not supported"
       Goto 1000
    End If

    !**************************************************************************
    !** Get global mesh parameters ********************************************
    Call Search_branch("Domain "//trim(nn_char), root, db  , success)
    
    desc = ''
    Write(desc,'(A,I0)')'Mesh info of '//trim(project_name)//'_',ddc_nn
    Call search_branch(trim(desc), db, mesh, success)

    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    call log_tree(mesh,un_lf,.FALSE.)
    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    call pd_get(mesh, 'No of nodes in mesh',  no_nodes)
    call pd_get(mesh, 'No of elements in mesh',  no_dat)
   
    micro_elem_nodes = 20

    If (out_amount /= "PRODUCTION" ) then
       write(un_lf,FMT_MSG_AI0)"Set number of nodes per micro element to : " , &
            micro_elem_nodes
    End If
    
    call pd_get(mesh, 'No of cdofs in parts',  no_cnodes_pp)
    no_cnodes=sum(no_cnodes_pp)/3

    If (out_amount /= "PRODUCTION" ) then
       write(un_lf,*)
       write(un_lf,FMT_MSG_AI0)"Read no_nodes  from domain branch : " , no_nodes
       write(un_lf,FMT_MSG_AI0)"Read no_cnodes from domain branch : " , no_cnodes
    End If
    
    !**************************************************************************
    !** Allocate fields and read data *****************************************
    !**************************************************************************

    !** Copy Coordinate list nodes(3,no_nodes) ********************************
    Allocate (nodes(3,no_nodes),stat=alloc_stat)
    call pd_get(mesh, "Coordinates", node_leaf_pointer)
    nodes = reshape( node_leaf_pointer%p_real8 ,[3,no_nodes])

    !** Build Cross reference for constrained nodes from "Boundary_Ids" *******
    Allocate(cref_cnodes(no_cnodes),stat=alloc_stat)
    call alloc_err("cref_cnodes",alloc_stat)

    Do ii = 1, mesh%branches(2)%branches(1)%leaves(1)%dat_no, 9
       cref_cnodes((ii+8)/9) = mesh%branches(2)%branches(1)%leaves(1)%p_int8(ii)
    end Do

    !** Read displacement results ***********************************************
    Allocate(uu(3,no_nodes,no_lc),stat=alloc_stat)
    call alloc_err("uu",alloc_stat)

    call get_leaf_list("Displacements", mesh, num_leaves, leaf_list)
    If (out_amount /= "PRODUCTION" ) then
       write(un_lf,FMT_MSG_AI0)"Number of leaves with desc = Displacements     : ",&
            num_leaves
    End If
    
    Do ii = 1,num_leaves
       uu(:,:,ii) = reshape(leaf_list(ii)%p_real8,[3,no_nodes])
    End Do

    Deallocate(leaf_list)

    !** Read stress and strain results ******************************************
    Allocate(edat(15,no_lc,no_dat),stat=alloc_stat)
    call alloc_err("edat",alloc_stat)

    call get_leaf_list("Avg. Element Data", mesh, num_leaves, leaf_list)
    If (out_amount /= "PRODUCTION" ) then
       write(un_lf,FMT_MSG_AI0)"Number of leaves with desc = Avg. Element Data : ",&
            num_leaves
    End If
    
    Do ii = 1,num_leaves
       edat(:,ii,:) = reshape(leaf_list(ii)%p_real8,[15,no_dat])
    End Do

    Deallocate(leaf_list)

    !** Read Reaction Forces ****************************************************
    Allocate(rforces(3,no_nodes,no_lc),stat=alloc_stat)
    call alloc_err("rforces",alloc_stat)

    call get_leaf_list("Reaction Forces", mesh, num_leaves, leaf_list)
    If (out_amount /= "PRODUCTION" ) then
       write(un_lf,FMT_MSG_AI0)"Number of leaves with desc = Reaction Forces   : ",&
            num_leaves
    End If
    
    Do ii = 1,num_leaves
       rforces(:,:,ii) = reshape(leaf_list(ii)%p_real8, [3,no_nodes])
    End Do

    Deallocate(leaf_list)

    !****************************************************************************

    Allocate(calc_rforces(3,no_nodes,no_lc),stat=alloc_stat)
    call alloc_err("calc_rforces",alloc_stat)
    calc_rforces = 0._rk

    allocate(vv(no_lc,no_lc),stat=alloc_stat)
    call alloc_err("vv",alloc_stat)
    allocate(ff(no_lc,no_lc),stat=alloc_stat)
    call alloc_err("ff",alloc_stat)

    ff = 0._rk

    allocate(stiffness(no_lc,no_lc),stat=alloc_stat)
    call alloc_err("stiffness",alloc_stat)

    allocate(tmp_nn(no_elem_nodes),stat=alloc_stat)
    call alloc_err("tmp_nn",alloc_stat)

    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    If (out_amount == "DEBUG") THEN
       Write(un_lf,fmt_dbg_sep)
       write(un_lf,*)
       write(un_lf,*)'min uu      = ',minval(uu),'max uu      = ',maxval(uu)
       write(un_lf,*)'min rforces = ',minval(rforces),'max rforces = ',maxval(rforces)
       write(un_lf,*)
       Write(un_lf,fmt_dbg_sep)
    end if
    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    Select Case (timer_level)
    Case (3)
       call end_timer  ("  +-- Loading input data "//trim(nn_char))
       call start_timer("  +-- Calc material data "//trim(nn_char))
    Case (2)
       call end_timer  ("  +-- Loading input data "//trim(nn_char))
       call start_timer("  +-- Calc material data "//trim(nn_char))
    Case default
       continue
    End Select

    !****************************************************************************
    !** Calc material data ******************************************************
    !****************************************************************************

    !** Init C ******************************************************************
    E_div  = 1._rk / E_Modul
    nu_div = -nu * E_div
    G_div  = 2._rk * (1._rk+nu) * E_div

    CC = 0._rk

    CC(1,1:3) = (/  E_div, nu_div, nu_div /)
    CC(2,1:3) = (/ nu_div,  E_div, nu_div /)
    CC(3,1:3) = (/ nu_div, nu_div,  E_div /)
    CC(4,4) = G_div
    CC(5,5) = G_div
    CC(6,6) = G_div

    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    If (out_amount == "DEBUG") THEN
       Write(un_lf,fmt_dbg_sep)
       Call Write_real_matrix(un_lf,cc,6_ik,6_ik," Isotropic Compliance -- CC")
       Write(un_lf,fmt_dbg_sep)
    End if
    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    !** Reorder stresses and strains ********************************************
    Do jj = 1, no_dat
       Do ii = 1, no_lc

          ro_stress(1:4)  = edat(1:4,ii,jj)
          ro_stress(5)    = edat(6  ,ii,jj)
          ro_stress(6)    = edat(5  ,ii,jj)
          edat(1:6,ii,jj) = ro_stress

          ro_stress(1:4)   = edat(7:10,ii,jj)
          ro_stress(5)     = edat(12  ,ii,jj)
          ro_stress(6)     = edat(11  ,ii,jj)
          edat(7:12,ii,jj) = ro_stress 

       End Do
    End Do

    !** Init element and RVE volume *********************************************
    v_elem = delta(1)*delta(2)*delta(3)
    v_cube = x_D_phy(1)*x_D_phy(2)*x_D_phy(3)
    
    If (out_amount /= "PRODUCTION" ) then
       Write(un_lf,"('-- Volume per element -- ',E15.6)")v_elem
       Write(un_lf,"('-- MVE Volume         -- ',E15.6)")v_cube
       Write(un_lf,*)
    End If

    !** Loadcase init ***********************************************************
    call init_loadcase(rve_strain,vv)

    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    If (out_amount == "DEBUG") THEN
       Write(un_lf,fmt_dbg_sep)
       call write_real_matrix(un_lf,vv,no_lc,no_lc,&
            "Displacement matrix")
       Write(un_lf,fmt_dbg_sep)
    End if
    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    !** Calculate inverse of displacement matrix with LinPack *****************
    Call inverse(vv, no_lc, un_lf)

    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    If (out_amount == "DEBUG") THEN
       Write(un_lf,fmt_dbg_sep)
       call write_real_matrix(un_lf,vv,no_lc,no_lc,&
            "Inverted displacement matrix")
       Write(un_lf,fmt_dbg_sep)
    End if
    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    !****************************************************************************
    !** calc load matrix ********************************************************
    ff = 0._rk

    Do ii = 1, no_lc   !** Cycle through all load cases
       Do jj = 1, no_nodes      !** Cycle through all boundary nodes

          !** t_geom_xi transforms coordinates from geometry to xi space 
          !** Result(phi_nn) :  Real(kind=rk), dimension(8)
          tmp_nn = phi_nn(t_geom_xi(nodes(:,jj),min_c,max_c))

          Do kk = 1,3
             Do ll = 1,no_elem_nodes
                !** Setup of load matrix ****************************************
                !** 1st index - rows : Acumulate the reaction forces of the 
                !**                    micro model via the trail functions to 
                !**                    the macro element
                !** 2nd index - cols : Do first index step for all load cases
                ff((kk-1)*8+ll,ii) = ff((kk-1)*8+ll,ii) + rforces(kk,jj,ii)*tmp_nn(ll)
             End Do
          End Do

       End Do
    End Do

    If (out_amount /= "PRODUCTION" ) then
       call write_real_matrix(un_lf,ff,no_lc,no_lc, &
            "ff by summation of single force formula")
    End If
    
    !call add_leaf_to_branch(res_tree, "Domain forces", no_lc*no_lc, reshape(ff,[no_lc*no_lc]))
    !    Call MPI_FILE_WRITE_AT(FH_MPI(5), &
    !     Int(root%branches(4)%leaves(1)%lbound-1+(nn-1)*no_lc*no_lc, MPI_OFFSET_KIND), &
    !     reshape(ff,[no_lc*no_lc]), &
    !     Int(no_lc*no_lc,pd_mpi_ik), MPI_Real8, &
    !     status_mpi, ierr)
    
    !****************************************************************************
    !** Calc effective nummerical stiffness *************************************
    stiffness = matmul(ff,vv)
    !call add_leaf_to_branch(res_tree, "Effective numerical stiffness", no_lc*no_lc, &
    !     reshape(stiffness,[no_lc*no_lc]))
    !Call MPI_FILE_WRITE_AT(FH_MPI(5), &
    !     Int(root%branches(4)%leaves(2)%lbound-1+(nn-1)*no_lc*no_lc, MPI_OFFSET_KIND), &
    !     reshape(stiffness,[no_lc*no_lc]), &
    !     Int(no_lc*no_lc,pd_mpi_ik), MPI_Real8, &
    !     status_mpi, ierr)
    
    If (out_amount /= "PRODUCTION" ) then
       Call Write_real_matrix(un_lf,stiffness,no_lc,no_lc,"Stiffness")
    End If
    
    call check_sym(A=stiffness,name="Stiffness",sym_out=tmp_real_fd1(1))

    !Call add_leaf_to_branch(res_tree, &
    !    "Symmetry deviation - effective numerical stiffness", 1_pd_ik, &
    !    tmp_real_fd1)
    !Call MPI_FILE_WRITE_AT(FH_MPI(5), &
    !     Int(root%branches(4)%leaves(3)%lbound-1+(nn-1), MPI_OFFSET_KIND), &
    !     tmp_real_fd1, &
    !     Int(1,pd_mpi_ik), MPI_Real8, &
    !     status_mpi, ierr)
    
    !****************************************************************************
    !** Calc inverse of effective nummerical stiffness for cross check against **
    !** dislacement matrix                                                     **
    !Call inverse(stiffness, no_lc, un_lf)
    !Call Write_real_matrix(un_lf,stiffness,no_lc,no_lc,"Inverse Stiffness")
    !Call Write_real_matrix(un_lf,matmul(stiffness,ff(:,1)), no_lc, 1_ik, &
    !     "Inverse Stiffness o f(:,1)")

    !****************************************************************************
    !** Calc consistent force matrix
    Do ii = 1, 6
       fv(1,ii) = (ff( 2,ii) + ff( 3,ii) + ff( 6,ii) + ff( 7,ii)) / ( x_D_phy(2) * x_D_phy(3) )
       fv(2,ii) = (ff(11,ii) + ff(12,ii) + ff(15,ii) + ff(16,ii)) / ( x_D_phy(1) * x_D_phy(3) )
       fv(3,ii) = (ff(21,ii) + ff(22,ii) + ff(23,ii) + ff(24,ii)) / ( x_D_phy(1) * x_D_phy(2) )

       fv(4,ii) = 0.5_rk * &
            ( (ff( 3,ii) + ff( 4,ii) + ff( 7,ii) + ff( 8,ii)) / ( x_D_phy(1) * x_D_phy(3) ) + &
            (ff(10,ii) + ff(11,ii) + ff(14,ii) + ff(15,ii)) / ( x_D_phy(2) * x_D_phy(3) ) )
       fv(5,ii) = 0.5_rk * &
            ( (ff( 5,ii) + ff( 6,ii) + ff( 7,ii) + ff( 8,ii)) / ( x_D_phy(1) * x_D_phy(2) ) + &
            (ff(18,ii) + ff(19,ii) + ff(22,ii) + ff(23,ii)) / ( x_D_phy(2) * x_D_phy(3) ) )
       fv(6,ii) = 0.5_rk * &
            ( (ff(13,ii) + ff(14,ii) + ff(15,ii) + ff(16,ii)) / ( x_D_phy(1) * x_D_phy(2) ) + &
            (ff(19,ii) + ff(20,ii) + ff(23,ii) + ff(24,ii)) / ( x_D_phy(1) * x_D_phy(3) ) )
    End Do
    If (out_amount /= "PRODUCTION" ) then
       Call Write_real_matrix(un_lf, fv ,6_ik, 6_ik, "Konsistent force matrix")
    End If
    
    !****************************************************************************
    !** Calc integrated force matrix
    !Do jj = 1, 6
    !   Do ii = 1, 6
    !      fv(ii,jj) = sum(edat(ii,jj,:))
    !   End do
    !End Do

    !fv = fv * v_elem / v_cube
    !If (out_amount /= "PRODUCTION" ) then
    !   Call Write_real_matrix(un_lf, fv ,6_ik, 6_ik, "Integrated force matrix")
    !End If
    
    !****************************************************************************
    !** Calc averaged strains and stresses 
    !Do jj = 1, 24
    !   Do ii = 1,6
    !      int_stress(ii,jj) = sum(edat(ii,jj,:))
    !   End do
    !   Do ii = 7, 12
    !     int_strain(ii-6,jj) = sum(edat(ii,jj,:))
    !   End Do
    !End Do

    !int_strain = int_strain * v_elem / v_cube
    !int_stress = int_stress * v_elem / v_cube

    !If (out_amount /= "PRODUCTION" ) then
    !   write(un_lf,*)
    !   write(un_lf,"(150('='))")
    !   Write(un_lf,"(A)")"--------- Averaged strains of loadcases ---------"
    !   Write(un_lf,"(7(A20))")'-- ii --','-- E11 --','-- E22 --','-- E33 --','-- E12 --','-- E13 --','-- E23 --'
    !   Do ii = 1, 24
    !      write(un_lf,"(I20,6(E20.9))")ii,int_strain(:,ii)
    !   End Do

    !   write(un_lf,*)
    !   write(un_lf,"(150('='))")
    !   Write(un_lf,"(A)")"--------- Averaged stresses of loadcases ---------"
    !   Write(un_lf,"(7(A20))")'-- ii --','-- S11 --','-- S22 --','-- S33 --','-- S12 --','-- S13 --','-- S23 --'
    !   Do ii = 1, 24
    !      write(un_lf,"(I20,6(E20.9))")ii,int_stress(:,ii)
    !   End Do
    !End If
    
    !call add_leaf_to_branch(res_tree, "Averaged stresses", 6*no_lc, reshape(int_stress,[6*no_lc]))
!!$    Call MPI_FILE_WRITE_AT(FH_MPI(5), &
!!$         Int(root%branches(4)%leaves(4)%lbound-1+(nn-1)*6*no_lc, MPI_OFFSET_KIND), &
!!$         reshape(int_stress,[6*no_lc]), &
!!$         Int(6*no_lc,pd_mpi_ik), MPI_Real8, &
!!$         status_mpi, ierr)
!!$    !call add_leaf_to_branch(res_tree, "Averaged strains",  6*no_lc, reshape(int_strain,[6*no_lc]))
!!$    Call MPI_FILE_WRITE_AT(FH_MPI(5), &
!!$         Int(root%branches(4)%leaves(5)%lbound-1+(nn-1)*6*no_lc, MPI_OFFSET_KIND), &
!!$         reshape(int_strain,[6*no_lc]), &
!!$         Int(6*no_lc,pd_mpi_ik), MPI_Real8, &
!!$         status_mpi, ierr)
    
    !****************************************************************************
    !** Calc integrated strain matrix for first 6 loadcases
    !meps = int_strain(1:6,1:6)
    !If (out_amount /= "PRODUCTION" ) then
    !   Call Write_real_matrix(un_lf, meps ,6_ik, 6_ik, &
    !        "Integrated strain matrix of first 6 loadcases")
    !End If
    
    !****************************************************************************
    !** Calc theoretical effective stiffness
    !** Linear loadcases with constant strain fields ****************************
    meps(1,:) = [ 1.00E+06 , 0.00E+00 ,  0.00E+06 , 0.00E+06 , 0.00E+00, 0.00E+06]
    meps(2,:) = [ 0.00E+06 , 1.00E+06 ,  0.00E+06 , 0.00E+06 , 0.00E+00 ,0.00E+06]
    meps(3,:) = [ 0.00E+06 , 0.00E+00 ,  1.00E+06 , 0.00E+06 , 0.00E+06 ,0.00E+06]
    meps(4,:) = [ 0.00E+06 , 0.00E+00 ,  0.00E+06 ,-1.00E+06 , 0.00E+00 ,0.00E+06]
    meps(5,:) = [ 0.00E+06 , 0.00E+06 ,  0.00E+06 , 0.00E+06 ,-1.00E+06 ,0.00E+06]
    meps(6,:) = [ 0.00E+06 , 0.00E+00 ,  0.00E+06 , 0.00E+00 , 0.00E+06 ,1.00E+06]

    cc_mean = matmul(fv,meps)

    If (out_amount /= "PRODUCTION" ) then
       Call Write_real_matrix(un_lf, cc_mean ,6_ik, 6_ik, "Effective stiffness")
    End If
    
    call check_sym(A=cc_mean,name="Effective stiffness",sym_out=tmp_real_fd1(1))

    !call add_leaf_to_branch(res_tree, "Effective stiffness",  36_pd_ik, &
    !     reshape(cc_mean,[36_pd_ik]))
!!$    Call MPI_FILE_WRITE_AT(FH_MPI(5), &
!!$         Int(root%branches(4)%leaves(6)%lbound-1+(nn-1)*36, MPI_OFFSET_KIND), &
!!$         reshape(cc_mean,[36]), &
!!$         Int(36,pd_mpi_ik), MPI_Real8, &
!!$         status_mpi, ierr)
       
    !Call add_leaf_to_branch(res_tree, &
    !    "Symmetry deviation - effective stiffness",  1_pd_ik, tmp_real_fd1)
!!$    Call MPI_FILE_WRITE_AT(FH_MPI(5), &
!!$         Int(root%branches(4)%leaves(7)%lbound-1+(nn-1), MPI_OFFSET_KIND), &
!!$         tmp_real_fd1, &
!!$         Int(1,pd_mpi_ik), MPI_Real8, &
!!$         status_mpi, ierr)
    
    Select Case (timer_level)
    Case (3)
       call end_timer  ("  +-- Calc material data "//trim(nn_char))
       call start_timer("  +-- Back rotation of material matrix "//trim(nn_char))
    Case (2)
       call end_timer  ("  +-- Calc material data "//trim(nn_char))
       call start_timer("  +-- Back rotation of material matrix "//trim(nn_char))
    Case default
       continue
    End Select

    !****************************************************************************
    !** Back rotation of material matrix according to different criterias *******
    !****************************************************************************

    !EE = matmul(fv,transpose(meps))
    EE = (cc_mean + transpose(cc_mean)) / 2._rk

    If (out_amount /= "PRODUCTION" ) then
       Call Write_real_matrix(un_lf, ee ,6_ik, 6_ik, "Averaged Effective stiffness")
    End If
    call check_sym(A=ee,name="Averaged Effective stiffness",sym_out=tmp_real_fd1(1))

    !call add_leaf_to_branch(res_tree, "Averaged Effective stiffness",  36_pd_ik, &
    !     reshape(ee,[36_pd_ik]))
!!$    Call MPI_FILE_WRITE_AT(FH_MPI(5), &
!!$         Int(root%branches(4)%leaves(8)%lbound-1+(nn-1)*36, MPI_OFFSET_KIND), &
!!$         reshape(ee,[36]), &
!!$         Int(36,pd_mpi_ik), MPI_Real8, &
!!$         status_mpi, ierr)
    
    !call add_leaf_to_branch(res_tree, &
    !     "Symmetry deviation - Averaged effective stiffness",  1_pd_ik, tmp_real_fd1)
!!$    Call MPI_FILE_WRITE_AT(FH_MPI(5), &
!!$         Int(root%branches(4)%leaves(9)%lbound-1+(nn-1), MPI_OFFSET_KIND), &
!!$         tmp_real_fd1, &
!!$         Int(1,pd_mpi_ik), MPI_Real8, &
!!$         status_mpi, ierr)
    
    EE_Orig = EE

    !###############################################################################
    !###############################################################################

!!$  !==========================================
!!$
!!$  oc%E1 = 1._rk
!!$  oc%E2 = 2._rk
!!$  oc%E3 = 3._rk
!!$
!!$  oc%v12 = 2._rk/5._rk
!!$  oc%v13 = 1._rk/10._rk
!!$  oc%v23 = 1._rk/3._rk
!!$
!!$  oc%G12 = 1._rk
!!$  oc%G13 = 2._rk
!!$  oc%G23 = 3._rk
!!$
!!$  EE= Matrix_Ortho(oc)
!!$  Call inverse(EE, 6, un_lf)
!!$
!!$ EE=1._rk
!!$
!!$  desc="CG_2.4_784_c1_mono.raw"
!!$  open(unit=1234,file=trim(desc),action="write",status="replace",access="stream")
!!$  desc="CG_2.4_784_c2_ortho.raw"
!!$  open(unit=1235,file=trim(desc),action="write",status="replace",access="stream")
!!$  !==========================================

    !###############################################################################
    !###############################################################################

    kk_eta = 0
    kk_phi = 0
    kk = 0

    Do ii_eta = 0 , 90 , 1

       kk_phi = 0

       Do ii_phi = 0 , 180 , 1

          kk = 0

          Do ii = 0 , 180 , 1

             alpha = Real(ii,rk)     * pi_div_180
             phi   = Real(ii_phi,rk) * pi_div_180
             eta   = Real(ii_eta,rk) * pi_div_180

             n = [cos(phi)*sin(eta) , sin(phi)*sin(eta) , cos(eta) ]
             n = n / sqrt(sum(n*n))

             !aa = rot_alg(n,alpha)

             cos_alpha           = cos(alpha)
             sin_alpha           = sin(alpha)
             One_Minus_cos_alpha = 1._8 - cos_alpha
             n12                 = n(1)*n(2)
             n13                 = n(1)*n(3)                
             n23                 = n(2)*n(3)

             aa(1,1) = cos_alpha + n(1)*n(1)* One_Minus_cos_alpha
             aa(2,2) = cos_alpha + n(2)*n(2)* One_Minus_cos_alpha
             aa(3,3) = cos_alpha + n(3)*n(3)* One_Minus_cos_alpha 

             aa(1,2) = n12 * One_Minus_cos_alpha  - n(3) * sin_alpha
             aa(2,1) = n12 * One_Minus_cos_alpha  + n(3) * sin_alpha

             aa(1,3) = n13 * One_Minus_cos_alpha  + n(2) * sin_alpha
             aa(3,1) = n13 * One_Minus_cos_alpha  - n(2) * sin_alpha

             aa(2,3) = n23 * One_Minus_cos_alpha  - n(1) * sin_alpha
             aa(3,2) = n23 * One_Minus_cos_alpha  + n(1) * sin_alpha

             !BB = tra_R6(aa)

             BB(:,1) = [ aa(1,1)*aa(1,1) , aa(2,1)*aa(2,1) , aa(3,1)*aa(3,1) , &
                  sq2*aa(2,1)*aa(1,1) , sq2*aa(1,1)*aa(3,1) , sq2*aa(2,1)*aa(3,1) ]
             BB(:,2) = [ aa(1,2)*aa(1,2) , aa(2,2)*aa(2,2) , aa(3,2)*aa(3,2) , &
                  sq2*aa(2,2)*aa(1,2) , sq2*aa(1,2)*aa(3,2) , sq2*aa(2,2)*aa(3,2) ]
             BB(:,3) = [ aa(1,3)*aa(1,3) , aa(2,3)*aa(2,3) , aa(3,3)*aa(3,3) , &
                  sq2*aa(2,3)*aa(1,3) , sq2*aa(1,3)*aa(3,3) , sq2*aa(2,3)*aa(3,3) ]
             BB(:,4) = [ sq2*aa(1,1)*aa(1,2) , sq2*aa(2,1)*aa(2,2) , sq2*aa(3,1)*aa(3,2) , &
                  aa(2,1)*aa(1,2)+aa(2,2)*aa(1,1) , aa(1,1)*aa(3,2)+aa(1,2)*aa(3,1) , aa(2,1)*aa(3,2)+aa(2,2)*aa(3,1) ]
             BB(:,5) = [ sq2*aa(1,1)*aa(1,3) , sq2*aa(2,1)*aa(2,3) , sq2*aa(3,1)*aa(3,3) , &
                  aa(2,1)*aa(1,3)+aa(2,3)*aa(1,1) , aa(1,1)*aa(3,3)+aa(1,3)*aa(3,1) , aa(2,1)*aa(3,3)+aa(2,3)*aa(3,1) ]
             BB(:,6) = [ sq2*aa(1,2)*aa(1,3) , sq2*aa(2,2)*aa(2,3) , sq2*aa(3,2)*aa(3,3) , &
                  aa(2,2)*aa(1,3)+aa(2,3)*aa(1,2) , aa(1,2)*aa(3,3)+aa(1,3)*aa(3,2) , aa(2,2)*aa(3,3)+aa(2,3)*aa(3,2) ]

             !tmp_r6x6 = matmul(matmul(transpose(BB),EE),BB)

             tmp_r12(1) = &
                  BB(6,1) * &
                  (BB(6,4)*EE(6,6)+BB(5,4)*EE(6,5)+BB(4,4)*EE(6,4)+BB(3,4)*EE(6,3)+BB(2,4)*EE(6,2)+BB(1,4)*EE(6,1)) + &
                  BB(5,1) * &
                  (EE(5,6)*BB(6,4)+BB(5,4)*EE(5,5)+BB(4,4)*EE(5,4)+BB(3,4)*EE(5,3)+BB(2,4)*EE(5,2)+BB(1,4)*EE(5,1)) + &
                  BB(4,1) * &
                  (EE(4,6)*BB(6,4)+EE(4,5)*BB(5,4)+BB(4,4)*EE(4,4)+BB(3,4)*EE(4,3)+BB(2,4)*EE(4,2)+BB(1,4)*EE(4,1)) + &
                  BB(3,1) * &
                  (EE(3,6)*BB(6,4)+EE(3,5)*BB(5,4)+EE(3,4)*BB(4,4)+EE(3,3)*BB(3,4)+BB(2,4)*EE(3,2)+BB(1,4)*EE(3,1)) + &
                  BB(2,1) * &
                  (EE(2,6)*BB(6,4)+EE(2,5)*BB(5,4)+EE(2,4)*BB(4,4)+EE(2,3)*BB(3,4)+EE(2,2)*BB(2,4)+BB(1,4)*EE(2,1)) + &
                  BB(1,1) * &
                  (EE(1,6)*BB(6,4)+EE(1,5)*BB(5,4)+EE(1,4)*BB(4,4)+EE(1,3)*BB(3,4)+EE(1,2)*BB(2,4)+EE(1,1)*BB(1,4))
             tmp_r12(2) =  &
                  BB(6,1) * &
                  (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) + &
                  BB(5,1) * &
                  (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) + &
                  BB(4,1) * &
                  (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) + &
                  BB(3,1) * &
                  (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) + &
                  BB(2,1) * &
                  (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) + &
                  BB(1,1) * &
                  (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
             tmp_r12(3) = &
                  BB(6,1) * &
                  (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) + &
                  BB(5,1) * &
                  (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) + &
                  BB(4,1) * &
                  (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) + &
                  BB(3,1) * &
                  (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) + &
                  BB(2,1) * &
                  (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) + &
                  BB(1,1) * &
                  (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))
             tmp_r12(4) =  &
                  BB(6,2) * &
                  (BB(6,4)*EE(6,6)+BB(5,4)*EE(6,5)+BB(4,4)*EE(6,4)+BB(3,4)*EE(6,3)+BB(2,4)*EE(6,2)+BB(1,4)*EE(6,1)) + &
                  BB(5,2) * &
                  (EE(5,6)*BB(6,4)+BB(5,4)*EE(5,5)+BB(4,4)*EE(5,4)+BB(3,4)*EE(5,3)+BB(2,4)*EE(5,2)+BB(1,4)*EE(5,1)) + &
                  BB(4,2) * &
                  (EE(4,6)*BB(6,4)+EE(4,5)*BB(5,4)+BB(4,4)*EE(4,4)+BB(3,4)*EE(4,3)+BB(2,4)*EE(4,2)+BB(1,4)*EE(4,1)) + &
                  BB(3,2) * &
                  (EE(3,6)*BB(6,4)+EE(3,5)*BB(5,4)+EE(3,4)*BB(4,4)+EE(3,3)*BB(3,4)+BB(2,4)*EE(3,2)+BB(1,4)*EE(3,1)) + &
                  BB(2,2) * &
                  (EE(2,6)*BB(6,4)+EE(2,5)*BB(5,4)+EE(2,4)*BB(4,4)+EE(2,3)*BB(3,4)+EE(2,2)*BB(2,4)+BB(1,4)*EE(2,1)) + &
                  BB(1,2) * &
                  (EE(1,6)*BB(6,4)+EE(1,5)*BB(5,4)+EE(1,4)*BB(4,4)+EE(1,3)*BB(3,4)+EE(1,2)*BB(2,4)+EE(1,1)*BB(1,4))
             tmp_r12( 5) = &
                  BB(6,2) * &
                  (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) + &
                  BB(5,2) * &
                  (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) + &
                  BB(4,2) * &
                  (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) + &
                  BB(3,2) * &
                  (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) + &
                  BB(2,2) * &
                  (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) + &
                  BB(1,2) * &
                  (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
             tmp_r12( 6) = &
                  BB(6,2) * &
                  (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) + &
                  BB(5,2) * &
                  (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) + &
                  BB(4,2) * &
                  (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) + &
                  BB(3,2) * &
                  (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) + &
                  BB(2,2) * &
                  (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) + &
                  BB(1,2) * &
                  (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))
             tmp_r12( 7) = &
                  BB(6,3) * &
                  (BB(6,4)*EE(6,6)+BB(5,4)*EE(6,5)+BB(4,4)*EE(6,4)+BB(3,4)*EE(6,3)+BB(2,4)*EE(6,2)+BB(1,4)*EE(6,1)) + &
                  BB(5,3) * &
                  (EE(5,6)*BB(6,4)+BB(5,4)*EE(5,5)+BB(4,4)*EE(5,4)+BB(3,4)*EE(5,3)+BB(2,4)*EE(5,2)+BB(1,4)*EE(5,1)) + &
                  BB(4,3) * &
                  (EE(4,6)*BB(6,4)+EE(4,5)*BB(5,4)+BB(4,4)*EE(4,4)+BB(3,4)*EE(4,3)+BB(2,4)*EE(4,2)+BB(1,4)*EE(4,1)) + &
                  BB(3,3) * &
                  (EE(3,6)*BB(6,4)+EE(3,5)*BB(5,4)+EE(3,4)*BB(4,4)+EE(3,3)*BB(3,4)+BB(2,4)*EE(3,2)+BB(1,4)*EE(3,1)) + &
                  BB(2,3) * &
                  (EE(2,6)*BB(6,4)+EE(2,5)*BB(5,4)+EE(2,4)*BB(4,4)+EE(2,3)*BB(3,4)+EE(2,2)*BB(2,4)+BB(1,4)*EE(2,1)) + &
                  BB(1,3) * &
                  (EE(1,6)*BB(6,4)+EE(1,5)*BB(5,4)+EE(1,4)*BB(4,4)+EE(1,3)*BB(3,4)+EE(1,2)*BB(2,4)+EE(1,1)*BB(1,4))
             tmp_r12( 8) = &
                  BB(6,3) * &
                  (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) + &
                  BB(5,3) * &
                  (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) + &
                  BB(4,3) * &
                  (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) + &
                  BB(3,3) * &
                  (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) + &
                  BB(2,3) * &
                  (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) + &
                  BB(1,3) * &
                  (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
             tmp_r12( 9) = &
                  BB(6,3) * &
                  (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) + &
                  BB(5,3) * &
                  (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) + &
                  BB(4,3) * &
                  (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) + &
                  BB(3,3) * &
                  (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) + &
                  BB(2,3) * &
                  (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) + &
                  BB(1,3) * &
                  (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))
             tmp_r12(10) = &
                  BB(6,4) * &
                  (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) + &
                  BB(5,4) * &
                  (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) + &
                  BB(4,4) * &
                  (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) + &
                  BB(3,4) * &
                  (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) + &
                  BB(2,4) * &
                  (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) + &
                  BB(1,4) * &
                  (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
             tmp_r12(11) = &
                  BB(6,4) * &
                  (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) + &
                  BB(5,4) * &
                  (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) + &
                  BB(4,4) * &
                  (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) + &
                  BB(3,4) * &
                  (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) + &
                  BB(2,4) * &
                  (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) + &
                  BB(1,4) * &
                  (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))
             tmp_r12(12) = &
                  BB(6,5) * &
                  (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) + &
                  BB(5,5) * &
                  (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) + &
                  BB(4,5) * &
                  (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) + &
                  BB(3,5) * &
                  (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) + &
                  BB(2,5) * &
                  (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) + &
                  BB(1,5) * &
                  (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))

             ang(:,kk,kk_phi,kk_eta)  = [ii,ii_phi,ii_eta]

             !-- CR1 Monotropic ----------------------------------------
!!$             crit_1(kk,kk_phi,kk_eta) = (&
!!$                  sum((tmp_r6x6(1:4,5:6))*(tmp_r6x6(1:4,5:6))))
             ! Calculated in tmp_r6x6 :
             !  1  2  3  4  5  6
             !     7  8  9 10 11
             !       12 13 14 15
             !          16 17 18
             !             19 20
             !                21
             crit_1(kk,kk_phi,kk_eta) = ( &
                  tmp_r12( 2)*tmp_r12( 2) + tmp_r12( 3)*tmp_r12( 3) + &
                  tmp_r12( 5)*tmp_r12( 5) + tmp_r12( 6)*tmp_r12( 6) + &
                  tmp_r12( 8)*tmp_r12( 8) + tmp_r12( 9)*tmp_r12( 9) + &
                  tmp_r12(10)*tmp_r12(10) + tmp_r12(11)*tmp_r12(11)   &
                  )
             !-- CR2 Orthotropic ---------------------------------------
!!$             crit_2(kk,kk_phi,kk_eta) = (&
!!$                  sum(tmp_r6x6(1:3,4:6) * tmp_r6x6(1:3,4:6)) + &
!!$                  sum(tmp_r6x6( 4 ,5:6) * tmp_r6x6( 4 ,5:6)) + &
!!$                  tmp_r6x6( 5 , 6 ) * tmp_r6x6( 5 , 6 )     )
             crit_2(kk,kk_phi,kk_eta) = (&
                  tmp_r12( 1)*tmp_r12( 1) + tmp_r12( 2)*tmp_r12( 2) + tmp_r12( 3)*tmp_r12( 3) + &
                  tmp_r12( 4)*tmp_r12( 4) + tmp_r12( 5)*tmp_r12( 5) + tmp_r12( 6)*tmp_r12( 6) + &
                  tmp_r12( 7)*tmp_r12( 7) + tmp_r12( 8)*tmp_r12( 8) + tmp_r12( 9)*tmp_r12( 9) + &
                  tmp_r12(10)*tmp_r12(10) + tmp_r12(11)*tmp_r12(11) + &
                  tmp_r12(12)*tmp_r12(12) &
                  )
             kk = kk + 1
          End Do
          kk_phi = kk_phi + 1
       end Do
       kk_eta = kk_eta + 1
    end Do

    !=========================================================================
    !== Iteration of Crit_1 ==================================================

    crit_min    = 0._rk
    crit_min(0) = minval(crit_1)
    mlc         = minloc(crit_1)-1

    If (out_amount /= "PRODUCTION" ) then
       write(un_lf,FMT_MSG_A3I0)'Initial Minloc  CR_1 : ',mlc
       write(un_lf,FMT_MSG_AF0) 'Initial Minimum CR_1 : ',crit_min(0)
    End If
    
    jj = 1

    div_10_exp_jj = pi_div_180

    Do 

       mlc = minloc(crit_1(0:kk-1,0:kk_phi-1,0:kk_eta-1))-1

       s_loop = (ang(:,mlc(1),mlc(2),mlc(3))-1)*10
       e_loop = (ang(:,mlc(1),mlc(2),mlc(3))+1)*10

       kk_eta = 0
       kk_phi = 0
       kk = 0

       If (out_amount /= "PRODUCTION" ) then
          write(un_lf,FMT_MSG_AI0) 'Iteration            : ',jj
          write(un_lf,FMT_MSG_A3I0)'Loop start           : ',s_loop
          write(un_lf,FMT_MSG_A3I0)'Loop end             : ',e_loop
       End If
       
       div_10_exp_jj = div_10_exp_jj * 0.1_rk

       Do ii_eta = s_loop(3), e_loop(3)

          kk_phi = 0

          Do ii_phi = s_loop(2), e_loop(2)

             kk = 0

             Do ii = s_loop(1), e_loop(1)

                alpha = Real(ii,rk)     * div_10_exp_jj
                phi   = Real(ii_phi,rk) * div_10_exp_jj
                eta   = Real(ii_eta,rk) * div_10_exp_jj

                n = [cos(phi)*sin(eta) , sin(phi)*sin(eta) , cos(eta) ] 
                n = n / sqrt(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))

                !aa = rot_alg(n,alpha)

                cos_alpha           = cos(alpha)
                sin_alpha           = sin(alpha)
                One_Minus_cos_alpha = 1._8 - cos_alpha
                n12                 = n(1)*n(2)
                n13                 = n(1)*n(3)                
                n23                 = n(2)*n(3)

                aa(1,1) = cos_alpha + n(1)*n(1)* One_Minus_cos_alpha
                aa(2,2) = cos_alpha + n(2)*n(2)* One_Minus_cos_alpha
                aa(3,3) = cos_alpha + n(3)*n(3)* One_Minus_cos_alpha 

                aa(1,2) = n12 * One_Minus_cos_alpha  - n(3) * sin_alpha
                aa(2,1) = n12 * One_Minus_cos_alpha  + n(3) * sin_alpha

                aa(1,3) = n13 * One_Minus_cos_alpha  + n(2) * sin_alpha
                aa(3,1) = n13 * One_Minus_cos_alpha  - n(2) * sin_alpha

                aa(2,3) = n23 * One_Minus_cos_alpha  - n(1) * sin_alpha
                aa(3,2) = n23 * One_Minus_cos_alpha  + n(1) * sin_alpha

                !BB = tra_R6(aa)

                BB(:,1) = [ aa(1,1)*aa(1,1) , aa(2,1)*aa(2,1) , aa(3,1)*aa(3,1) , &
                     sq2*aa(2,1)*aa(1,1) , sq2*aa(1,1)*aa(3,1) , sq2*aa(2,1)*aa(3,1) ]
                BB(:,2) = [ aa(1,2)*aa(1,2) , aa(2,2)*aa(2,2) , aa(3,2)*aa(3,2) , &
                     sq2*aa(2,2)*aa(1,2) , sq2*aa(1,2)*aa(3,2) , sq2*aa(2,2)*aa(3,2) ]
                BB(:,3) = [ aa(1,3)*aa(1,3) , aa(2,3)*aa(2,3) , aa(3,3)*aa(3,3) , &
                     sq2*aa(2,3)*aa(1,3) , sq2*aa(1,3)*aa(3,3) , sq2*aa(2,3)*aa(3,3) ]
                BB(:,4) = [ sq2*aa(1,1)*aa(1,2) , sq2*aa(2,1)*aa(2,2) , sq2*aa(3,1)*aa(3,2) , &
                     aa(2,1)*aa(1,2)+aa(2,2)*aa(1,1) , aa(1,1)*aa(3,2)+aa(1,2)*aa(3,1) , aa(2,1)*aa(3,2)+aa(2,2)*aa(3,1) ]
                BB(:,5) = [ sq2*aa(1,1)*aa(1,3) , sq2*aa(2,1)*aa(2,3) , sq2*aa(3,1)*aa(3,3) , &
                     aa(2,1)*aa(1,3)+aa(2,3)*aa(1,1) , aa(1,1)*aa(3,3)+aa(1,3)*aa(3,1) , aa(2,1)*aa(3,3)+aa(2,3)*aa(3,1) ]
                BB(:,6) = [ sq2*aa(1,2)*aa(1,3) , sq2*aa(2,2)*aa(2,3) , sq2*aa(3,2)*aa(3,3) , &
                     aa(2,2)*aa(1,3)+aa(2,3)*aa(1,2) , aa(1,2)*aa(3,3)+aa(1,3)*aa(3,2) , aa(2,2)*aa(3,3)+aa(2,3)*aa(3,2) ]

                !tmp_r6x6 = matmul(matmul(transpose(BB),EE),BB)

                tmp_r8(1) = &
                     BB(6,1) * &
                     (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) &
                     +BB(5,1) * &
                     (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) &
                     +BB(4,1) * &
                     (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) &
                     +BB(3,1) * &
                     (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) &
                     +BB(2,1) * &
                     (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) &
                     +BB(1,1) * &
                     (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
                tmp_r8(2) = &
                     BB(6,1) * &
                     (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) &
                     +BB(5,1) * &
                     (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) &
                     +BB(4,1) * &
                     (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) &
                     +BB(3,1) * &
                     (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) &
                     +BB(2,1) * &
                     (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) &
                     +BB(1,1) * &
                     (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))
                tmp_r8( 3) = &
                     BB(6,2) * &
                     (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) &
                     +BB(5,2) * &
                     (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) &
                     +BB(4,2) * &
                     (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) &
                     +BB(3,2) * &
                     (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) &
                     +BB(2,2) * &
                     (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) &
                     +BB(1,2) * &
                     (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
                tmp_r8( 4) = &
                     BB(6,2) * &
                     (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) &
                     +BB(5,2) * &
                     (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) &
                     +BB(4,2) * &
                     (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) &
                     +BB(3,2) * &
                     (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) &
                     +BB(2,2) * &
                     (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) &
                     +BB(1,2) * &
                     (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))
                tmp_r8( 5) = &
                     BB(6,3) * &
                     (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) &
                     +BB(5,3) * &
                     (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) &
                     +BB(4,3) * &
                     (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) &
                     +BB(3,3) * &
                     (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) &
                     +BB(2,3) * &
                     (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) &
                     +BB(1,3) * &
                     (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
                tmp_r8( 6) = &
                     BB(6,3) * &
                     (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) &
                     +BB(5,3) * &
                     (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) &
                     +BB(4,3) * &
                     (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) &
                     +BB(3,3) * &
                     (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) &
                     +BB(2,3) * &
                     (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) &
                     +BB(1,3) * &
                     (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))
                tmp_r8( 7) = &
                     BB(6,4) * &
                     (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) &
                     +BB(5,4) * &
                     (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) &
                     +BB(4,4) * &
                     (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) &
                     +BB(3,4) * &
                     (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) &
                     +BB(2,4) * &
                     (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) &
                     +BB(1,4) * &
                     (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
                tmp_r8( 8) = &
                     BB(6,4) * &
                     (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) &
                     +BB(5,4) * &
                     (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) &
                     +BB(4,4) * &
                     (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) &
                     +BB(3,4) * &
                     (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) &
                     +BB(2,4) * &
                     (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) &
                     +BB(1,4) * &
                     (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))

                ang(:,kk,kk_phi,kk_eta)  = [ii,ii_phi,ii_eta]

                !-- CR1 Monotropic ----------------------------------------
!!$                crit_1(kk,kk_phi,kk_eta) = (&
!!$                     sum((tmp_r6x6(1:4,5:6))*(tmp_r6x6(1:4,5:6))))
                crit_1(kk,kk_phi,kk_eta) = ( &
                     tmp_r8( 1)*tmp_r8( 1) + tmp_r8( 2)*tmp_r8( 2) + &
                     tmp_r8( 3)*tmp_r8( 3) + tmp_r8( 4)*tmp_r8( 4) + &
                     tmp_r8( 5)*tmp_r8( 5) + tmp_r8( 6)*tmp_r8( 6) + &
                     tmp_r8( 7)*tmp_r8( 7) + tmp_r8( 8)*tmp_r8( 8)   &
                     )
                kk = kk + 1

             End Do
             kk_phi = kk_phi + 1
          end Do
          kk_eta = kk_eta + 1
       end Do

       crit_min(jj) = minval(crit_1(0:kk-1,0:kk_phi-1,0:kk_eta-1))

       If (out_amount /= "PRODUCTION" ) then
          write(un_lf,FMT_MSG_AF0)'Minimum CR_1         : ',crit_min(jj)
       End If
       
       If ( (abs(crit_min(jj-1)-crit_min(jj)) < num_zero) .OR. (jj >= 16)) Exit

       jj = jj + 1

    End Do

    !** Be aware that minloc starts off at field index 1 !!!
    mlc = minloc(crit_1(0:kk-1,0:kk_phi-1,0:kk_eta-1))-1

    alpha = Real( ang(1,mlc(1),mlc(2),mlc(3)),rk ) * pi / (180._rk*(10._rk**jj-1))
    phi   = Real( ang(2,mlc(1),mlc(2),mlc(3)),rk ) * pi / (180._rk*(10._rk**jj-1))
    eta   = Real( ang(3,mlc(1),mlc(2),mlc(3)),rk ) * pi / (180._rk*(10._rk**jj-1))

    n = [cos(phi)*sin(eta) , sin(phi)*sin(eta) , cos(eta) ] 
    n = n / sqrt(sum(n*n))

    If (out_amount /= "PRODUCTION" ) then
       write(un_lf,*)
       Write(un_lf,"(A,I0,A) ")" Solution converged after :",jj," iterations"
       Write(un_lf,"(A,F17.6,F17.6)")" With final citerion 1    :",&
            minval(crit_1(0:kk-1,0:kk_phi-1,0:kk_eta-1)),crit_1(mlc(1),mlc(2),mlc(3))
       Write(un_lf,"(A,F17.6)")" With final epsilon       :",&
            crit_min(jj-1)-crit_min(jj)
       Write(un_lf,"(A,F17.9,A)")" Final rotation angle  is : ",alpha
       Write(un_lf,"(A,3F17.9,A)")" Final rotation vector is : ",n
       Write(un_lf,*)
    End If
    
    tmp_real_fd1 = alpha 
    !Call add_leaf_to_branch(res_tree,"Rotation Angle CR_1" , 1_pd_ik, tmp_real_fd1)
!!$    Call MPI_FILE_WRITE_AT(FH_MPI(5), &
!!$         Int(root%branches(4)%leaves(10)%lbound-1+(nn-1), MPI_OFFSET_KIND), &
!!$         tmp_real_fd1, &
!!$         Int(1,pd_mpi_ik), MPI_Real8, &
!!$         status_mpi, ierr)
    
    !call add_leaf_to_branch(res_tree,"Rotation Vector CR_1", 3_pd_ik, n)
!!$    Call MPI_FILE_WRITE_AT(FH_MPI(5), &
!!$         Int(root%branches(4)%leaves(11)%lbound-1+(nn-1)*3, MPI_OFFSET_KIND), &
!!$         n, &
!!$         Int(3,pd_mpi_ik), MPI_Real8, &
!!$         status_mpi, ierr)
       
    !=========================================================================
    !== Inlining of EE =======================================================
    aa = rot_alg(n,alpha)
    BB = tra_R6(aa)
    EE = matmul(matmul(transpose(BB),EE_Orig),BB)

    If (out_amount /= "PRODUCTION" ) then
       call write_real_matrix(un_lf,EE,6_ik,6_ik,"Backrotated anisotropic stiffness CR_1")
    End If
    
    !=========================================================

    If ( (EE(1,1) < EE(2,2)) .AND.  &
         (EE(1,1) < EE(3,3)) .AND.  (EE(2,2) < EE(3,3))         ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"123"
       Continue
       
    Else If ( (EE(1,1) < EE(2,2)) .AND.  &
         (EE(1,1) < EE(3,3)) .AND.  (EE(2,2) > EE(3,3)) ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"132"

       !** 132 => 123 *****************
       n = aa(:,1)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

    Else If ( (EE(1,1) < EE(2,2)) .AND.  &
         (EE(1,1) > EE(3,3)) .AND.  (EE(2,2) > EE(3,3)) ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"231"

       !** 231 => 132 *****************
       n = aa(:,2)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)

       !** 132 => 123 *****************
       n = aa(:,1)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

    Else If ( (EE(1,1) > EE(2,2)) .AND.  &
         (EE(1,1) < EE(3,3)) .AND.  (EE(2,2) < EE(3,3)) ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"213"

       !** 213 => 123 *****************
       n = aa(:,3)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

    Else If ( (EE(1,1) > EE(2,2)) .AND.  &
         (EE(1,1) > EE(3,3)) .AND.  (EE(2,2) < EE(3,3)) ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"312"

       !** 312 => 132 *****************
       n = aa(:,3)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

       !** 132 => 123 *****************
       n = aa(:,1)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)  

    Else If ( (EE(1,1) > EE(2,2)) .AND.  &
         (EE(1,1) > EE(3,3)) .AND.  (EE(2,2) > EE(3,3)) ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"321"

       !** 321 => 123 *****************
       n = aa(:,2)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)  

    End If

    BB = tra_R6(aa)
    EE = matmul(matmul(transpose(BB),EE_Orig),BB)

    If (out_amount /= "PRODUCTION" ) then
       call write_real_matrix(un_lf, aa,3_ik,3_ik,"Final coordinate system CR_1")
       call write_real_matrix(un_lf, EE,6_ik,6_ik,"Inlined anisotropic stiffness CR_1")
    End If
    
    call check_sym(A=EE,name="Inlined anisotropic stiffness CR_1")

    !call add_leaf_to_branch(res_tree,"Final coordinate system CR_1", 9_pd_ik, &
    !     reshape(aa,[9_pd_ik]))
!!$    Call MPI_FILE_WRITE_AT(FH_MPI(5), &
!!$         Int(root%branches(4)%leaves(12)%lbound-1+(nn-1)*9, MPI_OFFSET_KIND), &
!!$         reshape(aa,[9_pd_ik]), &
!!$         Int(9,pd_mpi_ik), MPI_Real8, &
!!$         status_mpi, ierr)

    !call add_leaf_to_branch(res_tree,"Optimized Effective stiffness CR_1", 36_pd_ik, &
    !     reshape(EE,[36_pd_ik]))
!!$    Call MPI_FILE_WRITE_AT(FH_MPI(5), &
!!$         Int(root%branches(4)%leaves(13)%lbound-1+(nn-1)*36, MPI_OFFSET_KIND), &
!!$         reshape(EE,[36_pd_ik]), &
!!$         Int(36,pd_mpi_ik), MPI_Real8, &
!!$         status_mpi, ierr)
 
    !=========================================================================
    !== Iteration of Crit_2 ==================================================
    EE = EE_Orig

    kk_eta = 0
    kk_phi = 0
    kk = 0

    Do ii_eta = 0 , 90 , 1
       kk_phi = 0
       Do ii_phi = 0 , 180 , 1
          kk = 0
          Do ii = 0 , 180 , 1
             ang(:,kk,kk_phi,kk_eta)  = [ii,ii_phi,ii_eta]
             kk = kk + 1
          End Do
          kk_phi = kk_phi + 1
       End Do
       kk_eta = kk_eta + 1
    End Do

    crit_min = 0._rk
    crit_min(0) = minval(crit_2)

    mlc = minloc(crit_2)-1

    If (out_amount /= "PRODUCTION" ) then
       write(un_lf,FMT_MSG_A3I0)'Initial Minloc  CR_2 : ',mlc
       write(un_lf,FMT_MSG_AF0) 'Initial Minimum CR_2 : ',crit_min(0)
    End If
    
    jj = 1

    div_10_exp_jj = pi_div_180

    Do 

       mlc = minloc(crit_2(0:kk-1,0:kk_phi-1,0:kk_eta-1))-1

       s_loop = (ang(:,mlc(1),mlc(2),mlc(3))-1)*10
       e_loop = (ang(:,mlc(1),mlc(2),mlc(3))+1)*10

       kk_eta = 0
       kk_phi = 0
       kk = 0

       If (out_amount /= "PRODUCTION" ) then
          write(un_lf,FMT_MSG_AI0) 'Iteration            : ',jj
          write(un_lf,FMT_MSG_A3I0)'Loop start           : ',s_loop
          write(un_lf,FMT_MSG_A3I0)'Loop end             : ',e_loop
       End If
       
       div_10_exp_jj = div_10_exp_jj * 0.1_rk

       Do ii_eta = s_loop(3), e_loop(3)

          kk_phi = 0

          Do ii_phi = s_loop(2), e_loop(2)

             kk = 0

             Do ii = s_loop(1), e_loop(1)

                alpha = Real(ii,rk)     * div_10_exp_jj
                phi   = Real(ii_phi,rk) * div_10_exp_jj
                eta   = Real(ii_eta,rk) * div_10_exp_jj

                n = [cos(phi)*sin(eta) , sin(phi)*sin(eta) , cos(eta) ] 

                cos_alpha           = cos(alpha)
                sin_alpha           = sin(alpha)
                One_Minus_cos_alpha = 1._8 - cos_alpha
                n12                 = n(1)*n(2)
                n13                 = n(1)*n(3)                
                n23                 = n(2)*n(3)

                aa(1,1) = cos_alpha + n(1)*n(1)* One_Minus_cos_alpha
                aa(2,2) = cos_alpha + n(2)*n(2)* One_Minus_cos_alpha
                aa(3,3) = cos_alpha + n(3)*n(3)* One_Minus_cos_alpha 

                aa(1,2) = n12 * One_Minus_cos_alpha  - n(3) * sin_alpha
                aa(2,1) = n12 * One_Minus_cos_alpha  + n(3) * sin_alpha

                aa(1,3) = n13 * One_Minus_cos_alpha  + n(2) * sin_alpha
                aa(3,1) = n13 * One_Minus_cos_alpha  - n(2) * sin_alpha

                aa(2,3) = n23 * One_Minus_cos_alpha  - n(1) * sin_alpha
                aa(3,2) = n23 * One_Minus_cos_alpha  + n(1) * sin_alpha

                BB(:,1) = [ aa(1,1)*aa(1,1) , aa(2,1)*aa(2,1) , aa(3,1)*aa(3,1) , &
                     sq2*aa(2,1)*aa(1,1) , sq2*aa(1,1)*aa(3,1) , sq2*aa(2,1)*aa(3,1) ]
                BB(:,2) = [ aa(1,2)*aa(1,2) , aa(2,2)*aa(2,2) , aa(3,2)*aa(3,2) , &
                     sq2*aa(2,2)*aa(1,2) , sq2*aa(1,2)*aa(3,2) , sq2*aa(2,2)*aa(3,2) ]
                BB(:,3) = [ aa(1,3)*aa(1,3) , aa(2,3)*aa(2,3) , aa(3,3)*aa(3,3) , &
                     sq2*aa(2,3)*aa(1,3) , sq2*aa(1,3)*aa(3,3) , sq2*aa(2,3)*aa(3,3) ]
                BB(:,4) = [ sq2*aa(1,1)*aa(1,2) , sq2*aa(2,1)*aa(2,2) , sq2*aa(3,1)*aa(3,2) , &
                     aa(2,1)*aa(1,2)+aa(2,2)*aa(1,1) , aa(1,1)*aa(3,2)+aa(1,2)*aa(3,1) , &
                     aa(2,1)*aa(3,2)+aa(2,2)*aa(3,1) ]
                BB(:,5) = [ sq2*aa(1,1)*aa(1,3) , sq2*aa(2,1)*aa(2,3) , sq2*aa(3,1)*aa(3,3) , &
                     aa(2,1)*aa(1,3)+aa(2,3)*aa(1,1) , aa(1,1)*aa(3,3)+aa(1,3)*aa(3,1) , &
                     aa(2,1)*aa(3,3)+aa(2,3)*aa(3,1) ]
                BB(:,6) = [ sq2*aa(1,2)*aa(1,3) , sq2*aa(2,2)*aa(2,3) , sq2*aa(3,2)*aa(3,3) , &
                     aa(2,2)*aa(1,3)+aa(2,3)*aa(1,2) , aa(1,2)*aa(3,3)+aa(1,3)*aa(3,2) , &
                     aa(2,2)*aa(3,3)+aa(2,3)*aa(3,2) ]

                tmp_r12(1) = &
                     BB(6,1) * &
                     (BB(6,4)*EE(6,6)+BB(5,4)*EE(6,5)+BB(4,4)*EE(6,4)+BB(3,4)*EE(6,3)+BB(2,4)*EE(6,2)+BB(1,4)*EE(6,1)) + &
                     BB(5,1) * &
                     (EE(5,6)*BB(6,4)+BB(5,4)*EE(5,5)+BB(4,4)*EE(5,4)+BB(3,4)*EE(5,3)+BB(2,4)*EE(5,2)+BB(1,4)*EE(5,1)) + &
                     BB(4,1) * &
                     (EE(4,6)*BB(6,4)+EE(4,5)*BB(5,4)+BB(4,4)*EE(4,4)+BB(3,4)*EE(4,3)+BB(2,4)*EE(4,2)+BB(1,4)*EE(4,1)) + &
                     BB(3,1) * &
                     (EE(3,6)*BB(6,4)+EE(3,5)*BB(5,4)+EE(3,4)*BB(4,4)+EE(3,3)*BB(3,4)+BB(2,4)*EE(3,2)+BB(1,4)*EE(3,1)) + &
                     BB(2,1) * &
                     (EE(2,6)*BB(6,4)+EE(2,5)*BB(5,4)+EE(2,4)*BB(4,4)+EE(2,3)*BB(3,4)+EE(2,2)*BB(2,4)+BB(1,4)*EE(2,1)) + &
                     BB(1,1) * &
                     (EE(1,6)*BB(6,4)+EE(1,5)*BB(5,4)+EE(1,4)*BB(4,4)+EE(1,3)*BB(3,4)+EE(1,2)*BB(2,4)+EE(1,1)*BB(1,4))
                tmp_r12(2) =  &
                     BB(6,1) * &
                     (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) + &
                     BB(5,1) * &
                     (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) + &
                     BB(4,1) * &
                     (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) + &
                     BB(3,1) * &
                     (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) + &
                     BB(2,1) * &
                     (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) + &
                     BB(1,1) * &
                     (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
                tmp_r12(3) = &
                     BB(6,1) * &
                     (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) + &
                     BB(5,1) * &
                     (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) + &
                     BB(4,1) * &
                     (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) + &
                     BB(3,1) * &
                     (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) + &
                     BB(2,1) * &
                     (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) + &
                     BB(1,1) * &
                     (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))
                tmp_r12(4) =  &
                     BB(6,2) * &
                     (BB(6,4)*EE(6,6)+BB(5,4)*EE(6,5)+BB(4,4)*EE(6,4)+BB(3,4)*EE(6,3)+BB(2,4)*EE(6,2)+BB(1,4)*EE(6,1)) + &
                     BB(5,2) * &
                     (EE(5,6)*BB(6,4)+BB(5,4)*EE(5,5)+BB(4,4)*EE(5,4)+BB(3,4)*EE(5,3)+BB(2,4)*EE(5,2)+BB(1,4)*EE(5,1)) + &
                     BB(4,2) * &
                     (EE(4,6)*BB(6,4)+EE(4,5)*BB(5,4)+BB(4,4)*EE(4,4)+BB(3,4)*EE(4,3)+BB(2,4)*EE(4,2)+BB(1,4)*EE(4,1)) + &
                     BB(3,2) * &
                     (EE(3,6)*BB(6,4)+EE(3,5)*BB(5,4)+EE(3,4)*BB(4,4)+EE(3,3)*BB(3,4)+BB(2,4)*EE(3,2)+BB(1,4)*EE(3,1)) + &
                     BB(2,2) * &
                     (EE(2,6)*BB(6,4)+EE(2,5)*BB(5,4)+EE(2,4)*BB(4,4)+EE(2,3)*BB(3,4)+EE(2,2)*BB(2,4)+BB(1,4)*EE(2,1)) + &
                     BB(1,2) * &
                     (EE(1,6)*BB(6,4)+EE(1,5)*BB(5,4)+EE(1,4)*BB(4,4)+EE(1,3)*BB(3,4)+EE(1,2)*BB(2,4)+EE(1,1)*BB(1,4))
                tmp_r12( 5) = &
                     BB(6,2) * &
                     (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) + &
                     BB(5,2) * &
                     (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) + &
                     BB(4,2) * &
                     (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) + &
                     BB(3,2) * &
                     (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) + &
                     BB(2,2) * &
                     (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) + &
                     BB(1,2) * &
                     (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
                tmp_r12( 6) = &
                     BB(6,2) * &
                     (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) + &
                     BB(5,2) * &
                     (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) + &
                     BB(4,2) * &
                     (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) + &
                     BB(3,2) * &
                     (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) + &
                     BB(2,2) * &
                     (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) + &
                     BB(1,2) * &
                     (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))
                tmp_r12( 7) = &
                     BB(6,3) * &
                     (BB(6,4)*EE(6,6)+BB(5,4)*EE(6,5)+BB(4,4)*EE(6,4)+BB(3,4)*EE(6,3)+BB(2,4)*EE(6,2)+BB(1,4)*EE(6,1)) + &
                     BB(5,3) * &
                     (EE(5,6)*BB(6,4)+BB(5,4)*EE(5,5)+BB(4,4)*EE(5,4)+BB(3,4)*EE(5,3)+BB(2,4)*EE(5,2)+BB(1,4)*EE(5,1)) + &
                     BB(4,3) * &
                     (EE(4,6)*BB(6,4)+EE(4,5)*BB(5,4)+BB(4,4)*EE(4,4)+BB(3,4)*EE(4,3)+BB(2,4)*EE(4,2)+BB(1,4)*EE(4,1)) + &
                     BB(3,3) * &
                     (EE(3,6)*BB(6,4)+EE(3,5)*BB(5,4)+EE(3,4)*BB(4,4)+EE(3,3)*BB(3,4)+BB(2,4)*EE(3,2)+BB(1,4)*EE(3,1)) + &
                     BB(2,3) * &
                     (EE(2,6)*BB(6,4)+EE(2,5)*BB(5,4)+EE(2,4)*BB(4,4)+EE(2,3)*BB(3,4)+EE(2,2)*BB(2,4)+BB(1,4)*EE(2,1)) + &
                     BB(1,3) * &
                     (EE(1,6)*BB(6,4)+EE(1,5)*BB(5,4)+EE(1,4)*BB(4,4)+EE(1,3)*BB(3,4)+EE(1,2)*BB(2,4)+EE(1,1)*BB(1,4))
                tmp_r12( 8) = &
                     BB(6,3) * &
                     (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) + &
                     BB(5,3) * &
                     (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) + &
                     BB(4,3) * &
                     (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) + &
                     BB(3,3) * &
                     (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) + &
                     BB(2,3) * &
                     (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) + &
                     BB(1,3) * &
                     (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
                tmp_r12( 9) = &
                     BB(6,3) * &
                     (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) + &
                     BB(5,3) * &
                     (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) + &
                     BB(4,3) * &
                     (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) + &
                     BB(3,3) * &
                     (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) + &
                     BB(2,3) * &
                     (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) + &
                     BB(1,3) * &
                     (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))
                tmp_r12(10) = &
                     BB(6,4) * &
                     (BB(6,5)*EE(6,6)+BB(5,5)*EE(6,5)+BB(4,5)*EE(6,4)+BB(3,5)*EE(6,3)+BB(2,5)*EE(6,2)+BB(1,5)*EE(6,1)) + &
                     BB(5,4) * &
                     (EE(5,6)*BB(6,5)+BB(5,5)*EE(5,5)+BB(4,5)*EE(5,4)+BB(3,5)*EE(5,3)+BB(2,5)*EE(5,2)+BB(1,5)*EE(5,1)) + &
                     BB(4,4) * &
                     (EE(4,6)*BB(6,5)+EE(4,5)*BB(5,5)+EE(4,4)*BB(4,5)+BB(3,5)*EE(4,3)+BB(2,5)*EE(4,2)+BB(1,5)*EE(4,1)) + &
                     BB(3,4) * &
                     (EE(3,6)*BB(6,5)+EE(3,5)*BB(5,5)+EE(3,4)*BB(4,5)+EE(3,3)*BB(3,5)+BB(2,5)*EE(3,2)+BB(1,5)*EE(3,1)) + &
                     BB(2,4) * &
                     (EE(2,6)*BB(6,5)+EE(2,5)*BB(5,5)+EE(2,4)*BB(4,5)+EE(2,3)*BB(3,5)+EE(2,2)*BB(2,5)+BB(1,5)*EE(2,1)) + &
                     BB(1,4) * &
                     (EE(1,6)*BB(6,5)+EE(1,5)*BB(5,5)+EE(1,4)*BB(4,5)+EE(1,3)*BB(3,5)+EE(1,2)*BB(2,5)+EE(1,1)*BB(1,5))
                tmp_r12(11) = &
                     BB(6,4) * &
                     (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) + &
                     BB(5,4) * &
                     (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) + &
                     BB(4,4) * &
                     (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) + &
                     BB(3,4) * &
                     (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) + &
                     BB(2,4) * &
                     (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) + &
                     BB(1,4) * &
                     (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))
                tmp_r12(12) = &
                     BB(6,5) * &
                     (BB(6,6)*EE(6,6)+BB(5,6)*EE(6,5)+BB(4,6)*EE(6,4)+BB(3,6)*EE(6,3)+BB(2,6)*EE(6,2)+BB(1,6)*EE(6,1)) + &
                     BB(5,5) * &
                     (EE(5,6)*BB(6,6)+EE(5,5)*BB(5,6)+BB(4,6)*EE(5,4)+BB(3,6)*EE(5,3)+BB(2,6)*EE(5,2)+BB(1,6)*EE(5,1)) + &
                     BB(4,5) * &
                     (EE(4,6)*BB(6,6)+EE(4,5)*BB(5,6)+EE(4,4)*BB(4,6)+BB(3,6)*EE(4,3)+BB(2,6)*EE(4,2)+BB(1,6)*EE(4,1)) + &
                     BB(3,5) * &
                     (EE(3,6)*BB(6,6)+EE(3,5)*BB(5,6)+EE(3,4)*BB(4,6)+EE(3,3)*BB(3,6)+BB(2,6)*EE(3,2)+BB(1,6)*EE(3,1)) + &
                     BB(2,5) * &
                     (EE(2,6)*BB(6,6)+EE(2,5)*BB(5,6)+EE(2,4)*BB(4,6)+EE(2,3)*BB(3,6)+EE(2,2)*BB(2,6)+BB(1,6)*EE(2,1)) + &
                     BB(1,5) * &
                     (EE(1,6)*BB(6,6)+EE(1,5)*BB(5,6)+EE(1,4)*BB(4,6)+EE(1,3)*BB(3,6)+EE(1,2)*BB(2,6)+EE(1,1)*BB(1,6))

                ang(:,kk,kk_phi,kk_eta)  = [ii,ii_phi,ii_eta]

                crit_2(kk,kk_phi,kk_eta) = (&
                     tmp_r12( 1)*tmp_r12( 1) + tmp_r12( 2)*tmp_r12( 2) + tmp_r12( 3)*tmp_r12( 3) + &
                     tmp_r12( 4)*tmp_r12( 4) + tmp_r12( 5)*tmp_r12( 5) + tmp_r12( 6)*tmp_r12( 6) + &
                     tmp_r12( 7)*tmp_r12( 7) + tmp_r12( 8)*tmp_r12( 8) + tmp_r12( 9)*tmp_r12( 9) + &
                     tmp_r12(10)*tmp_r12(10) + tmp_r12(11)*tmp_r12(11) + &
                     tmp_r12(12)*tmp_r12(12) &
                     )
                kk = kk + 1

             End Do
             kk_phi = kk_phi + 1
          end Do
          kk_eta = kk_eta + 1
       end Do

       crit_min(jj) = minval(crit_2(0:kk-1,0:kk_phi-1,0:kk_eta-1))

       !write(un_lf,FMT_MSG_AF0)'Minimum CR_2         : ',crit_min(jj)
       If (out_amount /= "PRODUCTION" ) then
          write(un_lf,*)'Minimum CR_2         : ',crit_min(jj)
          write(un_lf,*)'Minloc  CR_2         : ',minloc(crit_2(0:kk-1,0:kk_phi-1,0:kk_eta-1))
          write(un_lf,*)'kk, kk_phi, kk_eta   : ',kk,kk_phi,kk_eta
       End If
       
       If ( (abs(crit_min(jj-1)-crit_min(jj)) < num_zero) .OR. (jj >= 16)) Exit

       jj = jj + 1

    End Do

    mlc = minloc(crit_2(0:kk-1,0:kk_phi-1,0:kk_eta-1))-1

    alpha = Real( ang(1,mlc(1),mlc(2),mlc(3)),rk ) * pi / (180._rk*(10._rk**jj-1))
    phi   = Real( ang(2,mlc(1),mlc(2),mlc(3)),rk ) * pi / (180._rk*(10._rk**jj-1))
    eta   = Real( ang(3,mlc(1),mlc(2),mlc(3)),rk ) * pi / (180._rk*(10._rk**jj-1))

    n = [cos(phi)*sin(eta) , sin(phi)*sin(eta) , cos(eta) ] 
    n = n / sqrt(sum(n*n))

    If (out_amount /= "PRODUCTION" ) then
       write(un_lf,*)
       Write(un_lf,"(A,I0,A) ")" Solution converged after :",jj," iterations"
       Write(un_lf,"(A,F17.6)")" With final citerion 2    :",&
            minval(crit_2(1:kk-2,1:kk_phi-2,1:kk_eta-2))
       Write(un_lf,"(A,F17.6)")" With final epsilon       :",&
            crit_min(jj-1)-crit_min(jj)
       Write(un_lf,"(A,3F17.9,A)")" Final rotation angle  is : ",alpha
       Write(un_lf,"(A,3F17.9,A)")" Final rotation vector is : ",n
       Write(un_lf,*)
    End If
    
    tmp_real_fd1 = alpha 
    !call add_leaf_to_branch(res_tree,"Rotation Angle CR_2" , 1_pd_ik, tmp_real_fd1)
!!$    Call MPI_FILE_WRITE_AT(FH_MPI(5), &
!!$         Int(root%branches(4)%leaves(14)%lbound-1+(nn-1), MPI_OFFSET_KIND), &
!!$         tmp_real_fd1, &
!!$         Int(1,pd_mpi_ik), MPI_Real8, &
!!$         status_mpi, ierr)
    
    !call add_leaf_to_branch(res_tree,"Rotation Vector CR_2", 3_pd_ik, n)
!!$    Call MPI_FILE_WRITE_AT(FH_MPI(5), &
!!$         Int(root%branches(4)%leaves(15)%lbound-1+(nn-1)*3, MPI_OFFSET_KIND), &
!!$         n, &
!!$         Int(3,pd_mpi_ik), MPI_Real8, &
!!$         status_mpi, ierr)
    
    !=========================================================================
    !== Inlining of EE =======================================================

    aa = rot_alg(n,alpha)
    BB = tra_R6(aa)
    EE = matmul(matmul(transpose(BB),EE),BB)

    If (out_amount /= "PRODUCTION" ) &
         call write_real_matrix(un_lf,EE,6_ik,6_ik,"Backrotated anisotropic stiffness CR_2")

    !=========================================================

    If ( (EE(1,1) < EE(2,2)) .AND.  &
         (EE(1,1) < EE(3,3)) .AND.  (EE(2,2) < EE(3,3))         ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"123"
       continue
       
    Else If ( (EE(1,1) < EE(2,2)) .AND.  &
         (EE(1,1) < EE(3,3)) .AND.  (EE(2,2) > EE(3,3)) ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"132"

       !** 132 => 123 *****************
       n = aa(:,1)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

    Else If ( (EE(1,1) < EE(2,2)) .AND.  &
         (EE(1,1) > EE(3,3)) .AND.  (EE(2,2) > EE(3,3)) ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"231"

       !** 231 => 132 *****************
       n = aa(:,2)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)

       !** 132 => 123 *****************
       n = aa(:,1)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

    Else If ( (EE(1,1) > EE(2,2)) .AND.  &
         (EE(1,1) < EE(3,3)) .AND.  (EE(2,2) < EE(3,3)) ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"213"

       !** 213 => 123 *****************
       n = aa(:,3)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

    Else If ( (EE(1,1) > EE(2,2)) .AND.  &
         (EE(1,1) > EE(3,3)) .AND.  (EE(2,2) < EE(3,3)) ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"312"

       !** 312 => 132 *****************
       n = aa(:,3)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

       !** 132 => 123 *****************
       n = aa(:,1)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)  

    Else If ( (EE(1,1) > EE(2,2)) .AND.  &
         (EE(1,1) > EE(3,3)) .AND.  (EE(2,2) > EE(3,3)) ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"321"

       !** 321 => 123 *****************
       n = aa(:,2)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)  

    End If

    BB = tra_R6(aa)
    EE = matmul(matmul(transpose(BB),EE_Orig),BB)

    If (out_amount /= "PRODUCTION" ) then
       call write_real_matrix(un_lf, aa,3_ik,3_ik,"Final coordinate system CR_2")
       call write_real_matrix(un_lf, EE,6_ik,6_ik,"Inlined anisotropic stiffness CR_2")
    End If
    
    call check_sym(A=EE,name="Inlined anisotropic stiffness CR_2")

    !Call add_leaf_to_branch(res_tree,"Final coordinate system CR_2", 9_pd_ik, &
    !    reshape(aa,[9_pd_ik]))
!!$    Call MPI_FILE_WRITE_AT(FH_MPI(5), &
!!$         Int(root%branches(4)%leaves(16)%lbound-1+(nn-1)*9, MPI_OFFSET_KIND), &
!!$         reshape(aa,[9_pd_ik]), &
!!$         Int(9,pd_mpi_ik), MPI_Real8, &
!!$         status_mpi, ierr)
    
    !call add_leaf_to_branch(res_tree,"Optimized Effective stiffness CR_2", 36_pd_ik, &
    !     reshape(EE,[36_pd_ik]))
!!$    Call MPI_FILE_WRITE_AT(FH_MPI(5), &
!!$         Int(root%branches(4)%leaves(17)%lbound-1+(nn-1)*36, MPI_OFFSET_KIND), &
!!$         reshape(EE,[36_pd_ik]), &
!!$         Int(36,pd_mpi_ik), MPI_Real8, &
!!$         status_mpi, ierr)
        
    Select Case (timer_level)
    Case (3)
       call end_timer("  +-- Back rotation of material matrix "//trim(nn_char))
    Case (2)
       call end_timer("  +-- Back rotation of material matrix "//trim(nn_char))
    Case default
       continue
    End Select

1000 Continue

    deallocate(nodes,uu,edat,cref_cnodes,rforces)

  End subroutine calc_effective_material_parameters

  !****************************************************************************
  subroutine init_loadcase(eps,vv)

    Real(Kind=rk), intent(in)     :: eps
    Real(Kind=rk), Dimension(:,:) :: vv

    Real(Kind=rk)                 :: eps2,eps3,eps4

    eps2 = eps*2._rk
    eps3 = eps*3._rk
    eps4 = eps*4._rk

    vv = 0._rk

    vv( 1,:) = [ 0._rk , 0._rk , 0._rk ,  eps  ,  eps  , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk]
    vv( 2,:) = [  eps  , 0._rk , 0._rk ,  eps  ,  eps  , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk]
    vv( 3,:) = [  eps  , 0._rk , 0._rk , 0._rk ,  eps  , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk]
    vv( 4,:) = [ 0._rk , 0._rk , 0._rk , 0._rk ,  eps  , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk ,  eps  , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk]
    vv( 5,:) = [ 0._rk , 0._rk , 0._rk ,  eps  , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk ,  eps  , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk]
    vv( 6,:) = [  eps  , 0._rk , 0._rk ,  eps  , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk ,  eps  , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk]
    vv( 7,:) = [  eps  , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 eps  , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk]
    vv( 8,:) = [ 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk ,  eps  , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk ,  eps  , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk]
    vv( 9,:) = [ 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk ,  eps  , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk ,  eps  , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk]
    vv(10,:) = [ 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 eps  , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk ,  eps  , 0._rk , 0._rk , 0._rk , 0._rk]
    vv(11,:) = [ 0._rk ,  eps  , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk ,  eps  , 0._rk , 0._rk , 0._rk]
    vv(12,:) = [ 0._rk ,  eps  , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk ,  eps  , 0._rk , 0._rk]
    vv(13,:) = [ 0._rk , 0._rk , 0._rk , 0._rk , 0._rk ,  eps  , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk ,  eps  , 0._rk]
    vv(14,:) = [ 0._rk , 0._rk , 0._rk , 0._rk , 0._rk ,  eps  , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk ,  eps ]
    vv(15,:) = [ 0._rk ,  eps  , 0._rk , 0._rk , 0._rk ,  eps  , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , -eps ]
    vv(16,:) = [ 0._rk ,  eps  , 0._rk , 0._rk , 0._rk ,  eps  , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk]
    vv(17,:) = [ 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk ,  eps  , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk]
    vv(18,:) = [ 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk ,  eps  , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk]
    vv(19,:) = [ 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk ,  eps  , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk]
    vv(20,:) = [ 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk ,  eps  , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk]
    vv(21,:) = [ 0._rk , 0._rk ,  eps  , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk]
    vv(22,:) = [ 0._rk , 0._rk ,  eps  , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , -eps2 , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk]
    vv(23,:) = [ 0._rk , 0._rk ,  eps  , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , -eps3 , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk]
    vv(24,:) = [ 0._rk , 0._rk ,  eps  , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , 0._rk , &
                 0._rk , 0._rk , 0._rk , -eps4 , 0._rk , 0._rk , 0._rk , 0._rk]

  end subroutine init_loadcase

End Module calcmat

  
