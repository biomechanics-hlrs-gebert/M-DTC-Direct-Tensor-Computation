!==============================================================================
!> This Program calculates the effective numerical stiffness of
!> microstructured Volume elements
!>
Module calcmat

USE tensors
USE decomp
USE mat_matrices
USE mechanical
USE chain_routines
USE auxiliaries
USE user_interaction
USE formatted_plain
USE linFE
USE mpi_system
USE mpi

implicit none

contains

subroutine calc_effective_material_parameters(root, comm_nn, ddc_nn, &
     fh_mpi_worker, size_mpi, comm_mpi, collected_logs)

Type(tBranch)     , Intent(InOut) :: root
INTEGER(mik) , Intent(In) :: size_mpi, comm_mpi
integer(ik)  , Intent(in) :: ddc_nn, comm_nn
Integer(mik), Dimension(no_streams), Intent(in) :: fh_mpi_worker

Real(rk) :: div_10_exp_jj, eff_density, n12, n13, n23, alpha, phi, eta
Real(rk) :: cos_alpha, sin_alpha, One_Minus_cos_alpha, sym

Real(rk), Dimension(:)    , allocatable :: tmp_nn, delta, x_D_phy
Real(rk), Dimension(:,:)  , allocatable :: nodes, vv, ff, stiffness
Real(rk), Dimension(:,:,:), allocatable :: calc_rforces, uu, rforces, edat, crit_1, crit_2
Real(rk), Dimension(1)    :: tmp_real_fd1
Real(rk), Dimension(3)    :: min_c, max_c, n
Real(rk), Dimension(6)    :: ro_stress
Real(rk), Dimension(8)    :: tmp_r8 
Real(rk), Dimension(12)   :: tmp_r12
Real(rk), Dimension(3,3)  :: aa
Real(rk), Dimension(6,6)  :: ee_orig, BB, CC, cc_mean, EE, fv,meps
Real(rk), Dimension(0:16) :: crit_min
Real(rk), Dimension(:,:), ALLOCATABLE :: int_strain, int_stress
Real(rk):: E_Modul, nu, rve_strain, v_elem, v_cube

Integer(mik), Dimension(MPI_STATUS_SIZE) :: status_mpi
Integer(mik) :: ierr, global_rank_mpi

integer(ik) :: ii, jj, kk, ll, no_elem_nodes, micro_elem_nodes, no_lc, num_leaves, alloc_stat, &
     no_elems, no_nodes, no_cnodes, macro_order, ii_phi, ii_eta, kk_phi, kk_eta, mem_global, status_global

Integer(ik), Dimension(:,:,:,:), Allocatable :: ang
Integer(ik), Dimension(:)      , Allocatable :: xa_n, xe_n, no_cnodes_pp, cref_cnodes
Integer(ik), Dimension(3)                    :: s_loop,e_loop, mlc
INTEGER(ik), DIMENSION(24) :: collected_logs ! timestamps, memory_usage, pid_returned

Logical :: success

Character(*), Parameter :: link_name="struct_calcmat_fmps"
Character(9)   :: nn_char
Character(mcl) :: desc

Type(tBranch), Pointer :: ddc, loc_ddc, meta_para, domain_branch, mesh_branch, result_branch

Type(tLeaf),  pointer :: node_leaf_pointer
Type(tLeaf), Allocatable, Dimension(:)   :: leaf_list

CALL MPI_COMM_RANK(MPI_COMM_WORLD, global_rank_mpi, ierr)
CALL print_err_stop(std_out, "MPI_COMM_RANK couldn't be retrieved", ierr)

write(nn_char,'(I0)') ddc_nn

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

! Load input meta_para
Call Search_branch("Input parameters", root, meta_para, success)

Call pd_get(meta_para, "Young_s modulus", E_Modul)
Call pd_get(meta_para, "Poisson_s ratio", nu)
Call pd_get(meta_para, "Average strain on RVE", rve_strain)
Call pd_get(meta_para, "Element order on macro scale", macro_order)

!------------------------------------------------------------------------------
! Load domain decomposition parameters
!------------------------------------------------------------------------------

! Get global DDC parameters from root 
Call Search_branch("Global domain decomposition", root, ddc, success)

call pd_get(ddc, "delta", delta)
call pd_get(ddc, "x_D_phy", x_D_phy)

! Get local DDC parameters from root
Write(desc,'(A,I0)') "Local domain Decomposition of domain no ", ddc_nn

Call Search_branch(trim(desc), root, loc_ddc, success)
call pd_get(loc_ddc, "xa_n", xa_n)
call pd_get(loc_ddc, "xe_n", xe_n)

min_c = Real(xa_n - 1,rk) * delta
max_c = Real(xe_n    ,rk) * delta

!------------------------------------------------------------------------------
! Set node number of macro element
!------------------------------------------------------------------------------
If (macro_order == 1) then
    no_elem_nodes = 8
    no_lc = 24
ELSE IF (macro_order == 2) THEN
     no_elem_nodes = 20
     no_lc = 60
Else
    CALL print_err_stop(std_out, "Element orders other than 1 or 2 are not supported", 1)
End If

IF (.NOT. ALLOCATED(int_strain)) ALLOCATE(int_strain(6,no_lc))
IF (.NOT. ALLOCATED(int_stress)) ALLOCATE(int_stress(6,no_lc))

!------------------------------------------------------------------------------
! Get global mesh_branch parameters
! project_name already is Rank_xxxxy
!------------------------------------------------------------------------------
Call Search_branch("Domain "//trim(nn_char), root, domain_branch, success)

desc = ''
Write(desc,'(A,I0)') 'Mesh info of '//trim(project_name)//'_', ddc_nn
Call search_branch(trim(desc), domain_branch, mesh_branch, success)

! DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Braucht man das noch?
! call log_tree(mesh_branch, un_lf, .FALSE.)
! DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

call pd_get(mesh_branch, 'No of nodes in mesh',  no_nodes)
call pd_get(mesh_branch, 'No of elements in mesh',  no_elems)

call pd_get(mesh_branch, 'No of cdofs in parts', no_cnodes_pp)
no_cnodes=sum(no_cnodes_pp)/3

If (out_amount /= "PRODUCTION") then
    write(un_lf,*)
    write(un_lf,FMT_MSG_xAI0)"Read no_nodes  from domain branch: ", no_nodes
    write(un_lf,FMT_MSG_xAI0)"Read no_cnodes from domain branch: ", no_cnodes
End If

!------------------------------------------------------------------------------
! Allocate fields and read data
!------------------------------------------------------------------------------

! Copy Coordinate list nodes(3,no_nodes)
Allocate (nodes(3,no_nodes),stat=alloc_stat)
call pd_get(mesh_branch, "Coordinates", node_leaf_pointer)
nodes = reshape( node_leaf_pointer%p_real8 ,[3,no_nodes])

!------------------------------------------------------------------------------
! Build Cross reference for constrained nodes from "Boundary_Ids"
!------------------------------------------------------------------------------
Allocate(cref_cnodes(no_cnodes),stat=alloc_stat)
call alloc_err( "cref_cnodes",alloc_stat)

Do ii = 1, mesh_branch%branches(2)%branches(1)%leaves(1)%dat_no, 9
    cref_cnodes((ii+8)/9) = mesh_branch%branches(2)%branches(1)%leaves(1)%p_int8(ii)
end Do

!------------------------------------------------------------------------------
! Read displacement results
!------------------------------------------------------------------------------
Allocate(uu(3, no_nodes, no_lc), stat=alloc_stat)
call alloc_err("uu", alloc_stat)

call get_leaf_list ("Displacements", mesh_branch, num_leaves, leaf_list)
If (out_amount /= "PRODUCTION" ) then
    write(un_lf,FMT_MSG_AxI0) "Number of leaves with desc = Displacements: ", num_leaves
End If

Do ii = 1, num_leaves
    uu(:,:,ii) = reshape(leaf_list(ii)%p_real8, [3, no_nodes])
End Do

Deallocate(leaf_list)

!------------------------------------------------------------------------------
! Read stress and strain results
!------------------------------------------------------------------------------
Allocate(edat(15, no_lc, no_elems), stat=alloc_stat)
call alloc_err( "edat",alloc_stat)

call get_leaf_list("Avg. Element Data", mesh_branch, num_leaves, leaf_list)
If (out_amount /= "PRODUCTION" ) then
    write(un_lf,FMT_MSG_xAI0)"Number of leaves with desc = Avg. Element Data: ", num_leaves
End If

Do ii = 1, num_leaves
    edat(:,ii,:) = reshape(leaf_list(ii)%p_real8, [15, no_elems])
End Do

Deallocate(leaf_list)

!------------------------------------------------------------------------------
! Read Reaction Forces
!------------------------------------------------------------------------------
Allocate(rforces(3,no_nodes,no_lc),stat=alloc_stat)
call alloc_err("rforces", alloc_stat)

call get_leaf_list("Reaction Forces", mesh_branch, num_leaves, leaf_list)
If (out_amount /= "PRODUCTION" ) then
    write(un_lf,FMT_MSG_AxI0) "Number of leaves with desc = Reaction Forces: ", num_leaves
End If

Do ii = 1, num_leaves
    rforces(:,:,ii) = reshape(leaf_list(ii)%p_real8, [3,no_nodes])
End Do

Deallocate(leaf_list)

Allocate(calc_rforces(3, no_nodes, no_lc), stat=alloc_stat)
call alloc_err("calc_rforces", alloc_stat)
calc_rforces = 0._rk

allocate(vv(no_lc, no_lc), stat=alloc_stat)
call alloc_err("vv", alloc_stat)
allocate(ff(no_lc, no_lc), stat=alloc_stat)
call alloc_err("ff", alloc_stat)

ff = 0._rk

allocate(stiffness(no_lc, no_lc), stat=alloc_stat)
call alloc_err("stiffness", alloc_stat)

allocate(tmp_nn(no_elem_nodes), stat=alloc_stat)
call alloc_err("tmp_nn", alloc_stat)

If (out_amount == "DEBUG") THEN
    WRITE(un_lf, FMT_DBG_SEP)
    write(un_lf, *)
    write(un_lf, *)'min uu      = ', minval(uu), 'max uu      = ', maxval(uu)
    write(un_lf, *)'min rforces = ', minval(rforces), 'max rforces = ', maxval(rforces)
    write(un_lf, *)
    WRITE(un_lf, FMT_DBG_SEP)
end if

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

!------------------------------------------------------------------------------
! Search effective results branch
!------------------------------------------------------------------------------
Call Search_branch("Results of domain "//nn_char, root, result_branch, success)

!------------------------------------------------------------------------------
! Init C
!------------------------------------------------------------------------------
CC = iso_compliance_voigt(E_Modul, nu)

If (out_amount == "DEBUG") THEN
    WRITE(un_lf,FMT_DBG_SEP)
    Call Write_matrix(un_lf, "Effective Isotropic Compliance -- CC", cc, fmti='std', unit='1/MPa')
    WRITE(un_lf,FMT_DBG_SEP)
End if

!------------------------------------------------------------------------------
! Reorder stresses and strains
!------------------------------------------------------------------------------
Do jj = 1, no_elems
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

!------------------------------------------------------------------------------
! Init element and RVE volume 
!------------------------------------------------------------------------------
v_elem = delta(1)   * delta(2)   * delta(3)
v_cube = x_D_phy(1) * x_D_phy(2) * x_D_phy(3)

If (out_amount /= "PRODUCTION" ) then
    Write(un_lf, "(A,E15.6)") 'Volume per element: ', v_elem
    Write(un_lf, "(A,E15.6)") 'MVE Volume:         ', v_cube
    Write(un_lf, *)
End If

!------------------------------------------------------------------------------
! Loadcase init
!------------------------------------------------------------------------------
If (macro_order == 1) then
     call init_loadcase_el_order_lin(rve_strain, vv)
     
ELSE IF (macro_order == 2) THEN
     call init_loadcase_el_order_quad(rve_strain, vv)

Else
     CALL print_err_stop(std_out, "Element orders other than 1 or 2 are not supported", 1)
 End If

If (out_amount == "DEBUG") THEN
    WRITE(un_lf, FMT_DBG_SEP)
    Call Write_matrix(un_lf, "Displacement matrix", vv, fmti='std', unit='mm')
    WRITE(un_lf, FMT_DBG_SEP)
End if

!------------------------------------------------------------------------------
! Calculate inverse of displacement matrix with LinPack
!------------------------------------------------------------------------------
Call inverse(vv, no_lc, un_lf)

If (out_amount == "DEBUG") THEN
    WRITE(un_lf, FMT_DBG_SEP)
    Call Write_matrix(un_lf, "Inverted displacement matrix", vv, fmti='std', unit='1/mm')
    WRITE(un_lf, FMT_DBG_SEP)
End if

!------------------------------------------------------------------------------
! calc load matrix
!------------------------------------------------------------------------------
ff = 0._rk

Do ii = 1, no_lc                         ! Cycle through all load cases
    Do jj = 1, no_nodes                  ! Cycle through all boundary nodes

          ! t_geom_xi transforms coordinates from geometry to xi space 
          ! Result(phi_nn) :  Real(rk), dimension(8)
          IF (macro_order == 1) THEN
               tmp_nn = phi_NN_hexe8(t_geom_xi(nodes(:,jj),min_c,max_c))
          ELSE IF (macro_order == 2) THEN
               tmp_nn = phi_NN_hexe20(t_geom_xi(nodes(:,jj),min_c,max_c))
          END IF 

          Do kk = 1,3                         ! Dof per node
                    Do ll = 1, no_elem_nodes  ! No of FE nodes per Macro element
                    ! Setup of load matrix
                    ! 1st index - rows : Acumulate the reaction forces of the 
                    !                    micro model via the trail functions to 
                    !                    the macro element
                    ! 2nd index - cols : Do first index step for all load cases
                    ff((kk-1)*no_elem_nodes+ll,ii) = ff((kk-1)*no_elem_nodes+ll,ii) + rforces(kk,jj,ii)*tmp_nn(ll)
                    End Do
          End Do

    End Do
End Do

If (out_amount /= "PRODUCTION" ) then
    Call Write_matrix(un_lf, "ff by summation of single force formula", ff, fmti='std', unit='N')
End If

!------------------------------------------------------------------------------
! Domain forces
!------------------------------------------------------------------------------
CALL add_leaf_to_branch(result_branch, "Domain forces", no_lc*no_lc, reshape(ff,[no_lc*no_lc]))
CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
     Int(root%branches(3)%leaves(7)%lbound-1+(comm_nn-1)*no_lc*no_lc, MPI_OFFSET_KIND), &
     reshape(ff,[no_lc*no_lc]), &
     Int(no_lc*no_lc,pd_mik), MPI_Real8, &
     status_mpi, ierr)

If (out_amount /= "PRODUCTION" ) then
     Call Write_matrix(un_lf, "Domain forces", ff, fmti='std')
End If 

!------------------------------------------------------------------------------
! Calc effective nummerical stiffness
!------------------------------------------------------------------------------
stiffness = matmul(ff,vv)

CALL add_leaf_to_branch(result_branch, "Effective numerical stiffness", no_lc*no_lc, &
    reshape(stiffness,[no_lc*no_lc]))
CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
    Int(root%branches(3)%leaves(8)%lbound-1+(comm_nn-1)*no_lc*no_lc, MPI_OFFSET_KIND), &
    reshape(stiffness, [no_lc*no_lc]), &
    Int(no_lc*no_lc,pd_mik), MPI_Real8, &
    status_mpi, ierr)

If (out_amount /= "PRODUCTION" ) then
    Call Write_matrix(std_out, "Stiffness", stiffness, fmti='std', unit='MPa')
End If


!------------------------------------------------------------------------------
! Calc Symmetry deviation - effective numerical stiffness
!------------------------------------------------------------------------------
sym = check_sym(stiffness)
tmp_real_fd1 = sym

CALL add_leaf_to_branch(result_branch, &
    "Symmetry deviation - effective numerical stiffness", 1_pd_ik, tmp_real_fd1)
CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
    Int(root%branches(3)%leaves(9)%lbound-1+(comm_nn-1), MPI_OFFSET_KIND), &
    tmp_real_fd1, &
    1_pd_mik, MPI_REAL8, status_mpi, ierr)

!------------------------------------------------------------------------------
! Calc inverse of effective nummerical stiffness for cross check against
! dislacement matrix
! Call inverse(stiffness, no_lc, un_lf)
! Call Write_real_matrix(un_lf,stiffness,no_lc,no_lc,"Inverse Stiffness")
! Call Write_real_matrix(un_lf,matmul(stiffness,ff(:,1)), no_lc, 1_ik, &
!     "Inverse Stiffness o f(:,1)")
!------------------------------------------------------------------------------

If (macro_order == 1) then
    !------------------------------------------------------------------------------
    ! Calc consistent force matrix
    ! 8x 8y 8z
    !------------------------------------------------------------------------------
    Do ii = 1, 6
        fv(1,ii) = (ff( 2,ii) + ff( 3,ii) + ff( 6,ii) + ff( 7,ii)) / ( x_D_phy(2) * x_D_phy(3) ) ! X to X | Node 2,3,6,7
        fv(2,ii) = (ff(11,ii) + ff(12,ii) + ff(15,ii) + ff(16,ii)) / ( x_D_phy(1) * x_D_phy(3) ) ! Y to Y | Node 3,4,7,8
        fv(3,ii) = (ff(21,ii) + ff(22,ii) + ff(23,ii) + ff(24,ii)) / ( x_D_phy(1) * x_D_phy(2) ) ! Z to Z | Node 5,6,7,8
    
        fv(4,ii) = 0.5_rk * &
            ( (ff( 3,ii) + ff( 4,ii) + ff( 7,ii) + ff( 8,ii)) / ( x_D_phy(1) * x_D_phy(3) ) + &  ! Y to X  | Node 3,4,7,8 
            (ff(10,ii) + ff(11,ii) + ff(14,ii) + ff(15,ii)) / ( x_D_phy(2) * x_D_phy(3) ) )      ! X to Y  | Node 2,3,6,7
        fv(5,ii) = 0.5_rk * &
            ( (ff( 5,ii) + ff( 6,ii) + ff( 7,ii) + ff( 8,ii)) / ( x_D_phy(1) * x_D_phy(2) ) + &  ! Z to X | Node 5,6,7,8
            (ff(18,ii) + ff(19,ii) + ff(22,ii) + ff(23,ii)) / ( x_D_phy(2) * x_D_phy(3) ) )      ! X to Z | Node 2,3,6,7 
        fv(6,ii) = 0.5_rk * &
            ( (ff(13,ii) + ff(14,ii) + ff(15,ii) + ff(16,ii)) / ( x_D_phy(1) * x_D_phy(2) ) + &  ! Z to Y | Node 5,6,7,8
            (ff(19,ii) + ff(20,ii) + ff(23,ii) + ff(24,ii)) / ( x_D_phy(1) * x_D_phy(3) ) )      ! Y to Z | Node 3,4,7,8
    End Do

    !  1  2  3  4  5  6  7  8   |   1  2  3  4  5  6  7  8  |   1  2  3  4  5  6  7  8
    !  1  2  3  4  5  6  7  8   |   9 10 11 12 13 14 15 16  |  17 18 19 20 21 22 23 24
    !  x  x  x  x  x  x  x  x   |   y  y  y  y  y  y  y  y  |   z  z  z  z  z  z  z  z       
ELSE IF (macro_order == 2) THEN
     !------------------------------------------------------------------------------
     ! Calc consistent force matrix
     !------------------------------------------------------------------------------
     Do ii = 1, 6
        fv(1,ii) = (ff( 2,ii) + ff( 3,ii) + ff( 6,ii) + ff( 7,ii) + &
                ff(10,ii) + ff(14,ii) + ff(18,ii) + ff(19,ii))  / ( x_D_phy(2) * x_D_phy(3) ) ! X to X | Node 2,3,6,7,10,14,18,19
                
        fv(2,ii) = (ff(23,ii) + ff(24,ii) + ff(27,ii) + ff(28,ii) + &
                ff(31,ii) + ff(35,ii) + ff(39,ii) + ff(40,ii))  / ( x_D_phy(1) * x_D_phy(3) ) ! Y to Y | Node 3,4,7,8,11,15,19,20
                
        fv(3,ii) = (ff(45,ii) + ff(46,ii) + ff(47,ii) + ff(48,ii) + &
                ff(53,ii) + ff(54,ii) + ff(55,ii) + ff(56,ii))  / ( x_D_phy(1) * x_D_phy(2) ) ! Z to Z | Node 5,6,7,8,13,14,15,16

        fv(4,ii) = 0.5_rk * &
            ( (ff( 3,ii) + ff( 4,ii) + ff( 7,ii) + ff( 8,ii) + ff(11,ii) + ff(15,ii) + ff(19,ii) + ff(20,ii)) / &
            ( x_D_phy(1) * x_D_phy(3)) + & ! Y to X | Node 3,4,7,8,11,15,19,20
            (ff(22,ii) + ff(23,ii) + ff(26,ii) + ff(27,ii) + ff(30,ii) + ff(34,ii) + ff(38,ii) + ff(39,ii)) / &
            ( x_D_phy(2) * x_D_phy(3) ))     ! X to Y | Node 2,3,6,7,10,14,18,19

        fv(5,ii) = 0.5_rk * &
            ( (ff( 5,ii) + ff( 6,ii) + ff( 7,ii) + ff( 8,ii) + ff(13,ii) + ff(14,ii) + ff(15,ii) + ff(16,ii)) / &
            ( x_D_phy(1) * x_D_phy(2)) + & ! Z to X | Node 5,6,7,8,13,14,15,16

            (ff(42,ii) + ff(43,ii) + ff(46,ii) + ff(47,ii) + ff(50,ii) + ff(54,ii) + ff(58,ii) + ff(59,ii)) / &
            ( x_D_phy(2) * x_D_phy(3) ))     ! X to Z | Node 2,3,6,7,10,14,18,19

        fv(6,ii) = 0.5_rk * &
            ( (ff(25,ii) + ff(26,ii) + ff(27,ii) + ff(28,ii) + ff(33,ii) + ff(34,ii) + ff(35,ii) + ff(36,ii)) / &
            ( x_D_phy(1) * x_D_phy(2)) + & ! Z to Y | Node 5,6,7,8,13,14,15,16
            (ff(43,ii) + ff(44,ii) + ff(47,ii) + ff(48,ii) + ff(51,ii) + ff(55,ii) + ff(59,ii) + ff(60,ii)) / &
            ( x_D_phy(1) * x_D_phy(3) ))     ! Y to Z | Node 3,4,7,8,11,15,19,20
     
    !  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20  |  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20  |  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
    !  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20  | 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40  | 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60
    !  x  x  x  x  x  x  x  x  x  x  x  x  x  x  x  x  x  x  x  x  |  y  y  y  y  y  y  y  y  y  y  y  y  y  y  y  y  y  y  y  y  |  z  z  z  z  z  z  z  z  z  z  z  z  z  z  z  z  z  z  z  z

    End Do

End If



If (out_amount /= "PRODUCTION" ) then
    Call Write_matrix(un_lf, "Konsistent force matrix", fv, fmti='std', unit='N')
End If

!------------------------------------------------------------------------------
! Calc integrated force matrix
! ------------------------------------------------------------------------------
Do jj = 1, 6
  Do ii = 1, 6
     fv(ii,jj) = sum(edat(ii,jj,:))
  End do
End Do

fv = fv * v_elem / v_cube
If (out_amount /= "PRODUCTION" ) then
  CALL write_matrix(un_lf, "Integrated force matrix", fv)
End If



!------------------------------------------------------------------------------
! Calc averaged strains and stresses 
!------------------------------------------------------------------------------
Do jj = 1, no_lc 
    Do ii = 1,6
       int_stress(ii,jj) = sum(edat(ii,jj,:))
    End do
    Do ii = 7, 12
      int_strain(ii-6,jj) = sum(edat(ii,jj,:))
    End Do
End Do

int_strain = int_strain * v_elem / v_cube
int_stress = int_stress * v_elem / v_cube

If (out_amount /= "PRODUCTION" ) then
  write(un_lf,*)
  write(un_lf,"(150('='))")
  Write(un_lf,"(A)")"--------- Averaged strains of loadcases ---------"
  Write(un_lf,"(7(A20))")'-- ii --','-- E11 --','-- E22 --','-- E33 --','-- E12 --','-- E13 --','-- E23 --'
  Do ii = 1, 24
     write(un_lf,"(I20,6(E20.9))")ii,int_strain(:,ii)
  End Do

  write(un_lf,*)
  write(un_lf,"(150('='))")
  Write(un_lf,"(A)")"--------- Averaged stresses of loadcases ---------"
  Write(un_lf,"(7(A20))")'-- ii --','-- S11 --','-- S22 --','-- S33 --','-- S12 --','-- S13 --','-- S23 --'
  Do ii = 1, 24
     write(un_lf,"(I20,6(E20.9))")ii,int_stress(:,ii)
  End Do
End If

!------------------------------------------------------------------------------
! Averaged stresses
!------------------------------------------------------------------------------
CALL add_leaf_to_branch(result_branch, "Averaged stresses", 6*no_lc, reshape(int_stress,[6*no_lc]))
CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
    Int(root%branches(3)%leaves(10)%lbound-1+(comm_nn-1)*6*no_lc, MPI_OFFSET_KIND), &
    reshape(int_stress, [6*no_lc]), &
    Int(6*no_lc, pd_mik), MPI_REAL8, status_mpi, ierr)

!------------------------------------------------------------------------------
! Averaged strains
!------------------------------------------------------------------------------
CALL add_leaf_to_branch(result_branch, "Averaged strains",  6*no_lc, reshape(int_strain,[6*no_lc]))
CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
    Int(root%branches(3)%leaves(11)%lbound-1+(comm_nn-1)*6*no_lc, MPI_OFFSET_KIND), &
    reshape(int_strain,[6*no_lc]), &
    Int(6*no_lc, pd_mik), MPI_REAL8, status_mpi, ierr)

!------------------------------------------------------------------------------
! Calc integrated strain matrix for first 6 loadcases
!meps = int_strain(1:6,1:6)
!If (out_amount /= "PRODUCTION" ) then
!   Call Write_real_matrix(un_lf, meps,6_ik, 6_ik, &
!        "Integrated strain matrix of first 6 loadcases")
!End If

!------------------------------------------------------------------------------
! Calc theoretical effective stiffness
! Linear loadcases with constant strain fields
!------------------------------------------------------------------------------
meps(1,:) = [ 1.00E+06,  0.00E+00,  0.00E+06,  0.00E+06,  0.00E+00,  0.00E+06 ]
meps(2,:) = [ 0.00E+06,  1.00E+06,  0.00E+06,  0.00E+06,  0.00E+00,  0.00E+06 ]
meps(3,:) = [ 0.00E+06,  0.00E+00,  1.00E+06,  0.00E+06,  0.00E+06,  0.00E+06 ]
meps(4,:) = [ 0.00E+06,  0.00E+00,  0.00E+06, -1.00E+06,  0.00E+00,  0.00E+06 ]
meps(5,:) = [ 0.00E+06,  0.00E+06,  0.00E+06,  0.00E+06, -1.00E+06,  0.00E+06 ]
meps(6,:) = [ 0.00E+06,  0.00E+00,  0.00E+06,  0.00E+00,  0.00E+06,  1.00E+06 ]

cc_mean = matmul(fv,meps)

If (out_amount /= "PRODUCTION" ) then
    Call Write_matrix(un_lf, "Effective stiffness", cc_mean, fmti='std', unit='MPa')
End If

!------------------------------------------------------------------------------
! Effective stiffness
!------------------------------------------------------------------------------
CALL add_leaf_to_branch(result_branch, "Effective stiffness", 36_pd_ik, reshape(cc_mean, [36_pd_ik]))
CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
    Int(root%branches(3)%leaves(12)%lbound-1+(comm_nn-1)*36, MPI_OFFSET_KIND), &
    reshape(cc_mean,[36]), &
    36_pd_mik, MPI_Real8, status_mpi, ierr)

!------------------------------------------------------------------------------
! Symmetry deviation - effective stiffness
!------------------------------------------------------------------------------
sym = check_sym(cc_mean)
tmp_real_fd1 = sym

CALL add_leaf_to_branch(result_branch, "Symmetry deviation - effective stiffness",  1_pd_ik, tmp_real_fd1)
CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
    Int(root%branches(3)%leaves(13)%lbound-1+(comm_nn-1), MPI_OFFSET_KIND), &
    tmp_real_fd1, &
    1_pd_mik, MPI_REAL8, status_mpi, ierr)

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

!------------------------------------------------------------------------------
! Back rotation of material matrix according to different criterias
!------------------------------------------------------------------------------

!EE = matmul(fv,transpose(meps))
EE = (cc_mean + transpose(cc_mean)) / 2._rk

If (out_amount /= "PRODUCTION" ) then
    Call Write_matrix(un_lf, "Averaged Effective stiffness", ee, fmti='std', unit='MPa')
End If

!------------------------------------------------------------------------------
! Symmetry deviation - effective stiffness
!------------------------------------------------------------------------------
CALL add_leaf_to_branch(result_branch, "Averaged Effective stiffness",  36_pd_ik, reshape(ee,[36_pd_ik]))
CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
    Int(root%branches(3)%leaves(14)%lbound-1+(comm_nn-1)*36, MPI_OFFSET_KIND), &
    reshape(ee,[36]), &
    36_pd_mik, MPI_REAL8, status_mpi, ierr)

!------------------------------------------------------------------------------
! Symmetry deviation - effective stiffness
!------------------------------------------------------------------------------
tmp_real_fd1 = check_sym(ee)

CALL add_leaf_to_branch(result_branch, &
    "Symmetry deviation - Averaged effective stiffness",  1_pd_ik, tmp_real_fd1)
CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
    Int(root%branches(3)%leaves(15)%lbound-1+(comm_nn-1), MPI_OFFSET_KIND), &
    tmp_real_fd1, &
    1_pd_mik, MPI_REAL8, status_mpi, ierr)

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
       write(un_lf,FMT_MSG_AxI0)'Initial Minloc  CR_1 : ',mlc
       write(un_lf,FMT_MSG_xAF0) 'Initial Minimum CR_1 : ',crit_min(0)
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
          write(un_lf,FMT_MSG_xAI0) 'Iteration            : ',jj
          write(un_lf,FMT_MSG_AxI0)'Loop start           : ',s_loop
          write(un_lf,FMT_MSG_AxI0)'Loop end             : ',e_loop
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
          write(un_lf, FMT_MSG_xAF0)'Minimum CR_1         : ',crit_min(jj)
       End If
       
       If ( (abs(crit_min(jj-1)-crit_min(jj)) < num_zero) .OR. (jj >= 16)) Exit

       jj = jj + 1

    End Do

    ! Be aware that minloc starts off at field index 1 !!!
    mlc = minloc(crit_1(0:kk-1,0:kk_phi-1,0:kk_eta-1))-1

    alpha = Real( ang(1,mlc(1),mlc(2),mlc(3)),rk ) * pi / (180._rk*(10._rk**jj-1))
    phi   = Real( ang(2,mlc(1),mlc(2),mlc(3)),rk ) * pi / (180._rk*(10._rk**jj-1))
    eta   = Real( ang(3,mlc(1),mlc(2),mlc(3)),rk ) * pi / (180._rk*(10._rk**jj-1))

    n = [cos(phi)*sin(eta) , sin(phi)*sin(eta) , cos(eta) ] 
    n = n / sqrt(sum(n*n))

    If (out_amount /= "PRODUCTION" ) then
       write(un_lf,*)
       Write(un_lf,FMT_MSG_xAI0) "Solution converged after : ",jj," iterations"
       Write(un_lf,FMT_MSG_AxF0) "With final citerion 1    : ",&
            minval(crit_1(0:kk-1,0:kk_phi-1,0:kk_eta-1)),crit_1(mlc(1),mlc(2),mlc(3))
       Write(un_lf,FMT_MSG_xAF0)  "With final epsilon       : ", crit_min(jj-1)-crit_min(jj)
       Write(un_lf,FMT_MSG_xAF0) "Final rotation angle  is : ", alpha
       Write(un_lf,FMT_MSG_AxF0) "Final rotation vector is : ", n
       Write(un_lf,*)
    End If
    
    !------------------------------------------------------------------------------
    ! Rotation Angle CR_1
    !------------------------------------------------------------------------------
    tmp_real_fd1 = alpha 

    CALL add_leaf_to_branch(result_branch, "Rotation Angle CR_1" , 1_pd_ik, tmp_real_fd1)
    CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
          Int(root%branches(3)%leaves(16)%lbound-1+(comm_nn-1), MPI_OFFSET_KIND), &
          tmp_real_fd1, &
          1_pd_mik, MPI_REAL8, status_mpi, ierr)

    !------------------------------------------------------------------------------
    ! Rotation Vector CR_1
    !------------------------------------------------------------------------------
    CALL add_leaf_to_branch(result_branch, "Rotation Vector CR_1", 3_pd_ik, n)
    CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
          Int(root%branches(3)%leaves(17)%lbound-1+(comm_nn-1)*3, MPI_OFFSET_KIND), &
          n, &
          3_pd_mik, MPI_REAL8, status_mpi, ierr)

       
    !=========================================================================
    !== Inlining of EE =======================================================
    aa = rot_alg(n,alpha)
    BB = tra_R6(aa)
    EE = matmul(matmul(transpose(BB),EE_Orig),BB)

    If (out_amount /= "PRODUCTION" ) then
       Call Write_matrix(un_lf, "Backrotated anisotropic stiffness CR_1", EE, fmti='std', unit='MPa')
    End If
    
    !=========================================================

    If ( (EE(1,1) < EE(2,2)) .AND.  &
         (EE(1,1) < EE(3,3)) .AND.  (EE(2,2) < EE(3,3))) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"123"
       Continue
       
    Else If ( (EE(1,1) < EE(2,2)) .AND.  &
         (EE(1,1) < EE(3,3)) .AND.  (EE(2,2) > EE(3,3))) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"132"

       ! 132 => 123 ********
       n = aa(:,1)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

    Else If ( (EE(1,1) < EE(2,2)) .AND.  &
         (EE(1,1) > EE(3,3)) .AND.  (EE(2,2) > EE(3,3))) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"231"

       ! 231 => 132 ********
       n = aa(:,2)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)

       ! 132 => 123 ********
       n = aa(:,1)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

    Else If ( (EE(1,1) > EE(2,2)) .AND.  &
         (EE(1,1) < EE(3,3)) .AND.  (EE(2,2) < EE(3,3)) ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"213"

       ! 213 => 123 ********
       n = aa(:,3)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

    Else If ( (EE(1,1) > EE(2,2)) .AND.  &
         (EE(1,1) > EE(3,3)) .AND.  (EE(2,2) < EE(3,3)) ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"312"

       ! 312 => 132 ********
       n = aa(:,3)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

       ! 132 => 123 ********
       n = aa(:,1)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)  

    Else If ( (EE(1,1) > EE(2,2)) .AND.  &
         (EE(1,1) > EE(3,3)) .AND.  (EE(2,2) > EE(3,3)) ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"321"

       ! 321 => 123 ********
       n = aa(:,2)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)  

    End If

    BB = tra_R6(aa)
    EE = matmul(matmul(transpose(BB),EE_Orig),BB)

    If (out_amount /= "PRODUCTION" ) then
       Call Write_matrix(un_lf, "Final coordinate system CR_1", aa, fmti='std')
       Call Write_matrix(un_lf, "Inlined anisotropic stiffness CR_1", EE, fmti='std', unit='MPa')
    End If
    
    !------------------------------------------------------------------------------
    ! Final coordinate system CR_1
    !------------------------------------------------------------------------------
    CALL add_leaf_to_branch(result_branch, "Final coordinate system CR_1", 9_pd_ik, reshape(aa,[9_pd_ik]))
    CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
          Int(root%branches(3)%leaves(18)%lbound-1+(comm_nn-1)*9, MPI_OFFSET_KIND), &
          reshape(aa,[9_pd_ik]), &
          9_pd_mik, MPI_REAL8, status_mpi, ierr)

    !------------------------------------------------------------------------------
    ! Optimized Effective stiffness CR_1
    !------------------------------------------------------------------------------
    CALL add_leaf_to_branch(result_branch, "Optimized Effective stiffness CR_1", 36_pd_ik, reshape(EE,[36_pd_ik]))
    CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
          Int(root%branches(3)%leaves(19)%lbound-1+(comm_nn-1)*36, MPI_OFFSET_KIND), &
          reshape(EE, [36_pd_ik]), &
          36_pd_mik, MPI_REAL8, status_mpi, ierr)

     If (out_amount /= "PRODUCTION" ) then
          Call Write_matrix(un_lf, "Optimized Effective stiffness CR_1", EE, fmti='std')
     End If
 

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
       write(un_lf,FMT_MSG_AxI0)'Initial Minloc  CR_2: ',mlc
       write(un_lf,FMT_MSG_xAF0) 'Initial Minimum CR_2: ',crit_min(0)
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
          write(un_lf,FMT_MSG_AxI0)'Iteration : ',jj
          write(un_lf,FMT_MSG_AxI0)'Loop start: ',s_loop
          write(un_lf,FMT_MSG_AxI0)'Loop end  : ',e_loop
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
          write(un_lf,FMT_MSG_AxF0)'Minimum CR_2         : ', crit_min(jj)
          write(un_lf,FMT_MSG_AxI0)'Minloc  CR_2         : ', minloc(crit_2(0:kk-1,0:kk_phi-1,0:kk_eta-1))
          write(un_lf,FMT_MSG_AxI0)'kk, kk_phi, kk_eta   : ', kk,kk_phi,kk_eta
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
       write(un_lf, *)
       Write(un_lf, FMT_MSG_xAI0) "Solution converged after : ", jj," iterations"
       Write(un_lf, FMT_MSG_AxF0) "With final citerion 2    : ", minval(crit_2(1:kk-2, 1:kk_phi-2, 1:kk_eta-2))
       Write(un_lf, FMT_MSG_AxF0) "With final epsilon       : ", crit_min(jj-1)-crit_min(jj)
       Write(un_lf, FMT_MSG_AxF0) "Final rotation angle  is : ", alpha
       Write(un_lf, FMT_MSG_AxF0) "Final rotation vector is : ", n
       Write(un_lf, *)
    End If
    
    !------------------------------------------------------------------------------
    ! Rotation Angle CR_2
    !------------------------------------------------------------------------------
    tmp_real_fd1 = alpha 

    CALL add_leaf_to_branch(result_branch, "Rotation Angle CR_2" , 1_pd_ik, tmp_real_fd1)
    CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
        Int(root%branches(3)%leaves(20)%lbound-1+(comm_nn-1), MPI_OFFSET_KIND), &
        tmp_real_fd1, &
        1_pd_mik, MPI_REAL8, status_mpi, ierr)
    
    !------------------------------------------------------------------------------
    ! Rotation Vector CR_2
    !------------------------------------------------------------------------------
    CALL add_leaf_to_branch(result_branch, "Rotation Vector CR_2", 3_pd_ik, n)
    CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
        Int(root%branches(3)%leaves(21)%lbound-1+(comm_nn-1)*3, MPI_OFFSET_KIND), &
        n, &
        3_pd_mik, MPI_REAL8, status_mpi, ierr)
    
    !------------------------------------------------------------------------------
    ! Inlining of EE
    !------------------------------------------------------------------------------
    aa = rot_alg(n,alpha)
    BB = tra_R6(aa)
    EE = matmul(matmul(transpose(BB),EE),BB)

    If (out_amount /= "PRODUCTION" ) &
         Call Write_matrix(un_lf, "Backrotated anisotropic stiffness CR_2", EE, fmti='std', unit='MPa')

    If ( (EE(1,1) < EE(2,2)) .AND.  &
         (EE(1,1) < EE(3,3)) .AND.  (EE(2,2) < EE(3,3))         ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"123"
       continue
       
    Else If ( (EE(1,1) < EE(2,2)) .AND.  &
         (EE(1,1) < EE(3,3)) .AND.  (EE(2,2) > EE(3,3)) ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"132"

       ! 132 => 123 ********
       n = aa(:,1)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

    Else If ( (EE(1,1) < EE(2,2)) .AND.  &
         (EE(1,1) > EE(3,3)) .AND.  (EE(2,2) > EE(3,3)) ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"231"

       ! 231 => 132 ********
       n = aa(:,2)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)

       ! 132 => 123 ********
       n = aa(:,1)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

    Else If ( (EE(1,1) > EE(2,2)) .AND.  &
         (EE(1,1) < EE(3,3)) .AND.  (EE(2,2) < EE(3,3)) ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"213"

       ! 213 => 123 ********
       n = aa(:,3)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

    Else If ( (EE(1,1) > EE(2,2)) .AND.  &
         (EE(1,1) > EE(3,3)) .AND.  (EE(2,2) < EE(3,3)) ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"312"

       ! 312 => 132 ********
       n = aa(:,3)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)   

       ! 132 => 123 ********
       n = aa(:,1)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)  

    Else If ( (EE(1,1) > EE(2,2)) .AND.  &
         (EE(1,1) > EE(3,3)) .AND.  (EE(2,2) > EE(3,3)) ) then

       If (out_amount /= "PRODUCTION" ) write(un_lf,*)"321"

       ! 321 => 123 ********
       n = aa(:,2)
       alpha = pi/2
       aa = matmul(rot_alg(n,alpha),aa)  

    End If

    BB = tra_R6(aa)
    EE = matmul(matmul(transpose(BB),EE_Orig),BB)

    If (out_amount /= "PRODUCTION" ) then
       Call Write_matrix(un_lf, "Final coordinate system CR_2", aa, fmti='std')
       Call Write_matrix(un_lf, "Inlined anisotropic stiffness CR_2", EE, fmti='std', unit='MPa')
    End If
    
    !------------------------------------------------------------------------------
    ! Final coordinate system CR_2
    !------------------------------------------------------------------------------
    CALL add_leaf_to_branch(result_branch,"Final coordinate system CR_2", 9_pd_ik, reshape(aa,[9_pd_ik]))
    CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
        Int(root%branches(3)%leaves(22)%lbound-1+(comm_nn-1)*9, MPI_OFFSET_KIND), &
        reshape(aa,[9_pd_ik]), &
        9_pd_mik, MPI_REAL8, status_mpi, ierr)
        
    !------------------------------------------------------------------------------
    ! Optimized Effective stiffness CR_2
    !------------------------------------------------------------------------------
    CALL add_leaf_to_branch(result_branch,"Optimized Effective stiffness CR_2", 36_pd_ik, reshape(EE,[36_pd_ik]))
    CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
        Int(root%branches(3)%leaves(23)%lbound-1+(comm_nn-1)*36, MPI_OFFSET_KIND), &
        reshape(EE,[36_pd_ik]), &
        36_pd_mik, MPI_REAL8, status_mpi, ierr)

    !------------------------------------------------------------------------------
    ! Domain number
    ! In some sense, the List of domain numbers represents a status file tailored
    ! to the PETSc sub comm.
    !------------------------------------------------------------------------------
    ! Last piece of information written to file. If it is the first entry, 
    ! data is more likely to get corrupted. For example, the domain number gets
    ! written to file while some computations within this module is still pending.
    !------------------------------------------------------------------------------
    CALL add_leaf_to_branch(result_branch, "Domain number", 1_ik, [ddc_nn])
    CALL MPI_FILE_WRITE_AT(fh_mpi_worker(4), &
        Int(root%branches(3)%leaves(1)%lbound-1+(comm_nn-1), MPI_OFFSET_KIND), &
        ddc_nn, 1_pd_mik, MPI_INTEGER8, status_mpi, ierr)


     !------------------------------------------------------------------------------
     ! Number of Elements
     !------------------------------------------------------------------------------
     CALL add_leaf_to_branch(result_branch, "Number of Elements", 1_pd_ik, [no_elems])
     CALL MPI_FILE_WRITE_AT(fh_mpi_worker(4), &
          Int(root%branches(3)%leaves(2)%lbound-1+(comm_nn-1), MPI_OFFSET_KIND), &
          no_elems, 1_pd_mik, MPI_INTEGER8, status_mpi, ierr)

     !------------------------------------------------------------------------------
     ! Number of Nodes
     !------------------------------------------------------------------------------
     CALL add_leaf_to_branch(result_branch, "Number of Nodes", 1_pd_ik, [no_nodes])
     CALL MPI_FILE_WRITE_AT(fh_mpi_worker(4), &
          Int(root%branches(3)%leaves(3)%lbound-1+(comm_nn-1), MPI_OFFSET_KIND), &
          no_nodes, 1_pd_mik, MPI_INTEGER8, status_mpi, ierr)

    If (out_amount /= "PRODUCTION" ) then
        Call Write_matrix(std_out, "Optimized Effective stiffness CR_2", EE, fmti='std')
    End If
    
    Select Case (timer_level)
    Case (3)
       call end_timer("  +-- Back rotation of material matrix "//trim(nn_char))
    Case (2)
       call end_timer("  +-- Back rotation of material matrix "//trim(nn_char))
    Case default
       continue
    End Select

    !------------------------------------------------------------------------------
    ! Effective density
    !------------------------------------------------------------------------------
    eff_density = 0._rk
    eff_density = REAL(no_elems, rk) / &
                REAL(ANINT(x_D_phy(1)/delta(1)) * &
                     ANINT(x_D_phy(2)/delta(2)) * &
                     ANINT(x_D_phy(3)/delta(3)), rk)

    CALL add_leaf_to_branch(result_branch, "Effective density", 1_pd_ik, [eff_density])
    CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
        Int(root%branches(3)%leaves(24)%lbound-1+(comm_nn-1), MPI_OFFSET_KIND), &
        eff_density, &
        1_pd_mik, MPI_REAL8, status_mpi, ierr)

     !------------------------------------------------------------------------------
     ! Write another memory log.
     !------------------------------------------------------------------------------
     IF (no_nodes /= 0) THEN
          collected_logs(22) = size_mpi
          collected_logs(23) = global_rank_mpi
          collected_logs(24) = SUM(collected_logs(15:21))

          CALL add_leaf_to_branch(result_branch, "Collected logs", 24_ik, [collected_logs])
          CALL MPI_FILE_WRITE_AT(fh_mpi_worker(4), &
               Int(root%branches(3)%leaves(4)%lbound-1+(comm_nn-1)*24, MPI_OFFSET_KIND), &
               collected_logs, 24_pd_mik, MPI_INTEGER8, status_mpi, ierr)

     END IF

    DEALLOCATE(tmp_nn, delta, x_D_phy)
    DEALLOCATE(nodes, vv, ff, stiffness)
    DEALLOCATE(calc_rforces, uu, rforces, edat, crit_1, crit_2)
    DEALLOCATE(ang)
    DEALLOCATE(no_cnodes_pp, cref_cnodes)

End subroutine calc_effective_material_parameters


!------------------------------------------------------------------------------
! SUBROUTINE: init_loadcase_el_order_lin
!------------------------------------------------------------------------------
!> @author Ralf Schneider - HLRS - NUM - schneider@hlrs.de
!
! @Brief:
!> 24 loadcases for lin macro elements.
!------------------------------------------------------------------------------
subroutine init_loadcase_el_order_lin(eps,vv)

    Real(rk), intent(in)     :: eps
    Real(rk), Dimension(:,:) :: vv

    Real(rk)                 :: eps2,eps3,eps4

    eps2 = eps*2._rk
    eps3 = eps*3._rk
    eps4 = eps*4._rk

    vv = 0._rk

    vv( 1,:) = [ 0._rk, 0._rk, 0._rk,   eps,   eps, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk]
    vv( 2,:) = [   eps, 0._rk, 0._rk,   eps,   eps, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk]
    vv( 3,:) = [   eps, 0._rk, 0._rk, 0._rk,   eps, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk]
    vv( 4,:) = [ 0._rk, 0._rk, 0._rk, 0._rk,   eps, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,   eps, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk]
    vv( 5,:) = [ 0._rk, 0._rk, 0._rk,   eps, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,   eps, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk]
    vv( 6,:) = [   eps, 0._rk, 0._rk,   eps, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,   eps, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk]
    vv( 7,:) = [   eps, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                  eps, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk]
    vv( 8,:) = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,   eps, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk,   eps, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk]
    vv( 9,:) = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,   eps, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk,   eps, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk]
    vv(10,:) = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                  eps, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk,   eps, 0._rk, 0._rk, 0._rk, 0._rk]
    vv(11,:) = [ 0._rk,   eps, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk,   eps, 0._rk, 0._rk, 0._rk]
    vv(12,:) = [ 0._rk,   eps, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,   eps, 0._rk, 0._rk]
    vv(13,:) = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,   eps, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,   eps, 0._rk]
    vv(14,:) = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,   eps, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,  eps ]
    vv(15,:) = [ 0._rk,   eps, 0._rk, 0._rk, 0._rk,   eps, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, -eps ]
    vv(16,:) = [ 0._rk,   eps, 0._rk, 0._rk, 0._rk,   eps, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk]
    vv(17,:) = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk,   eps, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk]
    vv(18,:) = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk,   eps, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk]
    vv(19,:) = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk,   eps, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk]
    vv(20,:) = [ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk,   eps, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk]
    vv(21,:) = [ 0._rk, 0._rk,   eps, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk]
    vv(22,:) = [ 0._rk, 0._rk,   eps, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, -eps2 , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk]
    vv(23,:) = [ 0._rk, 0._rk,   eps, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, -eps3 , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk]
    vv(24,:) = [ 0._rk, 0._rk,   eps, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, &
                 0._rk, 0._rk, 0._rk, -eps4, 0._rk, 0._rk, 0._rk, 0._rk]

  end subroutine init_loadcase_el_order_lin

!------------------------------------------------------------------------------
! SUBROUTINE: init_loadcase_el_order_quad
!------------------------------------------------------------------------------
!> @author Ralf Schneider - HLRS - NUM - schneider@hlrs.de
!
! @Brief:
!> 24 loadcases for lin macro elements.
!------------------------------------------------------------------------------
   subroutine init_loadcase_el_order_quad(e,vv)

     Real(rk), intent(in) :: e
     Real(rk), Dimension(:,:) :: vv
     REAL(rk), PARAMETER :: z = 0._rk, f = 0.99794_rk, oh=1.5_rk, h=0.5_rk, zz=2.0_rk
     vv = 0._rk
     
vv( 1,:) = [   z ,   z ,   z ,   e ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
    3._rk*e,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
       z ,   z,    z ,   z ,   z ,   z]
vv( 2,:) = [   e ,   z ,   z ,   e ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   &
      z ,   z,    z ,   z ,   z ,   z]
vv( 3,:) = [   e ,   z ,   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,-3._rk*e,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
         z ,   z,    z ,   z ,   z ,   z]
vv( 4,:) = [   z ,   z ,   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , 0.2_rk*e  ,   z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   &
      z ,-1.5 ,   z  ,  z  ,  z  ,  z ,   z]
vv( 5,:) = [   z ,   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   e ,   z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z , 5._rk*e,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
         z ,   z,    z ,   z ,   z ,   z]
vv( 6,:) = [   e ,   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,    5._rk*e,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , 0.1_rk*e,   z ,  &
    z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z &
    ,   z ,   z ,   z,    z ,   z ,   z ,   z]
vv( 7,:) = [   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,-4.  ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   &
      z ,   z,    z ,   z ,   z ,   z]
vv( 8,:) = [   z ,   z ,   z ,   z ,   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   &
      z ,   z,    z ,   z ,   z ,   z]
vv( 9,:) = [ 0.5_rk*e ,   z ,   z ,   e ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,    z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   e ,   z ,   z ,   &
      z ,   z ,   z,    z ,   z ,   z ,   z]
vv(10,:) = [   e ,   z ,   z , 0.5_rk*e ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,    z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,-0.9_rk*e ,   z ,   z ,   z ,   z&
       ,   z ,   z ,   z,    z ,   z ,   z ,   z]
vv(11,:) = [ 0.5_rk*e ,   z ,   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,    z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,-  e ,   &
      z ,   z ,   z,    z ,   z ,   z ,   z]
vv(12,:) = [   z ,   z ,   z , 0.5_rk*e ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , 0.5_rk*e   ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   &
      z ,   z ,   z ,   z,    z ,   z ,   z ,   z]
vv(13,:) = [ 0.5_rk*e ,   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , 0.5_rk*e   ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   &
      z ,   z ,   z ,   z,    z ,   z ,   z ,   z]
vv(14,:) = [   e ,   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,    z ,&
 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,-0.1_rk*e,   z ,   z ,   z ,   &
 z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,-0.7 ,   z ,   z ,   z ,&
    z ,   z ,   z ,   z,    z ,   z ,   z ,   z]
vv(15,:) = [ 0.5_rk*e ,   z ,   z ,   z ,   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,     z ,&
   z ,   z , 0.1_rk*e,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  &
    z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z &
    ,   z ,   z ,   z ,   z, -  e ,   z ,   z ,   z]
vv(16,:) = [   z ,   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,    z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   e ,   z ,   z ,   z ,   z ,   &
      z ,   z ,   z,    z ,-  e ,   z ,   z]
vv(17,:) = [   z ,   z ,   z ,   e , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,    z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,-  e ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   &
      z ,   z ,   z,    z ,   z ,   z ,   z]
vv(18,:) = [   e ,   z ,   z ,   e , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , z,0.2_rk*e ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,-  e ,   z ,   z ,   z ,   z ,   z ,   z ,   &
      z ,   z ,   z,    z ,   z ,   z ,   z]
vv(19,:) = [   e ,   z ,   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,    z ,&
 0.3_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , &
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,-  e ,   z ,   z ,   z ,   z ,   z ,   z&
    ,   z ,   z ,   z,    z ,   z ,   z ,   z]
vv(20,:) = [   z ,   z ,   z ,   z , 0.5_rk*e ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z , 5._rk*e,  z ,&
   z ,   z , 0.1_rk*e,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  &
    z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , 2.  ,   z ,   z ,   z ,   z ,   z &
    ,   z ,   z ,   z ,   z,    z ,   z ,   z ,   z]
vv(21,:) = [   z ,   z ,   z ,   z ,   z ,   z ,   z ,-2.  ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   &
      z ,   z,    z ,   z ,   z ,   z]
vv(22,:) = [   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,-  e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,-5._rk*e,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
         z ,   z,    z ,   z ,   z ,   z]
vv(23,:) = [   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z , 5._rk*e,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   &
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z,    z ,   z ,   z ,   z]
vv(24,:) = [   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,-2.  ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   &
      z ,   z,    z ,   z ,   z ,   z]
vv(25,:) = [   z ,   z ,   z ,   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   &
      z ,   z,    z ,   z ,   z ,   z]
vv(26,:) = [   z ,   z ,   z ,   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z , 3._rk*e,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
         z ,   z,    z ,   z ,   z ,   z]
vv(27,:) = [   z ,   e ,   z ,   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   &
      z ,   z,    z ,   z ,   z ,   z]
vv(28,:) = [   z ,   e ,   z ,   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z , 4.  ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   &
      z ,   z,    z ,   z ,   z ,   z]
vv(29,:) = [   z ,   z ,   z ,   z ,   z ,   z ,   z , 0.5_rk*e , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,    z ,&
   z ,   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  &
       z ,   z ,   z ,   z ,   z,    z ,   z ,   z ,   z]
vv(30,:) = [   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , 3._rk*e&
      ,   z ,   z ,   z ,   z,    z ,   z ,   z ,   z]
vv(31,:) = [   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  z ,&
   z ,   z ,   z ,   z , 5._rk*e,   z ,   z ,   z ,   z ,-0.1_rk*e,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z&
    ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
       z ,   z ,   z,    z ,   z ,   z ,   z]
vv(32,:) = [   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  &
       z , 5._rk*e,   z ,   z,    z ,   z ,   z ,   z]
vv(33,:) = [   z ,   z ,   z ,   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  z ,&
   z ,   z ,   z ,   z ,   z ,   z , 5._rk*e,   z ,   z ,   z ,   z , 0.1_rk*e,   z ,   z ,   z ,   z ,   z ,   z&
    ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
       z ,   z ,   z,    z ,   z ,   z ,   z]
vv(34,:) = [   z , 0.5_rk*e ,   z ,   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z &
   ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , &
     z ,   z ,-5. ,    z ,   z ,   z ,   z]
vv(35,:) = [   z ,   e ,   z ,   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,-0.1_rk*e,   z ,   z ,   z ,   z , &
     z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   e ,   z ,   z ,   z ,   &
     z ,   z ,   z,    z ,   z ,   z ,   z]
vv(36,:) = [   z , 0.5_rk*e ,   z ,   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z , 0.1_rk*e,   z ,   z , 2.  ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , &
     z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   &
     z ,   z ,   z ,   z,    z ,   z ,   z ,   z]
vv(37,:) = [   z ,   z ,   z ,   z ,   z , 0.5_rk*e ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,    z ,&
   z ,   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , 2.  ,   z ,   z ,   z ,   z ,   z ,   z ,  &
       z ,   z ,   z ,   z ,   z,    z ,   z ,   z ,   z]
vv(38,:) = [   z ,   z ,   z ,   z ,   z , 0.5_rk*e ,   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,    z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,-4.  &
   ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , &
     z ,   z ,   z ,   z,    z ,   z ,   z ,   z]
vv(39,:) = [   z ,   e ,   z ,   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z ,   z , 0.5_rk*e ,   z ,   z ,-0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , &
     z ,-  e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   &
     z ,   z ,   z ,   z ,   z,    z ,   z ,   z ,   z]
vv(40,:) = [   z ,   e ,   z ,   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z ,   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,&
      z , 2.  ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  &
       z ,   z ,   z ,   z,    z ,   z ,   z ,   z]
vv(41,:) = [   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,-  e ,   z ,   z ,   z ,   z ,   z ,   z ,  z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z &
   ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , &
     z ,   z,    z ,   z ,   z ,   z]
vv(42,:) = [   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,-2.  ,   z ,   z ,   z ,   z ,   z ,  z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z &
   ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , &
     z ,   z,    z ,   z ,   z ,   z]
vv(43,:) = [   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , 5._rk*e,   z ,   z ,   z ,   z ,  &
    z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z&
     ,   z ,   z,    z ,   z ,   z ,   z]
vv(44,:) = [   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z &
   ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,-2.  ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , &
     z ,   z,    z ,   z ,   z ,   z]
vv(45,:) = [   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z &
   ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , &
     z ,   z,    z ,   z ,   z ,   z]
vv(46,:) = [   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  z ,&
-0.7 ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , 2.  ,   z ,   z ,   z ,   z ,   z ,   z &
,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , &
  z ,   z,    z ,   z ,   z ,   z]
vv(47,:) = [   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  z ,&
   z ,-3._rk*e,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  &
    z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z&
     ,   z ,   z,    z ,   z ,   z ,   z]
vv(48,:) = [   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  z ,&
   z ,   z ,-3._rk*e,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  &
    z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z&
     ,   z ,   z,    z ,   z ,   z ,   z]
vv(59,:) = [   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , 2.  ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z &
   ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , &
     z ,   z ,   z,    z ,   z ,   z ,   z]
vv(50,:) = [   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , 5._rk*e,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  &
    z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z&
     ,   z ,   z ,   z,    z ,   z ,   z ,   z]
vv(51,:) = [   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , 0.1_rk*e,   z ,   z ,   z ,    z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z &
   ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , &
     z ,   z ,   z,    z ,   z ,   z ,   z]
vv(52,:) = [   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , 0.1_rk*e,   z ,   z ,   z ,   z ,   z ,    z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , 4.  ,   z ,   z ,   z ,   z ,   z ,   z ,   z &
   ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , &
     z ,   z ,   z,    z ,   z ,   z ,   z]
vv(53,:) = [   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  z ,&
-0.1_rk*e,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,-4.  ,   z ,   z ,   z ,   z ,   z ,   z , &
  z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   &
  z ,   z ,   z,    z ,   z ,   z ,   z]
vv(54,:) = [   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  z ,&
-  e ,-1.5 ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z &
,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,-1.5 ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , &
  z ,   z , z,  z,  z ,   z ,   z]
vv(55,:) = [   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  z ,&
   z ,   z ,-0.1_rk*e,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , &
     z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,-  e ,   z ,   z ,   z ,   z ,   z ,   &
     z ,   z ,   z,    z ,   z ,   z ,   z]
vv(56,:) = [   z ,   z ,   e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,-0.5_rk*e ,   z ,   z ,   z ,&
      z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,  &
       z ,   z ,   z,    z ,   z ,   z ,   z]
vv(57,:) = [   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,    z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z &
   ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , &
     z ,   z ,   z ,   z,    z ,   z ,   e ,   z]
vv(58,:) = [   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,    z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z &
   ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , &
     z ,   z ,   z ,   z,    z ,   z ,   z ,-  e]
vv(59,:) = [   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , 0.9_rk*e ,   z ,   z ,    z ,&
   z ,   z ,   z ,-1.5 ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z &
   ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , &
     z ,   z ,   z ,   z,    z ,   z ,   z ,   z]
vv(60,:) = [   z ,   z , 0.5_rk*e ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , 5._rk*e,   z ,      z ,&
   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z &
   ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z ,   z , &
     z ,   z ,   z ,   z,    z ,   z ,   z ,   z]

         end subroutine init_loadcase_el_order_quad
     
     End Module calcmat
      