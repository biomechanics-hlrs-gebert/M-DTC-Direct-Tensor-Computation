Module write_deck

  use linFE
  use mesh_partitioning
  use strings
  
  implicit none

  Real(Kind=rk), Parameter     :: delta_b = 1.E-9_rk

  Character(Len=*) , Parameter :: fmt_filename = "(A,I0,A,I0,A)"

contains

  !****************************************************************************
  !**
  !** Write FMPS Loadcase Input Decks
  !** 
  Subroutine write_fmps_decks(job_dir, ddc_nn, meta_para, db)

    !-- Input Parameters ------------------------------------------------------
    Character(LEN=*), Intent(in)        :: job_dir
    Integer(Kind=ik), Intent(in)        :: ddc_nn
    Type(tBranch)   , Intent(In)        :: meta_para
    Type(tBranch)   , Intent(InOut)     :: db

    Type(tBranch)   , Pointer           :: fb, ifb

    Integer                             :: elo_macro
    Real(Kind=rk)                       :: e_modul,nu

    Character(Len=8)                    :: elt_micro
    Integer                             :: no_solver,pscratch, parts

    !-- Parameters ------------------------------------------------------------
    Character(len=*), Parameter         :: inpsep = "('#',79('='))"

    !--------------------------------------------------------------------------

    Character(Len=mcl)                 :: desc, line

    integer                            :: no_elem_nodes
    integer                            :: no_nodes_macro

    integer(Kind=ik)                                 :: ii, jj
    integer(Kind=ik), Dimension(128)                 :: crp
    Character, Dimension(:), Allocatable             :: char_arr
    Character(Len=:), Allocatable                    :: alloc_str

    !--------------------------------------------------------------------------

    call pd_get(meta_para,'Element type  on micro scale',char_arr)
    elt_micro = char_to_str(char_arr)
    deallocate(char_arr)

    Call pd_get(meta_para, 'No of mesh parts'     , parts)
    Call pd_get(meta_para, "Young_s modulus"      , e_modul)
    Call pd_get(meta_para, "Poisson_s ratio"      , nu)

    Call pd_get(meta_para, "Element order on macro scale", elo_macro)
    Call pd_get(meta_para, 'FMPS Solver to use'          , no_solver)
    Call pd_get(meta_para, 'FMPS Parser Scratch'         , pscratch)

    IF (elo_macro == 1) no_nodes_macro = 8

    IF (elt_micro == "HEX20") then
       no_elem_nodes = 20
    Else
       no_elem_nodes = 8
    End IF

    Call add_branch_to_branch(db,ifb)
    Call raise_branch("Input files", 0_ik, 0_ik, ifb)

    Do ii = 1, no_nodes_macro*3 

       jj        = 1
       crp       = 0
       alloc_str = ""

       Call add_branch_to_branch(ifb,fb)
       Write(desc,'(A,I0)')"Input file LC ",ii
       Call raise_branch(trim(desc), 0_ik, 0_ik, fb)

       if (out_amount == "PRODUCTION") then
          alloc_str = alloc_str//"PARAMETER"//char(10)
          jj = jj + 1
          crp(jj) = len(alloc_str)

          alloc_str = alloc_str//"I IURP 0"//char(10)
          jj = jj + 1
          crp(jj) = len(alloc_str)
       End if

       alloc_str = alloc_str//"ANALYSIS"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)     

       alloc_str = alloc_str//"PARTS"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"STATIC 'ELASTIC'"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"OPERATIONS"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"SOLVE 'ELASTIC'"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"END ANALYSIS"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"PARAMETER"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"I INCE 1"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       Write(line,'(A,I0)')"I NSCR ",pscratch
       alloc_str = alloc_str//trim(line)//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       alloc_str = alloc_str//"I MOPR 1"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       Write(line,'(A,I0)')"I MEQS ",no_solver
       alloc_str = alloc_str//trim(line)//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       alloc_str = alloc_str//"I IPDC 1     # Append to existing PureDat project"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"I ISRP 1     # Suppress generation of Puredat report file"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"I MTOL 1     # Output of Integral / Monitor Data to Log-File IUPR"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"I NITE 10"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"R EPS1 1.E-06"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"R EPS2 1.E-06"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"I IBCO 1"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"I IBTO 1"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       Write(line,'(A,I0)')"I MPRT ",parts
       alloc_str = alloc_str//trim(line)//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       alloc_str = alloc_str//"PART 'ELASTIC'"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"KNOP"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       If (no_elem_nodes == 20) then
          alloc_str = alloc_str//"'HEX20E'"//char(10)
          jj = jj + 1
          crp(jj) = len(alloc_str)
       else if (no_elem_nodes == 8) then
          alloc_str = alloc_str//"'HEX08E'"//char(10)
          jj = jj + 1
          crp(jj) = len(alloc_str)
       End If

       alloc_str = alloc_str//'COORDINATES'//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//'PUREDAT'//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//trim(job_dir)//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       Write(line,'(A,I0)')trim(project_name)//'_',ddc_nn
       alloc_str = alloc_str//trim(line)//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       Write(line,'(A,I0)')"Node Numbers"
       alloc_str = alloc_str//trim(line)//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       Write(line,'(A,I0)')"Coordinates"
       alloc_str = alloc_str//trim(line)//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       alloc_str = alloc_str//'TOPOLOGY'//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       If (no_elem_nodes == 20) then
          alloc_str = alloc_str//"'HEX20E'  PROPERTY=1"//char(10)
          jj = jj + 1
          crp(jj) = len(alloc_str)
       else if (no_elem_nodes == 8) then
          alloc_str = alloc_str//"'HEX08E'  PROPERTY=1"//char(10)
          jj = jj + 1
          crp(jj) = len(alloc_str)
       End If

       alloc_str = alloc_str//'PUREDAT'//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//trim(job_dir)//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       Write(line,'(A,I0)')trim(project_name)//'_',ddc_nn
       alloc_str = alloc_str//trim(line)//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       alloc_str = alloc_str//"Element Numbers"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       alloc_str = alloc_str//"Topology"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       alloc_str = alloc_str//"MATERIAL"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"1     &  PID"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       Write(line,'(F15.4,A)')e_modul," &  1 modulus of elasticity"
       alloc_str = alloc_str//trim(line)//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       Write(line,'(F15.4,A)')nu     ," &  2 Poissons number"
       alloc_str = alloc_str//trim(line)//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       alloc_str = alloc_str//"0.0       &  3 density"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"0.0       &  4 coefficient of thermal expansion"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"0.0       &  5 heat capacity"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"0.0       &  6 conversion factor of inelastic energy into heat"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"0.0       &  7 heat conduction "//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"0.0       &  8 dynamic viscosity"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"0.0       &  9 penalty coefficient"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"0.0       & 10 material index"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"0.0       & 11 solidus temperature"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"0.0       & 12 liquidus temperature"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"0.0       & 13 latent heat"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"0.0       & 14 convective heat transfer coefficient"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"0.0       & 15 emission number for radiation"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"0.0       & 16 radiation coefficient"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"0.0       & 17 ..."//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"0.0       ! 18 ..."//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"PRESCRIBE"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//'PUREDAT'//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//trim(job_dir)//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       Write(line,'(A,I0)')trim(project_name)//'_',ddc_nn
       alloc_str = alloc_str//trim(line)//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       Write(line,'(A,I0,A,I0)')"Boundaries"//'_',ddc_nn,'_',ii
       alloc_str = alloc_str//trim(line)//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       if (parts > 1) then

          alloc_str = alloc_str//"METIS"//char(10)
          jj = jj + 1
          crp(jj) = len(alloc_str)
          alloc_str = alloc_str//'PUREDAT'//char(10)
          jj = jj + 1
          crp(jj) = len(alloc_str)
          Write(line,'(A,I0,A)')trim(project_name)//'_',ddc_nn,'_Mesh'
          alloc_str = alloc_str//trim(line)//char(10)
          jj = jj + 1
          crp(jj) = len(alloc_str)

       End if

       alloc_str = alloc_str//"OUTPUT"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"   MEDIA      PUREDAT"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"   SELECTION  DISC EDAT FORC"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)
       alloc_str = alloc_str//"   AT_INCR    FROM 1 TO 100"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       alloc_str = alloc_str//"FIN"//char(10)
       jj = jj + 1
       crp(jj) = len(alloc_str)

       Call add_leaf_to_branch(fb,"No lines"  , 1_pd_ik       , [jj])
       Call add_leaf_to_branch(fb,"CR pointer", jj            , crp(1:jj) )
       Call add_leaf_to_branch(fb,"Input"     , len(alloc_str), str_to_char(alloc_str))

       alloc_str = ''
       jj        = 0

    End Do

  End Subroutine write_fmps_decks

  !****************************************************************************
  !**
  !** Write FMPS Loadcase Input Decks
  !** 
  Subroutine write_fmps_pdTree(job_dir, ddc,loc_ddc, PMesh, elo_macro)

    !-- Parameters ------------------------------------------------------------
    Character(LEN=*), Intent(in)               :: job_dir
    Type(tBranch),    Intent(In)               :: ddc, loc_ddc

    Type(T_PMesh),    Intent(In), Dimension(:),Allocatable :: PMesh

    Integer         , Intent(In)               :: elo_macro

    !--------------------------------------------------------------------------

    Integer(Kind=ik)                                      :: ddc_nn

    Real(Kind=rk)     , Dimension(:), Allocatable         :: dim_c, delta
    Integer(Kind=ik)  , Dimension(:), Allocatable         :: xa_n, xe_n
    Real(Kind=rk)     , Dimension(3)                      :: min_c, max_c

    Integer(kind=ik)                                      :: parts

    Character(Len=mcl)                                    :: desc

    integer                            :: no_elem_nodes

    integer                            :: no_nodes_macro

    integer(Kind=ik)                   :: ii, jj, no_bnodes, b_items


    Character(len=*), Parameter :: inpsep = "('#',79('='))"

    Type(tBranch)                                    :: loc_mesh
    Type(tBranch), Pointer                           :: part_b, bounds_b, conn_b
    Integer(Kind=ik), Dimension(:), Allocatable      :: numbers

    Integer(Kind=ik), Dimension(:), Allocatable      :: no_nodes_all
    Integer(Kind=ik), Dimension(:), Allocatable      :: no_elems_all
    Integer(Kind=ik), Dimension(:), Allocatable      :: no_cdofs_all

    Integer(Kind=ik), Dimension(:,:), Allocatable      :: bnode_ids
    Real(Kind=rk)   , Dimension(:,:), Allocatable      :: bnode_vals

    Integer(kind=pd_ik), Dimension(no_streams)         :: removed_data

    !===========================================================================
    !== Code ===================================================================
    !===========================================================================

    !** Get Parameters of domain decomposition *******
    call pd_get(loc_ddc, "nn",      ddc_nn)
    call pd_get(ddc, "x_D_phy", dim_c )
    call pd_get(loc_ddc, "xa_n",    xa_n  )
    call pd_get(loc_ddc, "xe_n",    xe_n  )
    call pd_get(ddc, "delta",   delta )

    min_c = Real(xa_n - 1,rk) * delta
    max_c = Real(xe_n    ,rk) * delta

    Write(un_lf,FMT_MSG_A3F0)'Minimum Coordinates',min_c
    Write(un_lf,FMT_MSG_A3F0)'Maximum Coordinates',max_c
    Write(un_lf,*)

    Write(un_lf,FMT_MSG_A3F0)'Cube dimensions',dim_c
    Write(un_lf,*)

    parts         = size(PMesh)
    no_elem_nodes = size(PMesh(1)%EIND(:,1))
    IF (elo_macro == 1) no_nodes_macro = 8

    !** init mesh pd-project **************************************************
    pro_path = trim(job_dir)
    Write(pro_name,'(A,I0,A,I0)')trim(project_name)//"_",ddc_nn,"_Mesh"

    Write(desc,'(A,1X,I0)')'Nodes, Elements and boundaries of '//trim(project_name),ddc_nn

    call raise_tree(trim(desc),loc_mesh)
    call raise_branch(trim(desc), Int(parts,pd_ik), 0_pd_ik, loc_mesh)

    Call open_stream_files(loc_mesh,"write","replace")

    Allocate(no_nodes_all(parts))
    Allocate(no_elems_all(parts))
    Allocate(no_cdofs_all(parts))

    Do jj = 1, parts

       no_nodes_all(jj) = PMesh(jj)%nnodes       
       no_elems_all(jj) = PMesh(jj)%nelems

       part_b => loc_mesh%branches(jj)

       Write(desc,'(A,I0)')'Part_',jj
       call raise_branch(trim(desc), 0, 0, part_b)

       !***********************************************************************
       !** Nodes **************************************************************

       call add_leaf_to_branch(part_b, "Node Numbers", 4_1, PMesh(jj)%nnodes)

       Allocate(numbers(PMesh(jj)%nnodes))
       Do ii = 1, PMesh(jj)%nnodes
          numbers(ii) = ii
       End Do

       call pd_store(loc_mesh%streams, part_b, "Node Numbers", numbers)

       deallocate(numbers)

       call add_leaf_to_branch(part_b, "Coordinates", 5_1, INT(PMesh(jj)%nnodes,pd_ik)*3_pd_ik)
       call pd_store(loc_mesh%streams, part_b, "Coordinates", PMesh(jj)%coor)

       call add_leaf_to_branch(part_b, "Global Node Numbers", 4_1, PMesh(jj)%nnodes)
       call pd_store(loc_mesh%streams, part_b, "Global Node Numbers", &
            PMesh(jj)%nn(1:PMesh(jj)%nnodes))

       !***********************************************************************    
       !** Elements ***********************************************************

       call add_leaf_to_branch(part_b, "Element Numbers", 4_1, PMesh(jj)%nelems)

       Allocate(numbers(PMesh(jj)%nelems))
       Do ii = 1, PMesh(jj)%nelems
          numbers(ii) = ii
       End Do

       call pd_store(loc_mesh%streams, part_b, "Element Numbers", numbers)

       deallocate(numbers)

       !** HEXE 20 ************************************************************
       If (no_elem_nodes == 20) then

          call add_leaf_to_branch(part_b, "Topology", 4_1, INT(PMesh(jj)%nelems,pd_ik)*20_pd_ik)

          !** HEXE 08 ************************************************************
       else if (no_elem_nodes == 8) then

          call add_leaf_to_branch(part_b, "Topology", 4_1, INT(PMesh(jj)%nelems,pd_ik)*8_pd_ik)

       End If

       call pd_store(loc_mesh%streams, part_b, "Topology", PMesh(jj)%eind)

       !***********************************************************************
       !** Connections to other parts *****************************************

       call add_branch_to_branch(part_b,conn_b)
       call raise_branch("Connections", 0, 0, conn_b)

       call add_leaf_to_branch(conn_b, "Outer Nodes", 4_1, PMesh(jj)%nouter_nds)
       call pd_store(loc_mesh%streams, conn_b, "Outer Nodes", PMesh(jj)%bnodes)

       call add_leaf_to_branch(conn_b, "Neighbours", 4_1, size(PMesh(jj)%neigh))
       call pd_store(loc_mesh%streams, conn_b, "Neighbours", PMesh(jj)%neigh)

       !***********************************************************************
       !** Boundaries *********************************************************

       !** Determine number of boundary nodes
       no_bnodes = 0

       Do ii = 1, PMesh(jj)%nnodes

          !** Nodes in facet 1 ************************************************
          if ( (PMesh(jj)%coor(3,ii) - min_c(3)) <= (dim_c(3) * delta_b) ) then
             no_bnodes = no_bnodes + 1

             !** Nodes in facet 6 ***************************************************
          Else if ( (max_c(3) - PMesh(jj)%coor(3,ii)) <= (dim_c(3) * delta_b) ) then
             no_bnodes = no_bnodes + 1

             !** Nodes in facet 2 ***************************************************
          Else if ( (PMesh(jj)%coor(2,ii) - min_c(2)) <= (dim_c(2) * delta_b) ) then
             no_bnodes = no_bnodes + 1

             !** Nodes in facet 4 ***************************************************
          Else if ( (max_c(2) - PMesh(jj)%coor(2,ii)) <= (dim_c(2) * delta_b) ) then
             no_bnodes = no_bnodes + 1

             !** Nodes in facet 5 ***************************************************
          Else if ( (PMesh(jj)%coor(1,ii) - min_c(1)) <= (dim_c(1) * delta_b) ) then
             no_bnodes = no_bnodes + 1

             !** Nodes in facet 3 ***************************************************
          Else if ( (max_c(1) - PMesh(jj)%coor(1,ii)) <= (dim_c(1) * delta_b) ) then
             no_bnodes = no_bnodes + 1

          End if

       End Do

       Write(un_lf,FMT_MSG_AI0)'Number of constrained nodes in Part ',jj," : ",no_bnodes

       call add_branch_to_branch(part_b,bounds_b)
       call raise_branch("Boundaries", no_nodes_macro*3, 0, bounds_b)

       Do ii = 1, no_nodes_macro*3
          Write(bounds_b%branches(ii)%desc,'(A,I0,A,I0)') "Boundaries"//'_',ddc_nn,'_',ii

          call add_leaf_to_branch(bounds_b%branches(ii), &
               "Boundary_Ids", 4_1, no_bnodes*3*3)
          call add_leaf_to_branch(bounds_b%branches(ii), &
               "Boundary_Values", 5_1, no_bnodes*3)
       End Do

       Allocate(bnode_ids(3,no_bnodes*3))
       Allocate(bnode_vals(no_nodes_macro*3,no_bnodes*3))

       !***********************************************************************
       !** Boundary application     
       b_items = 0

       Do ii = 1, PMesh(jj)%nnodes

          !** Nodes in facet 1 ************************************************
          if ( (PMesh(jj)%coor(3,ii) - min_c(3)) <= (dim_c(3) * delta_b) ) then
             b_items = b_items + 1
             bnode_ids(:,b_items) = [ ii, 1_pd_ik, 1_pd_ik ]
             b_items = b_items + 1
             bnode_ids(:,b_items) = [ ii, 2_pd_ik, 1_pd_ik ]
             b_items = b_items + 1
             bnode_ids(:,b_items) = [ ii, 3_pd_ik, 1_pd_ik ]

             !** Nodes in facet 6 ************************************************
          Else if ( (max_c(3) - PMesh(jj)%coor(3,ii)) <= (dim_c(3) * delta_b) ) then
             b_items = b_items + 1
             bnode_ids(:,b_items) = [ ii, 1_pd_ik, 1_pd_ik ]
             b_items = b_items + 1
             bnode_ids(:,b_items) = [ ii, 2_pd_ik, 1_pd_ik ]
             b_items = b_items + 1
             bnode_ids(:,b_items) = [ ii, 3_pd_ik, 1_pd_ik ]

             !** Nodes in facet 2 ************************************************
          Else if ( (PMesh(jj)%coor(2,ii) - min_c(2)) <= (dim_c(2) * delta_b) ) then
             b_items = b_items + 1
             bnode_ids(:,b_items) = [ ii, 1_pd_ik, 1_pd_ik ]
             b_items = b_items + 1
             bnode_ids(:,b_items) = [ ii, 2_pd_ik, 1_pd_ik ]
             b_items = b_items + 1
             bnode_ids(:,b_items) = [ ii, 3_pd_ik, 1_pd_ik ]

             !** Nodes in facet 4 ************************************************
          Else if ( (max_c(2) - PMesh(jj)%coor(2,ii)) <= (dim_c(2) * delta_b) ) then
             b_items = b_items + 1
             bnode_ids(:,b_items) = [ ii, 1_pd_ik, 1_pd_ik ]
             b_items = b_items + 1
             bnode_ids(:,b_items) = [ ii, 2_pd_ik, 1_pd_ik ]
             b_items = b_items + 1
             bnode_ids(:,b_items) = [ ii, 3_pd_ik, 1_pd_ik ]

             !** Nodes in facet 5 ************************************************
          Else if ( (PMesh(jj)%coor(1,ii) - min_c(1)) <= (dim_c(1) * delta_b) ) then
             b_items = b_items + 1
             bnode_ids(:,b_items) = [ ii, 1_pd_ik, 1_pd_ik ]
             b_items = b_items + 1
             bnode_ids(:,b_items) = [ ii, 2_pd_ik, 1_pd_ik ]
             b_items = b_items + 1
             bnode_ids(:,b_items) = [ ii, 3_pd_ik, 1_pd_ik ]

             !** Nodes in facet 3 ************************************************
          Else if ( (max_c(1) - PMesh(jj)%coor(1,ii)) <= (dim_c(1) * delta_b) ) then
             b_items = b_items + 1
             bnode_ids(:,b_items) = [ ii, 1_pd_ik, 1_pd_ik ]
             b_items = b_items + 1
             bnode_ids(:,b_items) = [ ii, 2_pd_ik, 1_pd_ik ]
             b_items = b_items + 1
             bnode_ids(:,b_items) = [ ii, 3_pd_ik, 1_pd_ik ]

          End if

       End Do

       b_items = 0

       Do ii = 1, PMesh(jj)%nnodes

          !** Nodes in facet 1 ************************************************
          if ( (PMesh(jj)%coor(3,ii) - min_c(3)) <= (dim_c(3) * delta_b) ) then
             Call Determine_prescribed_displ(PMesh(jj)%coor(:,ii), min_c, max_c, &
                  no_nodes_macro, bnode_vals(:,b_items+1:b_items+3))
             b_items = b_items+3

             !** Nodes in facet 6 ************************************************
          Else if ( (max_c(3) - PMesh(jj)%coor(3,ii)) <= (dim_c(3) * delta_b) ) then
             Call Determine_prescribed_displ(PMesh(jj)%coor(:,ii), min_c, max_c, &
                  no_nodes_macro, bnode_vals(:,b_items+1:b_items+3))
             b_items = b_items+3

             !** Nodes in facet 2 ************************************************
          Else if ( (PMesh(jj)%coor(2,ii) - min_c(2)) <= (dim_c(2) * delta_b) ) then
             Call Determine_prescribed_displ(PMesh(jj)%coor(:,ii), min_c, max_c, &
                  no_nodes_macro, bnode_vals(:,b_items+1:b_items+3))
             b_items = b_items+3

             !** Nodes in facet 4 ************************************************
          Else if ( (max_c(2) - PMesh(jj)%coor(2,ii)) <= (dim_c(2) * delta_b) ) then
             Call Determine_prescribed_displ(PMesh(jj)%coor(:,ii), min_c, max_c, &
                  no_nodes_macro, bnode_vals(:,b_items+1:b_items+3))
             b_items = b_items+3

             !** Nodes in facet 5 ************************************************
          Else if ( (PMesh(jj)%coor(1,ii) - min_c(1)) <= (dim_c(1) * delta_b) ) then
             Call Determine_prescribed_displ(PMesh(jj)%coor(:,ii), min_c, max_c, &
                  no_nodes_macro, bnode_vals(:,b_items+1:b_items+3))
             b_items = b_items+3

             !** Nodes in facet 3 ************************************************
          Else if ( (max_c(1) - PMesh(jj)%coor(1,ii)) <= (dim_c(1) * delta_b) ) then
             Call Determine_prescribed_displ(PMesh(jj)%coor(:,ii), min_c, max_c, &
                  no_nodes_macro, bnode_vals(:,b_items+1:b_items+3))
             b_items = b_items+3

          End if

       End Do

       Write(un_lf,FMT_MSG_AI0)'Number of constrained DOF',b_items

       no_cdofs_all(jj) = b_items

       Do ii = 1, no_nodes_macro*3

          call pd_store(loc_mesh%streams,bounds_b%branches(ii), &
               "Boundary_Ids", bnode_ids)
          call pd_store(loc_mesh%streams,bounds_b%branches(ii), &
               "Boundary_Values", bnode_vals(ii,:))

       End Do

       DeAllocate(bnode_ids)
       DeAllocate(bnode_vals)

    End Do

    call add_leaf_to_branch(loc_mesh, "No of nodes in parts", 4_1, parts)
    call add_leaf_to_branch(loc_mesh, "No of elems in parts", 4_1, parts)
    call add_leaf_to_branch(loc_mesh, "No of cdofs in parts", 4_1, parts)

    call pd_store(loc_mesh%streams,loc_mesh, "No of nodes in parts", no_nodes_all)
    call pd_store(loc_mesh%streams,loc_mesh, "No of elems in parts", no_elems_all)
    call pd_store(loc_mesh%streams,loc_mesh, "No of cdofs in parts", no_cdofs_all)

    call write_tree(loc_mesh)

    Call close_stream_files(loc_mesh,.TRUE.)

    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if (out_amount == "DEBUG") then
      WRITE(un_lf,FMT_DBG_SEP)
       call log_tree(loc_mesh,un_lf)
      WRITE(un_lf,FMT_DBG_SEP)
    End if
    !** DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    Call destroy_tree(loc_mesh, removed_data)

  End Subroutine write_fmps_pdTree

  !****************************************************************************
  !**
  !** Determine prescribed displacements
  !**
  Subroutine determine_prescribed_displ(node, xa, xe, elnodes, bnode_vals)

    Real(Kind=rk)   , Intent(In), Dimension(3)  :: node, xa, xe
    Integer         , Intent(in)                :: elnodes

    Real(Kind=rk), dimension(elnodes,3,elnodes*3) :: uu 

    integer                                       :: ii

    Real(Kind=rk)   , Dimension(:,:), Intent(Out) :: bnode_vals

    !**************************************************************************

    !** Loadcase init *********************************************************  
    uu = init_displ(elnodes)

    Do ii = 1, elnodes*3

       bnode_vals(ii,1) = sum(phi_nn(t_geom_xi(node,xa,xe))*uu(:,1,ii))
       bnode_vals(ii,2) = sum(phi_nn(t_geom_xi(node,xa,xe))*uu(:,2,ii))
       bnode_vals(ii,3) = sum(phi_nn(t_geom_xi(node,xa,xe))*uu(:,3,ii))

    End Do

  End Subroutine determine_prescribed_displ

  !****************************************************************************
  !** Initialisation of loadcases                                 
  !**
  Function init_displ(elnodes) Result(uu)

    Integer      , Intent(in)                     :: elnodes

    Real(Kind=rk), dimension(elnodes,3,elnodes*3) :: uu 
    Real(Kind=rk)                                 :: eps=1.E-6_rk

    uu = 0._rk

    uu(:,1,1)   = (/ 0._rk, eps, eps, 0._rk, 0._rk, eps, eps, 0._rk /)
    uu(:,2,2)   = (/ 0._rk, 0._rk, eps, eps, 0._rk, 0._rk, eps, eps /)
    uu(:,3,3)   = (/ 0._rk, 0._rk, 0._rk, 0._rk, eps, eps, eps, eps /)
    uu(:,1,4)   = (/ eps, eps, 0._rk, 0._rk, eps, eps, 0._rk, 0._rk /)
    uu(:,1,5)   = (/ eps, eps, eps, eps, 0._rk, 0._rk, 0._rk, 0._rk /)
    uu(:,2,6)   = (/ 0._rk, 0._rk, 0._rk, 0._rk, eps, eps, eps, eps /)

    uu(:,1, 7)   = (/ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, eps   /)
    uu(:,2, 8)   = (/ eps  , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk /)
    uu(:,2, 9)   = (/ 0._rk, eps  , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk /)
    uu(:,3,10)   = (/ eps  , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk /)
    uu(:,3,11)   = (/ 0._rk, eps  , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk /)
    uu(:,3,12)   = (/ 0._rk, 0._rk, eps  , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk /)
    uu(:,3,13)   = (/ 0._rk, 0._rk, 0._rk, eps  , 0._rk, 0._rk, 0._rk, 0._rk /)

    uu(:,1,14)   = (/ 0._rk, 0._rk, 0._rk, eps  , 0._rk, 0._rk, 0._rk, 0._rk /)
    uu(:,1,15)   = (/ 0._rk, 0._rk, 0._rk, 0._rk, eps  , 0._rk, 0._rk, 0._rk /)
    uu(:,1,16)   = (/ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, eps  , 0._rk, 0._rk /)
    uu(:,1,17)   = (/ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, eps  , 0._rk /)

    uu(:,1,18)   = (/ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, eps   /)
    uu(:,3,18)   = (/ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk,-2*eps, 0._rk, 0._rk /)

    uu(:,2,19)   = (/ eps  , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk /)
    uu(:,3,19)   = (/ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, -3*eps,0._rk /)

    uu(:,2,20)   = (/ 0._rk, eps  , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk /)
    uu(:,3,20)   = (/ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, -4*eps/)

    uu(:,2,21)   = (/ 0._rk, 0._rk, eps  , 0._rk, 0._rk, 0._rk, 0._rk, 0._rk /)
    uu(:,2,22)   = (/ 0._rk, 0._rk, 0._rk, eps  , 0._rk, 0._rk, 0._rk, 0._rk /)
    uu(:,2,23)   = (/ 0._rk, 0._rk, 0._rk, 0._rk, eps  , 0._rk, 0._rk, 0._rk /)

    uu(:,2,24)   = (/ 0._rk, 0._rk, 0._rk, 0._rk, 0._rk, eps  , -eps , 0._rk /)

  End Function init_displ

End Module write_deck
