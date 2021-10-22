  !>\page current Current Implementation Status
  !>The flowchart below shows the current implementation status.
  !> <hr>
  !>
  !>\dot
  !!digraph example {
  !!    compound=true;
  !!    fontname=Helvetica;
  !!    fontsize=10;
  !!    node [fontname=Helvetica, fontsize=10];
  !!
  !!    node [shape=component, color=royalblue];
  !!    a_ranks [ label="MPI_COMM_WORLD" ];
  !!
  !!    subgraph cluster_0 {
  !!        color=royalblue;
  !!        style=rounded;
  !!        label = "MPI_COMM_WORLD - Rank 0 - 'Master'";
  !!
  !!        node [shape=ellipse, color=black];
  !!        read_input [ label="Read Input" ];
  !!        create_outdir [ label="Create output directory" ];
  !!        case_restart [ label="Init new root or init restart" ];
  !!        init_resb [ label="Init result branch in root" ];
  !!        serialize_root [ label = "Serialize root" ];
  !!        broadcast_root_0 [ label = "Broadcast root" ];
  !!        split_comm_world_0 [ label=
  !!"Split MPI_COMM_WORLD
  !! to worker communicators
  !! with Rank 0 beeing not part of
  !! any worker com." ];
  !!        read_input -> create_outdir -> case_restart -> init_resb;
  !!        init_resb -> serialize_root -> broadcast_root_0 -> split_comm_world_0;
  !!    }
  !!
  !!    a_ranks -> read_input [lhead=cluster_0];
  !!
  !!    subgraph cluster_1 {
  !!        color=royalblue;
  !!        style=rounded;
  !!        label = "MPI_COMM_WORLD - Ranks > 0 - 'Workers'";
  !!
  !!        node [shape=ellipse, color=black];
  !!        broadcast_root_bt0 [ label = "Broadcast root" ];
  !!        split_comm_world_bt0 [ label=
  !!"Split MPI_COMM_WORLD \n to worker_comm by
  !! Int((rank_mpi-1)/parts)" ];
  !!        set_petsc_comm_world [label="PETSC_COMM_WORLD = worker_comm"];
  !!        node [shape=component, color=royalblue];
  !!        worker_comm_1st [label="1st domain worker_comm"];
  !!        worker_comm_2nd [label="2nd domain worker_comm"];
  !!        worker_comm_3rd [label="... domain worker_comm"];
  !!        node [shape=ellipse, color=black];
  !!        petsc_init_1st [ label="1st worker_comm \n PETSc Init " ];
  !!        petsc_init_2nd [ label="2nd worker_comm \n PETSc Init" ];
  !!        petsc_init_3rd [ label="... worker_comm \n PETSc Init" ];
  !!        broadcast_root_bt0 -> split_comm_world_bt0;
  !!        split_comm_world_bt0 -> set_petsc_comm_world;
  !!        set_petsc_comm_world -> { worker_comm_1st  worker_comm_2nd  worker_comm_3rd};
  !!        worker_comm_1st -> petsc_init_1st;
  !!        worker_comm_2nd -> petsc_init_2nd;
  !!        worker_comm_3rd -> petsc_init_3rd;
  !!    }
  !!
  !!    a_ranks ->  broadcast_root_bt0 [lhead=cluster_1];
  !!
  !!    broadcast_root_0     -> broadcast_root_bt0
  !!    [ dir="both", color="blue", style="dashed"];
  !!
  !!    split_comm_world_0   -> split_comm_world_bt0
  !!    [ dir="both", color="blue", style="dashed"];
  !!
  !!    node [shape=ellipse, color=black];
  !!    open_stream_files [ label="Open stream files"];
  !!    {petsc_init_1st petsc_init_2nd petsc_init_3rd split_comm_world_0} -> open_stream_files;
  !!
  !!    subgraph cluster_2 {
  !!        color=royalblue;
  !!        style=rounded;
  !!        label = "MPI_COMM_WORLD - Rank 0 - 'Master'";
  !!
  !!        node [shape=ellipse, color=black];
  !!        send_activity_0 [label="Send Activity(ii)"]
  !!
  !!        subgraph cluster_3 {
  !!           color=seagreen3;
  !!           style=dashed;
  !!           label = "Worker loop";
  !!
  !!           node [shape=ellipse, color=black];
  !!           Irecv_activity_0 [label="Irecv Activity(ii)"];
  !!           waitany_activity_0 [label="WaitAny for IRecv Activity(ii)"];
  !!           send_activity_0 -> Irecv_activity_0 -> waitany_activity_0;
  !!           waitany_activity_0 ->  send_activity_0 [ color="darkgreen", style="dashed"];
  !!        }
  !!        waitall_activity_0 [label="WaitAll for IRecv Activity(ii)"];
  !!        waitany_activity_0 ->  waitall_activity_0;
  !!    }
  !!
  !!    open_stream_files -> send_activity_0 [lhead=cluster_2];
  !!
  !!    subgraph cluster_4 {
  !!        color=royalblue;
  !!        style=rounded;
  !!        label = "MPI_COMM_WORLD - Ranks > 0 - 'Workers'";
  !!
  !!        node [shape=ellipse, color=black];
  !!        subgraph cluster_3 {
  !!           color=seagreen3;
  !!           style=dashed;
  !!           label = "Worker loop";
  !!
  !!           node [shape=ellipse, color=black];
  !!           recv_activity_bt0 [label="Recv Active"]; 
  !!           send_activity_bt0 [label="Send Active"];
  !!
  !!           node [shape=component, color=royalblue];
  !!           worker_comm_1st_wl [label="1st domain worker_comm"];
  !!           worker_comm_2nd_wl [label="2nd domain worker_comm"];
  !!           worker_comm_3rd_wl [label="... domain worker_comm"];
  !!
  !!           recv_activity_bt0 -> {worker_comm_1st_wl worker_comm_2nd_wl worker_comm_3rd_wl};
  !!
  !!           node [shape=ellipse, color=black];
  !!           execute_single_d_2nd [label="Execute Single Domain"]
  !!           execute_single_d_3rd [label="Execute Single Domain"]
  !!           worker_comm_2nd_wl -> execute_single_d_2nd; 
  !!           worker_comm_3rd_wl -> execute_single_d_3rd;
  !!
  !!           subgraph cluster_5 {
  !!              color=black;
  !!              style=rounded;
  !!              label = "Execute Single Domain";
  !!              subgraph cluster_6 {
  !!                 color=royalblue;
  !!                 style=rounded;
  !!                 label = "WORKER_COMM - Rank 0";
  !!
  !!                 node [shape=ellipse, color=black];
  !!                 generate_geom [label="Generate Geometry"];
  !!                 serialize_pb  [label="Serialize Part-Branches"];
  !!                 send_pb       [label="Send Part-Branches"];
  !!                 bc_matsize_0  [label="Send - Broadcast Matrix Size"];
  !!                 generate_geom -> serialize_pb -> send_pb -> bc_matsize_0;
  !!              }
  !!
  !!              subgraph cluster_7 {
  !!                 color=royalblue;
  !!                 style=rounded;
  !!                 label = "WORKER_COMM - Ranks > 0";
  !!
  !!                 node [shape=ellipse, color=black];
  !!                 recv_pb        [label="Recieve Part-Branch"];
  !!                 deserialize_pb [label="Deserialize Part-Branch"];
  !!                 bc_matsize_bt0 [label="Recieve - Broadcast Matrix Size"];
  !!                 recv_pb -> deserialize_pb -> bc_matsize_bt0;
  !!              }
  !!
  !!              node [shape=ellipse, color=black];
  !!              petsc_mat_create_A [label = "PETSc - Create Stiffness A \n
  !!MatCreate \n MatSetSizes \n MatSetFromOptions \n MatSetUp(A) \n
  !!MatSetValues(A,...) \n  MatAssemblyBegin \n MatAssemblyEnd"];
  !!              bc_matsize_0   -> petsc_mat_create_A [ltail=cluster_6];
  !!              bc_matsize_bt0 -> petsc_mat_create_A [ltail=cluster_7];
  !!
  !!
  !!              subgraph cluster_8 {
  !!                 color=seagreen3;
  !!                 style=dashed;
  !!                 label = "Do ii = 1, 24";
  !!                 c_load_vec [label = "Create load vector F \n
  !!VecCreate \n VecSetSizes \n VecSetFromOptions \n VecSet(F, 0._rk,...)"];
  !!                 c_load_vec -> VecAssemblyBeginF;
  !!                 VecAssemblyBeginF -> c_load_vec [color="darkgreen", style="dashed"];
  !!              }
  !!              petsc_mat_create_A -> c_load_vec;
  !!
  !!              subgraph cluster_9 {
  !!                 color=seagreen3;
  !!                 style=dashed;
  !!                 label = "Do ii = 1, 24";
  !!                 VecAssemblyEndF;
  !!              }
  !!              VecAssemblyBeginF -> VecAssemblyEndF;
  !!
  !!              search_bound_branch [label= "Search Boundary-Branch of LC1"];
  !!              VecAssemblyEndF -> search_bound_branch;
  !!
  !!              allocate_gnid_zeros [label= "Allocate global node-id
  !!cross reference and \n zero values for dof elimination"];
  !!              search_bound_branch -> allocate_gnid_zeros;
  !!
  !!              setup_gnid [label= "Setup global node-id cross reference"];
  !!              allocate_gnid_zeros -> setup_gnid;
  !!
  !!              c_sol_vec [label = "Create solution vector X \n
  !!VecCreate \n VecSetSizes \n VecSetFromOptions \n VecSet(X, 0._rk,...)"];
  !!
  !!              VecSetValues [label="VecSetValues(X,....)
  !!Set Vector Values from \n LC1 Boundary-Branch Values"]
  !!
  !!              setup_gnid -> c_sol_vec -> VecSetValues;
  !!              VecSetValues -> VecAssemblyBeginX -> VecAssemblyEndX;
  !!
  !!              apply_boundaries [label = "Apply Dirichlet Boundaries by
  !!MatZeroRowsColumns with \n Boundary-Branch node-ids \n applied to AX=F"];
  !!              VecAssemblyEndX -> apply_boundaries;
  !!
  !!              Create_lin_solver [label="Create linear solver context
  !!KSPCreate \n KSPSetOperators(ksp, A, A,...) \n KSPSetTolerances \n KSPSetFromOptions"]
  !!              apply_boundaries ->  Create_lin_solver;
  !!
  !!              solve_lin_sys [label="Solve the linear system
  !!KSPSolve(ksp, F, X,...)"];
  !!              Create_lin_solver -> solve_lin_sys;
  !! 
  !!           }
  !!
  !!           worker_comm_1st_wl -> generate_geom [lhead=cluster_6];
  !!           worker_comm_1st_wl -> recv_pb [lhead=cluster_7];
  !!           solve_lin_sys -> send_activity_bt0  [ltail=cluster_5];
  !!
  !!           {execute_single_d_2nd execute_single_d_3rd} -> send_activity_bt0;
  !!           send_activity_bt0 ->  recv_activity_bt0 [ color="darkgreen", style="dashed"];
  !!        }
  !!        node [shape=ellipse];
  !!        petsc_finalize [label="PETSc finalize"];
  !!        send_activity_bt0 -> petsc_finalize;
  !!    }
  !!
  !!    open_stream_files -> recv_activity_bt0 [lhead=cluster_4];
  !!
  !!    send_activity_0 -> recv_activity_bt0
  !!    [ dir="both", constraint=false, color="blue", style="dashed"];
  !!
  !!    send_activity_bt0 -> {waitany_activity_0 waitall_activity_0}
  !!    [ dir="both", constraint=false, color="blue", style="dashed"];
  !!
  !!    node [shape=ellipse, color=black];
  !!    mpi_finalize [label="MPI finalize"];
  !!    {petsc_finalize waitall_activity_0} -> mpi_finalize;
  !!
  !!
  !!}
  !!\enddot
  !> <hr>
