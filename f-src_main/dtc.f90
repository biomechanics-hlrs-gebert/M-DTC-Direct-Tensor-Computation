!==============================================================================
!> \mainpage HLRS - Material structure process chain - Code Documentation
!> 
!> <hr>
!> \section desc Description
!>
!> A program collection organised in a tool-chain for the automatic 
!> determination of anisotropic elastic constants of microstructured 
!> materials from micro computer tomography.
!> The Flowchart of how it works can be seen here: \ref current.
!> <hr>
!> \section developers Developers 
!> 
!> Ralf Schneider, Xin Wang, Dmitry Khabi, Johannes Gebert
!>
!> FMPS is developed by Hans Wüstenberg
!>
!> for and with support from
!>
!> <b>High Performance Computing Center Stuttgart &mdash; HLRS</b><BR>
!> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Johannes Gebert<BR>
!> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Uwe Küster<BR>
!> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Dmitry Khabi<BR>
!> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Michael M. Resch
!>
!> <b>LASSO Ingenieurgesellschaft mbH</b><BR>
!> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Ulrich Hindenlang<BR>
!> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Artur Kurz<BR>
!>
!> <hr>
!> \section version Version
!>
!> 4.0 alpha
!>
!> <hr>
!> \section svn SVN 
!>
!> Versioning control is provided by FusionForge
!> https://gforge.hlrs.de/projects/struct-process
!>
!> <hr>
!> \section install Installation
!>
!> Installation is done by
!> <ol>
!>   <li> source setenv.sh <system>
!>   <li> make
!>
!> </ol>
!>
!> Where <system> can currently take the values "zeus, "julius" or "hawk".
!> Please not that in order to install the material structure process chain working copies of metis and petsc have to be installed on the executing system. Installation instructions for the dependencies can be found at this place: 
!> \ref prerequisites.
!>
!> <hr>
!> \section testing Testing
!>
!> Description of the provided regression and correctness checks can be found at this place: 
!> \ref regression.
!>
!> <hr>
!>
!==============================================================================

!------------------------------------------------------------------------------
!> Struct Process main programm
!------------------------------------------------------------------------------
!>  \section written Written by:
!>  Johannes Gebert
!>  Ralf Schneider
!------------------------------------------------------------------------------
Program dtc

USE iso_c_binding
USE global_std
USE user_interaction
USE mpi_user_interaction
USE formatted_plain
USE puredat 
USE meta
USE meta_puredat_interface
USE auxiliaries
USE chain_routines
USE MPI
USE decomp 
USE dtc_main_subroutines
USE PETSC
USE petsc_opt
USE system

Implicit None

INTEGER(ik), PARAMETER :: debug = 2   ! Choose an even integer!!

TYPE(materialcard) :: bone

INTEGER(mik) :: ierr, rank_mpi, size_mpi, petsc_ierr, &
    worker_rank_mpi, worker_size_mpi, status_un, comms_un, parts_un, &
    worker_comm, comms_dmn_un, runtime_un, &
    Active, request, finished = -1

INTEGER(mik), Dimension(no_streams) :: fh_mpi_worker
INTEGER(mik), Dimension(MPI_STATUS_SIZE)  :: status_mpi

INTEGER(c_int) :: stat_c_int

Type(tBranch) :: root, phi_tree
Type(tBranch), pointer :: ddc, meta_para, result_branch

CHARACTER, DIMENSION(4*mcl) :: c_char_array
CHARACTER, DIMENSION(:), ALLOCATABLE :: char_arr
CHARACTER(4*mcl), DIMENSION(:), ALLOCATABLE :: domain_path
CHARACTER(mcl)  , DIMENSION(:), ALLOCATABLE :: m_rry      

CHARACTER(4*mcl) :: job_dir
CHARACTER(mcl)   :: cmd_arg_history='', link_name = 'struct process', comms_dmn_file, &
    muCT_pd_path, muCT_pd_name, binary, status_file, groups_file, cores_file, &
    desc="", memlog_file="", type_raw="", restart='N', restart_cmd_arg='U', &
    ios="", map_sgmnttn=""
CHARACTER(10) :: probably_failed_char
CHARACTER(8) :: elt_micro, output
CHARACTER(3) :: file_status
CHARACTER(1) :: bin_sgmnttn=""

REAL(rk) :: strain, t_start, t_end, t_duration, time_wasted, now, comm_fin_time, &
    final_time, mem_critical=1._rk

INTEGER(ik), DIMENSION(:), ALLOCATABLE :: Domains, nn_D, Domain_status, &
    parts_list, comm_groups, comm_ranges, dmns_per_comm
    
INTEGER(ik), DIMENSION(:,:), ALLOCATABLE :: comms_array
INTEGER(ik), DIMENSION(3) :: xa_d=0, xe_d=0

INTEGER(ik) :: nn, ii, jj, kk, ll, dc, computed_domains = 0, comm_nn = 1, max_domains_per_comm, &
    No_of_domains, No_of_domains_files, path_count, activity_size=0, No_of_comm_groups, comms_feedback, &
    alloc_stat, fh_cluster_log, free_file_handle, stat, No_of_cores_requested, my_global_rank, &
    no_lc=0, nl=0, Domain, llimit, parts, elo_macro, vdim(3), groups_size, parts_size, &
    global_color, rank, no_of_comms, comm_floor, comm_ceiling, comm_color, round_robin_skip, mesh_p_per_dmn, &
    max_skip, comm_counter, dmn_status, add_leaf_pntr, probably_failed, prts, comms_array_ptr, corecount
INTEGER(pd_ik), DIMENSION(:), ALLOCATABLE :: serial_root
INTEGER(pd_ik), DIMENSION(no_streams) :: dsize

INTEGER(pd_ik) :: serial_root_size, add_leaves

LOGICAL :: success, status_exists, groups_exists, parts_exists, comms_dmn_exists, found_proper_prt_count, &
    heaxist, abrt = .FALSE., already_finished=.FALSE., skip_active = .TRUE.,  &
    create_new_header = .FALSE., fex=.TRUE., no_restart_required = .FALSE., runtime_exists=.FALSE.

!----------------------------------------------------------------------------
CALL mpi_init(ierr)
CALL print_err_stop(std_out, "MPI_INIT didn't succeed", ierr)

CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank_mpi, ierr)
CALL print_err_stop(std_out, "MPI_COMM_RANK couldn't be retrieved", ierr)

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size_mpi, ierr)
CALL print_err_stop(std_out, "MPI_COMM_SIZE couldn't be retrieved", ierr)

!------------------------------------------------------------------------------
! Rank 0 -- Init (Master) Process and broadcast init parameters 
!------------------------------------------------------------------------------
If (rank_mpi == 0) THEN

    If (size_mpi < 2) THEN
        mssg = "We need at least 2 MPI processes to execute this program."
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF 

    CALL Start_Timer("Init Process")


    !------------------------------------------------------------------------------
    ! Parse the command arguments
    !------------------------------------------------------------------------------
    CALL get_cmd_args(binary, in%full, restart_cmd_arg, cmd_arg_history, prts)

    IF (in%full=='') THEN
        CALL usage(binary)    

        !------------------------------------------------------------------------------
        ! On std_out since file of std_out is not spawned
        !------------------------------------------------------------------------------
        mssg = "No input file given"
        CALL print_err_stop_slaves(mssg); GOTO 1000
    ELSE
        IF(in%full(LEN_TRIM(in%full)-4 : LEN_TRIM(in%full)) /= ".meta") THEN
            mssg = "No meta file given."
            CALL print_err_stop_slaves(mssg); GOTO 1000
        END IF         
    END IF

    !------------------------------------------------------------------------------
    ! Check and open the input file; Modify the Meta-Filename / Basename
    ! Define the new application name first
    !------------------------------------------------------------------------------
    global_meta_prgrm_mstr_app = 'dtc' 
    global_meta_program_keyword = 'TENSOR_COMPUTATION'
    CALL meta_append(m_rry, size_mpi, binary, ios)
        
    IF(ios /= "") CALL print_err_stop(6, ios, 1)

    !------------------------------------------------------------------------------
    ! Redirect std_out into a file in case std_out is not useful by environment.
    ! Place these lines before handle_lock_file :-)
    !------------------------------------------------------------------------------
    CALL determine_std_fh(std_out, std_err)
     
    !------------------------------------------------------------------------------
    ! Spawn standard out/err after(!) the basename is known
    !------------------------------------------------------------------------------
    IF(std_out/=6) CALL meta_start_ascii(std_out, '.std_out')
    IF(std_err/=0) CALL meta_start_ascii(std_err, '.std_err')

    CALL show_title(["Dr.-Ing. Ralf Schneider (HLRS, NUM)", &
        "Johannes Gebert, M.Sc. (HLRS, NUM) "])

    IF(debug >=0) WRITE(std_out, FMT_TXT) "Post mortem info probably in ./datasets/temporary.std_out"
    WRITE(std_out, FMT_TXT) "Program invocation:"//TRIM(cmd_arg_history)          
    WRITE(std_out, FMT_TXT_SEP)

    !------------------------------------------------------------------------------
    ! Set input paths
    !------------------------------------------------------------------------------
    ! It is strongly recommended not to play around with these paths carelessly.
    ! Some of the dependencies are easily overlooked.
    !------------------------------------------------------------------------------
    muCT_pd_path = TRIM(in%path)
    muCT_pd_name = TRIM(in%bsnm)

    !------------------------------------------------------------------------------
    ! Output directory and 
    ! Implicitly creates a subdirectory.
    !------------------------------------------------------------------------------    
    outpath = TRIM(out%path)//TRIM(out%bsnm)//"/"
    project_name = "results"

    pro_path = outpath
    pro_name = project_name

    CALL meta_read('MICRO_ELMNT_TYPE' , m_rry, elt_micro, ios); CALL mest(ios, abrt)
    CALL meta_read('OUT_FMT'          , m_rry, output, ios); CALL mest(ios, abrt)
    CALL meta_read('RESTART'          , m_rry, restart, ios); CALL mest(ios, abrt)
    CALL meta_read('SIZE_DOMAIN'      , m_rry, bone%phdsize, ios); CALL mest(ios, abrt)
    CALL meta_read('SPACING'          , m_rry, bone%delta, ios); CALL mest(ios, abrt)
    CALL meta_read('DIMENSIONS'       , m_rry, vdim, ios); CALL mest(ios, abrt)
    CALL meta_read('LO_BNDS_DMN_RANGE', m_rry, xa_d, ios); CALL mest(ios, abrt)
    CALL meta_read('UP_BNDS_DMN_RANGE', m_rry, xe_d, ios); CALL mest(ios, abrt)
    CALL meta_read('BINARIZE_LO'      , m_rry, llimit, ios); CALL mest(ios, abrt)
    CALL meta_read('RVE_STRAIN'       , m_rry, strain, ios); CALL mest(ios, abrt)
    CALL meta_read('YOUNG_MODULUS'    , m_rry, bone%E, ios); CALL mest(ios, abrt)
    CALL meta_read('POISSON_RATIO'    , m_rry, bone%nu, ios); CALL mest(ios, abrt)
    CALL meta_read('MACRO_ELMNT_ORDER', m_rry, elo_macro, ios); CALL mest(ios, abrt)
    CALL meta_read('TYPE_RAW'         , m_rry, type_raw, ios); CALL mest(ios, abrt)
    CALL meta_read('BINARY_SEGMENTATION', m_rry, bin_sgmnttn, ios); CALL mest(ios, abrt)
    CALL meta_read('SEGMENTATION_MAP'   , m_rry, map_sgmnttn, ios); CALL mest(ios, abrt)

    IF (elo_macro == 1) THEN
        no_lc = 24_ik
    ELSE IF (elo_macro == 2) THEN
        no_lc = 60_ik
    END IF 

    IF ((bin_sgmnttn /= "Y") .AND. (bin_sgmnttn /= "y") .AND. &
        (bin_sgmnttn /= "N") .AND. (bin_sgmnttn /= "n")) THEN
        CALL print_err_stop(std_out, "Keyword BINARY_SEGMENTATION not recognized.", 1)
    END IF 

    IF ((restart == "Y") .OR. (restart == "y")) restart = "Y"
    IF ((restart == "N") .OR. (restart == "n")) restart = "N"

    IF (bin_sgmnttn == "y") bin_sgmnttn = "Y"

    IF (map_sgmnttn /= "Hounsfield") THEN
        CALL print_err_stop(std_out, "Keyword SEGMENTATION_MAP not recognized. &
            &Currently only 'Hounsfield' is implemented.", 1)
    END IF 

    IF(abrt) CALL print_err_stop(std_out, "A keyword error occured.", 1)

    !------------------------------------------------------------------------------
    ! Restart handling
    !------------------------------------------------------------------------------
    ! Standard lock file of MeRaDat not used, because it would interfere with the 
    ! status file of the struct process. 
    !------------------------------------------------------------------------------
    CALL meta_compare_restart(restart, restart_cmd_arg)

    IF (std_out/=6) CALL meta_start_ascii(std_out, '.std_out')
    IF (std_out/=0) CALL meta_start_ascii(std_err, '.std_err')

    !------------------------------------------------------------------------------
    ! Warning / Error handling
    !------------------------------------------------------------------------------
    IF ( (bone%phdsize(1) /= bone%phdsize(2)) .OR. (bone%phdsize(1) /= bone%phdsize(3)) ) THEN
        mssg = 'Currently, all 3 dimensions of the physical domain size must be equal!'
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF
    
    IF ( (xa_d(1) > xe_d(1)) .OR. (xa_d(2) > xe_d(2)) .or. (xa_d(3) > xe_d(3)) ) THEN
        mssg = 'Input parameter error: Start value of domain range larger than end value.'
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF

    !------------------------------------------------------------------------------
    ! Program breaks if the phdsize is not taking the boundary nodes into account (!).
    ! Therefore, the boundaries are calculated with + 2 Voxels
    !------------------------------------------------------------------------------
    IF ((bone%phdsize(1) > (vdim(1) + 2_ik)*bone%delta(1)) .OR. & 
        (bone%phdsize(2) > (vdim(2) + 2_ik)*bone%delta(2)) .OR. & 
        (bone%phdsize(3) > (vdim(3) + 2_ik)*bone%delta(3))) THEN
        mssg = 'The domains are larger than the field of view.'
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF

    !------------------------------------------------------------------------------
    ! Check the existence of the required files.
    !------------------------------------------------------------------------------
    ! Not in a function/subroutine, because it does not save a line. Filename and 
    ! error message have to be done here anyway.
    !------------------------------------------------------------------------------
    status_file = TRIM(out%p_n_bsnm)//".status"
    INQUIRE(FILE = TRIM(status_file), EXIST=status_exists, SIZE=activity_size)
    IF(.NOT. status_exists) THEN
        mssg = "No *.status file found. This version of the struct_process &
            &requires the morphometric evaluation (binary) and the topology &
            &aware scheduler (Python script). Those procedures create the status &
            &file. Only the round-robin struct_process can create a status file &
            &on its own."
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF 

    groups_file = TRIM(out%p_n_bsnm)//".groups"
    INQUIRE(FILE = TRIM(groups_file), EXIST=groups_exists, SIZE=groups_size)

    IF(.NOT. groups_exists) THEN
        CALL print_err_stop_slaves("No *.groups file found.")
    GOTO 1000
    END IF 

    ! groups_size is a 2D array with 3 rows 
    ! (1 for comm size, 1 for number of such comms and 1 for no of domains assigned to those comms)
    No_of_comm_groups = ANINT(REAL(groups_size, rk) / 8._rk / 3._rk)

    cores_file = TRIM(out%p_n_bsnm)//".cores"
    INQUIRE(FILE = TRIM(cores_file), EXIST=parts_exists, SIZE=parts_size)
    IF(.NOT. parts_exists) THEN
        CALL print_err_stop_slaves("No *.cores file found.")
        GOTO 1000
    END IF 

    ! comm_of_domain -> which communicator computed which domain
    comms_dmn_file = TRIM(out%p_n_bsnm)//".comms"
    INQUIRE(FILE = TRIM(comms_dmn_file), EXIST=comms_dmn_exists)
    IF((comms_dmn_exists) .AND. (restart == "N")) THEN
        CALL print_err_stop_slaves("*.comms file already exists, no restart requested.")
        GOTO 1000
    END IF 

    !------------------------------------------------------------------------------
    ! Check the condition of the required files.
    !------------------------------------------------------------------------------
    IF(activity_size .NE. parts_size) THEN
        mssg = 'The size of the *.parts and the *.status files do not match.'
        CALL print_err_stop_slaves(mssg)
        GOTO 1000
    END IF 

    No_of_domains = (xe_d(1)-xa_d(1)+1) * (xe_d(2)-xa_d(2)+1) * (xe_d(3)-xa_d(3)+1)
    No_of_domains_files = ANINT(REAL(activity_size, rk) / 8._rk)

    IF(No_of_domains .NE. No_of_domains_files) THEN
        mssg = 'The size of the auxiliary files and the requestd domain range do not match.'
        CALL print_err_stop_slaves(mssg)
        GOTO 1000
    END IF 

END IF 

CALL MPI_BCAST(status_file, INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(groups_file,  INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(cores_file,  INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(comms_dmn_file,  INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)

CALL MPI_BCAST(No_of_domains,     1_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(No_of_comm_groups, 1_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(prts,              1_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)

CALL MPI_BCAST(type_raw, INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)

!------------------------------------------------------------------------------
! Allocate and init field for selected domain range
! Required for the worker_is_active tracking
! Required for steering the allocation of memory for writing data into 
! the stream files (raise leaves etc.)
!------------------------------------------------------------------------------
ALLOCATE(Domains(No_of_domains), stat=alloc_stat)
CALL alloc_err("Domains", alloc_stat)
Domains = 0_ik


!------------------------------------------------------------------------------
! Initial setup of the domain status tracker
!------------------------------------------------------------------------------
! Tracker -100000000-Domain_No   | Domain is planned, but not started computing
! Tracker           -Domain_No   | Domain is started computing, but not finished
! Tracker            Domain_No   | Domain is finished
!
! If the position within the tracking file is implicitely connected to the 
! domain number, changes to the domain range before restart will crash the 
! computations (!)
! This is even more dangerous for topology aware approaches
!------------------------------------------------------------------------------
ALLOCATE(Domain_status(No_of_domains), stat=alloc_stat)
CALL alloc_err("Domain_status", alloc_stat)

CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(status_file), MPI_MODE_RDWR, MPI_INFO_NULL, status_un, ierr)
CALL MPI_FILE_READ(status_un, Domain_status, INT(SIZE(Domain_status), mik), MPI_INTEGER8, MPI_STATUS_IGNORE, ierr)
CALL MPI_FILE_CLOSE(status_un, ierr)

!------------------------------------------------------------------------------
! Read the communicators for use with the dataset
!------------------------------------------------------------------------------
ALLOCATE(comms_array(3_ik,No_of_comm_groups), stat=alloc_stat)
CALL alloc_err("comms_array", alloc_stat)

ALLOCATE(comm_groups(No_of_comm_groups), stat=alloc_stat)
CALL alloc_err("comm_groups", alloc_stat)

ALLOCATE(comm_ranges(No_of_comm_groups), stat=alloc_stat)
CALL alloc_err("comm_ranges", alloc_stat)

ALLOCATE(dmns_per_comm(No_of_comm_groups), stat=alloc_stat)
CALL alloc_err("dmns_per_comm", alloc_stat)

CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(groups_file), MPI_MODE_RDONLY, MPI_INFO_NULL, comms_un, ierr)
CALL MPI_FILE_READ(comms_un, comms_array, INT(SIZE(comms_array), mik), MPI_INTEGER8, MPI_STATUS_IGNORE, ierr)
CALL MPI_FILE_CLOSE(comms_un, ierr)


!------------------------------------------------------------------------------
! Read the number of parts assigned to each domain (ppd). The ppd must fit
! to the communicators read from the *.comms file.
!------------------------------------------------------------------------------
ALLOCATE(parts_list(No_of_domains), stat=alloc_stat)
CALL alloc_err("parts_list", alloc_stat)
parts_list = 0_ik

CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(cores_file), MPI_MODE_RDONLY, MPI_INFO_NULL, parts_un, ierr)
CALL MPI_FILE_READ(parts_un, parts_list, INT(SIZE(parts_list), mik), MPI_INTEGER8, MPI_STATUS_IGNORE, ierr)
CALL MPI_FILE_CLOSE(parts_un, ierr)

ALLOCATE(domain_path(0:No_of_domains))
domain_path = ''

CALL MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(comms_dmn_file), MPI_MODE_RDWR, MPI_INFO_NULL,  comms_dmn_un, ierr)

!------------------------------------------------------------------------------
! Check if the computation is already done or not.
!------------------------------------------------------------------------------
If (rank_mpi == 0) THEN

    comms_array_ptr = 0_ik

    IF (prts == 0_ik) THEN

        ! +1 to account for the master process of the monolihtic topology aware approach
        No_of_cores_requested = SUM(comms_array(1,:)*comms_array(2,:)) + 1_ik

        ! Worst case assumption -> Max number of domains that may occur.
        max_domains_per_comm = MAXVAL(comms_array(3,:))


    ELSE
        No_of_cores_requested = size_mpi - 1_ik

        DO ii=1, SIZE(comms_array(1,:))
            IF (comms_array(1,ii ) == prts) comms_array_ptr = ii
            
        END DO

        ! Worst case assumption -> Max number of domains that may occur.
        max_domains_per_comm = comms_array(3,comms_array_ptr)
    END IF 
        
    WRITE(std_out, FMT_TXT_AxI0) "Groups of cores per domain (cpd):", comms_array(1,:)
    WRITE(std_out, FMT_TXT_AxI0) "Number of such communicators:    ", comms_array(2,:)
    WRITE(std_out, FMT_TXT_AxI0) "Number of domains of these cpds: ", comms_array(3,:)
    WRITE(std_out, FMT_TXT_AxI0) "Number of cores requested:       ", No_of_cores_requested
    WRITE(std_out, FMT_TXT_AxI0) "prts:                            ", prts
    
    IF ((No_of_cores_requested .NE. size_mpi) .AND. (prts == 0_ik)) THEN
        mssg = "The number of cores requested /= size_mpi and parts == 0_ik."
        CALL print_err_stop_slaves(mssg)
        GOTO 1000
    ELSE
        found_proper_prt_count = .FALSE.
        DO ii=1,17
            corecount = 2**ii
            IF (corecount == prts) found_proper_prt_count = .TRUE.
        END DO

        IF (.NOT. found_proper_prt_count) THEN
            mssg = "The number of cores given is not of a power of 2."
            CALL print_err_stop_slaves(mssg)
            GOTO 1000
        END IF 

    END IF 

    DO ii=1, SIZE(Domain_status)
        IF (Domain_status(ii) >= 0) computed_domains = computed_domains + 1 
    END DO

    IF (No_of_domains == computed_domains) THEN 
        already_finished=.TRUE.
        
        mssg = "Job is already finished. No restart required."
        CALL print_err_stop_slaves(mssg, "message")
        GOTO 1000
    END IF 

    !------------------------------------------------------------------------------
    ! Check if INPUT header exists. If not --> create one.
    !------------------------------------------------------------------------------
    INQUIRE(FILE = TRIM(in%p_n_bsnm)//".head", EXIST=heaxist)
    IF(.NOT. heaxist) THEN
        free_file_handle = give_new_unit()
        CALL convert_meta_to_puredat(free_file_handle, m_rry)
    END IF

    !------------------------------------------------------------------------------
    ! Check if OUTPUT header exists. If not --> Internally reset to restart='N'.
    !------------------------------------------------------------------------------
    INQUIRE(FILE = TRIM(pro_path)//TRIM(pro_name)//".head", EXIST=heaxist)


    !------------------------------------------------------------------------------
    ! Raise and build meta_para tree
    ! Hardcoded, implicitly given order of the leafs. 
    ! DO NOT CHANGE ORDER WITHOUT MODIFYING ALL OTHER INDICES REGARDING »meta_para«
    !------------------------------------------------------------------------------
    Allocate(meta_para)
    CALL raise_tree("Input parameters", meta_para)

    CALL add_leaf_to_branch(meta_para, "muCT puredat pro_name" , mcl , str_to_char(muCT_pd_name)) 
    CALL add_leaf_to_branch(meta_para, "muCT puredat pro_path" , mcl , str_to_char(muCT_pd_path)) 
    CALL add_leaf_to_branch(meta_para, "Physical domain size"                 , 3_ik, bone%phdsize) 
    CALL add_leaf_to_branch(meta_para, "Lower bounds of selected domain range", 3_ik, xa_d) 
    CALL add_leaf_to_branch(meta_para, "Upper bounds of selected domain range", 3_ik, xe_d)      
    
    CALL add_leaf_to_branch(meta_para, "Grid spacings"                 , 3_ik, bone%delta) 
    CALL add_leaf_to_branch(meta_para, "Lower limit of iso value"      , 1_ik, [llimit])      
    CALL add_leaf_to_branch(meta_para, "Element type  on micro scale"  , len(elt_micro) , str_to_char(elt_micro))      
    CALL add_leaf_to_branch(meta_para, "Output Format"         , len(output)    , str_to_char(output)) 
    CALL add_leaf_to_branch(meta_para, "Average strain on RVE" , 1_ik, [strain])    

    CALL add_leaf_to_branch(meta_para, "Young_s modulus"               , 1_ik, [bone%E]) 
    CALL add_leaf_to_branch(meta_para, "Poisson_s ratio"               , 1_ik, [bone%nu]) 
    CALL add_leaf_to_branch(meta_para, "Element order on macro scale"  , 1_ik, [elo_macro]) 
    CALL add_leaf_to_branch(meta_para, "Output amount"                 , len(out_amount), str_to_char(out_amount)) 
    add_leaf_pntr = 14_ik
    
    CALL add_leaf_to_branch(meta_para, "Restart", 1_ik, str_to_char(restart(1:1))) 

    CALL add_leaf_to_branch(meta_para, "Number of voxels per direction", 3_ik , vdim) 
    CALL add_leaf_to_branch(meta_para, "Domains per communicator", 1_ik, [max_domains_per_comm]) 

    CALL add_leaf_to_branch(meta_para, "Binary segmentation" , LEN(bin_sgmnttn) , str_to_char(bin_sgmnttn))      
    CALL add_leaf_to_branch(meta_para, "Segmentation map"    , LEN(map_sgmnttn) , str_to_char(map_sgmnttn))      

    !------------------------------------------------------------------------------
    ! Prepare output directory via Calling the c function.
    ! Required, because INQUIRE only acts on files, not on directories.
    ! File exists if stat_c_int = 0 
    !------------------------------------------------------------------------------
    c_char_array(1:LEN(TRIM(outpath)//CHAR(0))) = str_to_char(TRIM(outpath)//CHAR(0))
    CALL Stat_Dir(c_char_array, stat_c_int)

    IF(stat_c_int /= 0) THEN

        CALL exec_cmd_line("mkdir -p "//TRIM(outpath), stat, 20)

        IF(stat /= 0) THEN
            mssg = 'Could not execute syscall »mkdir -p '//trim(outpath)//'«.'
            CALL print_err_stop_slaves(mssg); GOTO 1000
        END IF
        
        CALL Stat_Dir(c_char_array, stat_c_int)

        IF(stat_c_int /= 0) THEN
            mssg = 'Could not create the output directory »'//TRIM(outpath)//'«.'
            CALL print_err_stop_slaves(mssg); GOTO 1000
        END IF
    ELSE 
        WRITE(std_out, FMT_MSG) "Reusing the output directory"
        WRITE(std_out, FMT_MSG) TRIM(outpath)
        WRITE(std_out, FMT_SEP)
    END IF

    CALL link_start(link_name, .FALSE., .FALSE., success)

    IF (.NOT. success) THEN
        mssg = "Something went wrong during link_start"
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF

    IF((restart == 'Y') .AND. (.NOT. heaxist)) create_new_header=.TRUE.
    !------------------------------------------------------------------------------
    ! The Output name normally is different than the input name.
    ! Therefore, an existing header implies a restart.
    !------------------------------------------------------------------------------
    ! Basically a completely new computation
    !------------------------------------------------------------------------------
    IF ((restart == 'N') .OR. (create_new_header)) THEN 

        !------------------------------------------------------------------------------
        ! project_name --> out%p_n_bsnm/bsnm --> subdirectory with file name = bsnm.suf
        !------------------------------------------------------------------------------
        IF ((status_exists) .AND. (.NOT. create_new_header)) THEN
            mssg = "No restart requested, but the file '"&
                //TRIM(out%p_n_bsnm)//".status' already exists."
            CALL print_err_stop_slaves(mssg, "warning"); GOTO 1000
        END IF 

        !------------------------------------------------------------------------------
        ! project_name --> out%p_n_bsnm/bsnm --> subdirectory with file name = bsnm.suf
        !------------------------------------------------------------------------------
        ! Tree which is fed back = root. A collection of stream paths and Null pointers
        !------------------------------------------------------------------------------
        CALL raise_tree(Trim(project_name), root)

        !------------------------------------------------------------------------------
        ! Source branch / target branch
        !------------------------------------------------------------------------------
        CALL include_branch_into_branch(s_b=meta_para, t_b=root, blind=.TRUE.)

        !------------------------------------------------------------------------------
        ! Read an existing input tree (with microfocus ct data).
        !------------------------------------------------------------------------------
        ! Load puredat tree of micro-CT data and calculate the global
        ! parameters of the domain decomposition
        !------------------------------------------------------------------------------
        pro_path = muCT_pd_path
        pro_name = muCT_pd_name

        phi_tree = read_tree()

        !------------------------------------------------------------------------------
        ! Set project name and path of global domain decomposition     
        !------------------------------------------------------------------------------
        pro_path = outpath
        pro_name = project_name

        allocate(ddc)
        ddc = calc_general_ddc_params(bone%phdsize, meta_para)
        
        CALL include_branch_into_branch(s_b=ddc, t_b=root, blind=.TRUE.)
        
    ELSE ! restart = 'Y'         

        !------------------------------------------------------------------------------
        ! Read an existing output tree (with microfocus ct data).
        ! Characteristics of PureDat and MeRaDat require this chicken dance...
        !------------------------------------------------------------------------------
        root = read_tree()

        If (out_amount == "DEBUG") THEN 
            WRITE(fh_log, fmt_dbg_sep)
            Write(fh_log,'(A)') "root right after restart read"
            CALL log_tree(root, fh_log,.FALSE.)
            WRITE(fh_log, fmt_dbg_sep)
            flush(fh_log)
        END If
       
        !------------------------------------------------------------------------------
        ! This Calling sequence is only valid since "Averaged Material 
        ! Properties" only contains r8 data added at the end of the 
        ! r8-stream. More correct would be a routine that ensures data
        ! integrity and compresses potentially missing stream data in
        ! an efficient way.
        !------------------------------------------------------------------------------
        CALL get_stream_size(root, dsize)
        root%streams%dim_st = dsize
        root%streams%ii_st  = dsize + 1

        CALL read_streams(root)

        CALL connect_pointers(root%streams, root)

        CALL search_branch("Global domain decomposition", root, ddc, success, out_amount)

        IF (.NOT. success) THEN
            mssg = "No branch named 'Global domain decomposition', however a restart was requested."
            CALL print_err_stop_slaves(mssg); GOTO 1000
        END IF

        CALL search_branch("Input parameters", root, meta_para, success, out_amount)

        IF (.NOT. success) then
            mssg = "No branch named 'Input parameters', however a restart was requested."
            CALL print_err_stop_slaves(mssg); GOTO 1000
        END IF

        !------------------------------------------------------------------------------
        ! Reset Output amount and Restart in loaded param branch
        ! Hardcoded, implicitly given order of the leaves. 
        ! DO NOT CHANGE INDICES WITHOUT MODYFING THE »add_leaf_to_branch« SEQUENCES.
        !------------------------------------------------------------------------------
        meta_para%leaves(add_leaf_pntr)%p_char = str_to_char(out_amount)
        meta_para%leaves(add_leaf_pntr+1_ik)%p_char = "Y"

        If (out_amount == "DEBUG") THEN 
            Write(fh_log,fmt_dbg_sep)
            Write(fh_log,'(A)') "root right after restart and deletion of avg mat props branch"
            CALL log_tree(root, fh_log,.FALSE.)
            Write(fh_log, fmt_dbg_sep)
            flush(fh_log)
        END If

    END IF ! restart == Yes/No

    !------------------------------------------------------------------------------
    ! Ensure to update the number of leaves and the indices of these in every
    ! line of code! Also update dat_ty and dat_no in "CALL raise_leaves"
    !------------------------------------------------------------------------------
    add_leaves = 24_pd_ik

    !------------------------------------------------------------------------------
    ! The name of the branche is a bit outdated. In the meantime, it contains 
    ! more data to accomodate the need for the doctoral project.
    ! It was done this way because it is the quickest and easiest way for 
    ! combining domain specfic output data.
    !------------------------------------------------------------------------------
    CALL add_branch_to_branch(root, result_branch)
    CALL raise_branch("Averaged Material Properties", 0_pd_ik, add_leaves, result_branch)


    ! for better formatting :-)
    nl = no_lc
    CALL raise_leaves(no_leaves = add_leaves, &
        desc = [ & ! DO NOT CHANGE THE LENGTH OF THE STRINGS      ! Leaf x bytes
        "Domain number                                     " , &  !  1 x  1
        "Number of Elements                                " , &  !  2 x  1
        "Number of Nodes                                   " , &  !  3 x  1
        "Collected logs                                    " , &  !  4 x  24
        "Start Time                                        " , &  !  5 x  1
        "Duration                                          " , &  !  6 x  1
        "Domain forces                                     " , &  !  7 x 24*24
        "Effective numerical stiffness                     " , &  !  8 x 24*24
        "Symmetry deviation - effective numerical stiffness" , &  !  9 x  1
        "Averaged stresses                                 " , &  ! 10 x  6*24
        "Averaged strains                                  " , &  ! 11 x  6*24
        "Effective stiffness                               " , &  ! 12 x  6* 6  
        "Symmetry deviation - effective stiffness          " , &  ! 13 x  1
        "Averaged Effective stiffness                      " , &  ! 14 x  6* 6
        "Symmetry deviation - Averaged effective stiffness " , &  ! 15 x  1
        "Rotation Angle CR_1                               " , &  ! 16 x  1
        "Rotation Vector CR_1                              " , &  ! 17 x  3
        "Final coordinate system CR_1                      " , &  ! 18 x  9
        "Optimized Effective stiffness CR_1                " , &  ! 19 x  6* 6
        "Rotation Angle CR_2                               " , &  ! 20 x  1
        "Rotation Vector CR_2                              " , &  ! 21 x  3
        "Final coordinate system CR_2                      " , &  ! 22 x  9
        "Optimized Effective stiffness CR_2                " , &  ! 23 x  6* 6
        "Effective density                                 "], &  ! 24 x  1
        dat_ty = [(4_1, ii=1, 4), (5_1, ii=5, add_leaves)], &
        dat_no = [ max_domains_per_comm, &
        max_domains_per_comm,         max_domains_per_comm, &
        max_domains_per_comm * 24   ,                                                    &
        max_domains_per_comm,         max_domains_per_comm, &                                  
        max_domains_per_comm * nl*nl, max_domains_per_comm * nl*nl, max_domains_per_comm       , &
        max_domains_per_comm *  6*nl, max_domains_per_comm *  6*nl, max_domains_per_comm *  6*6, &
        max_domains_per_comm        , max_domains_per_comm *  6* 6, max_domains_per_comm       , &
        max_domains_per_comm        , max_domains_per_comm *     3, max_domains_per_comm *    9, &
        max_domains_per_comm *  6* 6, max_domains_per_comm        , max_domains_per_comm *    3, &
        max_domains_per_comm *     9, max_domains_per_comm *  6* 6, max_domains_per_comm      ], &
        branch = result_branch)

    result_branch%leaves(:)%pstat = -1
    
    CALL set_bounds_in_branch(root, root%streams)
     
    If (out_amount == "DEBUG") THEN 
        Write(fh_log, FMT_DBG_SEP)
        Write(fh_log, '(A)') "root right before serialisation"
        CALL log_tree(root,fh_log,.True.)
        Write(fh_log, FMT_DBG_SEP)
        flush(fh_log)
    END If

    !------------------------------------------------------------------------------
    ! serialize root branch
    !------------------------------------------------------------------------------
    CALL serialize_branch(root, serial_root, serial_root_size, .TRUE.)

    !------------------------------------------------------------------------------
    ! Get number of domains on each axis to broadcast the information to all 
    ! workers. This is done globally to avoid overly complex MPI send/receive calls
    ! for targeting the exact rank of the worker masters, which is not transparent 
    ! due to the inhomogeneous comm sizes.
    !------------------------------------------------------------------------------
    CALL pd_get(ddc,"nn_D", nn_D)

    CALL End_Timer("Init Process")

    !------------------------------------------------------------------------------
    ! Start Workers
    !------------------------------------------------------------------------------
    CALL Start_Timer("Broadcast Init meta_para")


END IF ! (rank_mpi == 0)

CALL MPI_BCAST(outpath      , INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(project_name , INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)

CALL MPI_BCAST(serial_root_size, 1_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(max_domains_per_comm, 1_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)

CALL MPI_BCAST(dmns_per_comm, INT(SIZE(dmns_per_comm),mik), MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)

IF(rank_mpi == 0) THEN

    CALL mpi_bcast(serial_root , INT(serial_root_size, mik), MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)

    deallocate(serial_root)
    CALL End_Timer("Broadcast Init meta_para")

    !------------------------------------------------------------------------------
    ! Execute collective mpi_comm_split. Since mpi_comm_world rank 0 is
    ! the head master, the worker_comm is not needed. It mustn't be part of 
    ! any worker communicator. With MPI_UNDEFINED passed as
    ! color, the worker_comm gets the value MPI_COMM_NULL
    !------------------------------------------------------------------------------
    CALL MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, rank_mpi, worker_comm, ierr)
    IF (ierr /=0) THEN
        mssg = "MPI_COMM_SPLIT couldn't split MPI_COMM_WORLD"
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF 

!------------------------------------------------------------------------------
! Ranks > 0 -- Worker slaves
!------------------------------------------------------------------------------
ELSE
    CALL Start_Timer("Broadcast Init meta_para")
     
    pro_path = outpath
    pro_name = project_name

    !------------------------------------------------------------------------------
    ! Serial_root_size == -1 ==> Signal that No_of_domains < size_mpi-1
    !------------------------------------------------------------------------------
    If (serial_root_size == -1) then
        CALL End_Timer("Broadcast Init meta_para")
        Goto 1000
    End If
    
    Allocate(serial_root(serial_root_size))
    
    CALL mpi_bcast(serial_root, INT(serial_root_size, mik), MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)

    CALL End_Timer("Broadcast Init meta_para")

    !------------------------------------------------------------------------------
    ! Deserialize root branch
    !------------------------------------------------------------------------------
    CALL Start_Timer("Deserialize root branch")

    CALL deserialize_branch(root, serial_root, .TRUE.)
    CALL assign_pd_root(root)
    CALL set_bounds_in_branch(root, root%streams)
    deallocate(serial_root)

    !------------------------------------------------------------------------------
    ! Extract input parameters
    !------------------------------------------------------------------------------
    Call Search_branch("Input parameters", root, meta_para, success)
    CALL search_branch("Global domain decomposition", root, ddc, success)

    !------------------------------------------------------------------------------
    ! Commented out since out_amount is a parametrized global variable 
    ! CALL pd_get(root%branches(1),"Output amount", char_arr)
    ! out_amount = char_to_str(char_arr)
    ! deallocate(char_arr)
    !------------------------------------------------------------------------------
    CALL pd_get(meta_para,"Restart", char_arr)
    restart = char_to_str(char_arr)
    deallocate(char_arr)

    CALL End_Timer("Deserialize root branch")


    !------------------------------------------------------------------------------
    ! Init Domain Cross Reference
    !------------------------------------------------------------------------------
    CALL pd_get(meta_para, "Lower bounds of selected domain range", xa_d, 3)
    CALL pd_get(meta_para, "Upper bounds of selected domain range", xe_d, 3)

    CALL pd_get(ddc, "nn_D", nn_D)
END IF 

!------------------------------------------------------------------------------
! Init Domain Cross Reference and domain paths.
!------------------------------------------------------------------------------
dc = 0
domain_path(0) = Trim(pro_path)//Trim(pro_name)//"_domains" ! "_results"

nn = 1
Do kk = xa_d(3), xe_d(3)
Do jj = xa_d(2), xe_d(2)
Do ii = xa_d(1), xe_d(1)
    dc         = dc + 1_mik
    path_count = dc / 2_mik

    Write(domain_path(dc),'(A,"/",I0)') Trim(domain_path(path_count)),dc

    Domains(nn) = ii + jj * (nn_D(1)+1) + kk * (nn_D(1)+1) * (nn_D(2)+1)

    nn = nn + 1_ik
End Do
End Do
End Do

!------------------------------------------------------------------------------
! During the computation of the domains by themselves, the main rank is not 
! engaged anymore. But before computing the domains, MPI_COMM_WORLD gets split
! into the subcomms, requested by the morphometric evaluation and the 
! topology aware scheduling.
!------------------------------------------------------------------------------
IF(rank_mpi /= 0) THEN  

    CALL pd_get(root%branches(1), "Domains per communicator", max_domains_per_comm)

    !------------------------------------------------------------------------------
    ! Get the communicator, the rank is part of.
    ! Color of MPI_COMM_SPLIT --> number of the comm
    !------------------------------------------------------------------------------
    my_global_rank = rank_mpi

    jj = 1_ik
    kk = 1_ik
    ll = 1_ik
    rank = 0_ik

    IF (parts == 0_ik) THEN
        comm_ceiling = comms_array(1,kk)

        no_of_comms = SUM(comms_array(2,:))

        DO ii = 1, no_of_comms ! ii -> color 1 to no_of_comms

            ! comms_array(1,kk) --> Size of the communicator (parts per domain)
            rank = rank + comms_array(1,kk)

            comm_floor = comm_ceiling
            comm_ceiling = rank

            ! comms_array(2,kk) --> Numbers of corresponding communicaotrs
            IF(jj==comms_array(2,kk)) THEN
                jj = 0_ik
                kk = kk + 1_ik
                ll = 1_ik
            END IF 

            jj = jj + 1_ik

            IF ((rank_mpi > comm_floor) .AND. (rank_mpi <= comm_ceiling)) THEN

                ! Color for MPI_Comm_Split
                global_color = ii

                ! Color within the groups of communicators with the same size to avoid race conditions.
                comm_color = ll

                EXIT
            END IF

            ll = ll + 1_ik
            
        END DO

        ! Makeshift workaround to prevent race conditions
        CALL SLEEP(comm_color)
    ELSE
        no_of_comms = INT((size_mpi-1_ik) / prts)

        DO ii = 1, no_of_comms ! ii -> color 1 to no_of_comms
            comm_floor = ii * prts
            comm_ceiling = ((ii+1_ik) * prts) - 1_ik

            IF ((rank_mpi >= comm_floor) .AND. (rank_mpi <= comm_ceiling)) THEN

                ! Color for MPI_Comm_Split
                global_color = ii

                ! Color within the groups of communicators with the same size to avoid race conditions.
                comm_color = rank_mpi
            END IF 
        END DO

    END IF 

    !------------------------------------------------------------------------------
    ! All Worker Ranks -- Init worker Communicators
    !------------------------------------------------------------------------------
    CALL MPI_Comm_split(MPI_COMM_WORLD, INT(global_color,mik), rank_mpi, worker_comm, ierr)
    CALL print_err_stop(std_out, "MPI_COMM_SPLIT couldn't split MPI_COMM_WORLD", ierr)
    
    CALL MPI_COMM_RANK(WORKER_COMM, worker_rank_mpi, ierr)
    CALL print_err_stop(std_out, "MPI_COMM_RANK couldn't retrieve worker_rank_mpi", ierr)

    CALL MPI_COMM_SIZE(WORKER_COMM, worker_size_mpi, ierr)
    CALL print_err_stop(std_out, "MPI_COMM_SIZE couldn't retrieve worker_size_mpi", ierr)

    !------------------------------------------------------------------------------
    ! This sets the options for PETSc in-core. To alter the options
    ! add them in Set_PETSc_Options in Module pets_opt in file
    ! f-src/mod_parameters.f90
    !------------------------------------------------------------------------------
    CALL Set_PETSc_Options()

    PETSC_COMM_WORLD = worker_comm

    CALL PetscInitialize(PETSC_NULL_CHARACTER, petsc_ierr)
End If

!------------------------------------------------------------------------------
! Global ranks > 0 -- Workers
!------------------------------------------------------------------------------
! In the struct process, the first round of distributing domains to 
! MPI communicators is done manually. In this f-src, it is incorporated 
! into the general routine for assigning domains to communicators.
!------------------------------------------------------------------------------
IF (rank_mpi /= 0) THEN

    IF (worker_rank_mpi == 0) THEN

        CALL get_stream_size(root, dsize)

        !------------------------------------------------------------------------------
        ! Extend project_name and outpath with rank - only for domain specific master!
        !------------------------------------------------------------------------------
        WRITE(outpath,'(A,A,I7.7,A)')    TRIM(outpath), "Rank_", rank_mpi, "/"
        WRITE(project_name,'(A,A,I7.7)') TRIM(project_name), "_", rank_mpi

        !------------------------------------------------------------------------------
        ! Prepare output directory via CALLing the c function.
        ! File exists if stat_c_int = 0 
        !------------------------------------------------------------------------------
        c_char_array(1:LEN(TRIM(outpath)//CHAR(0))) = str_to_char(TRIM(outpath)//CHAR(0))
        CALL Stat_Dir(c_char_array, stat_c_int)

        IF(stat_c_int /= 0) THEN

            CALL exec_cmd_line("mkdir -p "//TRIM(outpath), stat)

            IF(stat /= 0) CALL print_err_stop(std_out, &
                'Could not execute syscall »mkdir -p '//trim(outpath)//'«.', 1)

            CALL Stat_Dir(c_char_array, stat_c_int)

            IF(stat_c_int /= 0) THEN
                mssg = 'Could not create the output directory »'//TRIM(outpath)//'«.'
                CALL print_err_stop(std_out, mssg, 1)
            END IF
        END IF

        CALL link_start(link_name, .FALSE., .FALSE.)

    END IF 

    CALL MPI_BCAST(outpath,      INT(mcl, mik), MPI_CHAR, 0_mik, WORKER_COMM, ierr)
    CALL MPI_BCAST(project_name, INT(mcl, mik), MPI_CHAR, 0_mik, WORKER_COMM, ierr)

    !------------------------------------------------------------------------------
    ! Open Stream files to write data to master-worker rank directories
    !------------------------------------------------------------------------------
    pro_path = outpath
    pro_name = project_name    

    CALL set_stream_filenames(root%streams)
    CALL Open_Stream_Files(root%streams, "write", "new", fh_mpi_worker, WORKER_COMM)

    CALL MPI_FILE_OPEN(WORKER_COMM, TRIM(status_file), MPI_MODE_RDWR, MPI_INFO_NULL, status_un, ierr)
    CALL MPI_FILE_READ(status_un, Domain_status, INT(SIZE(Domain_status), mik), MPI_INTEGER8, MPI_STATUS_IGNORE, ierr)

    IF ((worker_rank_mpi == 0) .AND. (restart == 'Y')) THEN
        !------------------------------------------------------------------------------
        ! Read and check the Local domain number list in case of a restart
        ! restart =  'Y' is given by user and may be given without any domain computed
        ! so far. Therefore, only increment comm_nn if this entry no 0
        !------------------------------------------------------------------------------
        CALL MPI_FILE_READ_AT(fh_mpi_worker(4), & 
            Int(root%branches(3)%leaves(1)%lbound-1+(max_domains_per_comm-1-1), MPI_OFFSET_KIND), &
            comm_nn, 1_mik, MPI_INTEGER8, MPI_STATUS_IGNORE, ierr)
    END IF 

    CALL MPI_BCAST(comm_nn, 1_mik , MPI_INTEGER8, 0_mik, WORKER_COMM, ierr)
    
    parts = worker_size_mpi   

    !------------------------------------------------------------------------------
    ! Distribute the first half of the batch of domains per size of communicator in a 
    ! round-robin manner. The second half is implemented on a first-come first-serve basis.
    ! This may reduce race conditions and improve the distribution of the domains to 
    ! communicators that fit their demand of processors in an efficient way. If we distribute
    ! the domains in a round-robin manner all the time, then many communicators may wait for 
    ! the last to complete which can scrap the whole gain of performance.
    !------------------------------------------------------------------------------
    max_skip = FLOOR(REAL(max_domains_per_comm)/2._rk)

    nn = 0_ik
    comm_counter = -1_ik
    round_robin_skip = 0_ik
    DO  WHILE (nn < No_of_domains) 

        IF (worker_rank_mpi == 0) THEN
            
            ! Iterate domain numbers
            nn = nn + 1_mik

            ! Only »accept« domains if they need as many parts as the communicator has processors
            IF (parts_list(nn) == worker_size_mpi) THEN
                comm_counter = comm_counter + 1_ik
            ELSE
                CYCLE
            END IF 

            ! The first batch of domains within a specific comm group are assigned to their 
            ! corresponding communicator in the order as they are mentioned in the *.parts files.
            ! This helps avoiding race conditions
             IF ((comm_color < comm_counter) .AND. (skip_active)) CYCLE
            
            ! Reset the comm_counter as long the round_robin skip is active. 
            ! This ensures distributing the domains in a controlled way without race conditions.
            !
            ! This if/else MUST come *after* checking parts_list(nn) == worker_size_mpi. Otherwise, 
            ! too many "foreign" domains are computed on this communicator.
            IF (round_robin_skip <= max_skip) THEN
                round_robin_skip = round_robin_skip + 1_ik
                comm_counter = -1_ik
            ELSE
                skip_active = .FALSE.
            END IF 

            CALL MPI_FILE_READ_AT(status_un, Int((nn-1)*ik, MPI_OFFSET_KIND), dmn_status, &
                1_mik, MPI_INTEGER8, status_mpi, ierr)
        
            IF (dmn_status > -100000000_ik) CYCLE

            dmn_status = dmn_status + 100000000

            ! Write to file which comm is responsible for the domain
            CALL MPI_FILE_WRITE_AT(comms_dmn_un, Int((nn-1)*mik, MPI_OFFSET_KIND), rank_mpi, &
                1_mik, MPI_INTEGER4, status_mpi, ierr)

            ! Mark the beginning of computing the Domain
            CALL MPI_FILE_WRITE_AT(status_un, Int((nn-1)*ik, MPI_OFFSET_KIND), dmn_status, &
                1_mik, MPI_INTEGER8, status_mpi, ierr)

            IF ((dmn_status*(-1) /= Domains(nn)) .AND. (out_amount == "DEBUG")) THEN
                WRITE(std_out,FMT_ERR_AxI0) "dmn_status*(-1)", dmn_status*(-1) 
                WRITE(std_out,FMT_ERR_AxI0) "Domains(nn)", Domains(nn)          
                CALL print_err_stop(std_out, "dmn_status*(-1) /= Domains(nn)", 1)
            END IF 

            Domain = Domains(nn)

        END IF 

        !------------------------------------------------------------------------------
        ! Basically "mark end of file"
        !------------------------------------------------------------------------------
        ! Write last segments to streams. 
        ! Required if a sub-comm will not write as many numbers into the stream as 
        ! Were requested by raise_leave. Not finalizing this procedure will end up in 
        ! "end of file errors" during reading with standard procedures that do not 
        ! parse the file size.
        !------------------------------------------------------------------------------
        ! First leaf (Integer 8)
        ! Last leaf (Real 8)
        !------------------------------------------------------------------------------
        ! @domain numbers --> The last domain number and its no_nodes/no_elems
        !  computed by this rank/comm!
        ! -2 to undo the last increment and to account for the first domain "0"
        !------------------------------------------------------------------------------
        ! leaves 3 -> 3 INTEGER 8 leaves
        ! leaves 22 --> Last leaf, contains density

        !------------------------------------------------------------------------------
        ! (max_domains_per_comm*24)-1 because the last int8 entry has 24 ints.
        !------------------------------------------------------------------------------
        CALL MPI_FILE_WRITE_AT(fh_mpi_worker(4), &
            Int(root%branches(3)%leaves(4)%lbound-1+((max_domains_per_comm*24)-1), MPI_OFFSET_KIND), &
            INT(Domain, KIND=ik), 1_pd_mik, MPI_INTEGER8, status_mpi, ierr)

        CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
            Int(root%branches(3)%leaves(24)%lbound-1+(max_domains_per_comm-1), MPI_OFFSET_KIND), &
            1_rk, 1_pd_mik, MPI_REAL8, status_mpi, ierr)
        
        !------------------------------------------------------------------------------
        ! Finish the communicator if all domains are done.
        ! Otherwise, receive a proper domain number (other than -94)
        !------------------------------------------------------------------------------
        CALL MPI_BCAST(Domain,  1_mik, MPI_INTEGER8, 0_mik, WORKER_COMM, ierr)
        IF(Domain == -94_ik) GOTO 1000

        CALL MPI_BCAST(comm_nn, 1_mik, MPI_INTEGER8, 0_mik, WORKER_COMM, ierr)

        IF (out_amount == "DEBUG") THEN
            job_dir = domain_path(nn)
        ELSE
            IF(worker_rank_mpi == 0)  job_dir = outpath
        END IF
        
        CALL MPI_BCAST(job_dir, INT(4*mcl,mik), MPI_CHAR, 0_mik, WORKER_COMM, ierr)

        IF (job_dir(len(job_dir):len(job_dir)) /= "/") job_dir = trim(job_dir)//"/"


        !------------------------------------------------------------------------------
        ! Track the start time of the computation of a domain.
        !------------------------------------------------------------------------------
        CALL CPU_TIME(t_start)
        
        !==============================================================================
        ! Compute a domain
        !==============================================================================
        CALL exec_single_domain(root, comm_nn, Domain, type_raw, job_dir, &
        Active, fh_mpi_worker, worker_comm, mem_critical)
        !==============================================================================

        !------------------------------------------------------------------------------
        ! Track the duration of the computation of a domain.
        !------------------------------------------------------------------------------
        CALL CPU_TIME(t_end)
        t_duration = t_end - t_start

        !------------------------------------------------------------------------------
        ! Organize Results
        !------------------------------------------------------------------------------
        IF (worker_rank_mpi==0) THEN

            !------------------------------------------------------------------------------
            ! Write the start time to file
            !------------------------------------------------------------------------------
            ! CALL add_leaf_to_branch(result_branch, "Start Time", 1_pd_ik, [t_start])
            CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
                Int(root%branches(3)%leaves(2)%lbound-1+(comm_nn-1), MPI_OFFSET_KIND), &
                t_start, 1_pd_mik, MPI_INTEGER8, status_mpi, ierr)

            !------------------------------------------------------------------------------
            ! Write the duration to file
            !------------------------------------------------------------------------------
            ! CALL add_leaf_to_branch(result_branch, "Duration", 1_pd_ik, [t_duration])
            CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
                Int(root%branches(3)%leaves(3)%lbound-1+(comm_nn-1), MPI_OFFSET_KIND), &
                t_duration, 1_pd_mik, MPI_INTEGER8, status_mpi, ierr)

            CALL Start_Timer("Write Worker Root Branch")

            !------------------------------------------------------------------------------
            ! Read this data by:
            ! ../../../bin/pd_dump_leaf_x86_64 $PWD/ results_0000001 7
            !------------------------------------------------------------------------------
            CALL Write_Tree(root%branches(3)) ! Branch with 'Averaged Material Properties'
            
            
            CALL End_Timer("Write Worker Root Branch")

            !------------------------------------------------------------------------------
            ! Store activity information
            !------------------------------------------------------------------------------
            Call MPI_FILE_WRITE_AT(status_un, INT(((nn-1) * ik), MPI_OFFSET_KIND), &
                Domain, 1_mik, MPI_INTEGER8, status_mpi, ierr)

            !------------------------------------------------------------------------------
            ! Deallocate results of this domain. Potential memory leak.
            !------------------------------------------------------------------------------
            WRITE(desc,'(A,I0)')"Domain ", Domain
            CALL delete_branch_from_branch(TRIM(desc), root, dsize)

        END IF

        !------------------------------------------------------------------------------
        ! Mark the current position within the stream.
        !------------------------------------------------------------------------------
        CALL MPI_FILE_WRITE_AT(fh_mpi_worker(4), &
            INT(root%branches(3)%leaves(1)%lbound-1+(max_domains_per_comm-1-1), MPI_OFFSET_KIND), &
            INT(comm_nn, KIND=ik), 1_mik, MPI_INTEGER8, status_mpi, ierr)

        !------------------------------------------------------------------------------
        ! Increment linear domain counter of the PETSc sub_comm instances
        ! It is important (!) to increment this variable after the check
        ! IF (Domain_status(nn) <= 0) THEN. Otherwise, the program will write into 
        ! slots that are used by previous domains.
        !------------------------------------------------------------------------------
        comm_nn = comm_nn + 1_ik

        !------------------------------------------------------------------------------
        ! Finish the domain
        !------------------------------------------------------------------------------
        IF (worker_rank_mpi == 0)  THEN

            dmn_status = dmn_status * (-1_ik)

            CALL MPI_FILE_WRITE_AT(status_un, Int((nn-1)*ik, MPI_OFFSET_KIND), dmn_status, &
                1_mik, MPI_INTEGER8, status_mpi, ierr)

        END IF 

        !------------------------------------------------------------------------------
        ! Stop (not very gracefully) if memory is exhausted.
        !------------------------------------------------------------------------------
        IF ((worker_rank_mpi==0) .AND. (mem_critical < 0._rk))  THEN
            WRITE(std_out, FMT_WRN_xAI0) "OOM on ranks", rank_mpi, "to", rank_mpi+parts
            WRITE(std_out, FMT_WRN_xAF0) "Memory usage normalized to a compute node of hawk:", -mem_critical, "GB"

            CALL MPI_ABORT(MPI_COMM_WORLD, 0_mik, ierr)

        END IF 

    End Do

    CALL PetscFinalize(petsc_ierr) 

    CALL CPU_TIME(final_time)

    !------------------------------------------------------------------------------
    ! Tell the global main rank that this communicator is done computing stuff.
    ! Tell the worker ranks by sending Domain=-94 to stop via GOTO 1000
    !------------------------------------------------------------------------------
    CALL MPI_SEND(final_time, 1_mik, MPI_DOUBLE_PRECISION, 0_mik, rank_mpi, MPI_COMM_WORLD, ierr)
    CALL print_err_stop(std_out, "MPI_SEND of final_time didn't succeed", ierr)

    CALL MPI_BCAST(-94_ik,  1_mik, MPI_INTEGER8, 0_mik, WORKER_COMM, ierr)

ELSE ! --> rank_mpi == 0

    !------------------------------------------------------------------------------
    ! Ask all the different communicator whether they are done computing stuff.
    ! Can be done with a trivial send/recv even if it is a blocking communication.
    ! All the threads are "finished" anyway.
    !------------------------------------------------------------------------------
    jj = 1_ik
    kk = 1_ik
    rank = 1_ik
    time_wasted = 0._rk

    CALL MPI_RECV(comm_fin_time, 1_mik, MPI_DOUBLE_PRECISION, 1_mik, 1_mik, MPI_COMM_WORLD, status_mpi, ierr)
    CALL print_err_stop(std_out, "MPI_RECV on active didn't succseed", ierr)

    !------------------------------------------------------------------------------
    ! Get the amount of time wasted
    !------------------------------------------------------------------------------
    CALL CPU_TIME(now)

    time_wasted = time_wasted + ((now - comm_fin_time) * mesh_p_per_dmn)

    DO ii = 1, no_of_comms-1_ik

        mesh_p_per_dmn = comms_array(1,kk)

        ! comms_array(1,kk) --> Size of the communicator (parts per domain)
        rank = rank + mesh_p_per_dmn

        ! comms_array(2,kk) --> Numbers of corresponding communicaotrs
        IF(jj==comms_array(2,kk)) THEN
            jj = 0_ik
            kk = kk + 1_ik
        END IF 

        jj = jj + 1_ik

        CALL MPI_RECV(comm_fin_time, 1_mik, MPI_DOUBLE_PRECISION, INT(rank,mik), INT(rank,mik), MPI_COMM_WORLD, status_mpi, ierr)
        CALL print_err_stop(std_out, "MPI_RECV on active didn't succseed", ierr)

        !------------------------------------------------------------------------------
        ! Get the amount of time wasted
        !------------------------------------------------------------------------------
        CALL CPU_TIME(now)

        time_wasted = time_wasted + ((now - comm_fin_time) * mesh_p_per_dmn)

    END DO

    !------------------------------------------------------------------------------
    ! Check the comms file whether there was a domain not computed
    !------------------------------------------------------------------------------
    DO nn = 1, no_of_domains

        CALL MPI_FILE_READ_AT(comms_dmn_un, Int((nn-1)*ik, MPI_OFFSET_KIND), comms_feedback, &
            1_mik, MPI_INTEGER8, status_mpi, ierr)

        IF(comms_feedback == 0_ik) probably_failed = probably_failed + 1_ik

    END DO

    IF (probably_failed > 0_ik) THEN
        WRITE(probably_failed_char, '(I10)') probably_failed

        mssg = TRIM(ADJUSTL(probably_failed_char))//" domains probably were not computed."
        CALL print_trimmed(std_out, TRIM(mssg), FMT_WRN)
    END IF 

END IF ! IF (rank_mpi /= 0) THEN

! May help to gracefully stop the program in case of an error.
1000 CONTINUE

CALL MPI_FILE_CLOSE(status_un, ierr)
CALL MPI_FILE_CLOSE(comms_dmn_un, ierr)

CALL MPI_FINALIZE(ierr)
CALL print_err_stop(std_out, "MPI_FINALIZE didn't succeed", ierr)

IF(rank_mpi==0) THEN

    !------------------------------------------------------------------------------
    ! Write the "JOB_FINISHED" keyword only if the job was finished during 
    ! this job (re)run.
    !------------------------------------------------------------------------------
    IF ((Domain_status(No_of_domains) == No_of_domains-1) .OR. (already_finished)) THEN ! counts from 0

        no_restart_required = .TRUE.
    END IF 

    CALL CPU_TIME(final_time)

    time_wasted = (time_wasted + (final_time - now) * size_mpi) / 3600._rk

    CALL meta_write('CORE_HOURS_WASTED', "(h)", time_wasted)

    CALL meta_close(m_rry)

    If (out_amount == "DEBUG") THEN
        Write(fh_log, fmt_dbg_sep)
        Write(fh_log, fmt_MSG_xAI0) "Final Root pointer proc", rank_mpi
        CALL log_tree(root, fh_log, .True.)
        Write(fh_log, fmt_dbg_sep)
    END If

    CALL link_end(link_name,.True.)
   
    IF (mem_critical > 0._rk)THEN
        WRITE(std_out, FMT_TXT_SEP)
        WRITE(std_out, FMT_TXT) "HLRS Direct Tensor Computation finished successfully."
        WRITE(std_out, FMT_TXT_SEP)

        CALL meta_close(m_rry)

    ELSE
        WRITE(std_out, FMT_WRN) "Stopping due to an OOM event."
    END IF 

    IF(std_err/=6) CALL meta_stop_ascii(std_out, '.std_out')
    IF(std_err/=0) CALL meta_start_ascii(std_err, '.std_err')

END IF ! (rank_mpi == 0)

End Program dtc
