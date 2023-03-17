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

INTEGER(mik) :: ierr, rank_mpi, size_mpi, petsc_ierr, mii, mjj, &
    worker_rank_mpi, worker_size_mpi, aun, par_domains, &
    Active, request, finished = -1, worker_comm

INTEGER(mik), Dimension(no_streams)       :: fh_mpi_worker
INTEGER(mik), Dimension(MPI_STATUS_SIZE)  :: status_mpi
INTEGER(mik), Dimension(:,:), Allocatable :: statuses_mpi
INTEGER(mik), Dimension(:)  , Allocatable :: worker_is_active, req_list

INTEGER(c_int) :: stat_c_int

Type(tBranch)          :: root, phi_tree
Type(tBranch), pointer :: ddc, meta_para, result_branch

CHARACTER, DIMENSION(4*mcl)          :: c_char_array
CHARACTER, DIMENSION(:), ALLOCATABLE :: char_arr
CHARACTER(4*mcl), DIMENSION(:), ALLOCATABLE :: domain_path
CHARACTER(mcl)  , DIMENSION(:), ALLOCATABLE :: m_rry      

CHARACTER(4*mcl) :: job_dir
CHARACTER(mcl)   :: cmd_arg_history='', link_name = 'struct process', &
    muCT_pd_path, muCT_pd_name, binary, status_file, desc="", memlog_file="", &
    typeraw="", restart='N', restart_cmd_arg='U',ios="", map_sgmnttn=""
CHARACTER(8)   :: elt_micro, output
CHARACTER(3)   :: file_status
CHARACTER(1) :: bin_sgmnttn=""

REAL(rk) :: strain, t_start, t_end, t_duration

INTEGER(ik), DIMENSION(:), ALLOCATABLE :: Domains, nn_D, Domain_status, comms_list, parts_list
INTEGER(ik), DIMENSION(3) :: xa_d=0, xe_d=0

INTEGER(ik) :: nn, ii, jj, kk, dc, computed_domains = 0, comm_nn = 1, &
    No_of_domains, No_of_domains_files, path_count, activity_size=0, No_of_comms, &
    alloc_stat, fh_cluster_log, free_file_handle, domains_per_comm, stat, &
    no_lc=0, nl=0, Domain, llimit, parts, elo_macro, vdim(3), comms_size, parts_size

INTEGER(pd_ik), DIMENSION(:), ALLOCATABLE :: serial_root
INTEGER(pd_ik), DIMENSION(no_streams) :: dsize

INTEGER(pd_ik) :: serial_root_size, add_leaves

LOGICAL :: success, status_exists, comms_exists, parts_exists, &
    heaxist, abrt = .FALSE., already_finished=.FALSE., &
    create_new_header = .FALSE., fex=.TRUE., no_restart_required = .FALSE.

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
    ! Batch feedback
    !------------------------------------------------------------------------------
    CALL execute_command_line ("echo 'JOB_STARTED' > BATCH_RUN")

    !------------------------------------------------------------------------------
    ! Parse the command arguments
    !------------------------------------------------------------------------------
    CALL get_cmd_args(binary, in%full, restart_cmd_arg, cmd_arg_history)
    
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
    project_name = "results" ! TRIM(out%bsnm)//

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
    CALL meta_read('MESH_PER_SUB_DMN' , m_rry, parts, ios); CALL mest(ios, abrt)
    CALL meta_read('RVE_STRAIN'       , m_rry, strain, ios); CALL mest(ios, abrt)
    CALL meta_read('YOUNG_MODULUS'    , m_rry, bone%E, ios); CALL mest(ios, abrt)
    CALL meta_read('POISSON_RATIO'    , m_rry, bone%nu, ios); CALL mest(ios, abrt)
    CALL meta_read('MACRO_ELMNT_ORDER', m_rry, elo_macro, ios); CALL mest(ios, abrt)
    CALL meta_read('TYPE_RAW'         , m_rry, typeraw, ios); CALL mest(ios, abrt)
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

    !------------------------------------------------------------------------------
    ! Spawn a results file
    !------------------------------------------------------------------------------
    ! This log file may collide with the original log file (!)
    ! The regular struct_process log file contains still has the "old" basename!
    !------------------------------------------------------------------------------
    CALL meta_start_ascii(fh_mon, mon_suf)

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
            &file. Only the round-robin struct_process can create a stauts file &
            &on its own."
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF 

    comms_file = TRIM(out%p_n_bsnm)//".comms"
    INQUIRE(FILE = TRIM(comms_file), EXIST=comms_exists, SIZE=comms_size)
    IF(.NOT. comms_exists) CALL print_err_stop_slaves("No *.comms file found."); GOTO 1000

    parts_file = TRIM(out%p_n_bsnm)//".parts"
    INQUIRE(FILE = TRIM(parts_file), EXIST=parts_exists, SIZE=parts_size)
    IF(.NOT. parts_exists) CALL print_err_stop_slaves("No *.parts file found."); GOTO 1000


    !------------------------------------------------------------------------------
    ! Check the condition of the required files.
    !------------------------------------------------------------------------------
    IF(activity_size .NE. parts_size) THEN
        mssg = 'The size of the *.parts and the *.status files do not match.'
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF 

    No_of_domains = (xe_d(1)-xa_d(1)+1) * (xe_d(2)-xa_d(2)+1) * (xe_d(3)-xa_d(3)+1)
    No_of_domains_files = ANINT(REAL(activity_size, rk) / 8._rk)

    IF(No_of_domains .NE. No_of_domains_files) THEN
        mssg = 'The size of the auxiliary files and the requestd domain range do not match.'
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF 

    No_of_comms = ANINT(REAL(comms_size, rk) / 8._rk)


    ! IF ((No_of_domains*ik /= activity_size ) .AND. (status_exists)) THEN
    !     mssg = 'Size of the status file and the number of requested domains do &
    !     &not match. It is very likely that the domain ranges were modified and a &
    !     &restart requested. This approach is not supported, as the stream sizes and &
    !     &boundaries will get corrupted. Please set up a completely new job.'
    !     CALL print_err_stop_slaves(mssg, "warning"); GOTO 1000
    ! END IF


    !------------------------------------------------------------------------------
    ! Raise and build meta_para tree
    ! Hardcoded, implicitly given order of the leafs. 
    ! DO NOT CHANGE ORDER WITHOUT MODIFYING ALL OTHER INDICES REGARDING »meta_para«
    !------------------------------------------------------------------------------
    Allocate(meta_para)
    CALL raise_tree("Input parameters", meta_para)

    CALL add_leaf_to_branch(meta_para, "muCT puredat pro_path"                , mcl , str_to_char(muCT_pd_path)) 
    CALL add_leaf_to_branch(meta_para, "muCT puredat pro_name"                , mcl , str_to_char(muCT_pd_name)) 
    CALL add_leaf_to_branch(meta_para, "Physical domain size"                 , 3_ik, bone%phdsize) 
    CALL add_leaf_to_branch(meta_para, "Lower bounds of selected domain range", 3_ik, xa_d) 
    CALL add_leaf_to_branch(meta_para, "Upper bounds of selected domain range", 3_ik, xe_d)      
    
    CALL add_leaf_to_branch(meta_para, "Grid spacings"                 , 3_ik, bone%delta) 
    CALL add_leaf_to_branch(meta_para, "Lower limit of iso value"      , 1_ik, [llimit])      
    CALL add_leaf_to_branch(meta_para, "Element type  on micro scale"  , len(elt_micro) , str_to_char(elt_micro))      
    CALL add_leaf_to_branch(meta_para, "No of mesh parts per subdomain", 1_ik           , [parts]) 
    CALL add_leaf_to_branch(meta_para, "Output Format"                 , len(output)    , str_to_char(output)) 
    
    CALL add_leaf_to_branch(meta_para, "Average strain on RVE"         , 1_ik, [strain])    
    CALL add_leaf_to_branch(meta_para, "Young_s modulus"               , 1_ik, [bone%E]) 
    CALL add_leaf_to_branch(meta_para, "Poisson_s ratio"               , 1_ik, [bone%nu]) 
    CALL add_leaf_to_branch(meta_para, "Element order on macro scale"  , 1_ik, [elo_macro]) 
    CALL add_leaf_to_branch(meta_para, "Output amount"                 , len(out_amount), str_to_char(out_amount)) 
    
    CALL add_leaf_to_branch(meta_para, "Restart", 1_ik, str_to_char(restart(1:1))) 
    CALL add_leaf_to_branch(meta_para, "Number of voxels per direction", 3_ik , vdim) 
    CALL add_leaf_to_branch(meta_para, "Domains per communicator", 1_ik, [domains_per_comm]) 

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
        WRITE(fh_mon, FMT_MSG) "Reusing the output directory"
        WRITE(fh_mon, FMT_MSG) TRIM(outpath)
    END IF


    CALL link_start(link_name, .TRUE., .FALSE., success)
    IF (.NOT. success) THEN
        mssg = "Something went wrong during link_start"
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF

END IF 

CALL MPI_BCAST(status_file, INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(comms_file,  INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(parts_file,  INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)


CALL MPI_BCAST(No_of_domains, 1_mik, MPI_INTEGER8   , 0_mik, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(No_of_comms,   1_mik, MPI_INTEGER8   , 0_mik, MPI_COMM_WORLD, ierr)

CALL MPI_BCAST(typeraw, INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)

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

Call MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(status_file), MPI_MODE_RDWR, MPI_INFO_NULL, aun, ierr)
CALL MPI_FILE_READ(aun, Domain_status, INT(SIZE(Domain_status), mik), MPI_INTEGER8, MPI_STATUS_IGNORE, ierr)


!------------------------------------------------------------------------------
! Read the communicators for use with the dataset
!------------------------------------------------------------------------------
ALLOCATE(comms_list(No_of_comms), stat=alloc_stat)
CALL alloc_err("comms_list", alloc_stat)

Call MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(comms_file), MPI_MODE_RDWR, MPI_INFO_NULL, aun, ierr)
CALL MPI_FILE_READ(aun, comms_list, INT(SIZE(comms_list), mik), MPI_INTEGER8, MPI_STATUS_IGNORE, ierr)


!------------------------------------------------------------------------------
! Read the number of parts assigned to each domain (ppd). The ppd must fit
! to the communicators read from the *.comms file.
!------------------------------------------------------------------------------
ALLOCATE(parts_list(No_of_domains), stat=alloc_stat)
CALL alloc_err("parts_list", alloc_stat)

Call MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(parts_file), MPI_MODE_RDWR, MPI_INFO_NULL, aun, ierr)
CALL MPI_FILE_READ(aun, parts_list, INT(SIZE(parts_list), mik), MPI_INTEGER8, MPI_STATUS_IGNORE, ierr)



ALLOCATE(domain_path(0:No_of_domains))
domain_path = ''




If (rank_mpi == 0) THEN

    IF (status_exists) THEN
        !------------------------------------------------------------------------------
        ! Read the Domain stats list
        !------------------------------------------------------------------------------
        CALL MPI_FILE_READ(aun, Domain_status, INT(SIZE(Domain_status), mik), &
            MPI_INTEGER8, MPI_STATUS_IGNORE, ierr)

        !------------------------------------------------------------------------------
        ! Check whether computation will use resources properly.
        !------------------------------------------------------------------------------
        DO ii=1, SIZE(Domain_status)
            IF (Domain_status(ii) >= 0) computed_domains = computed_domains + 1 
        END DO

        IF (No_of_domains == computed_domains) THEN 
            already_finished=.TRUE.
            mssg = "Job is already finished. No restart required."

            CALL print_err_stop_slaves(mssg, "message"); GOTO 1000
        END IF 

    ELSE
        CALL MPI_FILE_WRITE(aun, Domain_status, INT(SIZE(Domain_status), mik), &
            MPI_INTEGER8, MPI_STATUS_IGNORE, ierr)
    END IF

    !------------------------------------------------------------------------------
    ! Check if INPUT header exists. If not --> create one.
    !------------------------------------------------------------------------------
    INQUIRE(FILE = TRIM( in%p_n_bsnm)//".head", EXIST=heaxist)
    IF(.NOT. heaxist) THEN
        free_file_handle = give_new_unit()
        CALL convert_meta_to_puredat(free_file_handle, m_rry)
    END IF

    !------------------------------------------------------------------------------
    ! Check if OUTPUT header exists. If not --> Internally reset to restart='N'.
    !------------------------------------------------------------------------------
    INQUIRE(FILE = TRIM(pro_path)//TRIM(pro_name)//".head", EXIST=heaxist)


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
            WRITE(un_lf, fmt_dbg_sep)
            Write(un_lf,'(A)') "root right after restart read"
            CALL log_tree(root, un_lf,.FALSE.)
            WRITE(un_lf, fmt_dbg_sep)
            flush(un_lf)
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
        meta_para%leaves(15)%p_char = str_to_char(out_amount)
        meta_para%leaves(16)%p_char = "Y"

        If (out_amount == "DEBUG") THEN 
            Write(un_lf,fmt_dbg_sep)
            Write(un_lf,'(A)') "root right after restart and deletion of avg mat props branch"
            CALL log_tree(root, un_lf,.FALSE.)
            Write(un_lf, fmt_dbg_sep)
            flush(un_lf)
        END If

    END IF ! restart == Yes/No

    !------------------------------------------------------------------------------
    ! Check whether computation will use resources properly.
    !------------------------------------------------------------------------------
    ! NOT NECESSARY? PLEASE OBSERVE...
    !------------------------------------------------------------------------------
    ! IF ( MOD((No_of_domains - computed_domains) * parts, size_mpi - 1) /= 0 ) THEN
    !     write(std_out, FMT_DBG_AI0xAI0) "No_of_domains    = ", No_of_domains
    !     write(std_out, FMT_DBG_AI0xAI0) "computed_domains = ", computed_domains
    !     write(std_out, FMT_DBG_AI0xAI0) "parts            = ", parts
    !     write(std_out, FMT_DBG_AI0xAI0) "size_mpi - 1     = ", size_mpi - 1
    !     mssg = "MOD((No_of_domains - computed_domains) * parts, size_mpi - 1) /= 0"
    !     CALL print_err_stop_slaves(mssg, "warning"); GOTO 1000
    ! END IF
    
    IF ( parts > size_mpi - 1 ) THEN
        mssg = "parts > size_mpi - 1"
        CALL print_err_stop_slaves(mssg, "warning"); GOTO 1000
    END IF

    !------------------------------------------------------------------------------
    ! Ensure to update the number of leaves and the indices of these in every
    ! line of code! Also update dat_ty and dat_no in "CALL raise_leaves"
    !------------------------------------------------------------------------------
    add_leaves = 24_pd_ik

    !------------------------------------------------------------------------------
    ! Number of loadcases
    !------------------------------------------------------------------------------
    ! no_lc = 24 
    
    !------------------------------------------------------------------------------
    ! The name of the branche is a bit outdated. In the meantime, it contains 
    ! more data to accomodate the need for the doctoral project.
    ! It was done this way because it is the quickest and easiest way for 
    ! combining domain specfic output data.
    !------------------------------------------------------------------------------
    CALL add_branch_to_branch(root, result_branch)
    CALL raise_branch("Averaged Material Properties", 0_pd_ik, add_leaves, result_branch)

    !------------------------------------------------------------------------------
    ! Former *.memlog file
    !------------------------------------------------------------------------------
    ! Operation, Domain, Nodes, Elems, Preallo, Mem_comm, Pids_returned, Size_mpi, time
    ! Before PETSc preallocation              , 3395     , 1672328, 1442071, 546, 37198100, 252, 252, 1671690714
    ! After PETSc preallocation               , 3395     , 1672328, 1442071, 546, 209263436, 252, 252, 1671690714
    ! Before matrix assembly                  , 3395     , 1672328, 1442071, 546, 269684188, 252, 252, 1671690714
    ! After matrix assembly                   , 3395     , 1672328, 1442071, 546, 242114948, 252, 252, 1671690715
    ! Before solving                          , 3395     , 1672328, 1442071, 546, 242514604, 252, 252, 1671690716
    ! After solving                           , 3395     , 1672328, 1442071, 546, 244625464, 252, 252, 1671690736
    ! End of domain                           , 3395     , 1672328, 1442071, 546, 73615236, 252, 252, 1671690747
    !------------------------------------------------------------------------------

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
        dat_no = [ domains_per_comm, &
        domains_per_comm,         domains_per_comm, &
        domains_per_comm * 24   ,                                                    &
        domains_per_comm,         domains_per_comm, &                                  
        domains_per_comm * nl*nl, domains_per_comm * nl*nl, domains_per_comm       , &
        domains_per_comm *  6*nl, domains_per_comm *  6*nl, domains_per_comm *  6*6, &
        domains_per_comm        , domains_per_comm *  6* 6, domains_per_comm       , &
        domains_per_comm        , domains_per_comm *     3, domains_per_comm *    9, &
        domains_per_comm *  6* 6, domains_per_comm        , domains_per_comm *    3, &
        domains_per_comm *     9, domains_per_comm *  6* 6, domains_per_comm      ], &
        branch = result_branch)

    result_branch%leaves(:)%pstat = -1
    
    CALL set_bounds_in_branch(root, root%streams)
    
    If (out_amount == "DEBUG") THEN 
        Write(un_lf, fmt_dbg_sep)
        Write(un_lf, '(A)') "root right before serialisation"
        CALL log_tree(root,un_lf,.True.)
        Write(un_lf, fmt_dbg_sep)
        flush(un_lf)
    END If

    !------------------------------------------------------------------------------
    ! serialize root branch
    !------------------------------------------------------------------------------
    CALL serialize_branch(root, serial_root, serial_root_size, .TRUE.)

    !------------------------------------------------------------------------------
    ! Get number of domains on each axis
    !------------------------------------------------------------------------------
    CALL pd_get(ddc,"nn_D", nn_D)

    !------------------------------------------------------------------------------
    ! Init Domain Cross Reference and domain paths
    !------------------------------------------------------------------------------
    dc = 0
    domain_path(0) = Trim(pro_path)//Trim(pro_name)//"_domains" ! "_results"

    nn = 1
    Do kk = xa_d(3), xe_d(3)
    Do jj = xa_d(2), xe_d(2)
    Do ii = xa_d(1), xe_d(1)
        dc         = dc + 1_mik
        path_count = dc / 2_mik

        Write(domain_path(dc),'(A,"/",I0)')Trim(domain_path(path_count)),dc

        Domains(nn) = ii + jj * nn_D(1) + kk * nn_D(1)*nn_D(2)

        nn = nn + 1_ik
    End Do
    End Do
    End Do

    !------------------------------------------------------------------------------
    ! Generate worker_is_active_List
    !------------------------------------------------------------------------------
    Allocate(worker_is_active(size_mpi-1), stat=alloc_stat)
    CALL alloc_err("worker_is_active_List", alloc_stat)

    worker_is_active=0

    CALL End_Timer("Init Process")

    !------------------------------------------------------------------------------
    ! Start Workers
    !------------------------------------------------------------------------------
    CALL Start_Timer("Broadcast Init meta_para")

    CALL mpi_bcast(outpath,       INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(project_name , INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
    
    CALL mpi_bcast(serial_root_size, 1_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(serial_root , INT(serial_root_size, mik), MPI_INTEGER8, 0_mik,&
        MPI_COMM_WORLD, ierr)

    deallocate(serial_root)
    CALL End_Timer("Broadcast Init meta_para")

    !------------------------------------------------------------------------------
    ! Execute collective mpi_comm_split. Since mpi_comm_world rank 0 is
    ! the head master worker_comm is not needed and it should not be in
    ! any worker group and communicator. With MPI_UNDEFINED passed as
    ! color worker_comm gets the value MPI_COMM_NULL
    !------------------------------------------------------------------------------
    CALL MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, rank_mpi, worker_comm, ierr)
    IF (ierr /=0) THEN
        mssg = "MPI_COMM_SPLIT couldn't split MPI_COMM_WORLD"
        CALL print_err_stop_slaves(mssg); GOTO 1000
    END IF 


!------------------------------------------------------------------------------
! Ranks > 0 -- Worker slaves
!------------------------------------------------------------------------------
Else
    CALL Start_Timer("Broadcast Init meta_para")
     
    !------------------------------------------------------------------------------
    ! Broadcast recieve init parameters
    !------------------------------------------------------------------------------
    CALL MPI_BCAST(outpath      , INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(project_name , INT(mcl,mik), MPI_CHAR, 0_mik, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(serial_root_size, 1_mik, MPI_INTEGER8, 0_mik, MPI_COMM_WORLD, ierr)
 
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
    
    CALL mpi_bcast(serial_root, INT(serial_root_size, mik), MPI_INTEGER8, 0_mik, &
        MPI_COMM_WORLD, ierr)

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
    
    !------------------------------------------------------------------------------
    ! Linear vs ddc (3D) domain numbers. nn is an internal variable.
    !------------------------------------------------------------------------------
    nn = 1_ik

    Do kk = xa_d(3), xe_d(3)
    Do jj = xa_d(2), xe_d(2)
    Do ii = xa_d(1), xe_d(1)
        Domains(nn) = ii + jj * nn_D(1) + kk * nn_D(1)*nn_D(2)
        nn = nn + 1_ik
    End Do
    End Do
    End Do

    CALL pd_get(root%branches(1), "No of mesh parts per subdomain", parts)
    CALL pd_get(root%branches(1), "Domains per communicator", domains_per_comm)

    !------------------------------------------------------------------------------
    ! All Worker Ranks -- Init worker Communicators
    !------------------------------------------------------------------------------
    CALL MPI_Comm_split(MPI_COMM_WORLD, Int((rank_mpi-1)/parts, mik), rank_mpi, worker_comm, ierr)
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

! WRITE_ROOT_COMM = MPI_COMM_WORLD

!------------------------------------------------------------------------------
! Open stream files
!------------------------------------------------------------------------------
! Call Open_Stream_Files(root%streams, "write", "new", fh_mpi_root, WRITE_ROOT_COMM)

!------------------------------------------------------------------------------
! All Ranks -- Init MPI request and status lists
!------------------------------------------------------------------------------
Allocate(req_list(size_mpi-1), stat=alloc_stat)
CALL alloc_err("req_list", alloc_stat)
req_list=0

Allocate(statuses_mpi(MPI_STATUS_SIZE, size_mpi-1), stat=alloc_stat)
CALL alloc_err("statuses_mpi", alloc_stat)

!------------------------------------------------------------------------------
! Global Rank 0 -- Process master Start working process
!------------------------------------------------------------------------------
If (rank_mpi==0) Then

    !------------------------------------------------------------------------------
    ! Only global rank 0 -- Init MPI IO-System
    !------------------------------------------------------------------------------
    CALL get_stream_size(root, dsize)

    If (out_amount == "DEBUG") &
        write(fh_mon,FMT_MSG_AxI0) "On rank zero, stream sizes: ", dsize
    
    !------------------------------------------------------------------------------
    ! Supply all worker masters  with their first work package
    ! ii is incremented by ii = ii + parts_per_subdomain
    !------------------------------------------------------------------------------
    nn = 1_ik; mii = 1_mik
    Do While (mii <= (size_mpi-1_mik))

        if (nn > No_of_domains) exit

        Do mjj = mii, mii + parts-1
            !------------------------------------------------------------------------------
            ! Start worker
            !------------------------------------------------------------------------------
            CALL mpi_send(worker_is_active(mjj), 1_mik, MPI_INTEGER4, mjj, mjj, MPI_COMM_WORLD, ierr)
            CALL print_err_stop(std_out, "MPI_SEND of worker_is_active didn't succeed", ierr)

            !------------------------------------------------------------------------------
            ! Send Domain number
            !------------------------------------------------------------------------------
            CALL mpi_send(nn, 1_mik, MPI_INTEGER8, mjj, mjj, MPI_COMM_WORLD,ierr)
            CALL print_err_stop(std_out, "MPI_SEND of Domain number didn't succeed", ierr)
            
            if (out_amount /= "PRODUCTION") then
                CALL mpi_send(domain_path(nn), Int(4*mcl,mik), MPI_CHARACTER, mjj, mjj, MPI_COMM_WORLD, ierr)
                CALL print_err_stop(std_out, "MPI_SEND of Domain path didn't succeed", ierr)
            End if
        End Do

        !------------------------------------------------------------------------------
        ! Log to Monitor file
        !------------------------------------------------------------------------------
        WRITE(fh_mon, FMT_MSG_xAI0) "Domain ", Domains(nn), " at Ranks ", mii, " to ", mii + parts-1
        flush(fh_mon)

        nn = nn + 1_ik
        
        Call MPI_IRECV(worker_is_active(mii), 1_mik, MPI_INTEGER4, mii, mii, MPI_COMM_WORLD, req_list(mii), ierr)
        CALL print_err_stop(std_out, "MPI_IRECV of worker_is_active(mii) didn't succeed", INT(ierr, ik))

        mii = mii + Int(parts,mik)

    End Do

    Call MPI_WAITANY(size_mpi-1_mik, req_list, finished, status_mpi, ierr)
    CALL print_err_stop(std_out, &
        "MPI_WAITANY on req_list for IRECV of Activity(mii) didn't succeed", INT(ierr, ik))

    IF(finished /= MPI_UNDEFINED) mii = finished

    !------------------------------------------------------------------------------
    ! mii is incremented by mii = mii + parts
    ! mii --> Worker master
    !------------------------------------------------------------------------------
    nn = 1_ik
    DO  WHILE (nn <= No_of_domains) 

        !------------------------------------------------------------------------------
        ! Send information for all parts to the corresponding ranks.
        ! From global "first worker rank of specific domain" to "global last worker 
        ! rank of specific domain"
        !------------------------------------------------------------------------------
        ! jj --> Worker
        !------------------------------------------------------------------------------
        Do mjj = mii, mii + parts-1

            worker_is_active(mjj) = 1_mik
            !------------------------------------------------------------------------------
            ! Start worker
            !------------------------------------------------------------------------------
            CALL mpi_send(worker_is_active(mjj), 1_mik, MPI_INTEGER4, mjj, mjj, MPI_COMM_WORLD, ierr)
            CALL print_err_stop(std_out, "MPI_SEND of worker_is_active didn't succeed", ierr)

            !------------------------------------------------------------------------------
            ! Send Domain number
            !------------------------------------------------------------------------------
            CALL mpi_send(nn, 1_mik, MPI_INTEGER8, mjj, mjj, MPI_COMM_WORLD,ierr)
            CALL print_err_stop(std_out, "MPI_SEND of Domain number didn't succeed", ierr)
            
            if (out_amount /= "PRODUCTION") then
                CALL mpi_send(domain_path(nn), Int(4*mcl,mik), MPI_CHARACTER, mjj, mjj, MPI_COMM_WORLD,ierr)
                CALL print_err_stop(std_out, "MPI_SEND of Domain path didn't succeed", ierr)
            End if
        End Do
        
        !------------------------------------------------------------------------------
        ! Log to Monitor file
        !------------------------------------------------------------------------------
        WRITE(fh_mon, FMT_MSG_xAI0) "Domain ", Domains(nn), " at Ranks ", mii, " to ", mii + parts-1
        flush(fh_mon)

        !------------------------------------------------------------------------------
        ! Iterate over domain
        !------------------------------------------------------------------------------
        nn = nn + 1_mik

        !------------------------------------------------------------------------------
        ! Worker has finished
        !------------------------------------------------------------------------------
        CALL MPI_IRECV(worker_is_active(mii), 1_mik, MPI_INTEGER4, mii, mii, MPI_COMM_WORLD, req_list(mii), ierr)
        CALL print_err_stop(std_out, "MPI_IRECV of worker_is_active(mii) didn't succeed", ierr)

        !------------------------------------------------------------------------------
        ! Only the master worker is allowed to send a 'finished' flag. 
        ! Otherwise, index mii will screw up.
        !------------------------------------------------------------------------------
        Call MPI_WAITANY(size_mpi-1_mik, req_list, finished, status_mpi, ierr)
        CALL print_err_stop(std_out, &
        "MPI_WAITANY on req_list for IRECV of worker_is_active(mii) didn't succeed", ierr)

        IF(finished /= MPI_UNDEFINED) mii = finished

    End Do
    
    !------------------------------------------------------------------------------
    ! store_parallel_branch(root, fh_mpi_root) will not print anything?!
    !------------------------------------------------------------------------------
    ! Storing/writing the full root tree results in unnecessarily written data.
    ! Storing/writing a mix of root/branches(1..2) will result in inconsistens
    ! header/stream files. The binaries to dump data after production runs will 
    ! subsequently crash. With a proper implementation, the bondaries of the 
    ! streaming files are recognized. 
    !------------------------------------------------------------------------------
    ! OPTIONAL. 
    ! Turned off at the moment, because the data is stored via the meta file 
    ! format, which is much more accesible by human beings.
    !------------------------------------------------------------------------------
    ! CALL Start_Timer("Write Root Branch")

    ! CALL store_branch(root%branches(1), root%streams, .TRUE. )
    ! CALL store_branch(root%branches(2), root%streams, .TRUE. )
    ! CALL Write_Tree(root%branches(1))
    ! CALL Write_Tree(root%branches(2))

    ! CALL End_Timer("Write Root Branch")
    !------------------------------------------------------------------------------
    ! OPTIONAL. 
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Wait for all workers to return
    !------------------------------------------------------------------------------
    CALL MPI_WAITALL(size_mpi-1_mik, req_list, statuses_mpi, ierr)
    CALL print_err_stop(std_out, &
        "MPI_WAITALL on req_list for IRECV of worker_is_active(ii) didn't succeed", ierr)

    !------------------------------------------------------------------------------
    ! Stop workers
    !------------------------------------------------------------------------------
    worker_is_active = -1_mik
    
    Do mii = 1_mik, size_mpi-1_mik
        CALL mpi_send(worker_is_active(mii), 1_mik, MPI_INTEGER4, mii, mii, MPI_COMM_WORLD, ierr)
        CALL print_err_stop(std_out, "MPI_SEND of worker_is_active didn't succeed", ierr)
    End Do

!------------------------------------------------------------------------------
! Global ranks > 0 -- Workers
!------------------------------------------------------------------------------
Else

    IF (worker_rank_mpi == 0) THEN
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

        !------------------------------------------------------------------------------
        ! Start the memory logging
        !------------------------------------------------------------------------------
                ! INQUIRE(file="./datasets/memlog.sh", exist=fex)

                ! IF(fex) THEN
                !     CALL EXECUTE_COMMAND_LINE (&
                !         './datasets/memlog.sh '&
                !         //TRIM(outpath)//TRIM(project_name)//'.memlog', CMDSTAT=stat)   

                !     IF(stat /= 0) WRITE(std_err, FMT_WRN_xAI0) &
                !         "Could not start memory logging! Rank: ", rank_mpi
                ! ELSE
                !     WRITE(std_err, FMT_WRN_xAI0) &
                !         "File for memory logging not found! Rank: ", rank_mpi
                ! END IF


        memlog_file=TRIM(outpath)//TRIM(project_name)//'.memlog'
        INQUIRE (FILE = memlog_file, EXIST = fex)

        file_status="NEW"
        IF(fex) file_status="OLD"

        fh_cluster_log = give_new_unit()

        OPEN(UNIT=fh_cluster_log, FILE=TRIM(memlog_file), ACTION='WRITE', &
            ACCESS="SEQUENTIAL", STATUS=file_status)

        WRITE(fh_cluster_log, '(A)') &
            "Operation, Domain, Nodes, Elems, Preallo, "//&
            "Mem_comm, Pids_returned, Size_mpi, time"
    
        ! IF(stat /= 0) WRITE(std_err, FMT_WRN_xAI0) &
        !     "Could not start memory logging! Rank: ", rank_mpi


        CALL link_start(link_name, .True., .True.)

    ELSE ! Worker threads of the sub communicator
        !------------------------------------------------------------------------------
        ! Set fh_cluster_log = -1 for clarifying the state of the worker threads
        !------------------------------------------------------------------------------
        fh_cluster_log = -1_ik
    END IF 

        
    !------------------------------------------------------------------------------
    ! Send new directory/filename
    ! It is the Rank_ and not the root directory! Do not delete this sequence!
    ! Root dir was send via broadcast, because the worker masters are not 
    ! determined at the beginning of the program.
    !------------------------------------------------------------------------------
    CALL MPI_BCAST(outpath,      INT(mcl, mik), MPI_CHAR, 0_mik, WORKER_COMM, ierr)
    CALL MPI_BCAST(project_name, INT(mcl, mik), MPI_CHAR, 0_mik, WORKER_COMM, ierr)
    CALL MPI_BCAST(project_name, INT(mcl, mik), MPI_CHAR, 0_mik, WORKER_COMM, ierr)

    !------------------------------------------------------------------------------
    ! Open Stream files to write data to master-worker rank directories
    !------------------------------------------------------------------------------
    pro_path = outpath
    pro_name = project_name     

    CALL set_stream_filenames(root%streams)
    CALL Open_Stream_Files(root%streams, "write", "new", fh_mpi_worker, WORKER_COMM)

    IF ((worker_rank_mpi == 0) .AND. (restart == 'Y')) THEN
        !------------------------------------------------------------------------------
        ! Read and check the Local domain number list in case of a restart
        ! restart =  'Y' is given by user and may be given without any domain computed
        ! so far. Therefore, only increment comm_nn if this entry no 0
        !------------------------------------------------------------------------------
        CALL MPI_FILE_READ_AT(fh_mpi_worker(4), & 
            Int(root%branches(3)%leaves(1)%lbound-1+(domains_per_comm-1-1), MPI_OFFSET_KIND), &
            comm_nn, 1_mik, & 
            MPI_INTEGER8, MPI_STATUS_IGNORE, ierr)
        
        !!!!!------------------------------------------------------------------------------
        !!!!! Experimental
        !!!!! May help computing domains which may be skipped
        !!!!! Some domains are skipped while restarting. Considered a minor issue
        !!!!!------------------------------------------------------------------------------
        !!!!!!!!!!!!! IF(comm_nn > 1) comm_nn = comm_nn - 1_ik

    END IF 

    CALL MPI_BCAST(comm_nn, 1_mik , MPI_INTEGER8, 0_mik, WORKER_COMM, ierr)

    !------------------------------------------------------------------------------
    ! Worker Loop
    !------------------------------------------------------------------------------
    Do
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
        ! (domains_per_comm*24)-1 because the last int8 entry has 24 ints.
        !------------------------------------------------------------------------------
        CALL MPI_FILE_WRITE_AT(fh_mpi_worker(4), &
            Int(root%branches(3)%leaves(4)%lbound-1+((domains_per_comm*24)-1), MPI_OFFSET_KIND), &
            INT(Domain, KIND=ik), 1_pd_mik, MPI_INTEGER8, status_mpi, ierr)

        CALL MPI_FILE_WRITE_AT(fh_mpi_worker(5), &
            Int(root%branches(3)%leaves(24)%lbound-1+(domains_per_comm-1), MPI_OFFSET_KIND), &
            1_rk, 1_pd_mik, MPI_REAL8, status_mpi, ierr)


        !------------------------------------------------------------------------------
        ! Start workers
        !------------------------------------------------------------------------------
        CALL MPI_RECV(active, 1_mik, MPI_INTEGER4, 0_mik, rank_mpi, MPI_COMM_WORLD, status_mpi, ierr)
        CALL print_err_stop(std_out, "MPI_RECV on Active didn't succseed", ierr)

        !------------------------------------------------------------------------------
        ! Stop workers
        !------------------------------------------------------------------------------
        IF (active == -1) EXIT

        !------------------------------------------------------------------------------
        ! Receive domain numbers
        !------------------------------------------------------------------------------
        CALL mpi_recv(nn, 1_mik, MPI_INTEGER8, 0_mik, rank_mpi, MPI_COMM_WORLD, status_mpi, ierr)
        CALL print_err_stop(std_out, "MPI_RECV on Domain didn't succeed", ierr)

        Domain = Domains(nn)

        !------------------------------------------------------------------------------
        ! Only compute the domain if it was not yet (entry in status file = -999)
        !------------------------------------------------------------------------------
        IF (Domain_status(nn) < 0) THEN

            !------------------------------------------------------------------------------
            ! Outpath depends on whether basic or detailed information are written to fs
            !------------------------------------------------------------------------------
            IF (out_amount /= "PRODUCTION") THEN
                !------------------------------------------------------------------------------
                ! Receive Job_Dir
                !------------------------------------------------------------------------------
                CALL mpi_recv(job_dir, Int(4*mcl,mik), mpi_character, 0_mik, rank_mpi, MPI_COMM_WORLD, status_mpi, ierr)
                CALL print_err_stop(std_out, "MPI_RECV on Domain path didn't succeed", ierr)
            ELSE
                job_dir = outpath
            End if
            
            IF (job_dir(len(job_dir):len(job_dir)) /= "/") job_dir = trim(job_dir)//"/"

            !------------------------------------------------------------------------------
            ! Track the start time of the computation of a domain.
            !------------------------------------------------------------------------------
            CALL CPU_TIME(t_start)

            !==============================================================================
            ! Compute a domain
            !==============================================================================
            CALL exec_single_domain(root, comm_nn, Domain, typeraw, job_dir, fh_cluster_log, &
                Active, fh_mpi_worker, worker_rank_mpi, worker_size_mpi, worker_comm)
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
                ! Set these variables again...
                !------------------------------------------------------------------------------
                pro_path = outpath
                pro_name = project_name     
                                    
                !------------------------------------------------------------------------------
                ! Store header file
                !------------------------------------------------------------------------------
                CALL Write_Tree(root%branches(3)) ! Branch with 'Averaged Material Properties'
                ! ../../../bin/pd_dump_leaf_x86_64 $PWD/ results_0000001 7
            
            !    IF (out_amount == "ALEXANDRIA") THEN
            !        DO ii = 1, SIZE(root%branches)FH01-2_mu_Dev_dtc_Tensors
            !            write(*, '(A, I0, 2A, T80, I0, T84, A)') &
            !                "Branch(", ii, ") of the tree: ", &
            !                TRIM(root%branches(ii)%desc), &
            !                SIZE(root%branches(ii)%leaves), " leaves."
            !        END DO
            !    END IF

                CALL End_Timer("Write Worker Root Branch")

                !------------------------------------------------------------------------------
                ! Store activity information
                !------------------------------------------------------------------------------
                Call MPI_FILE_WRITE_AT(aun, INT(((nn-1) * ik), MPI_OFFSET_KIND), &
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
                INT(root%branches(3)%leaves(1)%lbound-1+(domains_per_comm-1-1), MPI_OFFSET_KIND), &
                INT(comm_nn, KIND=ik), 1_mik, MPI_INTEGER8, status_mpi, ierr)

            !------------------------------------------------------------------------------
            ! Increment linear domain counter of the PETSc sub_comm instances
            ! It is important (!) to increment this variable after the check
            ! IF (Domain_status(nn) <= 0) THEN. Otherwise, the program will write into 
            ! slots that are used by previous domains.
            !------------------------------------------------------------------------------
            comm_nn = comm_nn + 1_ik

        END IF 

        !------------------------------------------------------------------------------
        ! "Deactivate" communicator
        !------------------------------------------------------------------------------
        active = 0_mik

        CALL MPI_ISEND(active, 1_mik, MPI_INTEGER4, 0_mik, rank_mpi, MPI_COMM_WORLD, request, ierr)
        CALL print_err_stop(std_out, "MPI_ISEND on Active didn't succeed", ierr)

        CALL MPI_WAIT(request, status_mpi, ierr)
        CALL print_err_stop(std_out, "MPI_WAIT on request for ISEND Active didn't succeed", ierr)

    End Do

    CALL PetscFinalize(petsc_ierr) 

End If

IF(rank_mpi == 0) THEN

    If (out_amount == "DEBUG") THEN
        Write(un_lf, fmt_dbg_sep)
        Write(un_lf, fmt_MSG_xAI0) "Final Root pointer proc", rank_mpi
        CALL log_tree(root, un_lf, .True.)
        Write(un_lf, fmt_dbg_sep)

        !------------------------------------------------------------------------------
        ! Give information about the tree structure
        !------------------------------------------------------------------------------
        DO ii = 1, SIZE(root%branches)
            write(un_lf, '(A, I0, 2A, T80, I0, T84, A)') &
                "Branch(", ii, ") of the tree: ", &
                TRIM(root%branches(ii)%desc), &
                SIZE(root%branches(ii)%leaves), " leaves."
        END DO
    END If

    CALL link_end(link_name,.True.)

    WRITE(std_out, FMT_TXT_SEP)
    WRITE(std_out, FMT_TXT) "HLRS Direct Tensor Computation finished successfully."
    WRITE(std_out, FMT_TXT_SEP)

END IF ! (rank_mpi == 0)

IF(rank_mpi==0) THEN

    !------------------------------------------------------------------------------
    ! Write the "JOB_FINISHED" keyword only if the job was finished during 
    ! this job (re)run.
    !------------------------------------------------------------------------------
    IF ((Domain_status(No_of_domains) == No_of_domains-1) .OR. (already_finished)) THEN ! counts from 0

        no_restart_required = .TRUE.
        CALL execute_command_line ("echo 'JOB_FINISHED' > BATCH_RUN")
    END IF 

    CALL meta_close(m_rry)

    CALL meta_stop_ascii(fh_mon, mon_suf)

    IF(std_err/=6) CALL meta_stop_ascii(std_out, '.std_out')
    IF(std_err/=0) CALL meta_start_ascii(std_err, '.std_err')

END IF 

1000 Continue

CALL MPI_FILE_CLOSE(aun, ierr)
CALL MPI_FINALIZE(ierr)
CALL print_err_stop(std_out, "MPI_FINALIZE didn't succeed", ierr)

End Program dtc
