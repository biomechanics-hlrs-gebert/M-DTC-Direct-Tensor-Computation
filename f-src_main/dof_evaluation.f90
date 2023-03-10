!------------------------------------------------------------------------------
!> dof_evaluation - Check number of numerical degrees of freedom per domain
!
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> Date:    09.03.2023
!> LastMod: 09.03.2023
!
!> @brief:
!> Evaluate number of numerical degrees of freedom per domain. 
!> Write the information to file in a useful manner.
!------------------------------------------------------------------------------
PROGRAM dof_evaluation

    USE global_std
    USE puredat
    USE meta
    USE user_interaction
    USE formatted_plain
    USE mechanical
    USE ser_binary
    
    IMPLICIT NONE
    
    CHARACTER(*), PARAMETER :: doev_dbg_lvl =  "PRODUCTION" ! "DEBUG"
    
    TYPE(tLeaf), POINTER :: domain_no, eff_num_stiffness, density, tensors, leaf_pos, &
        no_elems, no_nodes, t_start, t_duration, collected_logs
    TYPE(tBranch)        :: rank_data
    TYPE(domain_data), DIMENSION(3) :: tensor
    
    INTEGER(INT16), DIMENSION(:,:,:), ALLOCATABLE :: rry_ik2
    INTEGER(INT32), DIMENSION(:,:,:), ALLOCATABLE :: rry_ik4
    INTEGER(ik), DIMENSION(:,:), ALLOCATABLE :: Domain_status
    INTEGER(ik), DIMENSION(3) :: xa_d=0, xe_d=0, dims=0, grid=0
    
    INTEGER(ik) :: alloc_stat, parts, fh_covo, fh_cr1, fh_cr2, fh_tens, &
        domains_crawled = 0_ik, nn_comm, par_domains, dof_size, size_mpi, fh_covo_num, &
        ii, kk, ll, mm, tt, xx, par_dmn_number, jj, current_domain, pntr, no_lc, ma_el_order, &
        dun, domains_per_comm = 0, rank_mpi, No_of_domains, local_domain_no, last_domain_rank
    
    REAL(rk) :: local_domain_density, sym, start, end
    REAL(rk), DIMENSION(3)   :: local_domain_opt_pos, spcng, dmn_size
    REAL(rk), DIMENSION(6,6) :: local_domain_tensor
    REAL(rk), DIMENSION(:,:), ALLOCATABLE :: local_num_tensor
    
    INTEGER(8), DIMENSION(:), ALLOCATABLE :: dat_domains, dat_no_elems, dat_no_nodes, dat_collected_logs
    REAL(8), DIMENSION(:), ALLOCATABLE :: dat_densities, dat_eff_num_stiffnesses, dat_tensors, &
        dat_pos, dat_t_start, dat_t_duration
    
    CHARACTER(mcl), DIMENSION(:), ALLOCATABLE :: m_rry      
    CHARACTER(mcl) :: cmd_arg_history='', target_path, file_head, binary, dof_file, stat=""
    CHARACTER(scl) :: mi_el_type, type
    
    CHARACTER(100*mcl) :: string
    
    CHARACTER(20)  :: dmn_char
    CHARACTER(1)   :: tt_char, restart_cmd_arg='U' ! U = 'undefined'
    
    LOGICAL :: success=.FALSE., dof_exists=.FALSE., opened=.FALSE., last_domain=.FALSE.
    
    CALL CPU_TIME(start)
    std_out = 6 

    !------------------------------------------------------------------------------
    ! Parse the command arguments
    !------------------------------------------------------------------------------
    CALL get_cmd_args(binary, in%full, restart_cmd_arg, cmd_arg_history)
    
    IF (in%full=='') THEN
        ! CALL usage(binary, [ &
        !     "Give the input meta File of the dataset.                 "])    
    
        CALL print_err_stop(6, "No input file given", 1)
    END IF
    
    !------------------------------------------------------------------------------
    ! Open the given meta file and parse its basename
    !------------------------------------------------------------------------------
    CALL meta_invoke(m_rry)
    ! Normally, one would call meta_append which wraps meta_invoke and meta_continue.
    ! However, meta_continue provides the assignment out=in which is required to set 
    ! the proper basename for meta_start_ascii.
    out = in
    
    !------------------------------------------------------------------------------
    ! Spawn standard out after(!) the basename is known
    !------------------------------------------------------------------------------
    IF(std_out/=6) CALL meta_start_ascii(std_out, '.std_out')
    
    CALL show_title(["Johannes Gebert, M.Sc. (HLRS, NUM) "])
    
    !------------------------------------------------------------------------------
    ! Parse input
    !------------------------------------------------------------------------------
    CALL meta_read('LO_BNDS_DMN_RANGE' , m_rry, xa_d, stat); IF(stat/="") STOP
    CALL meta_read('UP_BNDS_DMN_RANGE' , m_rry, xe_d, stat); IF(stat/="") STOP
    CALL meta_read('MESH_PER_SUB_DMN'  , m_rry, parts, stat); IF(stat/="") STOP
    
    CALL meta_read('MACRO_ELMNT_ORDER' , m_rry, ma_el_order, stat); IF(stat/="") STOP
    CALL meta_read('MICRO_ELMNT_TYPE'  , m_rry, mi_el_type, stat); IF(stat/="") STOP
    
    CALL meta_read('PROCESSORS' , m_rry, size_mpi, stat); IF(stat/="") STOP
    CALL meta_read('SIZE_DOMAIN', m_rry, dmn_size, stat); IF(stat/="") STOP
    CALL meta_read('SPACING'    , m_rry, spcng, stat); IF(stat/="") STOP
    CALL meta_read('DIMENSIONS' , m_rry, dims, stat); IF(stat/="") STOP
    CALL meta_read('TYPE_RAW'   , m_rry, type, stat); IF(stat/="") STOP

    
    IF (ma_el_order == 1) THEN
        no_lc = 24
    ELSE IF (ma_el_order == 2) THEN
        no_lc = 60
    ELSE
        CALL print_err_stop(std_out, "Macro element order supported.", 1)
    END IF 
    
    ALLOCATE(local_num_tensor(no_lc, no_lc))
    
    INQUIRE(UNIT=fhmei, OPENED=opened)
    IF(opened) CLOSE (fhmei)
    
    !------------------------------------------------------------------------------
    ! IMPORTANT: This calculations must be the exact same like in the 
    ! struct_process.f90 main program! Otherwise, the data may be garbage.
    !------------------------------------------------------------------------------
    No_of_domains = (xe_d(1)-xa_d(1)+1) * (xe_d(2)-xa_d(2)+1) * (xe_d(3)-xa_d(3)+1)
    
    !------------------------------------------------------------------------------
    ! Open *.dof file
    !------------------------------------------------------------------------------
    dun = give_new_unit()
    
    dof_file = TRIM(in%p_n_bsnm)//".dof"
    INQUIRE(FILE = TRIM(dof_file), EXIST=dof_exists, SIZE=dof_size)
    
    IF(dof_exists) CALL print_err_stop(std_out, "A dof file already exists.", 1_ik) 
    
    OPEN (UNIT=dun, FILE=TRIM(dof_file), ACCESS="STREAM")
    

    !------------------------------------------------------------------------------
    ! Open raw file
    !------------------------------------------------------------------------------
    SELECT CASE(type)
        CASE('ik2') 
            ALLOCATE(rry_ik2(dims(1), dims(2), dims(3)))
            CALL ser_read_binary(51, TRIM(in%p_n_bsnm)//raw_suf, rry_ik2)

        CASE('ik4') 
            ALLOCATE(rry_ik4(dims(1), dims(2), dims(3)))
            CALL ser_read_binary(51, TRIM(in%p_n_bsnm)//raw_suf, rry_ik4)

    END SELECT



    !------------------------------------------------------------------------------
    ! User feedback
    !------------------------------------------------------------------------------
    WRITE(std_out, FMT_MSG_AxI0) "size_mpi:         ", size_mpi
    WRITE(std_out, FMT_MSG_AxI0) "parts:            ", parts
    WRITE(std_out, FMT_MSG_AxI0) "No_of_domains:    ", No_of_domains
    
    CALL CPU_TIME(end)
    
    WRITE(std_out, FMT_TXT_SEP)
    WRITE(std_out, FMT_TXT_xAF0) "Program finished successfully in ", end-start, " seconds."
    WRITE(std_out, FMT_TXT_SEP)
    
    IF(std_out /= 6) CALL meta_stop_ascii(fh=std_out, suf='.std_out')
    
    END PROGRAM dof_evaluation
    
    