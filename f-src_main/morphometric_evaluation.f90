!------------------------------------------------------------------------------
!> morphometric_evaluation - Check number of numerical degrees of freedom per domain
!
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> Date:    09.03.2023
!> LastMod: 11.03.2023
!
!> @brief:
!> Evaluate number of numerical degrees of freedom per domain. 
!> Write the information to file in a useful manner.
!------------------------------------------------------------------------------
PROGRAM morphometric_evaluation

    USE global_std
    USE puredat
    USE meta
    USE user_interaction
    USE formatted_plain
    USE mechanical
    USE ser_binary
    
    IMPLICIT NONE
    
    CHARACTER(*), PARAMETER :: doev_dbg_lvl =  "PRODUCTION" ! "DEBUG"
    
    INTEGER(INT16), DIMENSION(:,:,:), ALLOCATABLE :: rry_ik2
    INTEGER(INT32), DIMENSION(:,:,:), ALLOCATABLE :: rry_ik4
    INTEGER(ik), DIMENSION(:), ALLOCATABLE :: vox_stats, Domains
    !
    ! You can write floats. But you will loose precision by doing that.
    ! INTEGER(ik), DIMENSION(:), ALLOCATABLE :: Domains
    ! REAL(rk), DIMENSION(:), ALLOCATABLE :: vox_stats
    INTEGER(ik), DIMENSION(3) :: xa_d=0, xe_d=0, dims=0, x_D_pos=0, x_D_end, nn_D, x_D
    
    INTEGER(ik) :: counter_BV, counter_TV, bin_hi, bin_lo, bytes, vox_dmn, &
        ii, jj, kk, ll, mm, nn, oo, vun, sun, No_of_domains, size_raw
    
    REAL(rk) :: start, end
    REAL(rk), DIMENSION(3) :: spcng(3), dmn_size, x_D_phy
    
    CHARACTER(mcl), DIMENSION(:), ALLOCATABLE :: m_rry      
    CHARACTER(mcl) :: cmd_arg_history='', binary, vox_file, sun_file, stat=""
    CHARACTER(scl) :: type, type_raw, datarep
    CHARACTER(1) :: bin_sgmnttn=""
    
    CHARACTER(1)  :: restart_cmd_arg='U' ! U = 'undefined'
    
    LOGICAL :: vox_exists=.FALSE., sun_exists=.FALSE., raw_exists=.FALSE., &
        opened=.FALSE., finished_successfully = .FALSE.
    
    CALL CPU_TIME(start)
    std_out = 6 

    !------------------------------------------------------------------------------
    ! Parse the command arguments
    !------------------------------------------------------------------------------
    CALL get_cmd_args(binary, in%full, restart_cmd_arg, cmd_arg_history)
    
    IF (in%full=='') CALL print_err_stop(6, "No input file given", 1)

    !------------------------------------------------------------------------------
    ! Open the given meta file and parse its basename
    !------------------------------------------------------------------------------
    CALL meta_invoke(m_rry)
    
    !------------------------------------------------------------------------------
    ! Spawn standard out after(!) the basename is known
    !------------------------------------------------------------------------------
    IF(std_out/=6) CALL meta_start_ascii(std_out, '.std_out')
    
    CALL show_title(["Morphometric evaluation            ", &
                     "                                   ", &
                     "Johannes Gebert, M.Sc. (HLRS, NUM) "])
    
    !------------------------------------------------------------------------------
    ! Parse input
    !------------------------------------------------------------------------------
    CALL meta_read('LO_BNDS_DMN_RANGE' , m_rry, xa_d, stat)
    IF(stat/="") WRITE(std_out, FMT_ERR) "Reading "//TRIM(stat)//" failed"

    CALL meta_read('UP_BNDS_DMN_RANGE' , m_rry, xe_d, stat)
    IF(stat/="") WRITE(std_out, FMT_ERR) "Reading "//TRIM(stat)//" failed"


    CALL meta_read('BINARIZE_HI', m_rry, bin_hi, stat)
    IF(stat/="") WRITE(std_out, FMT_ERR) "Reading "//TRIM(stat)//" failed"

    CALL meta_read('BINARIZE_LO', m_rry, bin_lo, stat)
    IF(stat/="") WRITE(std_out, FMT_ERR) "Reading "//TRIM(stat)//" failed"

    
    CALL meta_read('SIZE_DOMAIN', m_rry, dmn_size, stat)
    IF(stat/="") WRITE(std_out, FMT_ERR) "Reading "//TRIM(stat)//" failed"

    CALL meta_read('SPACING'    , m_rry, spcng, stat)
    IF(stat/="") WRITE(std_out, FMT_ERR) "Reading "//TRIM(stat)//" failed"

    CALL meta_read('DIMENSIONS' , m_rry, dims, stat)
    IF(stat/="") WRITE(std_out, FMT_ERR) "Reading "//TRIM(stat)//" failed"
    
    CALL meta_read('DATA_BYTE_ORDER', m_rry, type_raw, stat)
    IF(stat/="") WRITE(std_out, FMT_ERR) "Reading "//TRIM(stat)//" failed"
    CALL meta_read('TYPE_RAW'   , m_rry, type, stat)
    IF(stat/="") WRITE(std_out, FMT_ERR) "Reading "//TRIM(stat)//" failed"

    CALL meta_read('BINARY_SEGMENTATION', m_rry, bin_sgmnttn, stat)
    IF(stat/="") WRITE(std_out, FMT_ERR) "Reading "//TRIM(stat)//" failed"


    IF ((bin_sgmnttn /= "Y") .AND. (bin_sgmnttn /= "y") .AND. &
        (bin_sgmnttn /= "N") .AND. (bin_sgmnttn /= "n")) THEN
        CALL print_err_stop(std_out, "Keyword BINARY_SEGMENTATION not recognized.", 1)
    END IF 

    IF (bin_sgmnttn == "y") bin_sgmnttn = "Y"
    
    !------------------------------------------------------------------------------
    ! Warning / Error handling
    !------------------------------------------------------------------------------
    IF ( (dmn_size(1) /= dmn_size(2)) .OR. (dmn_size(1) /= dmn_size(3)) ) THEN
        WRITE(std_out, FMT_ERR) 'Currently, all 3 dimensions of the physical domain size must be equal!'
        STOP
    END IF
    
    ! IF ( (spcng(1) /= spcng(2)) .OR. (spcng(1) /= spcng(3)) ) THEN
    !     mssg = 'Currently, the spacings of all 3 dimensions must be equal!'
    !     STOP
    ! END IF

    IF ( (xa_d(1) > xe_d(1)) .OR. (xa_d(2) > xe_d(2)) .or. (xa_d(3) > xe_d(3)) ) THEN
        WRITE(std_out, FMT_ERR) 'Input parameter error: Start value of domain range larger than end value.'
        STOP
    END IF

    !------------------------------------------------------------------------------
    ! Program breaks if the phdsize is not taking the boundary nodes into account (!).
    ! Therefore, the boundaries are calculated with + 2 Voxels
    !------------------------------------------------------------------------------
    IF ((dmn_size(1) > (dims(1) + 2_ik)*spcng(1)) .OR. & 
        (dmn_size(2) > (dims(2) + 2_ik)*spcng(2)) .OR. & 
        (dmn_size(3) > (dims(3) + 2_ik)*spcng(3))) THEN
        WRITE(std_out, FMT_ERR) 'The domains are larger than the field of view.'
        STOP
    END IF

    
    INQUIRE(UNIT=fhmei, OPENED=opened)
    IF(opened) CLOSE (fhmei)
    
    !------------------------------------------------------------------------------
    ! IMPORTANT: This calculations must be the exact same like in the 
    ! struct_process.f90 main program! Otherwise, the data may be garbage.
    !------------------------------------------------------------------------------
    No_of_domains = (xe_d(1)-xa_d(1)+1) * (xe_d(2)-xa_d(2)+1) * (xe_d(3)-xa_d(3)+1)

    x_D = Nint(dmn_size/spcng)

    x_D_phy = REAL(x_D, rk)*spcng

    nn_D = INT((dims - (2*bpoints)) / x_D, ik)

    vox_dmn = PRODUCT((x_D+(bpoints*2))) 

    ALLOCATE(vox_stats(No_of_domains))
    ALLOCATE(Domains(No_of_domains))

    !------------------------------------------------------------------------------
    ! Open raw file
    !------------------------------------------------------------------------------
    IF(TRIM(type_raw) == "BigEndian") THEN
        datarep = "big_endian"
    ELSE
        datarep = "native"
    END IF 

    INQUIRE(FILE = TRIM(in%p_n_bsnm)//raw_suf, EXIST=raw_exists, SIZE=size_raw)

    IF (.NOT. raw_exists) THEN
        WRITE(std_out, FMT_ERR) "The input *.raw "//TRIM(in%p_n_bsnm)//raw_suf//" file was not found."
        GOTO 1000
    END IF 

    SELECT CASE(type)
        CASE('ik2'); bytes = 2_ik 
        CASE('ik4'); bytes = 4_ik 
    END SELECT

    WRITE(std_out,FMT_TXT_AxI0) "dims:                          ", dims
    WRITE(std_out,FMT_TXT_AxI0) "size:                          ", size_raw
    WRITE(std_out,FMT_TXT_AxI0) "PRODUCT(dims)*bytes:           ", PRODUCT(dims) * bytes
    WRITE(std_out,FMT_TXT_xAI0) "Voxels per domain:             ", vox_dmn
    WRITE(std_out,FMT_TXT_xAI0) "Number of requested domains:   ", No_of_domains
    WRITE(std_out,FMT_TXT_xAI0) "Number of voxels in the model: ", No_of_domains * vox_dmn

    IF (size_raw /= PRODUCT(dims)*bytes) THEN
        WRITE(std_out,FMT_ERR) "The size of the input *.raw file and the dimensions do not match."
        GOTO 1000              
    END IF 

    OPEN (UNIT=51, FILE=TRIM(in%p_n_bsnm)//raw_suf, CONVERT="big_endian", ACCESS="STREAM", STATUS='OLD')

    SELECT CASE(type)
        CASE('ik2') 
            ALLOCATE(rry_ik2(dims(1), dims(2), dims(3)))
            CALL ser_read_binary(51, TRIM(in%p_n_bsnm)//raw_suf, rry_ik2)

        CASE('ik4') 
            ALLOCATE(rry_ik4(dims(1), dims(2), dims(3)))
            CALL ser_read_binary(51, TRIM(in%p_n_bsnm)//raw_suf, rry_ik4)

    END SELECT
    CLOSE(51)
    !------------------------------------------------------------------------------
    ! No need to compute voxs if they are the same for all domains.
    ! Goto is more convenient than nested if/else branches.
    ! Only one GOTO/CONTINUE
    !------------------------------------------------------------------------------
    IF (bin_sgmnttn == "Y") vox_stats = PRODUCT((x_D+(bpoints*2))) 

    !------------------------------------------------------------------------------
    ! Decomposition
    ! Meta contains domain ranges 0-(n-1)
    !------------------------------------------------------------------------------
    oo = 0
    DO kk = xa_d(3), xe_d(3)
    DO jj = xa_d(2), xe_d(2)
    DO ii = xa_d(1), xe_d(1)
        oo = oo + 1_ik

        Domains(oo) = ii + jj * nn_D(1) + kk * nn_D(1)*nn_D(2)
        
        IF(bin_sgmnttn == "Y") CYCLE

        x_D_pos = ANINT(ii * dmn_size / spcng, ik)
        x_D_end = (x_D_pos + x_D)            

        counter_BV = 0_ik
        counter_TV = 0_ik
        Do ll = x_D_pos(3)+1, x_D_end(3)
        Do mm = x_D_pos(2)+1, x_D_end(2)
        Do nn = x_D_pos(1)+1, x_D_end(1)


            ! There are more elegant solutions, but it is a rather transparent approach.
            SELECT CASE(type)
            CASE('ik2') 
                IF ((rry_ik2(nn,mm,ll) >= bin_lo) .AND. (rry_ik2(nn,mm,ll) <= bin_hi)) THEN
                    counter_BV = counter_BV + 1_ik
                END IF 


            CASE('ik4') 
                IF ((rry_ik4(nn,mm,ll) >= bin_lo) .AND. (rry_ik4(nn,mm,ll) <= bin_hi)) THEN
                    counter_BV = counter_BV + 1_ik
                END IF 
    
            END SELECT

            counter_TV = counter_TV + 1_ik 

        END DO
        END DO
        END DO

        ! You can write floats. But you will loose precision by doing that.
        ! vox_stats(oo) = REAL(counter_BV, rk) / REAL(counter_TV, rk)
        vox_stats(oo) = counter_BV
        
        !------------------------------------------------------------------------------
        ! Count domain number
        !------------------------------------------------------------------------------

    END DO
    END DO
    END DO

    5000 CONTINUE
        
    !------------------------------------------------------------------------------
    ! Open *.vox file
    !------------------------------------------------------------------------------
    vun = give_new_unit()
    
    vox_file = TRIM(in%p_n_bsnm)//".vox"
    INQUIRE(FILE = TRIM(vox_file), EXIST=vox_exists)
    
    IF(vox_exists) CALL print_err_stop(std_out, "A vox file already exists.", 1_ik) 
    
    OPEN (UNIT=vun, FILE=TRIM(vox_file), ACCESS="STREAM")

    !------------------------------------------------------------------------------
    ! Open *.status file
    !------------------------------------------------------------------------------
    sun = give_new_unit()
    
    sun_file = TRIM(in%p_n_bsnm)//".status_preprocess" ! compare to DTC
    INQUIRE(FILE = TRIM(sun_file), EXIST=sun_exists)
    
    IF(sun_exists) CALL print_err_stop(std_out, "A status file already exists.", 1_ik) 
    
    OPEN (UNIT=sun, FILE=TRIM(sun_file), ACCESS="STREAM")

    !------------------------------------------------------------------------------
    ! Write everything to file once. Not with every updated domain
    !------------------------------------------------------------------------------
    WRITE(vun) vox_stats
    CLOSE(vun)

    ! Assuming, that no image will contain more than 10 Mio. domains (which is possible in the future!)
    WRITE(sun) -Domains-10000000
    CLOSE(sun)
    
    finished_successfully = .TRUE.    
    1000 CONTINUE
    CALL CPU_TIME(end)
    
    WRITE(std_out, FMT_TXT_SEP)
    
    IF(finished_successfully) THEN 
        WRITE(std_out, FMT_TXT_xAF0) "Program finished successfully in ", end-start, " seconds."
    ELSE
        WRITE(std_out, FMT_TXT_xAF0) "Program stopped without success."
    END IF 

    WRITE(std_out, FMT_TXT_SEP)


    IF(std_out /= 6) CALL meta_stop_ascii(fh=std_out, suf='.std_out')
    
END PROGRAM morphometric_evaluation
    
    