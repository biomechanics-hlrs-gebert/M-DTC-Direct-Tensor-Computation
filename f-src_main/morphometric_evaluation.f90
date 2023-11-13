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
! This program has to read with the provenance basename and to write with the
! actual basename!
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
    
    INTEGER(INT16), DIMENSION(:,:,:), ALLOCATABLE :: Aik2, Aik2ext
    INTEGER(INT32), DIMENSION(:,:,:), ALLOCATABLE :: Aik4, Aik4ext
    INTEGER(INT16), DIMENSION(:,:), ALLOCATABLE :: ConPlane1_ik2,ConPlane2_ik2,ConPlane3_ik2
    INTEGER(INT32), DIMENSION(:,:), ALLOCATABLE :: ConPlane1_ik4,ConPlane2_ik4,ConPlane3_ik4

    INTEGER(ik), DIMENSION(:), ALLOCATABLE :: vox_stats, Domains
    INTEGER(ik), DIMENSION(3), PARAMETER :: bpoints=[1_ik,1_ik,1_ik]
    !
    ! You can write floats. But you will loose precision by doing that.
    ! INTEGER(ik), DIMENSION(:), ALLOCATABLE :: Domains
    ! REAL(rk), DIMENSION(:), ALLOCATABLE :: vox_stats
    INTEGER(ik), DIMENSION(3) :: xa_d=0, xe_d=0, dims=0, x_D_pos=0, &
        x_D_end, nn_D, x_D, vox_min, vox_max
    
    INTEGER(ik) :: counter_BV, counter_TV, hi, lo, bytes=0, vox_dmn, &
        ii, jj, kk, ll, mm, nn, oo, vun, sun, sixfun, twsixfun, coplun, &
        No_of_domains, size_raw, cntr, tt = 0_ik, cx, cy, cz

    REAL(rk) :: start, end, sixf_tmp, twsixf_tmp, sixf_tmp_vox, twsixf_tmp_vox, &
        CoPl1Avg, CoPl2Avg, CoPl3Avg, max_copl, min_copl
    REAL(rk), DIMENSION(:), ALLOCATABLE :: sixf, twsixf, CoPl
    REAL(rk), DIMENSION(3) :: spcng(3), dmn_size, x_D_phy
    
    CHARACTER(mcl), DIMENSION(:), ALLOCATABLE :: m_rry      
    CHARACTER(mcl) :: cmd_arg_history='', binary, vox_file, sun_file, raw_file, &
        stat="", twsixf_file, sixf_file, copl_file
    CHARACTER(scl) :: type_raw, data_byte_order, datarep
    CHARACTER(1) :: bin_sgmnttn=""
    
    CHARACTER(1)  :: restart_cmd_arg='U' ! U = 'undefined'
    
    LOGICAL :: fex=.FALSE., vox_found = .FALSE., &
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
    ! Has to be given in the programs scope! And it refers to the tensor computation
    global_meta_prgrm_mstr_app = 'dtc' 
    global_meta_program_keyword = 'TENSOR_COMPUTATION'

    CALL meta_append(m_rry, 1_4, binary, stat)
    
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
    stat = ""
    CALL meta_read('LO_BNDS_DMN_RANGE' , m_rry, xa_d, stat)
    IF(stat/="") WRITE(std_out, FMT_ERR) "Reading "//TRIM(stat)//" failed"

    CALL meta_read('UP_BNDS_DMN_RANGE' , m_rry, xe_d, stat)
    IF(stat/="") WRITE(std_out, FMT_ERR) "Reading "//TRIM(stat)//" failed"


    CALL meta_read('BINARIZE_HI', m_rry, hi, stat)
    IF(stat/="") WRITE(std_out, FMT_ERR) "Reading "//TRIM(stat)//" failed"

    CALL meta_read('BINARIZE_LO', m_rry, lo, stat)
    IF(stat/="") WRITE(std_out, FMT_ERR) "Reading "//TRIM(stat)//" failed"

    
    CALL meta_read('SIZE_DOMAIN', m_rry, dmn_size, stat)
    IF(stat/="") WRITE(std_out, FMT_ERR) "Reading "//TRIM(stat)//" failed"

    CALL meta_read('SPACING'    , m_rry, spcng, stat)
    IF(stat/="") WRITE(std_out, FMT_ERR) "Reading "//TRIM(stat)//" failed"

    CALL meta_read('DIMENSIONS' , m_rry, dims, stat)
    IF(stat/="") WRITE(std_out, FMT_ERR) "Reading "//TRIM(stat)//" failed"
    
    CALL meta_read('DATA_BYTE_ORDER', m_rry, data_byte_order, stat)
    IF(stat/="") WRITE(std_out, FMT_ERR) "Reading "//TRIM(stat)//" failed"

    CALL meta_read('TYPE_RAW'   , m_rry, type_raw, stat)
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
    
    vun = 61_ik
    vox_file = TRIM(out%p_n_bsnm)//".vox"
    INQUIRE(FILE = TRIM(vox_file), EXIST=fex)
    IF(fex) CALL print_err_stop(std_out, "A vox file already exists.", 1_ik) 
    
    sun = 62_ik
    sun_file = TRIM(out%p_n_bsnm)//".status" ! compare to DTC
    INQUIRE(FILE = TRIM(sun_file), EXIST=fex)
    IF(fex) CALL print_err_stop(std_out, "A status file already exists.", 1_ik)    

    sixfun = 63_ik
    sixf_file = TRIM(out%p_n_bsnm)//".6fold" ! compare to DTC
    INQUIRE(FILE = TRIM(sixf_file), EXIST=fex)
    IF(fex) CALL print_err_stop(std_out, "A status file already exists.", 1_ik)    

    twsixfun = 64_ik
    twsixf_file = TRIM(out%p_n_bsnm)//".27fold" ! compare to DTC
    INQUIRE(FILE = TRIM(twsixf_file), EXIST=fex)
    IF(fex) CALL print_err_stop(std_out, "A status file already exists.", 1_ik) 

    CoPlun = 64_ik
    CoPl_file = TRIM(out%p_n_bsnm)//".CoPl" ! compare to DTC
    INQUIRE(FILE = TRIM(Copl_file), EXIST=fex)
    IF(fex) CALL print_err_stop(std_out, "A control plane file already exists.", 1_ik) 


    WRITE(std_out, FMT_TXT) "vox-file:      ", TRIM(vox_file)
    WRITE(std_out, FMT_TXT) "sun-file:      ", TRIM(sun_file)
    WRITE(std_out, FMT_TXT) "6-fold-file:   ", TRIM(sixf_file)
    WRITE(std_out, FMT_TXT) "27-fold-file:  ", TRIM(twsixf_file)
    WRITE(std_out, FMT_TXT) "CoPl-file:     ", TRIM(CoPl_file)

    !------------------------------------------------------------------------------
    ! IMPORTANT: This calculations must be the exact same like in the 
    ! struct_process.f90 main program! Otherwise, the data may be garbage.
    !------------------------------------------------------------------------------
    No_of_domains = (xe_d(1)-xa_d(1)+1) * (xe_d(2)-xa_d(2)+1) * (xe_d(3)-xa_d(3)+1)

    x_D = Nint(dmn_size/spcng)

    x_D_phy = REAL(x_D, rk)*spcng

    nn_D = INT((dims - (2*bpoints)) / x_D, ik)

    vox_dmn = PRODUCT((x_D+(bpoints*2))) 

    ALLOCATE(vox_stats(No_of_domains)); vox_stats = 0_ik
    ALLOCATE(Domains(No_of_domains)); Domains = 0_ik
    ALLOCATE(sixf(No_of_domains)); sixf = 0._rk
    ALLOCATE(twsixf(No_of_domains)); twsixf = 0._rk
    ALLOCATE(copl(No_of_domains)); copl = 0._rk

    !------------------------------------------------------------------------------
    ! Open raw file
    !------------------------------------------------------------------------------
    IF(TRIM(data_byte_order) == "BigEndian") THEN
        datarep = "big_endian"
    ELSE
        datarep = "native"
    END IF 

    raw_file = TRIM(in%p_n_bsnm)//raw_suf
    WRITE(std_out, FMT_TXT) "raw-file:      ", TRIM(raw_file)

    INQUIRE(FILE = TRIM(raw_file), EXIST=fex, SIZE=size_raw)

    IF (.NOT. fex) THEN
        WRITE(std_out, FMT_ERR) "The input *.raw "//TRIM(raw_file)//" file was not found."
        STOP
    END IF 

    SELECT CASE(type_raw)
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

    OPEN (UNIT=51, FILE=TRIM(in%p_n_bsnm)//raw_suf, CONVERT=datarep, ACCESS="STREAM", STATUS='OLD')

    SELECT CASE(type_raw)
        CASE('ik2') 
            ALLOCATE(Aik2(dims(1), dims(2), dims(3))); Aik2 = 0_2
            READ(UNIT=51) Aik2

            ! ALLOCATE(ConPlane1_ik2(x_D(1), x_D(2))); ConPlane1_ik2 = 0_2
            ! ALLOCATE(ConPlane2_ik2(x_D(1), x_D(3))); ConPlane2_ik2 = 0_2
            ! ALLOCATE(ConPlane3_ik2(x_D(2), x_D(3))); ConPlane3_ik2 = 0_2
            ALLOCATE(ConPlane1_ik2(x_D(1)+1_ik, x_D(2)+1_ik)); ConPlane1_ik2 = 0_2
            ALLOCATE(ConPlane2_ik2(x_D(1)+1_ik, x_D(3)+1_ik)); ConPlane2_ik2 = 0_2
            ALLOCATE(ConPlane3_ik2(x_D(2)+1_ik, x_D(3)+1_ik)); ConPlane3_ik2 = 0_2

        CASE('ik4') 
            ALLOCATE(Aik4(dims(1), dims(2), dims(3))); Aik4 = 0_4
            READ(UNIT=51) Aik4

            ! ALLOCATE(ConPlane1_ik4(x_D(1), x_D(2))); ConPlane1_ik4 = 0_2
            ! ALLOCATE(ConPlane2_ik4(x_D(1), x_D(3))); ConPlane2_ik4 = 0_2
            ! ALLOCATE(ConPlane3_ik4(x_D(2), x_D(3))); ConPlane3_ik4 = 0_2
            ALLOCATE(ConPlane1_ik4(x_D(1)+1_ik, x_D(2)+1_ik)); ConPlane1_ik4 = 0_2
            ALLOCATE(ConPlane2_ik4(x_D(1)+1_ik, x_D(3)+1_ik)); ConPlane2_ik4 = 0_2
            ALLOCATE(ConPlane3_ik4(x_D(2)+1_ik, x_D(3)+1_ik)); ConPlane3_ik4 = 0_2

        END SELECT

    CLOSE(51)
    !------------------------------------------------------------------------------
    ! No need to compute voxs if they are the same for all domains.
    ! Goto is more convenient than nested if/else branches.
    ! Only one GOTO/CONTINUE
    !------------------------------------------------------------------------------
    IF (bin_sgmnttn == "Y") vox_stats = PRODUCT((x_D+(bpoints*2))) 


    vox_min = INT(xa_d * dmn_size / spcng, ik) + 1_ik
    vox_max = INT(xe_d * dmn_size / spcng, ik) + 1_ik
    WRITE(std_out,FMT_TXT_AxI0) "vox_min: ", vox_min
    WRITE(std_out,FMT_TXT_AxI0) "vox_max: ", vox_max



    !------------------------------------------------------------------------------
    ! Decomposition
    ! Meta contains domain ranges 0-(n-1)
    !------------------------------------------------------------------------------

    WRITE(std_out,FMT_SEP)
    WRITE(std_out,"(A)") "-- Done with domain"
    WRITE(std_out,"(A)", ADVANCE='NO') "-- "

    oo = 1
    DO kk = xa_d(3), xe_d(3) 
    DO jj = xa_d(2), xe_d(2) 
    DO ii = xa_d(1), xe_d(1) 

        Domains(oo) = ii + jj * (nn_D(1)+1) + kk * (nn_D(1)+1) * (nn_D(2)+1)

        IF(bin_sgmnttn == "N") CYCLE

        x_D_pos = INT([ii,jj,kk] * dmn_size / spcng, ik) + 1_ik
        x_D_end = (x_D_pos + x_D)            

        counter_BV = 0_ik
        counter_TV = 0_ik

        ! Workaround
        IF(x_D_end(1) > dims(1)) x_D_end(1) = dims(1)
        IF(x_D_end(2) > dims(2)) x_D_end(2) = dims(2)
        IF(x_D_end(3) > dims(3)) x_D_end(3) = dims(3)

        Do ll = x_D_pos(3), x_D_end(3)
        Do mm = x_D_pos(2), x_D_end(2)
        Do nn = x_D_pos(1), x_D_end(1)

            ! There are more elegant solutions, but it is a rather transparent approach.
            SELECT CASE(type_raw)
            CASE('ik2') 
                IF ((Aik2(nn,mm,ll) >= lo) .AND. (Aik2(nn,mm,ll) <= hi)) THEN
                    counter_BV = counter_BV + 1_ik
                    Aik2(nn,mm,ll) = 1_2
                ELSE
                    Aik2(nn,mm,ll) = 0_2
                END IF 


            CASE('ik4') 
                IF ((Aik4(nn,mm,ll) >= lo) .AND. (Aik4(nn,mm,ll) <= hi)) THEN
                    counter_BV = counter_BV + 1_ik
                    Aik4(nn,mm,ll) = 1_4
                ELSE
                    Aik4(nn,mm,ll) = 0_4
                END IF 
    
            END SELECT

            counter_TV = counter_TV + 1_ik 

        END DO
        END DO
        END DO



        !------------------------------------------------------------------------------
        ! NHS
        !------------------------------------------------------------------------------

        ! There are more elegant solutions, but it is a rather transparent approach.
        SELECT CASE(type_raw)
        CASE('ik2') 
            cy=1_ik
            Do mm = x_D_pos(2), x_D_end(2)
                cx=1_ik
                Do nn = x_D_pos(1), x_D_end(1)
                    ConPlane1_ik2(cx, cy) = SUM(Aik2(nn,mm,:))
                    cx = cx + 1_ik
                END DO
                cy = cy + 1_ik
            END DO
            CoPl1Avg = REAL(SUM(ConPlane1_ik2)/(x_D(1)*x_D(2)),rk)

            cz=1_ik
            DO ll = x_D_pos(3), x_D_end(3)
                cx=1_ik
                Do nn = x_D_pos(1), x_D_end(1)
                    ConPlane2_ik2(cx, cz) = SUM(Aik2(nn,:,ll))
                    cx = cx + 1_ik
                END DO
                cz = cz + 1_ik
            END DO
            CoPl2Avg = REAL(SUM(ConPlane1_ik2)/(x_D(1)*x_D(3)),rk)
            
            cz=1_ik
            DO ll = x_D_pos(3), x_D_end(3)
                cy=1_ik
                Do mm = x_D_pos(2), x_D_end(2)
                    ConPlane3_ik2(cy, cz) = SUM(Aik2(:,mm,ll))
                    cy = cy + 1_ik
                END DO
                cz = cz + 1_ik
            END DO
            CoPl3Avg = REAL(SUM(ConPlane1_ik2)/(x_D(2)*x_D(3)),rk)

        CASE('ik4')
            cy=1_ik
            Do mm = x_D_pos(2), x_D_end(2)
                cx=1_ik
                Do nn = x_D_pos(1), x_D_end(1)
                    ConPlane1_ik4(cx, cy) = SUM(Aik4(nn,mm,:))
                    cx = cx + 1_ik
                END DO
                cy = cy + 1_ik
            END DO
            CoPl1Avg = REAL(SUM(ConPlane1_ik4)/(x_D(1)*x_D(2)),rk)

            cz=1_ik
            DO ll = x_D_pos(3), x_D_end(3)
                cx=1_ik
                Do nn = x_D_pos(1), x_D_end(1)
                    ConPlane2_ik4(cx, cz) = SUM(Aik4(nn,:,ll))
                    cx = cx + 1_ik
                END DO
                cz = cz + 1_ik
            END DO
            CoPl2Avg = REAL(SUM(ConPlane1_ik4)/(x_D(1)*x_D(3)),rk)
            
            cz=1_ik
            DO ll = x_D_pos(3), x_D_end(3)
                cy=1_ik
                Do mm = x_D_pos(2), x_D_end(2)
                    ConPlane3_ik4(cy, cz) = SUM(Aik4(:,mm,ll))
                    cy = cy + 1_ik
                END DO
                cx = cx + 1_ik
            END DO
            CoPl3Avg = REAL(SUM(ConPlane1_ik4)/(x_D(2)*x_D(3)),rk)

        END SELECT

        IF ((CoPl1Avg > CoPl2Avg) .AND. (CoPl1Avg > CoPl3Avg)) max_CoPl = CoPl1Avg
        IF ((CoPl2Avg > CoPl1Avg) .AND. (CoPl2Avg > CoPl3Avg)) max_CoPl = CoPl2Avg
        IF ((CoPl3Avg > CoPl1Avg) .AND. (CoPl3Avg > CoPl2Avg)) max_CoPl = CoPl3Avg

        IF ((CoPl1Avg < CoPl2Avg) .AND. (CoPl1Avg < CoPl3Avg)) min_CoPl = CoPl1Avg
        IF ((CoPl2Avg < CoPl1Avg) .AND. (CoPl2Avg < CoPl3Avg)) min_CoPl = CoPl2Avg
        IF ((CoPl3Avg < CoPl1Avg) .AND. (CoPl3Avg < CoPl2Avg)) min_CoPl = CoPl3Avg

        Copl(oo) = max_CoPl/min_CoPl


        sixf_tmp = 0._rk
        twsixf_tmp = 0._rk

        DO ll = x_D_pos(3)+1_ik, x_D_end(3)-1_ik
        DO mm = x_D_pos(2)+1_ik, x_D_end(2)-1_ik
        DO nn = x_D_pos(1)+1_ik, x_D_end(1)-1_ik

            ! There are more elegant solutions, but it is a rather transparent approach.
            SELECT CASE(type_raw)
            CASE('ik2') 
                IF (Aik2(nn,mm,ll) == 1_2) THEN
                    vox_found = .TRUE.

                    sixf_tmp = REAL((Aik2(nn-1_ik, mm     , ll     )) + &
                                    (Aik2(nn     , mm-1_ik, ll     )) + &
                                    (Aik2(nn     , mm     , ll-1_ik)) + &
                                    (Aik2(nn+1_ik, mm     , ll     )) + &
                                    (Aik2(nn     , mm+1_ik, ll     )) + &
                                    (Aik2(nn     , mm     , ll+1_ik)), rk)
                    
                    twsixf_tmp = (REAL(SUM(Aik2(nn-1_ik:nn+1_ik,mm-1_ik:mm+1_ik,ll-1_ik:ll+1_ik)), rk))
                    
                END IF 

            CASE('ik4') 
                IF (Aik4(nn,mm,ll)  == 1_4) THEN
                    vox_found = .TRUE.

                    sixf_tmp = REAL((Aik4(nn-1_ik, mm     , ll     )) + &
                                    (Aik4(nn     , mm-1_ik, ll     )) + &
                                    (Aik4(nn     , mm     , ll-1_ik)) + &
                                    (Aik4(nn+1_ik, mm     , ll     )) + &
                                    (Aik4(nn     , mm+1_ik, ll     )) + &
                                    (Aik4(nn     , mm     , ll+1_ik)),rk)
                        
                    twsixf_tmp = (REAL(SUM(Aik4(nn-1_ik:nn+1_ik,mm-1_ik:mm+1_ik,ll-1_ik:ll+1_ik)), rk))

                END IF 
    
            END SELECT

        END DO
        END DO
        END DO
        
        !------------------------------------------------------------------------------
        ! User Feedback
        !------------------------------------------------------------------------------
        WRITE(std_out,"(I0, A)", ADVANCE='NO') Domains(oo), " "
        tt = tt+1_ik
        if (tt == 20_ik) THEN
            tt = 0_ik
            WRITE(std_out,"(A)") ""
            WRITE(std_out,"(A)", ADVANCE='NO') "-- "
        END IF 

        !------------------------------------------------------------------------------
        ! Write morphometric quantities
        !------------------------------------------------------------------------------
        ! sixf(oo) = sixf_tmp
        ! twsixf(oo) = twsixf_tmp

        sixf(oo) = sixf_tmp
        twsixf(oo) = twsixf_tmp

        ! You can write floats. But you will loose precision by doing that.
        ! vox_stats(oo) = REAL(counter_BV, rk) / REAL(counter_TV, rk)
        vox_stats(oo) = counter_BV
        
        !------------------------------------------------------------------------------
        ! Count domain number
        !------------------------------------------------------------------------------
        oo = oo + 1_ik

    END DO
    END DO
    END DO
        
    OPEN (UNIT=vun, FILE=TRIM(vox_file), ACCESS="STREAM")
    OPEN (UNIT=sun, FILE=TRIM(sun_file), ACCESS="STREAM")
    OPEN (UNIT=sixfun, FILE=TRIM(sixf_file), ACCESS="STREAM")
    OPEN (UNIT=twsixfun, FILE=TRIM(twsixf_file), ACCESS="STREAM")
    OPEN (UNIT=coplun, FILE=TRIM(copl_file), ACCESS="STREAM")

    WRITE(vun) vox_stats
    ! Assuming, that no image will contain more than 10 Mio. domains (which is possible in the future!)
    WRITE(sun) (-Domains-100000000)

    WRITE(sixfun) (sixf)
    WRITE(twsixfun) (twsixf)
    WRITE(coplun) (copl)
    
    CLOSE(vun)
    CLOSE(sun)
    CLOSE(sixfun)
    CLOSE(twsixfun)
    CLOSE(coplun)
    
    finished_successfully = .TRUE.    
    1000 CONTINUE
    CALL CPU_TIME(end)
    
    WRITE(std_out,"(A)") ""
    WRITE(std_out, FMT_TXT_SEP)

    IF(finished_successfully) THEN 
        WRITE(std_out, FMT_TXT_xAF0) "Program finished successfully in ", end-start, " seconds."
    ELSE
        WRITE(std_out, FMT_TXT_xAF0) "Program stopped without success."
    END IF 

    WRITE(std_out, FMT_TXT_SEP)


    IF(std_out /= 6) CALL meta_stop_ascii(fh=std_out, suf='.std_out')
    
END PROGRAM morphometric_evaluation
    
    