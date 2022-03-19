!------------------------------------------------------------------------------
!> MeRaDat - Crawl Tensors from the Direct Tensor Computation
!
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> Date:    19.03.2022
!> LastMod: 19.03.2022
!------------------------------------------------------------------------------
PROGRAM mrd_crawl_tensors

USE global_std
USE puredat
USE meta
USE user_interaction
USE formatted_plain

IMPLICIT NONE

Real(Kind=pd_rk) :: gstart_time, gend_time

Integer          :: num_args, leaf_no, alloc_stat

CHARACTER(LEN=mcl) :: arg

TYPE(tBranch)        :: tree
TYPE(tLeaf), POINTER :: leaf
LOGICAL              :: success

INTEGER(KIND=1), DIMENSION(:), ALLOCATABLE :: dat_int1
INTEGER(KIND=2), DIMENSION(:), ALLOCATABLE :: dat_int2
INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: dat_int4
INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: dat_int8
REAL   (KIND=8), DIMENSION(:), ALLOCATABLE :: dat_real8
CHARACTER      , DIMENSION(:), ALLOCATABLE :: dat_char

CHARACTER(LEN=5) :: ten_suf = ".glco"


CALL CPU_TIME(global_start)

!------------------------------------------------------------------------------
! Parse the command arguments
!------------------------------------------------------------------------------
CALL get_cmd_args(binary, in%full, stp, restart_cmd_arg, cmd_arg_history)
IF(stp) GOTO 1000

IF (in%full=='') THEN
    CALL usage(binary)    

    !------------------------------------------------------------------------------
    ! On std_out since file of std_out is not spawned
    !------------------------------------------------------------------------------
    CALL print_err_stop(6, "No input file given", 1)
END IF

!------------------------------------------------------------------------------
! Open the given meta file and parse its basename
!------------------------------------------------------------------------------
CALL meta_invoke(m_rry)
CALL parse_basename(in%full, meta_suf)

!------------------------------------------------------------------------------
! Redirect std_out into a file in case std_out is not useful by environment.
!------------------------------------------------------------------------------
std_out = determine_stout()

!------------------------------------------------------------------------------
! Spawn standard out after(!) the basename is known
!------------------------------------------------------------------------------
IF(std_out/=6) CALL meta_start_ascii(std_out, '.std_out')

CALL show_title(["Johannes Gebert, M.Sc. (HLRS, NUM) "], "Crawl Tensors of DTC")

!------------------------------------------------------------------------------
! Open the file to store the tensors within their global coordinate system into
!------------------------------------------------------------------------------
CALL meta_start_ascii(fh_glco, ten_suf)

!------------------------------------------------------------------------------
! Parse input
!------------------------------------------------------------------------------
CALL meta_read('TENSOR-COMPUTATION', m_rry, dev_null)

CALL meta_read('LO_BNDS_DMN_RANGE' , m_rry, xa_d)
CALL meta_read('UP_BNDS_DMN_RANGE' , m_rry, xe_d)
CALL meta_read('MESH_PER_SUB_DMN'  , m_rry, parts)
CALL meta_read('PROCESSORS'        , m_rry, processors)

CALL meta_read('YOUNG_MODULUS'     , m_rry, bone%E)
CALL meta_read('POISSON_RATIO'     , m_rry, bone%nu)

!------------------------------------------------------------------------------
! Calculate crawling conditions
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Crawl Ranks
!------------------------------------------------------------------------------







!------------------------------------------------------------------------------
! Finish program
!------------------------------------------------------------------------------
CALL meta_stop_ascii(fh_glco, ten_suf)
CALL meta_stop_ascii(fhmei, meta_suf)

WRITE(std_out, FMT_TXT) 'Program finished successfully.'
WRITE(std_out, FMT_TXT_SEP)

IF(std_out /= 6) CALL meta_stop_ascii(fh=std_out, suf='.std_out')

END PROGRAM mrd_crawl_tensors




! CALL Cpu_time(gstart_time)


! num_args = command_argument_count()

! If (num_args /= 3) then
!     Write(*,FMT_TXT_SEP)
!     Write(*,FMT_TXT) "Usage:"
!     Write(*,FMT_TXT) "arg 1: Path of the branch."
!     Write(*,FMT_TXT) "arg 2: basename of the branch."
!     Write(*,FMT_TXT) "arg 3: Number of leaf to be dumped"
!     Write(*,FMT_TXT_SEP)
!     Stop
! End If

! CALL GET_COMMAND_ARGUMENT(1, pro_path)
! CALL GET_COMMAND_ARGUMENT(2, pro_name)
! CALL GET_COMMAND_ARGUMENT(3, arg)

! !   pro_name="/results"

! Read(arg,*)leaf_no

! write(*,PDF_SEP  )
! write(*,PDF_M_A  ) "Pro path: "//trim(pro_path)
! write(*,PDF_M_A  ) "Pro name: "//trim(pro_name)
! write(*,PDF_M_AI0) "Leaf no.: ",leaf_no
! write(*,PDF_M_A  )"=="

! tree = read_tree()
! CALL open_stream_files(tree, "read" , "old")
! CALL get_leaf_with_num(tree, leaf_no, leaf, success)

! If (.not. Success) then
!     Write(*,PDF_E_AI0)"There is no leaf with number ",leaf_no
! ELSE

!     write(*,PDF_M_A  ) "Description          : "//trim(leaf%desc)
!     write(*,PDF_M_AI0) "Number of data       :", leaf%dat_no
!     write(*,PDF_M_AI0) "Type of data         :", leaf%dat_ty
!     write(*,PDF_M_AI0) "Lower bound in stream:", leaf%lbound
!     write(*,PDF_M_AI0) "Upper bound in stream:", leaf%ubound
!     write(*,PDF_M_A  )"=="

!     Select case (leaf%dat_ty)
    
!     Case(1)

!         Allocate(dat_int1(leaf%dat_no),stat=alloc_stat)
!         If (alloc_stat /= 0) Then
!             WRITE(*,*)
!             WRITE(*,PDF_SEP)
!             WRITE(*,PDF_E_A)   'Allocation of the dat_int1 failed !!'
!             WRITE(*,PDF_E_AI0) 'The requested dimension was : ',leaf%dat_no
!             WRITE(*,PDF_E_STOP)
!             STOP
!         End If

!         CALL pd_read_leaf(tree%streams,leaf,dat_int1)

!         write(*,'(20I4)')dat_int1
!         write(*,*)
        
!     Case(2)

!         Allocate(dat_int2(leaf%dat_no),stat=alloc_stat)
!         If (alloc_stat /= 0) Then
!             WRITE(*,*)
!             WRITE(*,PDF_SEP)
!             WRITE(*,PDF_E_A)   'Allocation of the dat_int2 failed !!'
!             WRITE(*,PDF_E_AI0) 'The requested dimension was : ',leaf%dat_no
!             WRITE(*,PDF_E_STOP)
!             STOP
!         End If

!         CALL pd_read_leaf(tree%streams,leaf,dat_int2)

!         write(*,'(10I8)')dat_int2
!         write(*,*)
        
!     Case(3)

!         Allocate(dat_int4(leaf%dat_no),stat=alloc_stat)
!         If (alloc_stat /= 0) Then
!             WRITE(*,*)
!             WRITE(*,PDF_SEP)
!             WRITE(*,PDF_E_A)   'Allocation of the dat_int4 failed !!'
!             WRITE(*,PDF_E_AI0) 'The requested dimension was : ',leaf%dat_no
!             WRITE(*,PDF_E_STOP)
!             STOP
!         End If

!         CALL pd_read_leaf(tree%streams,leaf,dat_int4)

!         write(*,'(5I16)')dat_int4 
!         write(*,*)
        
!     Case(4)

!         Allocate(dat_int8(leaf%dat_no),stat=alloc_stat)
!         If (alloc_stat /= 0) Then
!             WRITE(*,*)
!             WRITE(*,PDF_SEP)
!             WRITE(*,PDF_E_A)   'Allocation of the dat_int8 failed !!'
!             WRITE(*,PDF_E_AI0) 'The requested dimension was : ',leaf%dat_no
!             WRITE(*,PDF_E_STOP)
!             STOP
!         End If

!         CALL pd_read_leaf(tree%streams,leaf,dat_int8)

!         write(*,'(5I16)')dat_int8
!         write(*,*)
        
!     Case(5)

!         Allocate(dat_real8(leaf%dat_no),stat=alloc_stat)
!         If (alloc_stat /= 0) Then
!             WRITE(*,*)
!             WRITE(*,PDF_SEP)
!             WRITE(*,PDF_E_A)   'Allocation of the dat_real8 failed !!'
!             WRITE(*,PDF_E_AI0) 'The requested dimension was : ',leaf%dat_no
!             WRITE(*,PDF_E_STOP)
!             STOP
!         End If

!         CALL pd_read_leaf(tree%streams,leaf,dat_real8)

!         write(*,'(4E18.9)')dat_real8
!         write(*,*)

!     Case(6)

!         Allocate(dat_char(leaf%dat_no),stat=alloc_stat)
!         If (alloc_stat /= 0) Then
!             WRITE(*,*)
!             WRITE(*,PDF_SEP)
!             WRITE(*,PDF_E_A)   'Allocation of the dat_cahr failed !!'
!             WRITE(*,PDF_E_AI0) 'The requested dimension was : ',leaf%dat_no
!             WRITE(*,PDF_E_STOP)
!             STOP
!         End If

!         CALL pd_read_leaf(tree%streams,leaf,dat_char)

!         write(*,'(*(A))')dat_char
!         write(*,*)
!     Case default
!         Write(*,PDF_E_AI0)"Unknown data type in leaf :",leaf%dat_ty
!     End Select

! End If

! CALL close_stream_files(tree)


! CALL Cpu_time(gend_time)
! Write(*,PDF_M_A )'=='
! Write(*,PDF_TIME)"User time :",gend_time-gstart_time
