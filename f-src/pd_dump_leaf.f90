!******************************************************************************
!>  Program for                                                              **
!**                                                                          **
!** ------------------------------------------------------------------------ **
!>  \section written Written by:
!>  Ralf Schneider
!>
!>  \section modified Last modified:
!>  by: Ralf Schneider \n
!>  on : 15.02.2012
!** ------------------------------------------------------------------------ **
Program pd_dump_leaf

  Use puredat         
  
  Implicit None

  !****************************************************************************
  !** Declarations ************************************************************

  !-- Chain Variables ---------------------------------------------------------
  Real(Kind=pd_rk)            :: gstart_time, gend_time

  Integer                     :: num_args, leaf_no, un, alloc_stat

  Character(len=pd_mcl)       :: outpath, outfile, arg

  Type(tBranch)               :: tree
  Type(tLeaf),   Pointer      :: leaf
  Logical                     :: success
  
  Integer(Kind=1)      , Dimension(:), Allocatable   :: dat_int1
  Integer(Kind=2)      , Dimension(:), Allocatable   :: dat_int2
  Integer(Kind=4)      , Dimension(:), Allocatable   :: dat_int4
  Integer(Kind=8)      , Dimension(:), Allocatable   :: dat_int8
  Real(Kind=8)         , Dimension(:), Allocatable   :: dat_real8
  Character            , Dimension(:), Allocatable   :: dat_char
  
  !== Code ====================================================================

  Call Cpu_time(gstart_time)

  !============================================================================

  num_args = command_argument_count()

  If (num_args < 3) then
     Write(*,'(80("="))')
     Write(*,'(A)')"== Usage:"
     Write(*,'(A)')"== arg 1   : Puredat project path"
     Write(*,'(A)')"== arg 2   : Puredat project name"
     Write(*,'(A)')"== arg 3   : Number of leaf to be dumped"
     Write(*,'(80("="))')
     Stop
  End If

  call get_command_argument(1, pro_path)
  call get_command_argument(2, pro_name)
  call get_command_argument(3, arg)
  Read(arg,*)leaf_no

  write(*,PDF_SEP  )
  write(*,PDF_M_A  ) "Pro path    : "//trim(pro_path)
  write(*,PDF_M_A  ) "Pro name    : "//trim(pro_name)
  write(*,PDF_M_AI0) "Leaf no.    : ",leaf_no
  write(*,PDF_M_A  )"=="

  tree = read_tree()
  call open_stream_files(tree, "read" , "old")
  call get_leaf_with_num(tree, leaf_no, leaf, success)

  If (.not. Success) then
     Write(*,PDF_E_AI0)"There is no leaf with number ",leaf_no
  ELSE

     write(*,PDF_M_A  ) "Description           :"//trim(leaf%desc)
     write(*,PDF_M_AI0) "Number of data        :",leaf%dat_no
     write(*,PDF_M_AI0) "Type of data          :",leaf%dat_ty
     write(*,PDF_M_AI0) "Lower bound in stream :",leaf%lbound
     write(*,PDF_M_AI0) "Upper bound in stream :",leaf%ubound
     write(*,PDF_M_A  )"=="

     Select case (leaf%dat_ty)
     
        Case(1)

           Allocate(dat_int1(leaf%dat_no),stat=alloc_stat)
           If (alloc_stat /= 0) Then
              WRITE(*,*)
              WRITE(*,PDF_SEP)
              WRITE(*,PDF_E_A)   'Allocation of the dat_int1 faild !!'
              WRITE(*,PDF_E_AI0) 'The requested dimension was : ',leaf%dat_no
              WRITE(*,PDF_E_STOP)
              STOP
           End If

           call pd_read_leaf(tree%streams,leaf,dat_int1)

           write(*,'(20I4)')dat_int1
           write(*,*)
           
        Case(2)

           Allocate(dat_int2(leaf%dat_no),stat=alloc_stat)
           If (alloc_stat /= 0) Then
              WRITE(*,*)
              WRITE(*,PDF_SEP)
              WRITE(*,PDF_E_A)   'Allocation of the dat_int2 faild !!'
              WRITE(*,PDF_E_AI0) 'The requested dimension was : ',leaf%dat_no
              WRITE(*,PDF_E_STOP)
              STOP
           End If

           call pd_read_leaf(tree%streams,leaf,dat_int2)

           write(*,'(10I8)')dat_int2
           write(*,*)
           
        Case(3)

           Allocate(dat_int4(leaf%dat_no),stat=alloc_stat)
           If (alloc_stat /= 0) Then
              WRITE(*,*)
              WRITE(*,PDF_SEP)
              WRITE(*,PDF_E_A)   'Allocation of the dat_int4 faild !!'
              WRITE(*,PDF_E_AI0) 'The requested dimension was : ',leaf%dat_no
              WRITE(*,PDF_E_STOP)
              STOP
           End If

           call pd_read_leaf(tree%streams,leaf,dat_int4)

           write(*,'(5I16)')dat_int4 
           write(*,*)
           
        Case(4)

           Allocate(dat_int8(leaf%dat_no),stat=alloc_stat)
           If (alloc_stat /= 0) Then
              WRITE(*,*)
              WRITE(*,PDF_SEP)
              WRITE(*,PDF_E_A)   'Allocation of the dat_int8 faild !!'
              WRITE(*,PDF_E_AI0) 'The requested dimension was : ',leaf%dat_no
              WRITE(*,PDF_E_STOP)
              STOP
           End If

           call pd_read_leaf(tree%streams,leaf,dat_int8)

           write(*,'(5I16)')dat_int8
           write(*,*)
           
        Case(5)

           Allocate(dat_real8(leaf%dat_no),stat=alloc_stat)
           If (alloc_stat /= 0) Then
              WRITE(*,*)
              WRITE(*,PDF_SEP)
              WRITE(*,PDF_E_A)   'Allocation of the dat_real8 faild !!'
              WRITE(*,PDF_E_AI0) 'The requested dimension was : ',leaf%dat_no
              WRITE(*,PDF_E_STOP)
              STOP
           End If

           call pd_read_leaf(tree%streams,leaf,dat_real8)

           write(*,'(4E18.9)')dat_real8
           write(*,*)

        Case(6)

           Allocate(dat_char(leaf%dat_no),stat=alloc_stat)
           If (alloc_stat /= 0) Then
              WRITE(*,*)
              WRITE(*,PDF_SEP)
              WRITE(*,PDF_E_A)   'Allocation of the dat_cahr faild !!'
              WRITE(*,PDF_E_AI0) 'The requested dimension was : ',leaf%dat_no
              WRITE(*,PDF_E_STOP)
              STOP
           End If

           call pd_read_leaf(tree%streams,leaf,dat_char)

           write(*,'(*(A))')dat_char
           write(*,*)
        Case default
           Write(*,PDF_E_AI0)"Unknown data type in leaf :",leaf%dat_ty
        End Select

  End If
  
  call close_stream_files(tree)

  !============================================================================

  Call Cpu_time(gend_time)
  Write(*,PDF_M_A )'=='
  Write(*,PDF_TIME)"User time :",gend_time-gstart_time
 
End Program pd_dump_leaf
