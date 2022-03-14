!******************************************************************************
!>  Program for ...                                                          **
!>                                                                           **
!** TODO: Refactoring to vector and matrix output.                           **
!**                                                                          **
!** ------------------------------------------------------------------------ **
!>  \section written Written by:
!>  Ralf Schneider
!>
!>  \section modified Last modified:
!>  by: Ralf Schneider \n
!>  on : 15.02.2012
!** ------------------------------------------------------------------------ **
Program pd_leaf_to_file 

  Use puredat         ! From libpuredat
  Use vtkio

  Implicit None

  !****************************************************************************
  !** Declarations ************************************************************

  !-- Chain Variables ---------------------------------------------------------
  Real(Kind=pd_rk)   :: gstart_time, gend_time

  Integer            :: num_args, un_st, un_out, un_log
  Character(len=256) :: pchars

  Character(len=pd_mcl), Allocatable,Dimension(:) :: leaf_desc
  Character(len=4*pd_mcl) :: long_arg  
  Character(len=pd_mcl)   :: branch_desc, arg
  Character(len=pd_mcl)   :: log_file_name, tmp_char
  Character(len=pd_mcl)   :: outfile, fmt_str="", data_format
  Character(len=pd_mcl)   :: vtk_location = "POINT_DATA"
  Character(len=12)       :: vtk_desc
  Character               :: lls
  
  Integer(Kind=pd_ik) :: tmp_leaf_num, leaf_num, leaf_desc_num
  Integer(Kind=pd_ik) :: ii, jj, kk
  Integer(Kind=pd_ik) :: char_pos, cl_pos, cl_tok, stride, char_len

  Type(tLeaf), Dimension(:), Allocatable :: leaf_list
  Type(tBranch), Pointer                 :: branch
  Type(tBranch)                          :: tree

  Logical :: log_data=.False., log_file=.False.
  Logical :: logtree =.False.
  Logical :: success, log_leaf_names=.FALSE.
  Logical :: silent=.False., vtk_data_head=.True.

  Integer(Kind=pd_ik), Dimension(3) :: extend
  Real(kind=pd_rk)   , Dimension(3) :: origin, spacing

  Type(tStreams) :: dat

  !** Init problematic chars in filenames ***
  Call init_pchars()

  !== Code ====================================================================

  Call Cpu_Time(gstart_time)

  !============================================================================

  num_args = command_argument_Count()

  !** We need at least -pro_name as argument ***
  If (num_args < 1) Then
     call print_help()
     Stop
  End If

  !** Set default arguments *******************************
  pro_name    = ""

  pro_path    = "./"
  branch_desc = "_MAIN_"
  outfile     = ""

  lls         = ","
  
  log_data    = .False.
  fmt_str     = ""

  log_file      = .False.
  log_file_name = ""
  data_format   = "raw"
  extend        = 0
  spacing       = 1._pd_rk
  origin        = 0._pd_rk

  silent         = .False.
  log_leaf_names = .FALSE.
  
  !** Parse arguments *************************************
  Do ii = 1, num_args

     arg=""
     Call get_command_Argument(ii, arg)

     char_pos=Scan(arg,"=")

     If (char_pos == 0) char_pos = len_Trim(arg)+1

     Select Case (arg(1:char_pos-1))

        !** Print help ************************************
     Case ("-help")
        Call print_help()
        stop
     Case ("--help")
        Call print_help()
        stop
     Case ("-h")
        Call print_help()
        stop
     Case ("--h")
        Call print_help()
        stop

        !** Project name and path *************************
     Case ("-pro_path")
        pro_path=Trim(arg(char_pos+1:len_Trim(arg)))

     Case ("-pro_name")
        pro_name=Trim(arg(char_pos+1:len_Trim(arg)))

        !** Leaf list seperator ***************************
     Case ("-sep_leaves")
        lls=Trim(arg(char_pos+1:len_Trim(arg)))
        
        !** Tree informations ***************************************
     Case ("-leaf")

        !** Get Argument string ***************************
        long_arg = Trim(arg(char_pos+1:len_Trim(arg)))

        !** Get length of argument **************
        char_len = len_trim(long_arg)

        !** Cuurent location of sep token *******
        cl_tok = Scan(long_arg(1:char_len),lls)
        
        !** No sep token so only one leaf description *****
        If ( cl_tok == 0 ) then

           leaf_desc_num = 1
           Allocate(leaf_desc(1))
           leaf_desc = trim(long_arg)

        !** More than one leaf description given **********
        Else

           leaf_desc_num = 1
           cl_pos = cl_tok + 1

           !** Count sep tokens ***************************
           do while ((cl_pos < char_len) .AND. (cl_tok /= 0))
              
              cl_tok = Scan(long_arg(cl_pos:char_len),",")
              cl_pos = cl_pos + cl_tok + 1

              leaf_desc_num = leaf_desc_num + 1
             
           End do

           !** Allocate filed for leaf descriptions *******
           Allocate(leaf_desc(leaf_desc_num))

           !** Get leaf descriptions from Argument ********
           cl_tok = Scan(long_arg(1:char_len),lls)

           leaf_desc_num = 1
           leaf_desc(leaf_desc_num) = long_arg(1:cl_tok-1)

           cl_pos = cl_tok + 1
           
           do while (cl_pos < char_len)

              cl_tok = Scan(long_arg(cl_pos:char_len),",")
              If (cl_tok == 0) then
                 cl_tok = len_trim(long_arg(cl_pos:char_len))+1
              End If

              leaf_desc_num = leaf_desc_num + 1
              leaf_desc(leaf_desc_num) = long_arg(cl_pos:cl_pos + cl_tok - 2)              
              cl_pos = cl_pos + cl_tok
              
           End do
           
        End If

     Case ("-branch")
        branch_desc=Trim(arg(char_pos+1:len_Trim(arg)))

        !** Where and what to output **********************
     Case ("-o")
        outfile=Trim(arg(char_pos+1:len_Trim(arg)))
        log_file = .True.

     Case ("-data_format")
        data_format=Trim(arg(char_pos+1:len_Trim(arg)))
        log_file = .True.

     Case ("-l")
        log_data = .True.

        If (len_Trim(arg) > 2) then
           log_file_name = Trim(arg(char_pos+1:len_Trim(arg)))
           un_log = pd_give_new_unit()
        Else
           un_log = std_out
        End If

     Case ("-lt")
        logtree = .True.

     Case ("-ll")
        log_leaf_names = .True.

     Case ("-ascii_format")
        fmt_str=Trim(arg(char_pos+1:len_Trim(arg)))

     Case ("-silent")
        silent=.True.

        !** VTK options ***********************************
     Case ("-extend")
        Read(arg(char_pos+1:len_Trim(arg)),*)extend

     Case ("-spacing")
        Read(arg(char_pos+1:len_Trim(arg)),*)spacing

     Case ("-origin")
        Read(arg(char_pos+1:len_Trim(arg)),*)origin

     Case default

        Write(*,'(80("!"))')
        Write(*,'(A,T78,"!!!")')"Found unknown argument"
        Write(*,'(A)')Trim(arg)
        Write(*,'(80("!"))')

     End Select

  End Do

  !** Set default leaf description ********************************************
  If (.NOT. Allocated(leaf_desc)) then
     leaf_desc_num = 1
     Allocate(leaf_desc(1))
     leaf_desc   = "_ALL_"
  End If
  
  !** We need a / at the end of pro_path **************************************
  If (pro_path(len_trim(pro_path):len_trim(pro_path)) /= "/") pro_path = trim(pro_path)//"/"

  !** Set default base output filename and path *******************************
  If (outfile == "") Then
     outfile = Trim(pro_path)//Trim(pro_name)//"."//Trim(data_format)
  End If

  !** Check whether there are only valid chars in the output filename *********
  outfile = check_valid_chars(outfile)

  !** Repeat input parameters to stdout ***************************************
  If (.Not. silent) Then

     Write(*,*)"================="
     Write(*,*)"Pro path            : ",Trim(pro_path)
     Write(*,*)"Pro name            : ",Trim(pro_name)
     Write(*,*)"-----------------"
     Write(*,*)"Branch desc         : ",Trim(branch_Desc)
     Write(*,*)
     Write(*,*)"Leaf descriptions   : ",Trim(leaf_desc(1))
     Do ii = 2, leaf_desc_num
        Write(*,*)"                    : ",Trim(leaf_desc(ii))
     End Do
     
     If (log_file) then

        Write(*,*)"-----------------"
        Write(*,*)"Out file            : ",Trim(outfile)
        Write(*,*)"Out format          : ",Trim(data_format)

        if (trim(data_format) == "vtk") then

           Write(*,*)"vtk extend          : ",extend
           Write(*,*)"vtk spacing         : ",spacing
           Write(*,*)"vtk origin          : ",origin

           !** Check if we have a valid extend ********************************
           If ( (extend(1)==0) .or. (extend(2)==0) .or. (extend(3)==0) ) then
              Write(*,*)"For VTK output a valid extend with all     !!!"
              Write(*,*)"three components not equal to 0 is needed  !!!"
              Write(*,*)"Program stopped !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
              stop
           End If
              
        End if

     End If

     If ( log_data ) then

        Write(*,*)"-----------------"
        Write(*,*)"Dump ascii data     : ",log_data
        Write(*,*)"Dump format         : ",Trim(fmt_str)
        If (trim(log_file_name) /= "" ) then
           Write(*,*)"Dump file name      : ",Trim(log_file_name)
        End If

     End If

     If ( logtree ) then

        Write(*,*)"-----------------"
        Write(*,*)"Dump tree to stdout : ",logtree

     End If

     Write(*,*)"-----------------"

  End If

  !** Read header *****************************************
  tree = read_tree()

  !** If the tree structure should be dumped to std out ***
  If (logtree) Call log_tree(tree,std_out)

  !** Assign branch from which to dump the leafs **********
  If (Trim(branch_desc) == "_MAIN_") Then

     Allocate(branch)
     Call assign_branches(branch,tree)

  Else

     Call search_branch(Trim(branch_desc),tree,branch,success)
     If (.Not. success) Then
        Write(*,*)"Couldn't find the branch ",Trim(branch_desc)
        Write(*,*)"given as argument 4      "
        Write(*,*)"in tree                  ",Trim(pro_name)
        Write(*,*)"with pro_path            ",Trim(pro_path)
        Write(*,*)"Program halted !!!!!!!!!!!"
        Stop
     End If

  End If

  leaf_dec_loop : Do kk = 1, leaf_desc_num
     
     !** Prepare output file *********************************
     If (log_file .AND. (Trim(data_format)=="vtk") ) Then
         
        !** If there is more than one leaf to dump we have to extend
        !** the vtk outfile name with the leaf descriptions
        char_pos=Scan(outfile,".",.TRUE.)
        If (char_pos == 0) char_pos = len_trim(outfile)
        
        Write(tmp_char,'(A,"-",A,A)')&
             outfile(1:char_pos-1),Trim(check_valid_chars(leaf_list(ii)%desc)),&
             outfile(char_pos:len_trim(outfile))

        Call write_vtk_structured_points(Trim(tmp_char), &
             extend, spacing, origin)
           
     End If
        
     If ( trim(leaf_desc(kk)) /= "_ALL_" ) then
        
        leaf_num = 0
        Do ii = 1, leaf_desc_num
           !** Look for number of leafs to dump **************
           tmp_leaf_num=0
           Call get_leaf_num(Trim(leaf_desc(1)), branch, tmp_leaf_num)
           leaf_num = leaf_num + tmp_leaf_num
        End Do
        
     Else
        
        leaf_num = branch%no_leaves
        
     End If

     If (.Not. silent) Then
        Write(*,*)"Leaf num            : ",leaf_num
        If ( (log_leaf_names) .AND. (trim(leaf_desc(1)) == "_ALL_") ) then
           Do ii = 1, leaf_num
              Write(*,"(I15,' - ',A)")ii,trim(branch%leaves(ii)%desc)
           End Do
        End If
        Write(*,'(60("="))')
     End If
  
     !** Leaf output *****************************************
     If (log_file .Or. log_data) Then

        !** Do we have something to write ? ******************
        number_of_leafs : If (leaf_num > 0) Then

           !** Get the leafs to write ************************
           If ( trim(leaf_desc(kk)) /= "_ALL_" ) then
              
              Call get_leaf_list(Trim(leaf_desc(kk)),branch,&
                   leaf_num,leaf_list)
           Else

              Allocate(leaf_list(leaf_num))
              leaf_list = branch%leaves

           End If

           !*******************************************************************
           !** Loop over all leafs in leaf_list *******************************
           Leaf_Loop : Do ii = 1, leaf_num

              !** Prepare ASCII log file *******************************
              If ( trim(log_file_name) /= "" ) then

                 !** If there is more than one leaf to dump we have to extend
                 !** the log_file_name with the leaf descriptions and numbers
                 char_pos=Scan(log_file_name,".",.TRUE.)
                 If (char_pos == 0) char_pos = len_trim(log_file_name)
                 
                 Write(tmp_char,'(A,"-",A,"-",I3.3,A)')&
                      log_file_name(1:char_pos-1),Trim(check_valid_chars(leaf_list(ii)%desc)),&
                      ii,log_file_name(char_pos:len_trim(log_file_name))
                 
                 Open(unit=un_log, file=Trim(tmp_char), &
                      status='replace', action='write')

              !** Prepare raw log file ******************************
              Else If ( (Trim(data_format) == "raw") .AND. log_file ) Then

                 !** If there is more than one leaf to dump we have to extend
                 !** the log_file_name with the leaf descriptions and numbers
                 char_pos=Scan(outfile,".",.TRUE.)
                 If (char_pos == 0) char_pos = len_trim(outfile)
                 
                 Write(tmp_char,'(A,"-",A,"-",I3.3,A)')&
                      outfile(1:char_pos-1),Trim(check_valid_chars(leaf_list(ii)%desc)),&
                      ii,outfile(char_pos:len_trim(outfile))

                 un_out=pd_give_new_unit()
                 Open(unit=un_out, file=Trim(tmp_char), &
                      status='replace', action='write', &
                      access='stream')

              End If
           
              If (.Not. silent) Write(*,*)"Leaf num    = ",ii

              !** If there are no data in the leaf cycle ***************
              If (leaf_list(ii)%dat_no <= 0) Then
                 Write(*,*)"No data in leaf No. ",ii
                 Cycle
              End If

              If  ( Trim(data_format)=="vtk" ) Then
              
                 !** Check whether we have an even distribution of leaf elements  ***
                 !** among the given vtk extend in case of multidimensional leaf  ***
                 !** elements                                                     ***
                 If ( mod(leaf_list(ii)%dat_no,(extend(1)*extend(2)*extend(3))) /= 0 ) then
                    
                    Write(*,*)"Found ",leaf_list(ii)%dat_no," elements in leaf"
                    write(*,*)trim(leaf_list(ii)%desc)
                    write(*,*)"Which leads to non zero mod against the product of the given vtk extend"
                    write(*,*)"extend(1)*extend(2)*extend(3) = ",extend(1)*extend(2)*extend(3)
                    write(*,*)"mod(leaf_list(ii)%dat_no,(extend(1)*extend(2)*extend(3))) = ",&
                         mod(leaf_list(ii)%dat_no,(extend(1)*extend(2)*extend(3)))
                    Cycle
                 
                 End If

              End If
           
              !** Open the stream files corresponding to the current leaf ********
              un_st = pd_give_new_unit()
              Open(unit=un_st, file=tree%streams%stream_files(leaf_list(ii)%dat_ty), &
                   status='old', action='read', access='stream')
              
              Select Case (leaf_list(ii)%dat_ty)
                 
              Case(1) !** 1Byte Integer Data ************************

                 Allocate(dat%int1_st(leaf_list(ii)%dat_no))
                 Read(un_st,pos=(leaf_list(ii)%lbound-1)*1+1)dat%int1_st
                 
                 !*************************************************************
                 !** Loop over the scalar leaf elements.    *******************
                 stride = leaf_list(ii)%dat_no / (extend(1)*extend(2)*extend(3))
                 If (.Not. silent) Write(*,*)"Leaf stride = ",stride
                 Leaf_Element_Loop_Int1 : Do jj = 1, stride

                    !** Binary output requested *****************
                    If (log_file) Then
                       
                       If (data_format == "raw") Then
                          Write(un_out)dat%int1_st
                          
                       Else If (data_format=="vtk") Then
                          
                          Write(vtk_desc,'(A4,2("-",I3.3))')&
                               Trim(check_valid_chars(leaf_list(ii)%desc)),&
                               ii,jj
                          
                          Call write_vtk_data_int4_scalar_1D (&
                               matrix=Int(dat%int1_st(jj:leaf_list(ii)%dat_no:stride),4), &
                               filename=Trim(tmp_char), desc=Trim(check_valid_chars(leaf_list(ii)%desc)), &
                               head=vtk_data_head, location=vtk_location)
                          vtk_data_head = .False.
                          
                       End If
                    End If

                    If (log_data) Then
                       If ( .Not. silent )Write(*,*)"=="
                       If (fmt_str=="") Then
                          Write(un_log,*)dat%int1_st
                       Else
                          Write(un_log,fmt_str)dat%int1_st
                       End If
                    End If
                    
                 end Do Leaf_Element_Loop_Int1
                 Deallocate(dat%int1_st)

              Case(2) !** 2Byte Integer Data ************************
                 
                 Allocate(dat%int2_st(leaf_list(ii)%dat_no))
                 Read(un_st,pos=(leaf_list(ii)%lbound-1)*2+1)dat%int2_st

                 stride = leaf_list(ii)%dat_no / (extend(1)*extend(2)*extend(3))
                 If (.Not. silent) Write(*,*)"Leaf stride = ",stride
                 Leaf_Element_Loop_Int2 : Do jj = 1, stride

                    If (log_file) Then
                       If (data_format=="raw") Then
                          Write(un_out)dat%int2_st

                       Else If (data_format=="vtk") Then

                          Write(vtk_desc,'(A4,2("-",I3.3))')&
                               Trim(check_valid_chars(leaf_list(ii)%desc)),&
                               ii,jj

                          Call write_vtk_data_int4_scalar_1D (&
                               matrix=Int(dat%int2_st(jj:leaf_list(ii)%dat_no:stride),4), &
                               filename=Trim(tmp_char), desc=Trim(check_valid_chars(leaf_list(ii)%desc)), &
                               head=vtk_data_head, location=vtk_location)
                          vtk_data_head = .False.

                       End If
                    End If

                    If (log_data) Then
                       If (.Not. silent)Write(*,*)"=="
                       If (fmt_str=="") Then
                          Write(un_log,*)dat%int2_st
                       Else
                          Write(un_log,fmt_str)dat%int2_st
                       End If
                    End If

                 End Do Leaf_Element_Loop_Int2
                 Deallocate(dat%int2_st)

              Case(3) !** 4Byte Integer Data ************************

                 Allocate(dat%int4_st(leaf_list(ii)%dat_no))
                 Read(un_st,pos=(leaf_list(ii)%lbound-1)*4+1)dat%int4_st

                 !*******************************************************************
                 !** Loop over the scalar leaf elements.    *************************
                 stride = leaf_list(ii)%dat_no / (extend(1)*extend(2)*extend(3))
                 If (.Not. silent) Write(*,*)"Leaf stride = ",stride
                 Leaf_Element_Loop_Int4 : Do jj = 1, stride

                    If (log_file) Then
                       If (data_format=="raw") Then
                          Write(un_out)dat%int4_st

                       Else If (data_format=="vtk") Then

                          Write(vtk_desc,'(A4,2("-",I3.3))')&
                               Trim(check_valid_chars(leaf_list(ii)%desc)),&
                               ii,jj

                          Call write_vtk_data_int4_scalar_1D (&
                               matrix=dat%int4_st(jj:leaf_list(ii)%dat_no:stride), &
                               filename=Trim(tmp_char), desc=Trim(check_valid_chars(leaf_list(ii)%desc)), &
                               head=vtk_data_head, location=vtk_location)
                          vtk_data_head = .False.

                       End If
                    End If

                    If (log_data) Then
                       If (.Not. silent)Write(*,*)"=="
                       If (fmt_str=="") Then
                          Write(un_log,*)dat%int4_st
                       Else
                          Write(un_log,fmt_str)dat%int4_st
                       End If
                    End If

                 End Do Leaf_Element_Loop_Int4
                 Deallocate(dat%int4_st)

              Case(4) !** 8Byte Integer Data ************************

                 Allocate(dat%int8_st(leaf_list(ii)%dat_no))
                 Read(un_st,pos=(leaf_list(ii)%lbound-1)*8+1)dat%int8_st

                 If (log_file) Then

                    If (data_format=="raw") Then
                       Write(un_out)dat%int8_st

                    Else If (data_format=="vtk") Then

                       !*******************************************************************
                       !** Loop over the scalar leaf elements.    *************************
                       stride = leaf_list(ii)%dat_no / (extend(1)*extend(2)*extend(3))
                       If (.Not. silent) Write(*,*)"Leaf stride = ",stride
                       Leaf_Element_Loop_Int8 : Do jj = 1, stride

                          Write(vtk_desc,'(A4,2("-",I3.3))')&
                               Trim(check_valid_chars(leaf_list(ii)%desc)),&
                               ii,jj

                          Call write_vtk_data_int8_scalar_1D (&
                               matrix=dat%int8_st(jj:leaf_list(ii)%dat_no:stride), &
                               filename=Trim(tmp_char), desc=vtk_desc, &
                               head=vtk_data_head, location=vtk_location)
                          vtk_data_head = .False.

                       End Do Leaf_Element_Loop_Int8

                    End If

                 End If

                 If (log_data) Then
                    If (.Not. silent)Write(*,*)"=="
                    If (fmt_str=="") Then
                       Write(un_log,*)dat%int8_st
                    Else
                       Write(un_log,fmt_str)dat%int8_st
                    End If
                 End If

                 Deallocate(dat%int8_st)

              Case(5) !** 8Byte Real Data ***************************

                 Allocate(dat%real8_st(leaf_list(ii)%dat_no))
                 Read(un_st,pos=(leaf_list(ii)%lbound-1)*8+1)dat%real8_st

                 If (log_file) Then

                    If (data_format=="raw") Then
                       Write(un_out)dat%real8_st

                    Else If (data_format=="vtk") Then

                       !*******************************************************************
                       !** Loop over the scalar leaf elements.    *************************
                       stride = leaf_list(ii)%dat_no / (extend(1)*extend(2)*extend(3))
                       If (.Not. silent) Write(*,*)"Leaf stride = ",stride
                       Leaf_Element_Loop_Real8 : Do jj = 1, stride

                          Write(vtk_desc,'(A4,2("-",I3.3))')&
                               Trim(check_valid_chars(leaf_list(ii)%desc)),&
                               ii,jj

                          Call write_vtk_data_real8_scalar_1D (&
                               matrix=dat%real8_st(jj:leaf_list(ii)%dat_no:stride), &
                               filename=Trim(tmp_char), desc=vtk_desc, &
                               head=vtk_data_head, location=vtk_location)
                          vtk_data_head = .False.

                       End Do Leaf_Element_Loop_Real8

                    End If

                 End If

                 If (log_data) Then
                    If (.Not. silent)Write(*,*)"=="
                    If (fmt_str=="") Then
                       Write(un_log,*)dat%real8_st
                    Else
                       Write(un_log,fmt_str)dat%real8_st
                    End If
                 End If

                 Deallocate(dat%real8_st)

              Case(6) !** Character Data **************************************

                 Allocate(dat%char_st(leaf_list(ii)%dat_no))
                 Read(un_st,pos=leaf_list(ii)%lbound)dat%char_st

                 If (log_file) Then
                    If (data_format=="raw") Then
                       Write(un_out)dat%char_st

                    Else If (data_format=="vtk") Then

                       Write(*,*)"Sorry character output in VTK format is not implemented"

                    End If
                 End If

                 If (log_data) Then
                    If (.Not. silent)Write(*,*)"=="
                    If (fmt_str=="") Then
                       Write(un_log,*)dat%char_st
                    Else
                       Write(un_log,fmt_str)dat%char_st
                    End If
                 End If

                 Deallocate(dat%char_st)

              Case default !** Default ******************************
                 Write(*,'(80("="))')
                 Write(*,'(A)')" == Something bad and unexpected happened !!!"
                 Write(*,'(A,I10,A)')" == Data type number ",leaf_list(ii)%dat_ty, &
                      " is not supported !!!"
                 Write(*,'(80("="))')
                 Stop
              End Select

              !** Close stream file *****************************
              Close(un_st)

              !** Close output file *****************************
              If (log_file  .And. Trim(data_format) == "raw" ) Then
                 Close(un_out)
              End If
              
           End Do Leaf_Loop

        Else

           Write(*,*)"Couldn't find the leaf ",Trim(leaf_desc(kk))
           Write(*,*)"given as argument 3    "
           Write(*,*)"in tree                ",Trim(pro_name)
           Write(*,*)"with pro_path          ",Trim(pro_path)

        End If number_of_leafs

     End If

     !** Close vtk output file ****************************
     If (log_file  .And. Trim(data_format) == "vtk" ) Then
        Close(un_out)
     End If

  ENd Do leaf_dec_loop
  !============================================================================

  Call Cpu_Time(gend_time)
  If (.Not. silent) Then
     Write(*,*)'=='
     Write(*,'(A,F9.3)')" Used CPU time [sec] : ",gend_time-gstart_time
  End If

contains

  !============================================================================
  Subroutine print_help()

    Write(*,'(80("="))')
    !              0        1         2         3         4         5         6
    !              123456789012345678901234567890123456789012345678901234567890
    Write(*,'(A)')"== Usage:"
    Write(*,'(A)')"== "
    Write(*,'(A)')"== pd_leaf_to_file -pro_name=pro_name [-leaf=leaf_desc]"
    Write(*,'(A)')"==                 [-pro_path=pro_path] [-branch=branch_desc]"
    Write(*,'(A)')"==                 [-o=filename] [-format=format_string]"
    Write(*,'(A)')"==                 [-l] [-silent]"
    Write(*,'(A)')"== "
    Write(*,'(A)')"== -ascii_format=format_string"
    Write(*,'(A)')"==      Format string for ascii dump to stdout."
    Write(*,'(A)')"== "
    Write(*,'(A)')"== -branch=branch_desc"
    Write(*,'(A)')"==      Description of branch to be searchred for leaf_desc."
    Write(*,'(A)')"==      Default is _MAIN_ which searches the main_branch of"
    Write(*,'(A)')"==      the tree given by pro_name."
    Write(*,'(A)')"== "
    Write(*,'(A)')"== -data_format=format"
    Write(*,'(A)')"==      Output format. Possible values are: raw and vtk."
    Write(*,'(A)')"==      Default is raw."
    Write(*,'(A)')"== "
    Write(*,'(A)')"== -extend=nx,ny,nz"
    Write(*,'(A)')"==      Extend of a regular grid in R3 to which the leaf"
    Write(*,'(A)')"==      data should be mapped when vtk output is performed."
    Write(*,'(A)')"==      The argument is mandatory for -data_format=vtk."
    Write(*,'(A)')"== "
    Write(*,'(A)')"== -help, --help, -h, --h"
    Write(*,'(A)')"==      Print a summary of the command-line usage and exit."
    Write(*,'(A)')"==      output."
    Write(*,'(A)')"== "
    Write(*,'(A)')"== -l   Whether data should be monitored in ascii format "
    Write(*,'(A)')"==      to stdout. Default is no monitoring."
    Write(*,'(A)')"== "    
    Write(*,'(A)')"== -ll  In case -leaf=_ALL_ is specified, whether the found"
    Write(*,'(A)')"==      leaf names should be dumped to stdout. Default is "
    Write(*,'(A)')"==      no dump of the leaf names."
    Write(*,'(A)')"== "
    Write(*,'(A)')"== -lt  Whether the tree structure should be dumped to"
    Write(*,'(A)')"==      stdout. Default is no dump of the tree structure."
    Write(*,'(A)')"== "
    Write(*,'(A)')"== -leaf=leaf_desc"
    Write(*,'(A)')"==      Description of leaf(s) to be searchred for. The search"
    Write(*,'(A)')"==      is performed recursively so returning all leafs in the"
    Write(*,'(A)')"==      Branch structure given by -branch with the same"
    Write(*,'(A)')"==      description given by leaf_desc. If -leaf is not given,"
    Write(*,'(A)')"==      the default is _ALL_ which returns all leafs of the "
    Write(*,'(A)')"==      branch given by parameter -branch. This operation is"
    Write(*,'(A)')"==      performed non-recursively."
    Write(*,'(A)')"== "
    Write(*,'(A)')"== -o=filename"
    Write(*,'(A)')"==      Path to output file with binary leaf data."
    Write(*,'(A)')"==      Default is pro_path/pro_name-leaf_desc.data_format"
    Write(*,'(A)')"== "
    Write(*,'(A)')"== -origin=Ox,Oy,Oz"
    Write(*,'(A)')"==      Origin of a regular grid in R3 to which the leaf"
    Write(*,'(A)')"==      data should be mapped when vtk output is performed."
    Write(*,'(A)')"==      Default is Ox = Oy = Oz = 0.0."
    Write(*,'(A)')"== "
    Write(*,'(A)')"== -pro_name=pro_name"
    Write(*,'(A)')"==      Puredat project name of tree to be searched for"
    Write(*,'(A)')"==      leaf_desc."
    Write(*,'(A)')"== "
    Write(*,'(A)')"== -pro_path=pro_path"
    Write(*,'(A)')"==      Puredat project path of tree to be searched for"
    Write(*,'(A)')"==      leaf_desc. Default is ./ ."
    Write(*,'(A)')"== "
    Write(*,'(A)')"== -silent"
    Write(*,'(A)')"==      Whether to do silent output. Default is verbose"
    Write(*,'(A)')"==      output."
    Write(*,'(A)')"== "
    Write(*,'(A)')"== -spacing=dx,dy,dz"
    Write(*,'(A)')"==      Spacing of a regular grid in R3 to which the leaf"
    Write(*,'(A)')"==      data should be mapped when vtk output is performed."
    Write(*,'(A)')"==      Default is dx = dy = dz = 1.0."
    Write(*,'(A)')"== "
    Write(*,'(80("="))')

  End Subroutine print_help

  !============================================================================
  Subroutine init_pchars()

    pchars = ""

    Do ii = 0, 31
       pchars(ii+1:ii+1)=achar(ii)
    End Do

    Do ii = 127, 256
       pchars(ii-94:ii-94)=achar(ii)
    End Do

    pchars=trim(pchars)//' <>|*?[]{}()$`´"!;\#^%&'//"'"

  End Subroutine init_pchars

  !============================================================================
  !** Check whether there are only valid chars in the output filename *********
  Function check_valid_chars(in_string) Result(string)

    Character(Len=*), Intent(in)   :: in_string
    Integer(Kind=8)                :: char_pos
    Character(Len=len(in_string))  :: string

    string = in_string

    char_pos = scan(trim(string),pchars)

    If (  char_pos > 0 ) then

       Do while ( char_pos > 0 )
          string(char_pos:char_pos) = "_"
          char_pos = scan(trim(string),pchars)
       End Do

    End If

  End Function check_valid_chars

End Program pd_leaf_to_file
