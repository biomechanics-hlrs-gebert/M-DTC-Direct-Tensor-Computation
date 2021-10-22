!==============================================================================
!> \file mod_chain.f90
!> Modules of the chain process library.
!>
!> The file holds modules which are necessary to use the chain process library
!>
!> \author Ralf Schneider
!> \date 07.05.2012

!==============================================================================
!> Character constants for unified log-file output
Module chain_constants

  use kinds

  Implicit None

  Integer         , Parameter :: timer_level = 3 ! 1 ! 2
  
  ! Character constants for nice output ---------------------------------------
  Character(Len=*), Parameter :: fmt_sep    = "('<',77('='),'>')"
  Character(LEN=*), Parameter :: fmt_inpsep = "('+',79('-'))"

  Character(Len=*), Parameter :: FMT_MSG     = "('MM ',A,T77,' MM')"
  Character(Len=*), Parameter :: FMT_MSG_BS  = "('MM ',A,T68,' ... ',$)"
  Character(Len=*), Parameter :: FMT_MSG_BE  = "('done MM')"

  Character(Len=*), Parameter :: FMT_WRN     = "('WW ',A,T77,' WW')"  
  Character(Len=*), Parameter :: FMT_ERR     = "('EE ',A,T77,' EE')"
  Character(Len=*), Parameter :: FMT_ERR_AI0 = "('EE ',*(A,I0))"  
  CHARACTER(Len=*), PARAMETER :: FMT_ERR_A   = "('EE ',A)"

  Character(Len=*), Parameter :: FMT_STOP    = "('EE PROGRAM STOPPED ..... ',&
                                                &T77,' EE',/,'<',77('='),'>')"

  Character(Len=*), Parameter :: FMT_TIME = "('MM ',A,1X,F0.6,' sec')"

  !** Warning fromats *********************************************************
  CHARACTER(Len=*), PARAMETER :: FMT_WRN_A    = "('WW ',A)"
  CHARACTER(Len=*), PARAMETER :: FMT_WRN_AI0  = "('WW ',A,1X,I0)"
  CHARACTER(Len=*), PARAMETER :: FMT_WRN_AI0A = "('WW ',A,1X,I0,1X,A)"
  CHARACTER(Len=*), PARAMETER :: FMT_WRN_AF0  = "('WW ',A,1X,F0.6)"

  !** Message formats *********************************************************
  CHARACTER(Len=*), PARAMETER :: FMT_MSG_AI0  = "('MM ',*(A,1X,I0,1X))"
  CHARACTER(Len=*), PARAMETER :: FMT_MSG_AI0A = "('MM ',A,1X,I0,1X,A)"
  CHARACTER(Len=*), PARAMETER :: FMT_MSG_A8I5 = "('MM ',A,1X,8(',',I5))"
  CHARACTER(Len=*), PARAMETER :: FMT_MSG_2AI0 = "('MM ',2(A,1X,I0,1X))"
  CHARACTER(Len=*), PARAMETER :: FMT_MSG_A3I0 = "('MM ',A,3(',',I0))"

  CHARACTER(Len=*), PARAMETER :: FMT_MSG_AF0  = "('MM ',A,F0.6)"
  CHARACTER(Len=*), PARAMETER :: FMT_MSG_AF0A = "('MM ',A,F0.6,A)"
  CHARACTER(Len=*), PARAMETER :: FMT_MSG_A2F0 = "('MM ',A,2(',',F0.6))"
  CHARACTER(Len=*), PARAMETER :: FMT_MSG_A3F0 = "('MM ',A,3(',',F0.6))"

  CHARACTER(Len=*), PARAMETER :: FMT_MSG_AL  = "('MM ',A,L1)"

  CHARACTER(Len=*), PARAMETER :: FMT_MSG_A   = "('MM ',A)"

  !** Error formats ***********************************************************

  !** Seperators **************************************************************
  CHARACTER(Len=*), PARAMETER :: FMT_HY_SEP  = "(80('-'))"
  CHARACTER(Len=*), PARAMETER :: FMT_EQ_SEP  = "(80('='))"
  CHARACTER(Len=*), PARAMETER :: FMT_DBG_SEP = "('#DBG#',75('='))"

End Module chain_constants

!==============================================================================
!> Global Variables for the chain process library
Module chain_variables
 
  Use ISO_FORTRAN_ENV
  use chain_constants
  
  Implicit None
 
  ! ---------------------------------------------------------------------------
  !> Logfile unit
  Integer                     :: un_lf   = 10000
  !> Monitor file unit (default = stdout)
  Integer                     :: un_mon  = OUTPUT_UNIT
  
  Character(len=mcl)          :: outpath = "./"
  Character(len=mcl)          :: inpath  = "./"
  Character(len=mcl)          :: project_name

  !-- Variables for reading input ---------------------------------------------
  Character(Len=mcl)            :: chp_char
  Integer(Kind=4)               :: chp_int4
  Integer(Kind=4), Dimension(3) :: chp_int4_3vector
  Integer(Kind=8)               :: chp_int8
  Integer(Kind=8), Dimension(3) :: chp_int8_3vector
  Real(Kind=8)                  :: chp_real
  Real(Kind=8)   , Dimension(3) :: chp_real_3vector

  !> Amount of monitor an log data written. Recognized keywords are:
  !> PRODUCTION, DEBUG and FE_RESULTS
  Character(Len=mcl)            :: out_amount  = "PRODUCTION" ! "DEBUG" ! 
  
End Module chain_variables

!==============================================================================
!> Routines for unified link handling in a process chain
!>
!> Routines for unified handling of program startup, logging, shutdown and 
!> error handling
Module chain_routines

  Use chain_variables
  Use timer
  Use ISO_FORTRAN_ENV

  Implicit None

  !> Interface for unified parameter reading from std input or file
  Interface read_input

     Module Procedure read_char_in
     
     Module Procedure read_int4_in
     Module Procedure read_int4_3vector_in

     Module Procedure read_int8_in
     Module Procedure read_int8_3vector_in
     
     Module Procedure read_real_in
     Module Procedure read_real_3vector_in
     Module Procedure read_real_matrix_in
     
     Module Procedure read_log_in
  
  End Interface read_input
  
!!$  Interface
!!$     subroutine fsystem(command) bind (c)
!!$        use, intrinsic :: iso_c_binding
!!$        character(kind=c_char,len=*)  :: command
!!$     end subroutine
!!$  end interface
  ! ---------------------------------------------------------------------------
  
Contains

!!$  !============================================================================
!!$  !> Subroutine for executing system calls (Implemented because of missing
!!$  !> target of execute_command_line in CRAY fortran library)
!!$  Subroutine system(command)
!!$
!!$     character(kind=c_char,len=*),intent(in) :: command
!!$     call fsystem(command//c_null_char)
!!$
!!$  end subroutine
  
  !============================================================================
  !> Subroutine that initializes the output to std-out according to the
  !> environment variable CHAIN_STDOUT
  Subroutine init_std_out(success)

    Character(len=mcl)             :: env_var
    Logical                        :: opened, exist
    Logical, optional, intent(Out) :: success

    Integer                        :: io_stat
    
    Call get_environment_Variable("CHAIN_STDOUT", env_var)

    if (present(success)) success = .TRUE.

    !** If CHAIN_STDOUT is set to a value *************************************
    If (Len_Trim(env_var) > 0) then

       Inquire(file=trim(env_var), opened = opened)
       Inquire(file=trim(env_var), exist  = exist)

       IF (exist .AND. opened) Then
          
          Inquire(file=trim(env_var), number=un_mon)

       Else if (exist .AND. (.NOT.opened)) then

          un_mon = give_new_unit()
          Open(unit=un_mon,  file=trim(env_var), Action='Write', &
               status='old', position='Append' )

       Else if (.NOT. exist) then

          un_mon = give_new_unit()
          Open(unit=un_mon,  file=trim(env_var), Action='Write', &
               status='NEW', iostat = io_stat )
          If (io_stat /= 0) then
             if (present(success)) then
                success = .FALSE.
             Else
                un_mon = OUTPUT_UNIT
                Write(un_mon,fmt_sep)
                Write(un_mon,FMT_WRN )"In init_std_out it was not possible to open the monitor file"
                Write(un_mon,FMT_WRN_A)trim(env_var)
                Write(un_mon,FMT_WRN_A)"Please check the value of the en CHAIN_STDOUT env variable"
                Write(un_mon,FMT_WRN )"Monitoring will be redirected to std output"
                Write(un_mon,fmt_sep)
                
             End if
          End IF
       End IF
       
    End If
   
  End Subroutine init_std_out
  
  !============================================================================
  !> Subroutine for writing link start tag and opening the global log-file
  Subroutine link_start(link_name, init_lf, silent_stdio, success)

    Character(Len=*), Intent(in) :: link_name
    
    Logical, Intent(In),      Optional :: init_lf

    !> If .TRUE. no output to std_out is done (default = .FALSE.)
    Logical, Intent(In)     , optional :: silent_stdio
    Logical, Intent(Out)    , optional :: success
    
    Logical :: opened, exist, loc_init_lf, loc_stdio

    Integer :: io_stat = 0, ii

    Character(len=mcl) :: lf=''
    !--------------------------------------------------------------------------

    Call init_std_out()
      
    !** Check Outpath ********************************************
    If (outpath(len_trim(outpath):len_trim(outpath)) /= "/") then
       outpath = trim(outpath)//"/"
    End If
    
    call start_timer(trim(link_name))

    If (present(silent_stdio)) then
       loc_stdio = .NOT.silent_stdio
    Else
       loc_stdio = .TRUE.
    End If
    
    If (present(init_lf)) Then
       loc_init_lf = init_lf
    Else
       loc_init_lf = .FALSE.
    End If

    if (present(success)) success = .TRUE.
    
    If (loc_stdio) then
       Write(un_mon,*)
       Write(un_mon,fmt_sep)
       Write(un_mon,'(A,A)')'Starting chain link :',trim(link_name)
       Write(un_mon,*)
    End If

    Inquire(file=trim(outpath)//trim(project_name)//'.log', opened=opened)
    Inquire(file=trim(outpath)//trim(project_name)//'.log', exist=exist)

    IF (opened) Then

       Inquire(file=trim(outpath)//trim(project_name)//'.log', number=un_lf)

       !** Message to std out *************************************************
       If (loc_stdio) then
          
          Write(un_mon,FMT_MSG )"The log-file was already opened"

          Write(un_mon,FMT_MSG)'Reusing open and existing log-file :'
          Write(lf,'(A)')trim(outpath)//trim(project_name)//'.log'

          If (Len_trim(lf) > 72) Then
             Do ii = 1, Len_trim(lf), 72
                Write(un_mon,FMT_MSG)lf(ii:ii+71)
             End Do
          Else
             Write(un_mon,FMT_MSG)trim(outpath)//trim(project_name)//'.log'
          End If
          
       End If
       
    Else if (exist .AND. (.NOT.loc_init_lf)) then

       un_lf = give_new_unit()
       Open(unit=un_lf, file=trim(outpath)//trim(project_name)//'.log', &
            Action='Write', status='old', position='Append')

       !** Message to std out *************************************************
       If (loc_stdio) then
          
          Write(un_mon,FMT_MSG)'Opened existing log-file :'
          Write(lf,'(A)')trim(outpath)//trim(project_name)//'.log'
          
          If (Len_trim(lf) > 72) Then
             Do ii = 1, Len_trim(lf), 72
                Write(un_mon,FMT_MSG)lf(ii:ii+71)
             End Do
          Else
             Write(un_mon,FMT_MSG)trim(outpath)//trim(project_name)//'.log'
          End If
       End If
       
    Else if (.NOT.exist) Then

       un_lf = give_new_unit()
       Open(unit=un_lf, file=trim(outpath)//trim(project_name)//'.log', &
            Action='Write', status='new', iostat=io_stat)

       If (io_stat /= 0) then
          Write(un_mon,fmt_sep)
          Write(un_mon,FMT_ERR )"In link_start it was not possible to open the file"
          Write(un_mon,FMT_ERR_A)trim(outpath)//trim(project_name)//'.log'
          Write(un_mon,FMT_ERR )"Please check the path and file naming conventions"
          Write(un_mon,FMT_STOP)
          
          If (present(success)) then
             success = .FALSE.
             Goto 1000 
          Else
             STOP
          End If
          
       End If
       
       !** Message to std out *************************************************
       If (loc_stdio) then
          
          Write(un_mon,FMT_MSG)'Opened new log-file :'
          Write(lf,'(A)')trim(outpath)//trim(project_name)//'.log'
          
          If (Len_trim(lf) > 72) Then
             Do ii = 1, Len_trim(lf), 72
                Write(un_mon,FMT_MSG)lf(ii:ii+71)
             End Do
          Else
             Write(un_mon,FMT_MSG)trim(outpath)//trim(project_name)//'.log'
          End If
          
       End If
       
    Else if (loc_init_lf) then

       un_lf = give_new_unit()
       Open(unit=un_lf, file=trim(outpath)//trim(project_name)//'.log', &
            Action='Write', status='replace')
       
       !** Message to std out *************************************************
       If (loc_stdio) then
          
          Write(un_mon,FMT_MSG)'Opened and replaced existing log-file :'
          Write(lf,'(A)')trim(outpath)//trim(project_name)//'.log'
          
          If (Len_trim(lf) > 72) Then
             Do ii = 1, Len_trim(lf), 72
                Write(un_mon,FMT_MSG)lf(ii:ii+71)
             End Do
          Else
             Write(un_mon,FMT_MSG)trim(outpath)//trim(project_name)//'.log'
          End If
       End If
       
    End If
       
    If (loc_stdio) write(un_mon,*)

1000 Continue
    !Write(un_lf,fmt_sep)
    !Write(un_lf,FMT_MSG_A)'Starting chain link :'//link_name
    !Write(un_lf,*)

  End Subroutine link_start

  !============================================================================
  !> Subroutine for writing link end tag and closing the global log file
  Subroutine link_end(link_name, silent_stdio)

    Character(Len=*), Intent(in)       :: link_name

    !> If .TRUE. no output to std_out is done (default = .FALSE.)
    Logical, Intent(In)     , optional :: silent_stdio
    
    Logical                            :: opened, loc_stdio

    !--------------------------------------------------------------------------

    If (present(silent_stdio)) then
       loc_stdio = .NOT.silent_stdio
    Else
       loc_stdio = .TRUE.
    End If
    
    call end_timer(trim(link_name))

    if (loc_stdio) then
       Write(un_mon    ,FMT_EQ_SEP)
       Write(un_mon    ,FMT_MSG_A)'Program terminated correctly !'
       Write(un_mon    ,FMT_EQ_SEP)
    End if

    INQUIRE(unit=un_lf, opened=opened)

    if (.not.opened) then
       un_lf = give_new_unit()
       Open(unit=un_lf, file=trim(outpath)//trim(project_name)//'.log', &
            Action='Write', status='replace')
    End if

    call write_timelist(unit=un_lf)
       
    !Write(un_lf,FMT_MSG_A)'Program terminated correctly !' 
    !Write(un_lf,fmt_sep)
    !Write(un_lf,*)
    close(un_lf)
    
  End Subroutine link_end

  !============================================================================
  !> Subroutine for writing link error and stop message if something fails
  Subroutine link_stop(link_name,time,msg)

    Character(Len=*), Intent(in)      :: link_name, msg
    Real(kind=rk),Intent(in),optional :: time
    Logical                           :: opened

    !--------------------------------------------------------------------------

    call end_timer(trim(link_name))

    Write(un_mon    ,FMT_EQ_SEP)
    Write(un_mon    ,FMT_ERR_A)'Program was halted with ERROR !!'
    Write(un_mon    ,FMT_ERR_A)'Error message was:'
    Write(un_mon    ,FMT_ERR_A)msg
    Write(un_mon    ,FMT_EQ_SEP)

    INQUIRE(unit=un_lf, opened=opened)

    if (.not.opened) then
       un_lf = give_new_unit()
       Open(unit=un_lf, file=trim(outpath)//trim(project_name)//'.log', &
            Action='Write', status='replace')
    End if

    call write_timelist(unit=un_lf)
       
    Write(un_lf,FMT_ERR_A)'Program was halted with ERROR !!'
    Write(un_lf,FMT_ERR_A)'Error message was:'
    Write(un_lf,FMT_ERR_A)msg 
    Write(un_lf,fmt_sep)
    Write(un_lf,*)
    close(un_lf)
    
    STOP

  End Subroutine link_stop

  !============================================================================
  !> Subroutine which returns the elapsed real time in hh:mm:ss.sss format
  !>
  !> The subroutines calculates the elapsed real time from two parameter sets
  !> returned by the intinsic date_and_time(values)
  subroutine calc_real_time(tstart, tend, elapsed, echo)

    Integer, Dimension(8), intent(in) :: tstart, tend
    Character(len=12), intent(out)    :: elapsed
    Integer                           :: hh, mm, ss, msec
    Integer(Kind=8)                   :: msec_s, msec_e, msec_diff
    Logical, intent(in), optional     :: echo

    msec_s = 0
    msec_s = msec_s + tstart(8) + tstart(7) * 1000 + tstart(6) * 60 * 1000
    msec_s = msec_s + tstart(5) * 60 * 60 * 1000

    msec_e = 0
    msec_e = msec_e + tend(8) + tend(7) * 1000 + tend(6) * 60 * 1000
    msec_e = msec_e + tend(5) * 60 * 60 * 1000

    if (msec_e < msec_s) then
       msec_e = msec_e + 24*60*60*1000
    End if

    msec_diff = msec_e - msec_s

    hh = msec_diff/(60 * 60 * 1000)
    mm = (msec_diff - hh * (60 * 60 * 1000)) / (60 * 1000)
    ss = (msec_diff - hh * (60 * 60 * 1000) - mm * 60 * 1000) / 1000
    msec = (msec_diff - hh * (60 * 60 * 1000) - mm * 60 * 1000 - ss * 1000)

    write(elapsed,"(I2.2,':',I2.2,':',I2.2,'.',I3.3)") hh,mm,ss,msec         

    If (present(echo)) then
       If (echo) Write(un_mon    ,FMT_TIME)&
            'Elapsed time was : '//elapsed//' : ',Real(msec_diff)/1000.
    End If

  End subroutine calc_real_time
  
  !============================================================================
  !> Subroutine for reading strings when parsing user input data
  Subroutine read_char_in(var,name,un)

    Integer(Kind=ik), Intent(In), Optional :: un
    Character(Len=*), Intent(IN)  :: name
    Character(Len=*), Intent(OUT) :: var
    
    Integer                       :: io_stat=0

    !==========================================================================

    If (present(un)) then
       Read (un,*,iostat=io_stat) var
    else
       Read (INPUT_UNIT,*,iostat=io_stat) var
    end If
    

    If (io_stat /=0) Then

       Write(un_mon,fmt_sep)
       Write(un_mon,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       Write(un_mon,'(A)')'!! ERROR while reading user input with read_char_in(var,name) !!'
       Write(un_mon,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       Write(un_mon,*)
       Write(un_mon,'(A,A)')'For variable : ',Trim(name)
       Write(un_mon,'(A)')'Something went wrong during read                              !!'
       Write(un_mon,*)
       Write(un_mon,'(A,I0)')'Program halted with I/O read status = ',io_stat
       Write(un_mon,fmt_sep)
       Stop

    Else

       Write(un_mon,fmt_inpsep)

       If (Len_trim(name) > 75) Then
          Write(un_mon,"('+--: ',A75)") name
       Else
          Write(un_mon,"('+--: ',A)") Trim(name)
       End If
       
       If (Len_trim(var) > 75) Then
          Write(un_mon,"('+--> ',A75)") var(1:75)
          Write(un_mon,"('+    ',A75)") var(76:150)
       Else
          Write(un_mon,"('+--> ',A)") Trim(var)
       End If

    End If

  End Subroutine read_char_in

  !============================================================================
  !> Subroutine for reading Integer kind=4 values when parsing user input data
  Subroutine read_int4_in(var,name,un)

    Integer(Kind=ik), Intent(In), Optional :: un

    Character(Len=*)   , Intent(IN)  :: name
    Integer(kind=4)    , Intent(OUT) :: var

    Integer                       :: io_stat

    !==========================================================================

    If (present(un)) then
       Read (un,*,iostat=io_stat) var
    else
       Read (INPUT_UNIT,*,iostat=io_stat) var
    end If

    If (io_stat /=0) Then

       Write(un_mon,fmt_sep)
       Write(un_mon,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       Write(un_mon,'(A)')'!! ERROR while reading user input !!'
       Write(un_mon,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       Write(un_mon,*)
       Write(un_mon,'(A,A)')'For variable : ',Trim(name)
       Write(un_mon,'(A)')'A INTEGER value is expected !!'
       Write(un_mon,*)
       Write(un_mon,'(A,I0)')'Program halted with I/O read status = ',io_stat
       Write(un_mon,fmt_sep)
       Stop

    Else

       Write(un_mon,fmt_inpsep)

       If (Len_trim(name) > 75) Then
          Write(un_mon,"('+--: ',A75)") name
       Else
          Write(un_mon,"('+--: ',A)") Trim(name)
       End If

       Write(un_mon,"('+--> ',I0)") var

    End If

  End Subroutine read_int4_in

  !============================================================================
  !> Subroutine for reading Integer Kind=8 values when parsing user input data
  Subroutine read_int8_in(var,name,un)

    Integer(Kind=ik), Intent(In), Optional :: un

    Character(Len=*)   , Intent(IN)  :: name
    Integer(Kind=8)   , Intent(OUT) :: var

    Integer                       :: io_stat

    !==========================================================================

    If (present(un)) then
       Read (un,*,iostat=io_stat) var
    else
       Read (INPUT_UNIT,*,iostat=io_stat) var
    end If

    If (io_stat /=0) Then

       Write(un_mon,fmt_sep)
       Write(un_mon,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       Write(un_mon,'(A)')'!! ERROR while reading user input !!'
       Write(un_mon,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       Write(un_mon,*)
       Write(un_mon,'(A,A)')'For variable : ',Trim(name)
       Write(un_mon,'(A)')'A INTEGER value is expected !!'
       Write(un_mon,*)
       Write(un_mon,'(A,I0)')'Program halted with I/O read status = ',io_stat
       Write(un_mon,fmt_sep)
       Stop

    Else

       Write(un_mon,fmt_inpsep)

       If (Len_trim(name) > 75) Then
          Write(un_mon,"('+--: ',A75)") name
       Else
          Write(un_mon,"('+--: ',A)") Trim(name)
       End If

       Write(un_mon,"('+--> ',I0)") var

    End If

  End Subroutine read_int8_in

  !============================================================================
  !> Subroutine for reading floating point numbers when parsing user input 
  !> data
  Subroutine read_real_in(var,name,un)

    Integer(Kind=ik), Intent(In), Optional :: un

    Character(Len=*)   , Intent(IN)  :: name
    Real(Kind=rk)   , Intent(OUT) :: var

    Integer                       :: io_stat

    !==========================================================================

    If (present(un)) then
       Read (un,*,iostat=io_stat) var
    else
       Read (INPUT_UNIT,*,iostat=io_stat) var
    end If

    If (io_stat /=0) Then

       Write(un_mon,fmt_sep)
       Write(un_mon,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       Write(un_mon,'(A)')'!! ERROR while reading user input !!'
       Write(un_mon,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       Write(un_mon,*)
       Write(un_mon,'(A,A)')'For variable : ',Trim(name)
       Write(un_mon,'(A)')'A FLOATING POINT VALUE value is expected !!'
       Write(un_mon,*)
       Write(un_mon,'(A,I0)')'Program halted with I/O read status = ',io_stat
       Write(un_mon,fmt_sep)
       Stop

    Else

       Write(un_mon,fmt_inpsep)

       If (Len_trim(name) > 75) Then
          Write(un_mon,"('+--: ',A75)") name
       Else
          Write(un_mon,"('+--: ',A)") Trim(name)
       End If

       Write(un_mon,"('+--> ',F0.3)") var

    End If

  End Subroutine read_real_in

  !============================================================================
  !> Subroutine for reading logical values when parsing user input data
  Subroutine read_log_in(var,name,un)

    Integer(Kind=ik), Intent(In), Optional :: un

    Character(Len=*)   , Intent(IN)  :: name
    logical            , Intent(OUT) :: var

    Integer                       :: io_stat

    !==========================================================================

    If (present(un)) then
       Read (un,*,iostat=io_stat) var
    else
       Read (INPUT_UNIT,*,iostat=io_stat) var
    end If

    If (io_stat /=0) Then

       Write(un_mon,fmt_sep)
       Write(un_mon,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       Write(un_mon,'(A)')'!! ERROR while reading user input !!'
       Write(un_mon,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       Write(un_mon,*)
       Write(un_mon,'(A,A)')'For variable : ',Trim(name)
       Write(un_mon,'(A)')'A Logical VALUE value is expected !!'
       Write(un_mon,*)
       Write(un_mon,'(A,I0)')'Program halted with I/O read status = ',io_stat
       Write(un_mon,fmt_sep)
       Stop

    Else

       Write(un_mon,fmt_inpsep)

       If (Len_trim(name) > 75) Then
          Write(un_mon,"('+--: ',A75)") name
       Else
          Write(un_mon,"('+--: ',A)") Trim(name)
       End If

       Write(un_mon,"('+--> ',L4)") var

    End If

  End Subroutine read_log_in

  !============================================================================
  !> Subroutine for reading floating point 3-vectors when parsing user input 
  !> data
  Subroutine read_real_3vector_in(var,name,un)

    Integer(Kind=ik), Intent(In), Optional :: un

    Character(Len=*)            , Intent(IN)  :: name
    Real(Kind=rk), Dimension(3) , Intent(OUT) :: var

    Integer                       :: io_stat

    !==========================================================================

    If (present(un)) then
       Read (un,*,iostat=io_stat) var
    else
       Read (INPUT_UNIT,*,iostat=io_stat) var
    end If

    If (io_stat /=0) Then

       Write(un_mon,fmt_sep)
       Write(un_mon,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       Write(un_mon,'(A)')'!! ERROR while reading user input !!'
       Write(un_mon,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       Write(un_mon,*)
       Write(un_mon,'(A,A)')'For variable : ',Trim(name)
       Write(un_mon,'(A)')'3 FLOATING POINT VALUES are expected !!'
       Write(un_mon,*)
       Write(un_mon,'(A,I0)')'Program halted with I/O read status = ',io_stat
       Write(un_mon,fmt_sep)
       Stop

    Else

       Write(un_mon,fmt_inpsep)

       If (Len_trim(name) > 75) Then
          Write(un_mon,"('+--: ',A75)") name
       Else
          Write(un_mon,"('+--: ',A)") Trim(name)
       End If

       Write(un_mon,"('+--> ',F0.3,2(' ,',F0.3))") var

    End If

  End Subroutine read_real_3vector_in

  !============================================================================
  !> Subroutine for reading integer 3-vectors of kind ik when parsing user 
  !> input data
  Subroutine read_int4_3vector_in(var,name,un)

    Integer(Kind=ik), Intent(In), Optional :: un

    Character(Len=*)               , Intent(IN)  :: name
    Integer(Kind=4),  Dimension(3) , Intent(OUT) :: var

    Integer                       :: io_stat

    !==========================================================================

    If (present(un)) then
       Read (un,*,iostat=io_stat) var
    else
       Read (INPUT_UNIT,*,iostat=io_stat) var
    end If

    If (io_stat /=0) Then

       Write(un_mon,fmt_sep)
       Write(un_mon,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       Write(un_mon,'(A)')'!! ERROR while reading user input !!'
       Write(un_mon,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       Write(un_mon,*)
       Write(un_mon,'(A,A)')'For variable : ',Trim(name)
       Write(un_mon,'(A)')'3 INTEGER VALUES are expected !!'
       Write(un_mon,*)
       Write(un_mon,'(A,I0)')'Program halted with I/O read status = ',io_stat
       Write(un_mon,fmt_sep)
       Stop

    Else

       Write(un_mon,fmt_inpsep)

       If (Len_trim(name) > 75) Then
          Write(un_mon,"('+--: ',A75)") name
       Else
          Write(un_mon,"('+--: ',A)") Trim(name)
       End If

       Write(un_mon,"('+--> ',I0,2(' ,',I0))") var

    End If

  End Subroutine read_int4_3vector_in

  !============================================================================
  !> Subroutine for reading integer 3-vectors of kind ik when parsing user 
  !> input data
  Subroutine read_int8_3vector_in(var,name,un)

    Integer(Kind=ik), Intent(In), Optional :: un

    Character(Len=*)               , Intent(IN)  :: name
    Integer(Kind=8),  Dimension(3) , Intent(OUT) :: var

    Integer                       :: io_stat

    !==========================================================================

    If (present(un)) then
       Read (un,*,iostat=io_stat) var
    else
       Read (INPUT_UNIT,*,iostat=io_stat) var
    end If

    If (io_stat /=0) Then

       Write(un_mon,fmt_sep)
       Write(un_mon,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       Write(un_mon,'(A)')'!! ERROR while reading user input !!'
       Write(un_mon,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       Write(un_mon,*)
       Write(un_mon,'(A,A)')'For variable : ',Trim(name)
       Write(un_mon,'(A)')'3 INTEGER VALUES are expected !!'
       Write(un_mon,*)
       Write(un_mon,'(A,I0)')'Program halted with I/O read status = ',io_stat
       Write(un_mon,fmt_sep)
       Stop

    Else

       Write(un_mon,fmt_inpsep)

       If (Len_trim(name) > 75) Then
          Write(un_mon,"('+--: ',A75)") name
       Else
          Write(un_mon,"('+--: ',A)") Trim(name)
       End If

       Write(un_mon,"('+--> ',I0,2(' ,',I0))") var

    End If

  End Subroutine read_int8_3vector_in

  !============================================================================
  !> Subroutine for reading floating point matrices when parsing user input
  !> data
  Subroutine read_real_matrix_in(var,name,un)

    Integer(Kind=ik), Intent(In), Optional :: un

    Character(Len=*)              , Intent(IN)  :: name
    Real(Kind=rk), Dimension(:,:) , Intent(OUT) :: var

    Integer                       :: io_stat
    Integer, Dimension(2)         :: sp_var
    Character(len=100)            :: fmt

    !==========================================================================

    If (present(un)) then
       Read (un,*,iostat=io_stat) var
    else
       Read (INPUT_UNIT,*,iostat=io_stat) var
    end If

    If (io_stat /=0) Then

       Write(un_mon,fmt_sep)
       Write(un_mon,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       Write(un_mon,'(A)')'!! ERROR while reading user input !!'
       Write(un_mon,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       Write(un_mon,*)
       Write(un_mon,'(A,A)')'For matrix variable : ',Trim(name)
       Write(un_mon,'(I0,A)')size(var),' FLOATING POINT VALUES are expected !!'
       Write(un_mon,*)
       Write(un_mon,'(A,I0)')'Program halted with I/O read status = ',io_stat
       Write(un_mon,fmt_sep)
       Stop

    Else

       Write(un_mon,fmt_inpsep)

       If (Len_trim(name) > 75) Then
          Write(un_mon,"('+--: ',A75)") name
       Else
          Write(un_mon,"('+--: ',A)") Trim(name)
       End If

       sp_var = shape(var)

       If ((sp_var(1) > 10) .Or. (sp_var(2) > 10)) then
          Write(un_mon,"('+--> ',A)")"Repetiton of input data is suppressed &
                                 &due to their amount"
       Else
          Write(un_mon,"('+--> +',T80,'+')")
          Write(fmt,"(A,I0,A,I0,A)")"(",sp_var(2),"('     |',",sp_var(1),&
               "('  ',E10.4),' |',/),'     +',T80,'+')"
          Write(un_mon,fmt) var

       End If

    End If

  End Subroutine read_real_matrix_in

  !============================================================================
  !> Subroutine for allocation error handling
  SUBROUTINE alloc_err(in_var,io_stat)

    INTEGER             :: io_stat
    CHARACTER (LEN=*)   :: in_var

    IF (io_stat /= 0) Then
       write(un_mon,*)
       WRITE(un_mon,fmt_sep)
       WRITE(un_mon,FMT_ERR)'Allocation of var :'       
       WRITE(un_mon,FMT_ERR) in_var
       WRITE(un_mon,FMT_ERR)'faild !!'
       WRITE(un_mon,FMT_ERR_AI0)'With Allocation Status ',io_stat
       WRITE(un_mon,FMT_STOP)
       STOP       
    End IF

  END SUBROUTINE alloc_err

  !============================================================================
  !> Subroutine for I/O error handling while operating on files
  SUBROUTINE file_err(in_file,io_stat)

    INTEGER             :: io_stat
    CHARACTER (LEN=*)   :: in_file

    IF (io_stat /= 0) Then
       WRITE(un_mon,fmt_sep)
       WRITE(un_mon,FMT_ERR)'Operation on file :'       
       WRITE(un_mon,FMT_ERR) in_file
       WRITE(un_mon,FMT_ERR)'faild !!'
       WRITE(un_mon,FMT_ERR_AI0)'With I/O Status ',io_stat
       WRITE(un_mon,FMT_STOP)
       STOP
    End IF

  END SUBROUTINE file_err

  !============================================================================
  !> Function which returns a new free unit
  function give_new_unit() result(new_unit)

    Integer :: new_unit

    Integer :: ii

    Logical :: unit_is_open

    Do ii = 300, huge(new_unit)-1

       inquire(unit=ii, opened=unit_is_open)

       if( .not.unit_is_open ) then
          new_unit = ii
          Exit
       end if

    End Do
    
    if ( unit_is_open ) then
       
       WRITE(un_mon,fmt_sep)
       WRITE(un_mon,FMT_ERR)'Something bad and unexpected happened during search for free unit'
       WRITE(un_mon,FMT_ERR)'Could not find a new unit between 100 and huge(Int(kind=4))'
       WRITE(un_mon,FMT_ERR)' '
       WRITE(un_mon,FMT_STOP)
       STOP
    END IF

  End function give_new_unit

  !============================================================================
  !> Subroutine for counting lines in an ascii file
  function count_lines(un) result(no_lines)

    Integer          ,Intent(in) :: un
    Integer(kind=ik)             :: no_lines

    Integer                      :: io_stat
    Character(len=2)             :: temp_char

    io_stat = 0
    no_lines=0

    Rewind(un)

    Do While (io_stat == 0)

       Read(un,'(A)', End=1000, iostat=io_stat) temp_char
       no_lines = no_lines + 1

    End Do

1000 Continue

    Rewind(un)

  End function count_lines

  !============================================================================
  !> Subroutine for formated matrix output of Real(kind=rk) values
  Subroutine Write_real_matrix(un,a,ii,jj,name,Fmt)

    Integer         , intent(in)                   :: un
    Integer(kind=ik), intent(in)                   :: ii, jj
    Real(kind=rk)   , intent(in), dimension(ii,jj) :: a
    Character(len=*), intent(in)                   :: name
    Character(len=8), Intent(in), Optional         :: fmt
    Character(len=8)                               :: fmt_loc
    Integer(kind=ik)                               :: kk

    Integer            :: prec , fw
    Character(len=mcl) :: fmt_a, sep

    If (present(fmt)) then
       fmt_loc = fmt
    Else
       fmt_loc = "default"
    End If

    prec = precision(a)
    fw   = prec+8

    If (trim(fmt_loc) == "wxmaxima") then
       
       write(fmt_a,'(5(A,I0),A)') &
            "(' [',",jj-1,"(E",fw,".",prec,"E2,','),E",fw,".",prec,"E2,'],' )"

       Write(un_lf,"(A,A)")trim(name),": matrix("

       Do kk = 1, ii-1
          Write(un_lf,fmt_a)a(kk,:)
       End Do

       write(fmt_a,'(5(A,I0),A)') &
            "(' [',",jj-1,"(E",fw,".",prec,"E2,','),E",fw,".",prec,"E2,']);' )"

       Write(un_lf,fmt_a)a(ii,:)

    Else

       !** generate formats *************************************
       write(fmt_a,'(3(A,I0),A)')'(',jj,'(E',fw,'.',prec,'E2))'
       write(sep  ,"(A,I0,A)")"(",fw*jj+10,"('='))"
       
       Write(un,sep)
       Write(un,"(6('-'),A,6('-'))")trim(name)
       write(un,fmt_a)transpose(a)

    End If

  End Subroutine Write_real_matrix

  !****************************************************************************
  !> Function which determins whether the input matrix is symmetric or not
  subroutine check_sym(A,name,un,sym_out,sup_log)

    Real(Kind=rk), Dimension(:,:), Intent(in)  :: A
    character(len=*), optional   , Intent(in)  :: name
    Integer         , optional   , Intent(in)  :: un   

    Real(kind=rk)   , optional   , Intent(out) :: sym_out

    Logical         , optional   , Intent(in)  :: sup_log

    Logical                                    :: spl_loc

    Real(kind=rk)                              :: sym

    Integer      , Dimension(2)   :: lb,ub
    Integer                       :: un_loc,ii,jj

    if (present(un)) then
       un_loc = un
    Else
       un_loc = un_lf
    End if

    if (present(sup_log)) then
       spl_loc = sup_log
    Else
       spl_loc = .FALSE.
    End if

    lb = lbound(A)
    ub = ubound(A)

    sym = 0._rk

    Do jj = lb(2), ub(2)
       DO ii = lb(1), ub(1)

          sym = sym + abs(A(ii,jj) -  A(jj,ii))

       End DO
    End Do
    
    If (.not.spl_loc) then
       Write(un_loc,*)
       If (present(name)) then
          Write(un_loc,*)"--------- check_sym ",trim(name)," ------------"
       Else
          Write(un_loc,*)"--------- check_sym ------------"
       End If
       write(un_loc,*)sym
    End If

    if (present(sym_out)) sym_out=sym

  End subroutine check_sym

End Module chain_routines
