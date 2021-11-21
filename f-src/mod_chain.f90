!==============================================================================
!> \file mod_chain.f90
!> Modules of the chain process library.
!>
!> The file holds modules which are necessary to use the chain process library
!>
!> \author Ralf Schneider
!> \date 07.05.2012

!==============================================================================
!> Global Variables for the chain process library
Module chain_variables
 
  USE global_std

  Implicit None
 
  ! ---------------------------------------------------------------------------
  !> Logfile unit
  Integer                     :: un_lf   = 10000
  !> Monitor file unit (default = stdout)
  Integer                     :: un_mon  = std_out
  
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
   
   USE error_handling
   USE chain_variables
   USE timer

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
                un_mon = std_out
                Write(un_mon,SEP_STD)
                Write(un_mon,FMT_WRN )"In init_std_out it was not possible to open the monitor file"
                Write(un_mon,FMT_WRN_A)trim(env_var)
                Write(un_mon,FMT_WRN_A)"Please check the value of the en CHAIN_STDOUT env variable"
                Write(un_mon,FMT_WRN )"Monitoring will be redirected to std output"
                Write(un_mon,SEP_STD)
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

   !  Call init_std_out()
      
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
       Write(un_mon,SEP_STD)
       Write(un_mon,'(A,A)')'Starting chain link: ',trim(link_name)
       Write(un_mon,*)
    End If

    Inquire(file=trim(outpath)//trim(project_name)//'.log', opened=opened)
    Inquire(file=trim(outpath)//trim(project_name)//'.log', exist=exist)

    IF (opened) Then

       Inquire(file=trim(outpath)//trim(project_name)//'.log', number=un_lf)

       !** Message to std out *************************************************
       If (loc_stdio) then
          
          Write(un_mon,FMT_MSG )"The log-file was already opened"

          Write(un_mon,FMT_MSG)'Reusing open and existing log-file: '
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
          Write(un_mon,SEP_STD)
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
          
          Write(un_mon,FMT_MSG)'Opened new log-file: '
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
    !Write(un_lf,SEP_STD)
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
       Write(un_mon, FMT_EQ_SEP)
       Write(un_mon, FMT_MSG_A) 'Program terminated correctly !'
       Write(un_mon, FMT_EQ_SEP)
    End if

    INQUIRE(unit=un_lf, opened=opened)

    IF (.not.opened) THEN
       un_lf = give_new_unit()
       OPEN(UNIT=un_lf, FILE=TRIM(outpath)//TRIM(project_name)//'.log', ACTION='WRITE', STATUS='REPLACE')
    END IF

    call write_timelist(unit=un_lf)
       
    !Write(un_lf,FMT_MSG_A)'Program terminated correctly !' 
    !Write(un_lf,SEP_STD)
    !Write(un_lf,*)
    close(un_lf)
    
  End Subroutine link_end

  !============================================================================
  !> Subroutine for writing link error and stop message if something fails
  Subroutine link_stop(link_name,msg)

   Character(Len=*), Intent(in)      :: link_name, msg
   Logical                           :: opened

   !--------------------------------------------------------------------------

   call end_timer(trim(link_name))

   CALL handle_err(un_mon, 'Program will be halted, error message was: '//msg, 0)

   INQUIRE(UNIT=un_lf, opened=opened)

   IF (.not.opened) THEN
      un_lf = give_new_unit()
      OPEN(UNIT=un_lf, FILE=TRIM(outpath)//TRIM(project_name)//'.log', ACTION='WRITE', STATUS='REPLACE')
   END IF

   call write_timelist(unit=un_lf)
      
   CALL handle_err(un_lf,   'Program will be halted, error message was: '//msg, 0)
   CALL handle_err(std_out, 'Program halted. View Log.', 1)

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
      Read (std_in,*,iostat=io_stat) var
   end If
    

   If (io_stat /=0) Then
      WRITE(mssg,'(3A)') "Reading character var ", TRIM(name), " failed."
      CALL handle_err(un_mon, mssg, io_stat)
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
      Read (std_in,*,iostat=io_stat) var
   end If

   If (io_stat /=0) Then
      WRITE(mssg,'(3A)') "Reading var ", TRIM(name), " failed. 1 integer expected."
      CALL handle_err(un_mon, mssg, io_stat)
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
      Read (std_in,*,iostat=io_stat) var
   end If

   If (io_stat /=0) Then
      WRITE(mssg,'(3A)') "Reading var ", TRIM(name), " failed. 1 integer expected."
      CALL handle_err(un_mon, mssg, io_stat)
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
      Read (std_in,*,iostat=io_stat) var
   end If

   If (io_stat /=0) Then
      WRITE(mssg,'(3A)') "Reading var ", TRIM(name), " failed. 1 floating point value expected."
      CALL handle_err(un_mon, mssg, io_stat)
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
      Read (std_in,*,iostat=io_stat) var
   end If

      If (io_stat /=0) Then
         WRITE(mssg,'(3A)') "Reading var ", TRIM(name), " failed. A logical value expected."
         CALL handle_err(un_mon, mssg, io_stat)
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
      Read (std_in,*,iostat=io_stat) var
   end If

      If (io_stat /=0) Then
         WRITE(mssg,'(3A)') "Reading var ", TRIM(name), " failed. 3 floating point values expected."
         CALL handle_err(un_mon, mssg, io_stat)
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
       Read (std_in,*,iostat=io_stat) var
    end If

    If (io_stat /=0) Then
      WRITE(mssg,'(3A)') "Reading var ", TRIM(name), " failed. 3 integers expected."
      CALL handle_err(un_mon, mssg, io_stat)
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
      Read (std_in,*,iostat=io_stat) var
   end If

   If (io_stat /=0) Then
      WRITE(mssg,'(3A)') "Reading var ", TRIM(name), " failed. 3 integers expected."
      CALL handle_err(un_mon, mssg, io_stat)
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
      Read (std_in,*,iostat=io_stat) var
   end If

   If (io_stat /=0) Then
      WRITE(mssg,'(3A,I15,A)') "Reading matrix var ", TRIM(name), " of size ", SIZE(var)," failed. Floating points expected."
      CALL handle_err(un_mon, mssg, io_stat)
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
      WRITE(mssg, '(A,I4,A)') "Allocation of var ", TRIM(in_var), " failed."
      CALL handle_err(un_mon, mssg, io_stat)
   End IF

   END SUBROUTINE alloc_err


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

       WRITE(un_mon, SEP_STD)
       WRITE(un_mon, FMT_ERR)'Something bad and unexpected happened during search for free unit'
       WRITE(un_mon, FMT_ERR)'Could not find a new unit between 100 and huge(Int(kind=4))'
       WRITE(un_mon, FMT_ERR)' '
       WRITE(un_mon, FMT_STOP)
       STOP
    END IF

  End function give_new_unit

End Module chain_routines
