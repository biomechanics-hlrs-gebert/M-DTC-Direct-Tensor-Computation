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
 
  USE error_handling
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
  
Contains

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


   !------------------------------------------------------------------------------
   ! SUBROUTINE: handle_err
   !------------------------------------------------------------------------------  
   !> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
   !
   !> @brief
   !> Wrapper for allocation error handling
   !
   !> @param[in] in_var Input variable
   !> @param[in] io_stat Status of allocation.
   !------------------------------------------------------------------------------  
   SUBROUTINE alloc_err(in_var, io_stat)

   INTEGER, INTENT(IN) :: io_stat
   CHARACTER(LEN=*), INTENT(IN) :: in_var

   IF (io_stat /= 0) Then
      WRITE(mssg, '(A,I4,A)') "Allocation of var ", TRIM(in_var), " failed."
      CALL handle_err(un_mon, mssg, io_stat)
   End IF

   END SUBROUTINE alloc_err

End Module chain_routines
