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
  USE messages_errors
  USE meta

  Implicit None
 
  ! ---------------------------------------------------------------------------
  !> Logfile unit
  Integer                     :: un_lf   = fhl
  !> Monitor file unit (default = stdout)
  Integer                     :: un_mon  = fhmon
  
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
   
   USE messages_errors
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

    Integer :: io_stat = 0
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
      CALL print_message(un_mon, 'Starting chain link: '//TRIM(link_name))
      CALL print_sep(un_mon)
    End If

    Inquire(file=trim(outpath)//trim(project_name)//'.log', opened=opened)
    Inquire(file=trim(outpath)//trim(project_name)//'.log', exist=exist)

    IF (opened) Then

       Inquire(file=trim(outpath)//trim(project_name)//'.log', number=un_lf)

       !** Message to std out *************************************************
       If (loc_stdio) then         
         CALL print_message(un_mon, "The log-file was already opened.")
         CALL print_message(un_mon, "Reusing open and existing log-file: "&
            //TRIM(outpath)//TRIM(project_name)//'.log')
         CALL print_sep(un_mon)
       End If
       
    Else if (exist .AND. (.NOT.loc_init_lf)) then

       un_lf = give_new_unit()
       Open(unit=un_lf, file=trim(outpath)//trim(project_name)//'.log', &
            Action='Write', status='old', position='Append')

       !** Message to std out *************************************************
       If (loc_stdio) then
         CALL print_message(un_mon, "Opened existing log-file: "//TRIM(outpath)//TRIM(project_name)//'.log')
       End If
       
    Else if (.NOT.exist) Then

       un_lf = give_new_unit()
       Open(unit=un_lf, file=trim(outpath)//trim(project_name)//'.log', &
            Action='Write', status='new', iostat=io_stat)

       If (io_stat /= 0) then
         mssg =  "In link_start it was not possible to open the file "//&
         trim(outpath)//trim(project_name)//'.log'//" Please check the path and file naming conventions"

         CALL print_err_stop(un_mon, mssg, io_stat)
          
          If (present(success)) then
             success = .FALSE.
             Goto 1000 
          Else
             STOP
          End If
          
       End If
       
       !** Message to std out *************************************************
       If (loc_stdio) then
         CALL print_message(un_mon, "Opened new log-file: "//TRIM(outpath)//TRIM(project_name)//'.log')          
       End If
       
    Else if (loc_init_lf) then

       un_lf = give_new_unit()
       Open(unit=un_lf, file=trim(outpath)//trim(project_name)//'.log', &
            Action='Write', status='replace')
       
       !** Message to std out *************************************************
       If (loc_stdio) then
         CALL print_message(un_mon, "Opened and replaced existing log-file: "&
            //TRIM(outpath)//TRIM(project_name)//'.log')          
       End If
       
    End If
       
    If (loc_stdio) CALL print_sep(un_mon)

1000 Continue

  End Subroutine link_start

  !============================================================================
  !> Subroutine for writing link end tag and closing the global log file
  Subroutine link_end(link_name, silent_stdio)

    Character(Len=*), Intent(in) :: link_name

    !> If .TRUE. no output to std_out is done (default = .FALSE.)
    Logical, Intent(In), optional :: silent_stdio
    
    Logical :: opened, loc_stdio

    !--------------------------------------------------------------------------

    If (present(silent_stdio)) then
       loc_stdio = .NOT.silent_stdio
    Else
       loc_stdio = .TRUE.
    End If
    
    call end_timer(trim(link_name))

    if (loc_stdio) CALL print_message(un_mon, "Program terminated correctly.")          

    INQUIRE(unit=un_lf, opened=opened)

    IF (.not.opened) THEN
       un_lf = give_new_unit()
       OPEN(UNIT=un_lf, FILE=TRIM(outpath)//TRIM(project_name)//'.log', ACTION='WRITE', STATUS='REPLACE')
    END IF

    call write_timelist(unit=un_lf)
       
    close(un_lf)
    
  End Subroutine link_end

  !============================================================================
  !> Subroutine for writing link error and stop message if something fails
  Subroutine link_stop(link_name,msg)

   Character(Len=*), Intent(in)      :: link_name, msg
   Logical                           :: opened

   !--------------------------------------------------------------------------

   call end_timer(trim(link_name))

   CALL print_err_stop(un_mon, 'Program will be halted, error message was: '//msg, 0)

   INQUIRE(UNIT=un_lf, opened=opened)

   IF (.not.opened) THEN
      un_lf = give_new_unit()
      OPEN(UNIT=un_lf, FILE=TRIM(outpath)//TRIM(project_name)//'.log', ACTION='WRITE', STATUS='REPLACE')
   END IF

   call write_timelist(unit=un_lf)
      
   CALL print_err_stop(un_lf,   'Program will be halted, error message was: '//msg, 0)
   CALL print_err_stop(std_out, 'Program halted. View Log.', 1)

  End Subroutine link_stop


   !------------------------------------------------------------------------------
   ! SUBROUTINE: alloc_err
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
      CALL print_err_stop(un_mon, mssg, io_stat)
   End IF

   END SUBROUTINE alloc_err

End Module chain_routines
