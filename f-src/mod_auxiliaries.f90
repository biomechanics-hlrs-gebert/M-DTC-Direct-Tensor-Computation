!------------------------------------------------------------------------------
! HLRS - Numerical Methods and Libraries, Nobelstraße 19, 70569 Stumgart
!------------------------------------------------------------------------------
! MODULE: auxiliaries
!
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
! DESCRIPTION: 
!> Module containing useful routines.
!
! REVISION HISTORY:
! 27 09 2021 - Initial Version
!------------------------------------------------------------------------------

MODULE auxiliaries

USE global_std
USE messages_errors
USE strings


IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: calc_real_time
!------------------------------------------------------------------------------  
!> @author Ralf Schneider, schneider@hlrs.de, HLRS/NUM
!
!> @brief
!> Returns the elapsed real time in hh:mm:ss.sss format
!
!> @description
!> The subroutines calculates the elapsed real time from two parameter sets
!> returned by the intinsic date_and_time(values)
!
!> @param[in] tstart Measured start time
!> @param[in] tend Measured end time
!> @param[in] elapsed Time elapsed
!> @param[in] echo Print elapsed time to a file handle
!------------------------------------------------------------------------------ 
subroutine calc_real_time(fh, tstart, tend, elapsed, echo)

Integer              , intent(in)  :: fh
Integer, Dimension(8), intent(in)  :: tstart
Integer, Dimension(8), intent(in)  :: tend
Character(len=12)    , intent(out) :: elapsed
Logical, optional    , intent(in)  :: echo

Integer                           :: hh, mm, ss, msec
Integer(Kind=8)                   :: msec_s, msec_e, msec_diff

msec_s = 0
msec_s = msec_s + tstart(8) + tstart(7) * 1000 + tstart(6) * 60 * 1000
msec_s = msec_s + tstart(5) * 60 * 60 * 1000

msec_e = 0
msec_e = msec_e + tend(8) + tend(7) * 1000 + tend(6) * 60 * 1000
msec_e = msec_e + tend(5) * 60 * 60 * 1000

if (msec_e < msec_s) msec_e = msec_e + 24*60*60*1000

msec_diff = msec_e - msec_s

hh = msec_diff/(60 * 60 * 1000)
mm = (msec_diff - hh * (60 * 60 * 1000)) / (60 * 1000)
ss = (msec_diff - hh * (60 * 60 * 1000) - mm * 60 * 1000) / 1000
msec = (msec_diff - hh * (60 * 60 * 1000) - mm * 60 * 1000 - ss * 1000)

write(elapsed,"(I2.2,':',I2.2,':',I2.2,'.',I3.3)") hh,mm,ss,msec         

If (present(echo)) then
    If (echo) Write(fh, "('MM ',A,1X,F0.6,' sec')")&
        'Elapsed time was: '//elapsed//' : ',Real(msec_diff)/1000.
End If

End subroutine calc_real_time

END MODULE auxiliaries
