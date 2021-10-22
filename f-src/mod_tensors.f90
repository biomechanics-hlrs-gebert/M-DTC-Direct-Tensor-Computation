module tensors

  implicit none

  Integer, Parameter, Private :: mt_rk = 8

  !============================================================================
  !== Interfaces

  Interface t_transform

     Module Procedure t_transform_1
     Module Procedure t_transform_2
     Module Procedure t_transform_4

  End Interface

contains

  !============================================================================
  !==
  !==
  function t_transform_1(a, base_i, base_m) Result(a_t)

    !** Input tensor**************************************
    Real(kind=mt_rk), Dimension(3),   Intent(in) :: a
    !** Original base ************************************
    Real(kind=mt_rk), Dimension(3,3), Intent(in) :: base_i
    !** Base to transform a to ***************************
    Real(kind=mt_rk), Dimension(3,3), Intent(in) :: base_m

    !** Output tensor ************************************
    Real(kind=mt_rk), Dimension(3)             :: a_t

    !*****************************************************

    !** Transformation tensor ****************************
    Real(kind=mt_rk), Dimension(3,3)           :: c

    Integer :: ii,jj
    !==========================================================================

    Do jj = 1, 3
       Do ii = 1, 3

          c(ii,jj) = sum(base_m(:,ii)*base_i(:,jj))

       End do
    End Do

    Do ii = 1, 3

       a_t(ii) = sum(c(ii,:)*a)

    End Do

  End function t_transform_1

  !============================================================================
  !==
  !==
  function t_transform_2(a, base_i, base_j, base_m, base_n) Result(a_t)

    !** Input tensor**************************************
    Real(kind=mt_rk), Dimension(3,3), Intent(in)   :: a
    !** Original bases ***********************************
    Real(kind=mt_rk), Dimension(3,3), Intent(in) :: base_i, base_j
    !** Base to transform a to ***************************
    Real(kind=mt_rk), Dimension(3,3), Intent(in) :: base_m, base_n

    !** Output tensor ************************************
    Real(kind=mt_rk), Dimension(3,3)             :: a_t

    !*****************************************************

    !** Transformation tensor ****************************
    Real(kind=mt_rk), Dimension(3,3)             :: c_mi, c_nj

    Integer :: ii,jj,mm,nn

    !==========================================================================
!!$    write(*,*)
!!$    write(*,*)' Input Parameters ============================================='
!!$    write(*,*)'base_i = '
!!$    write(*,'(3E15.5)')base_i
!!$    write(*,*)
!!$    write(*,*)'base_j = '
!!$    write(*,'(3E15.5)')base_j
!!$    write(*,*)
!!$    write(*,*)'base_m = '
!!$    write(*,'(3E15.5)')base_m
!!$    write(*,*)
!!$    write(*,*)'base_n = '
!!$    write(*,'(3E15.5)')base_n
!!$    write(*,*)
!!$    write(*,*)
!!$    write(*,*)'a = '
!!$    write(*,'(3E15.5)')a
!!$    write(*,*)

    Do jj = 1, 3
       Do ii = 1, 3

          c_mi(ii,jj) = sum(base_m(:,ii)*base_i(:,jj))
          c_nj(ii,jj) = sum(base_n(:,ii)*base_j(:,jj))

       End do
    End Do

!!$    write(*,*)
!!$    write(*,*)'c_mi = '
!!$    write(*,'(3E15.5)')c_mi
!!$    write(*,*)
!!$    write(*,*)
!!$    write(*,*)'c_nj = '
!!$    write(*,'(3E15.5)')c_nj
!!$    write(*,*)

    a_t = 0.

    Do nn = 1, 3
       Do mm = 1, 3

          Do jj = 1, 3
             Do ii = 1, 3
                
                a_t(mm,nn) =  a_t(mm,nn) + a(ii,jj)*c_mi(mm,ii)*c_nj(nn,jj)

             End Do
          End Do

       End Do
    End Do

!!$    write(*,*)'a_t = '
!!$    write(*,'(3E15.5)')a_t
!!$    write(*,*)

  End function t_transform_2

  !============================================================================
  !==
  !== Transforms the 4 bases of a 3D rank 4 tensor simultaneously from org_base 
  !== to target base
  !==
  function t_transform_4(a, org_base, target_base) Result(a_t)

    !** Input tensor
    Real(kind=mt_rk), Dimension(3,3,3,3), Intent(in)   :: a

    !** Original base 
    Real(kind=mt_rk), Dimension(3,3), Intent(in)       :: org_base

    !** Base to transform a to 
    Real(kind=mt_rk), Dimension(3,3), Intent(in)       :: target_base

    !** Output tensor
    Real(kind=mt_rk), Dimension(3,3,3,3)               :: a_t

    !**************************************************************************

    !** Transformation tensors
    Real(kind=mt_rk), Dimension(3,3)                   :: c_mi, c_nj, c_ok, c_pl

    Integer :: ii,jj,kk,ll,mm,nn,oo,pp

    !==========================================================================

    Do jj = 1, 3
       Do ii = 1, 3

          c_mi(ii,jj) = target_base(1,ii)*org_base(1,jj)+&
                        target_base(2,ii)*org_base(2,jj)+&
                        target_base(3,ii)*org_base(3,jj)
       End do
    End Do

    c_nj = c_mi
    c_ok = c_mi
    c_pl = c_mi

    a_t = 0._mt_rk

    Do pp = 1, 3
       Do oo = 1, 3
          Do nn = 1, 3
             Do mm = 1, 3

                Do ll = 1, 3
                   Do kk = 1, 3               
                      Do jj = 1, 3
                         Do ii = 1, 3
                
                            a_t(mm,nn,oo,pp) =  a_t(mm,nn,oo,pp) + &
                                 a(ii,jj,kk,ll)*c_pl(pp,ll)*c_ok(oo,kk)*&
                                                c_nj(nn,jj)*c_mi(mm,ii)

                         end Do
                      end Do
                   End Do
                End Do
                
             end Do
          end Do
       End Do
    End Do

  End function t_transform_4

  Function matrix66_to_tensor_rank4(EE) Result(EE_t)

    Real(Kind=mt_rk), Dimension(6,6)    , Intent(in) :: EE
    Real(Kind=mt_rk), Dimension(3,3,3,3)             :: EE_t

    EE_t = 0._mt_rk

    !-- 1,j ---------------
    EE_t(1,1,1,1) = EE(1,1)
    EE_t(1,1,2,2) = EE(1,2)
    EE_t(1,1,3,3) = EE(1,3)

    EE_t(1,1,1,2) = EE(1,4)
    EE_t(1,1,2,1) = EE(1,4)

    EE_t(1,1,1,3) = EE(1,5)
    EE_t(1,1,3,1) = EE(1,5)

    EE_t(1,1,2,3) = EE(1,6)
    EE_t(1,1,3,2) = EE(1,6)

    !-- 2,j ---------------
    EE_t(2,2,1,1) = EE(2,1)
    EE_t(2,2,2,2) = EE(2,2)
    EE_t(2,2,3,3) = EE(2,3)

    EE_t(2,2,1,2) = EE(2,4)
    EE_t(2,2,2,1) = EE(2,4)

    EE_t(2,2,1,3) = EE(2,5)
    EE_t(2,2,3,1) = EE(2,5)

    EE_t(2,2,2,3) = EE(2,6)
    EE_t(2,2,3,2) = EE(2,6)

    !-- 3,j ---------------
    EE_t(3,3,1,1) = EE(3,1)
    EE_t(3,3,2,2) = EE(3,2)
    EE_t(3,3,3,3) = EE(3,3)

    EE_t(3,3,1,2) = EE(3,4)
    EE_t(3,3,2,1) = EE(3,4)

    EE_t(3,3,1,3) = EE(3,5)
    EE_t(3,3,3,1) = EE(3,5)

    EE_t(3,3,2,3) = EE(3,6)
    EE_t(3,3,3,2) = EE(3,6)

    !-- 4,j ---------------
    EE_t(1,2,1,1) = EE(4,1)
    EE_t(2,1,1,1) = EE(4,1)

    EE_t(1,2,2,2) = EE(4,2)
    EE_t(2,1,2,2) = EE(4,2)

    EE_t(1,2,3,3) = EE(4,3)
    EE_t(2,1,3,3) = EE(4,3)

    EE_t(1,2,1,2) = EE(4,4)
    EE_t(2,1,1,2) = EE(4,4)
    EE_t(1,2,2,1) = EE(4,4)
    EE_t(2,1,2,1) = EE(4,4)

    EE_t(1,2,1,3) = EE(4,5)
    EE_t(2,1,1,3) = EE(4,5)
    EE_t(1,2,3,1) = EE(4,5)
    EE_t(2,1,3,1) = EE(4,5)

    EE_t(1,2,2,3) = EE(4,6)
    EE_t(2,1,2,3) = EE(4,6)
    EE_t(1,2,3,2) = EE(4,6)
    EE_t(2,1,3,2) = EE(4,6)

    !-- 5,j ---------------
    EE_t(1,3,1,1) = EE(5,1)
    EE_t(3,1,1,1) = EE(5,1)

    EE_t(1,3,2,2) = EE(5,2)
    EE_t(3,1,2,2) = EE(5,2)

    EE_t(1,3,3,3) = EE(5,3)
    EE_t(3,1,3,3) = EE(5,3)

    EE_t(1,3,1,2) = EE(5,4)
    EE_t(3,1,1,2) = EE(5,4)
    EE_t(1,3,2,1) = EE(5,4)
    EE_t(3,1,2,1) = EE(5,4)

    EE_t(1,3,1,3) = EE(5,5)
    EE_t(3,1,1,3) = EE(5,5)
    EE_t(1,3,3,1) = EE(5,5)
    EE_t(3,1,3,1) = EE(5,5)

    EE_t(1,3,2,3) = EE(5,6)
    EE_t(3,1,2,3) = EE(5,6)
    EE_t(1,3,3,2) = EE(5,6)
    EE_t(3,1,3,2) = EE(5,6)

    !-- 6,j ---------------
    EE_t(2,3,1,1) = EE(6,1)
    EE_t(3,2,1,1) = EE(6,1)

    EE_t(2,3,2,2) = EE(6,2)
    EE_t(3,2,2,2) = EE(6,2)

    EE_t(2,3,3,3) = EE(6,3)
    EE_t(3,2,3,3) = EE(6,3)

    EE_t(2,3,1,2) = EE(6,4)
    EE_t(3,2,1,2) = EE(6,4)
    EE_t(2,3,2,1) = EE(6,4)
    EE_t(3,2,2,1) = EE(6,4)

    EE_t(2,3,1,3) = EE(6,5)
    EE_t(3,2,1,3) = EE(6,5)
    EE_t(2,3,3,1) = EE(6,5)
    EE_t(3,2,3,1) = EE(6,5)

    EE_t(2,3,2,3) = EE(6,6)
    EE_t(3,2,2,3) = EE(6,6)
    EE_t(2,3,3,2) = EE(6,6)
    EE_t(3,2,3,2) = EE(6,6)

  End Function matrix66_to_tensor_rank4

  Function tensor_rank4_to_matrix66(EE_t) Result(EE)
  
    Real(Kind=mt_rk), Dimension(6,6)                 :: EE
    Real(Kind=mt_rk), Dimension(3,3,3,3), Intent(in) :: EE_t

    EE = 0._mt_rk

    !-- 1,j ---------------
     EE(1,1) = EE_t(1,1,1,1)
     EE(1,2) = EE_t(1,1,2,2)
     EE(1,3) = EE_t(1,1,3,3)
            
     EE(1,4) = EE_t(1,1,1,2)
     EE(1,4) = EE_t(1,1,2,1)
            
     EE(1,5) = EE_t(1,1,1,3)
     EE(1,5) = EE_t(1,1,3,1)
            
     EE(1,6) = EE_t(1,1,2,3)
     EE(1,6) = EE_t(1,1,3,2)

    !-- 2,j ---------------
     EE(2,1) = EE_t(2,2,1,1)
     EE(2,2) = EE_t(2,2,2,2)
     EE(2,3) = EE_t(2,2,3,3)
            
     EE(2,4) = EE_t(2,2,1,2)
     EE(2,4) = EE_t(2,2,2,1)
            
     EE(2,5) = EE_t(2,2,1,3)
     EE(2,5) = EE_t(2,2,3,1)
            
     EE(2,6) = EE_t(2,2,2,3)
     EE(2,6) = EE_t(2,2,3,2)

    !-- 3,j ---------------
     EE(3,1) = EE_t(3,3,1,1)
     EE(3,2) = EE_t(3,3,2,2)
     EE(3,3) = EE_t(3,3,3,3)
            
     EE(3,4) = EE_t(3,3,1,2)
     EE(3,4) = EE_t(3,3,2,1)
            
     EE(3,5) = EE_t(3,3,1,3)
     EE(3,5) = EE_t(3,3,3,1)
            
     EE(3,6) = EE_t(3,3,2,3)
     EE(3,6) = EE_t(3,3,3,2)

    !-- 4,j ---------------
     EE(4,1) = EE_t(1,2,1,1)
     EE(4,1) = EE_t(2,1,1,1)
            
     EE(4,2) = EE_t(1,2,2,2)
     EE(4,2) = EE_t(2,1,2,2)
            
     EE(4,3) = EE_t(1,2,3,3)
     EE(4,3) = EE_t(2,1,3,3)
            
     EE(4,4) = EE_t(1,2,1,2)
     EE(4,4) = EE_t(2,1,1,2)
     EE(4,4) = EE_t(1,2,2,1)
     EE(4,4) = EE_t(2,1,2,1)
            
     EE(4,5) = EE_t(1,2,1,3)
     EE(4,5) = EE_t(2,1,1,3)
     EE(4,5) = EE_t(1,2,3,1)
     EE(4,5) = EE_t(2,1,3,1)
            
     EE(4,6) = EE_t(1,2,2,3)
     EE(4,6) = EE_t(2,1,2,3)
     EE(4,6) = EE_t(1,2,3,2)
     EE(4,6) = EE_t(2,1,3,2)

    !-- 5,j ---------------
     EE(5,1) = EE_t(1,3,1,1)
     EE(5,1) = EE_t(3,1,1,1)
            
     EE(5,2) = EE_t(1,3,2,2)
     EE(5,2) = EE_t(3,1,2,2)
            
     EE(5,3) = EE_t(1,3,3,3)
     EE(5,3) = EE_t(3,1,3,3)
            
     EE(5,4) = EE_t(1,3,1,2)
     EE(5,4) = EE_t(3,1,1,2)
     EE(5,4) = EE_t(1,3,2,1)
     EE(5,4) = EE_t(3,1,2,1)
            
     EE(5,5) = EE_t(1,3,1,3)
     EE(5,5) = EE_t(3,1,1,3)
     EE(5,5) = EE_t(1,3,3,1)
     EE(5,5) = EE_t(3,1,3,1)
            
     EE(5,6) = EE_t(1,3,2,3)
     EE(5,6) = EE_t(3,1,2,3)
     EE(5,6) = EE_t(1,3,3,2)
     EE(5,6) = EE_t(3,1,3,2)

    !-- 6,j ---------------
     EE(6,1) = EE_t(2,3,1,1)
     EE(6,1) = EE_t(3,2,1,1)
            
     EE(6,2) = EE_t(2,3,2,2)
     EE(6,2) = EE_t(3,2,2,2)
            
     EE(6,3) = EE_t(2,3,3,3)
     EE(6,3) = EE_t(3,2,3,3)
            
     EE(6,4) = EE_t(2,3,1,2)
     EE(6,4) = EE_t(3,2,1,2)
     EE(6,4) = EE_t(2,3,2,1)
     EE(6,4) = EE_t(3,2,2,1)
            
     EE(6,5) = EE_t(2,3,1,3)
     EE(6,5) = EE_t(3,2,1,3)
     EE(6,5) = EE_t(2,3,3,1)
     EE(6,5) = EE_t(3,2,3,1)
            
     EE(6,6) = EE_t(2,3,2,3)
     EE(6,6) = EE_t(3,2,2,3)
     EE(6,6) = EE_t(2,3,3,2)
     EE(6,6) = EE_t(3,2,3,2)

   End Function tensor_rank4_to_matrix66

  Subroutine monitor_rank4_tensor(cc,tname,un_lf)

    Character(Len=*), Intent(In) :: tname

    real(kind=mt_rk), Dimension(3,3,3,3), intent(In) :: cc

    Integer , Intent(in) :: un_lf

    Integer :: ii,jj,kk,ll

    !----------------------------------------------------------------------------

    Write(un_lf,"('<',77('='),'>')")
    Write(un_lf,*)trim(tname)

    Write(un_lf,"(8X,3(32('-'),'+'))")
    Write(un_lf,"(' r/s =  ',3(10('-'),I6,6X,10('-'),'|'))")1,2,3

    Do kk = 1, 3

       Write(un_lf,"(I2,5X,'+',3(32('-'),'+'))")kk

       Write(un_lf,"('   k/l |',3(2X,3(2('-'),I3,2X,2('-'),1X),'|'))")1,2,3,1,2,3,1,2,3
       Do ii = 1, 3

          Write(un_lf,"(2X,I2,'   | ',$)")ii

          Do ll = 1, 3

             Do jj = 1, 3

                Write(un_lf,'(F10.2,$)') CC(ii,jj,kk,ll)

             End Do

             Write(un_lf,"(' | ',$)")

          End Do

          Write(un_lf,*)

       End Do

    End Do

    Write(un_lf,"(7X,'+'3(32('-'),'+'))")
    write(un_lf,*)

  End Subroutine monitor_rank4_tensor

End module tensors
