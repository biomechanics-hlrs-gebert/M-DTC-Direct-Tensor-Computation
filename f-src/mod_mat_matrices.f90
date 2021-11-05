!******************************************************************************
!*
!* Module to generate material matrices
!*
!******************************************************************************
Module mat_matrices

  Use linpack
  use math

  implicit none

  Type tPs_ortho

     Real(kind=rk)  ::     E1,E2,E3
     Real(kind=rk)  ::     v12,v13,v23
     Real(kind=rk)  ::     G12,G13,G23

  End Type tPs_ortho

  Type tMat_Tr

     !> Elasticity Matrix ---------------------------
     Real(Kind=rk), Dimension(6,6,8) :: EE
     !> Cos(Alpha) ----------------------------------
     Real(Kind=rk)                 :: cos_a
     !> Alpha ---------------------------------------
     Real(Kind=rk)                 :: alpha

     !> Rotation vector components squared ----------
     Real(Kind=rk), Dimension(3,2) :: nc2
     !> Rotation vector components ------------------
     Real(Kind=rk), Dimension(3,8) :: nc   
     !> Absolute value of rotation vector squared ---
     Real(Kind=rk), dimension(8)   :: sum_nc2

     !> switch for valid cos(alpha) -----------------
     Logical                       :: ang_valid
     !> switches for valid normal vectors -----------
     Logical, Dimension(8)         :: vec_valid

     Real(Kind=rk)                 :: cond
     Real(Kind=rk), Dimension(4,8) :: crit
     
  End type tMat_Tr

Contains 

  Subroutine write_eval_evec(un,eval,evali,evec)

    Integer                      , Intent(in) :: un
    Real(Kind=rk), Dimension(:)  , Intent(in) :: eval,evali
    Real(Kind=rk), Dimension(:,:), Intent(in) :: evec

    Character(Len=256) :: fmt

    Real(Kind=rk)                             :: sum_evi

!    Write(fmt,"(A,)")

    fmt="(6('|',1X,E13.6,1X),'|')"
    sum_evi=sum(abs(evali))

    Write(un,"('+',95('-') ,'+')")
    Write(un,"('+',1X,A,T97,'+')")"Eigenvalues with corresponding Eigenvectors"
    Write(un,"(6('+',15('-')),'+',$)")
    If (sum_evi > num_zero) then
       Write(un,"(6('+',15('-')),'+')")
    Else
       Write(un,*)
    End If

    Write(un,"(6('|',1X,E13.6,1X),'|',$)")eval

    If (sum_evi > num_zero) then
       Write(un,"(6(1X,E13.6,1X,'|'))")eval
    Else
       Write(un,*)
    End If

    Write(un,"(6('+',15('-')),'+',$)")
    If (sum_evi > num_zero) then
       Write(un,"(6('+',15('-')),'+')")
    Else
       Write(un,*)
    End If

!    Write(un,"(6('|',1X,E13.6,1X),'|')")transpose(evec)
    Write(un,"(6('[',1X,E13.6,5(' , ',E13.6),']',/))")transpose(evec)
!    Write(un,"(6('+',15('-')),'+')")
!    Write(un,*)

    If (sum_evi <= num_zero) then
       Write(un,*)"There are no eigenvalues with imaginary parts"
       Write(un,*)
    End If

  End Subroutine write_eval_evec

  Subroutine order_eval_evec(dim,wr,wi,zz)

    Integer(kind=ik),                     Intent(in)   :: dim
    Real(kind=rk)   , Dimension(dim)    , Intent(inOut):: wr, wi
    Real(kind=rk)   , Dimension(dim,dim), Intent(inOut):: zz

    !**************************************************************************
    !* Variables declaration    
    Logical                                  :: unsorted
    Integer                                  :: ii
    Real(Kind=rk)                            :: tmp_r
    Real(kind=rk), Dimension(:), allocatable :: tmp_rdim

    !==========================================================================

    allocate(tmp_rdim(dim))

    unsorted = .True.

    Do While (unsorted)

       unsorted = .false.

       Do ii = 2, dim
          
          If (wr(ii-1) < wr(ii)) Then

             tmp_r    = wr(ii-1)
             wr(ii-1) = wr(ii)
             wr(ii)   = tmp_r

             tmp_r    = wi(ii-1)
             wi(ii-1) = wi(ii)
             wi(ii)   = tmp_r

             tmp_rdim   = zz(:,ii-1)
             zz(:,ii-1) = zz(:,ii)
             zz(:,ii)   = tmp_rdim

             unsorted = .True.

          End IF

       End Do

    End Do

    deallocate(tmp_rdim)

  End Subroutine order_eval_evec

  Function write_disk_case(un_lf,disk,zzi,r) Result(prod)

    Real(kind=rk), dimension(3)  , intent(in) :: r
    Integer                      , intent(in) :: un_lf
    real(kind=rk), dimension(3)  , intent(in) :: disk
    real(kind=rk), dimension(6)  , intent(in) :: zzi
    
    Real(kind=rk)                             :: prod

    If ( ( disk(1) > num_zero ) .AND. ( disk(2) > num_zero ) .AND. &
         ( disk(3) > num_zero )) Then
       
       If ( abs(sq2**3*zzi(1)*zzi(2)*zzi(3) - &
                ((zzi(4) + sqrt(disk(1))*r(1)) * &
                 (zzi(5) + sqrt(disk(2))*r(2)) * &
                 (zzi(6) + sqrt(disk(3))*r(3))    )) > num_zero ) Then

          prod =  sq2**3*zzi(1)*zzi(2)*zzi(3) - &
               ((zzi(4) + sqrt(disk(1))*r(1)) * &
                (zzi(5) + sqrt(disk(2))*r(2)) * &
                (zzi(6) + sqrt(disk(3))*r(3))    )

          Write(un_lf,"(E17.9)",Advance='No')prod
              

       Else

          prod = 0._rk
          Write(un_lf,"(A17)",Advance='No')" 0.0"

       End If

    else if ( (disk(1) <= num_zero) .AND. ( disk(2) > num_zero ) .AND. &
         ( disk(3) > num_zero ) ) then

       If ( abs(sq2**3*zzi(1)*zzi(2)*zzi(3) - &
                ((zzi(4)                     ) * &
                 (zzi(5) + sqrt(disk(2))*r(2)) * &
                 (zzi(6) + sqrt(disk(3))*r(3))    )) > num_zero ) Then
 
          prod = sq2**3*zzi(1)*zzi(2)*zzi(3) - &
               ((zzi(4)                     ) * &
                (zzi(5) + sqrt(disk(2))*r(2)) * &
                (zzi(6) + sqrt(disk(3))*r(3))    )

          Write(un_lf,"(E17.9)",Advance='No')prod

       Else
          prod = 0._rk
          Write(un_lf,"(A17)",Advance='No')" 0.0"

       End If

    Else if ( (disk(1) > num_zero) .AND. ( disk(2) <= num_zero ) .AND. &
         ( disk(3) > num_zero ) ) then

       If ( abs(sq2**3*zzi(1)*zzi(2)*zzi(3) - &
                ((zzi(4) + sqrt(disk(1))*r(1)) * &
                 (zzi(5)                     ) * &
                 (zzi(6) + sqrt(disk(3))*r(3))    )) > num_zero ) Then
 
          prod = sq2**3*zzi(1)*zzi(2)*zzi(3) - &
               ((zzi(4) + sqrt(disk(1))*r(1)) * &
                (zzi(5)                     ) * &
                (zzi(6) + sqrt(disk(3))*r(3))    )

          Write(un_lf,"(E17.9)",Advance='No')prod

       Else

          prod = 0._rk
          Write(un_lf,"(A17)",Advance='No')" 0.0"

       End If

    Else if ( (disk(1) > num_zero) .AND. ( disk(2) > num_zero ) .AND. &
         ( disk(3) <= num_zero ) ) then

       If ( abs(sq2**3*zzi(1)*zzi(2)*zzi(3) - &
                ((zzi(4) + sqrt(disk(1))*r(1)) * &
                 (zzi(5) + sqrt(disk(2))*r(2)) * &
                 (zzi(6)                     )    )) > num_zero ) Then
 
          prod = sq2**3*zzi(1)*zzi(2)*zzi(3) - &
               ((zzi(4) + sqrt(disk(1))*r(1)) * &
                (zzi(5) + sqrt(disk(2))*r(2)) * &
                (zzi(6)                     )    )

          Write(un_lf,"(E17.9)",Advance='No')prod

       Else

          prod = 0._rk
          Write(un_lf,"(A17)",Advance='No')" 0.0"

       End If

    else If ( ( disk(1) <= num_zero ) .AND. ( disk(2) <= num_zero ) .AND. &
         ( disk(3) > num_zero )) Then
       
       If ( abs(sq2**3*zzi(1)*zzi(2)*zzi(3) - &
                ((zzi(4)                     ) * &
                 (zzi(5)                     ) * &
                 (zzi(6) + sqrt(disk(3))*r(3))    )) > num_zero ) Then
 
          prod = sq2**3*zzi(1)*zzi(2)*zzi(3) - &
               ((zzi(4)                     ) * &
                (zzi(5)                     ) * &
                (zzi(6) + sqrt(disk(3))*r(3))    )

          Write(un_lf,"(E17.9)",Advance='No')prod

       Else

          prod = 0._rk
          Write(un_lf,"(A17)",Advance='No')" 0.0"

       End If

    else If ( ( disk(1) <= num_zero ) .AND. ( disk(2) > num_zero ) .AND. &
         ( disk(3) <= num_zero )) Then
       
       If ( abs(sq2**3*zzi(1)*zzi(2)*zzi(3) - &
                ((zzi(4)                     ) * &
                 (zzi(5) + sqrt(disk(2))*r(2)) * &
                 (zzi(6)                     )    )) > num_zero ) Then
 
          prod = sq2**3*zzi(1)*zzi(2)*zzi(3) - &
               ((zzi(4)                     ) * &
                (zzi(5) + sqrt(disk(2))*r(2)) * &
                (zzi(6)                     )    )
          Write(un_lf,"(E17.9)",Advance='No')prod

       Else

          prod = 0._rk
          Write(un_lf,"(A17)",Advance='No')" 0.0"

       End If

    else If ( ( disk(1) > num_zero ) .AND. ( disk(2) <= num_zero ) .AND. &
         ( disk(3) <= num_zero )) Then
       
       If ( abs(sq2**3*zzi(1)*zzi(2)*zzi(3) - &
                ((zzi(4) + sqrt(disk(1))*r(1)) * &
                 (zzi(5)                     ) * &
                 (zzi(6)                     )    )) > num_zero ) Then
 
          prod = sq2**3*zzi(1)*zzi(2)*zzi(3) - &
               ((zzi(4) + sqrt(disk(1))*r(1)) * &
                (zzi(5)                     ) * &
                (zzi(6)                     )    )
          Write(un_lf,"(E17.9)",Advance='No')prod

       Else

          prod = 0._rk
          Write(un_lf,"(A17)",Advance='No')" 0.0"

       End If

    else If ( ( disk(1) <= num_zero ) .AND. ( disk(2) <= num_zero ) .AND. &
         ( disk(3) <= num_zero )) Then
       
       If ( abs(sq2**3*zzi(1)*zzi(2)*zzi(3) - &
                ((zzi(4)                     ) * &
                 (zzi(5)                     ) * &
                 (zzi(6)                     )    )) > num_zero ) Then
 
          prod = sq2**3*zzi(1)*zzi(2)*zzi(3) - &
               ((zzi(4)                     ) * &
                (zzi(5)                     ) * &
                (zzi(6)                     )    )
          Write(un_lf,"(E17.9)",Advance='No')prod

       Else

          prod = 0._rk
          Write(un_lf,"(A17)",Advance='No')" 0.0"

       End If

    End If
    
  end Function write_disk_case

  Function write_disk_case_zero(un_lf,disk,zzi,r) Result(prod)

    Real(kind=rk), dimension(3)  , intent(in) :: r
    Integer                      , intent(in) :: un_lf
    real(kind=rk), dimension(3)  , intent(in) :: disk
    real(kind=rk), dimension(6)  , intent(in) :: zzi
       
    Real(kind=rk)                             :: prod

    If ( ((zzi(4) + sqrt(disk(1))*r(1)) * &
          (zzi(5) + sqrt(disk(2))*r(2)) * &
          (zzi(6) + sqrt(disk(3))*r(3))    ) > num_zero ) Then
 
       prod = ((zzi(4) + sqrt(disk(1))*r(1)) * &
            (zzi(5) + sqrt(disk(2))*r(2)) * &
            (zzi(6) + sqrt(disk(3))*r(3))    )
       Write(un_lf,"(E17.9)",Advance='No')prod

    Else

       prod = 0._rk
       Write(un_lf,"(A17)",Advance='No')" 0.0"
       
    End If
    
  end Function write_disk_case_zero

End Module mat_matrices
