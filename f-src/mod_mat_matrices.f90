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

  !****************************************************************************
  !*
  !* Function to generate an isotropic stiffness matrix
  !*
  Function Matrix_Iso(E,v) Result (M_iso)
    
    Implicit None
    Real(kind=rk)                   ::     v, E
    Real(kind=rk), Dimension(6,6)   ::     M_iso
    
    M_Iso(:,1)=(/ 1._rk/E, -v/E, -v/E, .0_rk, .0_rk, .0_rk /)
    M_Iso(:,2)=(/ -v/E, 1._rk/E, -v/E, .0_rk, .0_rk, .0_rk /)
    M_Iso(:,3)=(/ -v/E, -v/E, 1._rk/E, .0_rk, .0_rk, .0_rk /)
    M_Iso(:,4)=(/ .0_rk, .0_rk, .0_rk, 2._rk*(1._rk+v)/E, .0_rk, .0_rk /)
    M_Iso(:,5)=(/ .0_rk, .0_rk, .0_rk, .0_rk, 2._rk*(1._rk+v)/E, .0_rk /)
    M_Iso(:,6)=(/ .0_rk, .0_rk, .0_rk, .0_rk, .0_rk, 2._rk*(1._rk+v)/E /)

  End Function Matrix_Iso

  !****************************************************************************
  !*
  !* Function to generate an isotropic stiffness matrix
  !*
  Function Matrix_Ortho(ps) Result (M_Ortho)

    Implicit None
    
    Type(tPs_ortho), Intent(in)     :: ps
    Real(kind=rk), Dimension(6,6)   :: M_Ortho
 
    M_Ortho(:,1)=(/   1._rk/ps%E1, -ps%v12/ps%E1, -ps%v13/ps%E1, &
         .0_rk, .0_rk, .0_rk /)
    M_Ortho(:,2)=(/ -ps%v12/ps%E1,   1._rk/ps%E2, -ps%v23/ps%E2, &
         .0_rk, .0_rk, .0_rk /)
    M_Ortho(:,3)=(/ -ps%v13/ps%E1, -ps%v23/ps%E2,   1._rk/ps%E3, &
         .0_rk, .0_rk, .0_rk /)
    M_Ortho(:,4)=(/ .0_rk, .0_rk, .0_rk, 1._rk/ps%G12, .0_rk, .0_rk /)
    M_Ortho(:,5)=(/ .0_rk, .0_rk, .0_rk, .0_rk, 1._rk/ps%G13, .0_rk /)
    M_Ortho(:,6)=(/ .0_rk, .0_rk, .0_rk, .0_rk, .0_rk, 1._rk/ps%G23 /)

  End Function Matrix_Ortho

  Subroutine write_matrix(un,aa,label)

    Integer                      , Intent(in) :: un
    Real(Kind=rk), Dimension(:,:), Intent(in) :: aa
    Character(Len=*), intent(in), optional    :: label

    Integer, dimension(2)                     :: lb
    Integer                                   :: fmt_mult
    Character(Len=256) :: fmt,fmt_sep


    If (present(label)) then

       fmt_mult = size(aa(:,lb(2)))*16
       Write(fmt_sep,"(A,I0,A)")"('+',",fmt_mult-1,"('-') ,'+')"
       Write(un,fmt_sep)
       Write(fmt_sep,"(A,I0,A)")"('+',1X,A,T",fmt_mult+1,",'+')"
       Write(un,fmt_sep)trim(label)

    end If

    lb = lbound(aa)

    Write(fmt_sep,"(A,I0,A)")"(",size(aa(:,lb(2))),"('+',15('-'))    ,'+')"
    Write(fmt,"(A,I0,A)")    "(",size(aa(:,lb(2))),"('|',1X,F13.3,1X),'|')"

    Write(un,trim(fmt_sep))
    Write(un,fmt)transpose(aa)
    Write(un,fmt_sep)
    Write(un,*)

  End Subroutine write_matrix

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

  Function write_disk_case(un_lf,disk,zzi,case,r) Result(prod)

    Real(kind=rk), dimension(3)  , intent(in) :: r
    Integer                      , intent(in) :: un_lf
    real(kind=rk), dimension(3)  , intent(in) :: disk
    real(kind=rk), dimension(6)  , intent(in) :: zzi
    Character(Len=*)             , intent(in) :: case
    
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
    
1000 Continue

  end Function write_disk_case

  Function write_disk_case_zero(un_lf,disk,zzi,case,r) Result(prod)

    Real(kind=rk), dimension(3)  , intent(in) :: r
    Integer                      , intent(in) :: un_lf
    real(kind=rk), dimension(3)  , intent(in) :: disk
    real(kind=rk), dimension(6)  , intent(in) :: zzi
    Character(Len=*)             , intent(in) :: case
       
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
