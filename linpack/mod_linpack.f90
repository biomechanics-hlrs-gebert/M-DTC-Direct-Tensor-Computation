Module LinPack

  USE global_std
  Use blas_1

  Implicit None

Contains

  Subroutine inverse(a, nn, un_lf)

    !**************************************************************************
    !* Parameter declaration  
    !*
    !* on entry
    !* -----------
    !*
    !*  a       Real(Kind=rk), Dimension = (nn,nn)
    !*          the matrix to be inverted.
    !*
    !*  nn      integer
    !*          the dimension of the matrix  a .
    !*
    !* on return
    !* -----------
    !*
    !*        a       inverse of original matrix
    !*
    !**************************************************************************

    Integer,                         Intent(In)    :: nn, un_lf
    Real(Kind=rk), Dimension(nn,nn), Intent(InOut) :: a

    !**************************************************************************
    !* Variables declaration      
    Integer,       Dimension(nn) :: ipvt
    Real(Kind=rk), Dimension(nn) :: z
    Real(Kind=rk), Dimension(2)  :: det
    Real(Kind=rk)                :: rcond
    

    !* Format Strings -------------------------------
    Character(Len=*), Parameter :: fmt_str_F    = '(A,E20.9,/)'
    !*-----------------------------------------------

    !write(    *      ,fmt_str_AdNo)'Factorisation of matrix .... '
    !write(un_lf,fmt_str_AdNo)'Factorisation of matrix .... '
    call dgeco(a, nn, nn, ipvt , rcond, z)
    !write(    *      ,fmt_str)' done !'
    !write(un_lf,fmt_str)' done !'
    
    !write(    *      ,fmt_str_F)'Estimated condition of  a = ',1._rk/rcond
    write(un_lf,fmt_str_F)'Estimated condition of  a = ',1._rk/rcond

    !write(    *      ,fmt_str_AdNo)'Calculation of inverse(a) .... '
    !write(un_lf,fmt_str_AdNo)'Calculation of inverse(a) .... '
    call dgedi(a, nn, nn, ipvt, det, z, 01)
    !write(    *      ,fmt_str)' done !'
    !write(un_lf,fmt_str)' done !'

  End Subroutine inverse

!******************************************************************************
      Subroutine dgedi(a,lda,n,ipvt,det,work,job)

      Integer :: lda, n, job
      Integer :: ipvt(n)

      Real(Kind=rk) :: a(lda,n)
      Real(Kind=rk) :: det(2)

      Real(Kind=rk) :: work(n)
!
!     dgedi computes the determinant and inverse of a matrix
!     using the factors computed by dgeco or dgefa.
!
!     on entry
!
!        a       double precision(lda, n)
!                the output from dgeco or dgefa.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!        ipvt    integer(n)
!                the pivot vector from dgeco or dgefa.
!
!        work    double precision(n)
!                work vector.  contents destroyed.
!
!        job     integer
!                = 11   both determinant and inverse.
!                = 01   inverse only.
!                = 10   determinant only.
!
!     on return
!
!        a       inverse of original matrix if requested.
!                otherwise unchanged.
!
!        det     double precision(2)
!                determinant of original matrix if requested.
!                otherwise not referenced.
!                determinant = det(1) * 10.0**det(2)
!                with  1.0 .le. dabs(det(1)) .lt. 10.0
!                or  det(1) .eq. 0.0 .
!
!     error condition
!
!        a division by zero will occur if the input factor contains
!        a zero on the diagonal and the inverse is requested.
!        it will not occur if the subroutines are called correctly
!        and if dgeco has set rcond .gt. 0.0 or dgefa has set
!        info .eq. 0 .
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,dscal,dswap
!     fortran dabs,mod
!
!     internal variables
!
      Real(Kind=rk) :: t
      Real(Kind=rk) :: ten
      Integer i,j,k,kb,kp1,l,nm1
!
!
!     compute determinant
!
      If (job/10 .Eq. 0) go to 70
         det(1) = 1.0d0
         det(2) = 0.0d0
         ten = 10.0d0
         Do 50 i = 1, n
            If (ipvt(i) .Ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
!        ...exit
            If (det(1) .Eq. 0.0d0) go to 60
   10       If (dabs(det(1)) .Ge. 1.0d0) go to 20
               det(1) = ten*det(1)
               det(2) = det(2) - 1.0d0
            go to 10
   20       Continue
   30       If (dabs(det(1)) .Lt. ten) go to 40
               det(1) = det(1)/ten
               det(2) = det(2) + 1.0d0
            go to 30
   40       Continue
   50    Continue
   60    Continue
   70 Continue
!
!     compute inverse(u)
!
      If (Mod(job,10) .Eq. 0) go to 150
         Do 100 k = 1, n
            a(k,k) = 1.0d0/a(k,k)
            t = -a(k,k)
            Call dscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            If (n .Lt. kp1) go to 90
            Do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0d0
               Call daxpy(k,t,a(1,k),1,a(1,j),1)
   80       Continue
   90       Continue
  100    Continue
!
!        form inverse(u)*inverse(l)
!
         nm1 = n - 1
         If (nm1 .Lt. 1) go to 140
         Do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            Do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = 0.0d0
  110       Continue
            Do 120 j = kp1, n
               t = work(j)
               Call daxpy(n,t,a(1,j),1,a(1,k),1)
  120       Continue
            l = ipvt(k)
            If (l .Ne. k) Call dswap(n,a(1,k),1,a(1,l),1)
  130    Continue
  140    Continue
  150 Continue
      Return
      End Subroutine dgedi

!***************************************************************************************************
      Subroutine dgeco(a,lda,n,ipvt,rcond,z)

      Integer       :: lda,n
      Integer       :: ipvt(n)
      Real(Kind=rk) :: a(lda,n),z(n)
      Real(Kind=rk) :: rcond
!
!     dgeco factors a double precision matrix by gaussian elimination
!     and estimates the condition of the matrix.
!
!     if  rcond  is not needed, dgefa is slightly faster.
!     to solve  a*x = b , follow dgeco by dgesl.
!     to compute  inverse(a)*c , follow dgeco by dgesl.
!     to compute  determinant(a) , follow dgeco by dgedi.
!     to compute  inverse(a) , follow dgeco by dgedi.
!
!     on entry
!
!        a       double precision(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        rcond   double precision
!                an estimate of the reciprocal condition of  a .
!                for the system  a*x = b , relative perturbations
!                in  a  and  b  of size  epsilon  may cause
!                relative perturbations in  x  of size  epsilon/rcond .
!                if  rcond  is so small that the logical expression
!                           1.0 + rcond .eq. 1.0
!                is true, then  a  may be singular to working
!                precision.  in particular,  rcond  is zero  if
!                exact singularity is detected or the estimate
!                underflows.
!
!        z       double precision(n)
!                a work vector whose contents are usually unimportant.
!                if  a  is close to a singular matrix, then  z  is
!                an approximate null vector in the sense that
!                norm(a*z) = rcond*norm(a)*norm(z) .
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     linpack dgefa
!     blas daxpy,ddot,dscal,dasum
!     fortran dabs,dmax1,dsign
!
!     internal variables
!
!Kuester      double precision ddot,ek,t,wk,wkm
      Real(Kind=rk) :: ek,t,wk,wkm
!Kuester      double precision anorm,s,dasum,sm,ynorm
      Real(Kind=rk) :: anorm,s,sm,ynorm
      Integer info,j,k,kb,kp1,l
!
!
!     compute 1-norm of a
!
      anorm = 0.0d0
      Do 10 j = 1, n
         anorm = dmax1(anorm,dasum(n,a(1,j),1))
   10 Continue
!
!     factor
!
      Call dgefa(a,lda,n,ipvt,info)
!
!     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
!     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
!     trans(a)  is the transpose of a .  the components of  e  are
!     chosen to cause maximum local growth in the elements of w  where
!     trans(u)*w = e .  the vectors are frequently rescaled to avoid
!     overflow.
!
!     solve trans(u)*w = e
!
      ek = 1.0d0
      Do 20 j = 1, n
         z(j) = 0.0d0
   20 Continue
      Do 100 k = 1, n
         If (z(k) .Ne. 0.0d0) ek = dsign(ek,-z(k))
         If (dabs(ek-z(k)) .Le. dabs(a(k,k))) go to 30
            s = dabs(a(k,k))/dabs(ek-z(k))
            Call dscal(n,s,z,1)
            ek = s*ek
   30    Continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = dabs(wk)
         sm = dabs(wkm)
         If (a(k,k) .Eq. 0.0d0) go to 40
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
         go to 50
   40    Continue
            wk = 1.0d0
            wkm = 1.0d0
   50    Continue
         kp1 = k + 1
         If (kp1 .Gt. n) go to 90
            Do 60 j = kp1, n
               sm = sm + dabs(z(j)+wkm*a(k,j))
               z(j) = z(j) + wk*a(k,j)
               s = s + dabs(z(j))
   60       Continue
            If (s .Ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               Do 70 j = kp1, n
                  z(j) = z(j) + t*a(k,j)
   70          Continue
   80       Continue
   90    Continue
         z(k) = wk
  100 Continue
      s = 1.0d0/dasum(n,z,1)
      Call dscal(n,s,z,1)
!
!     solve trans(l)*y = w
!
      Do 120 kb = 1, n
         k = n + 1 - kb
         If (k .Lt. n) z(k) = z(k) + ddot(n-k,a(k+1,k),1,z(k+1),1)
         If (dabs(z(k)) .Le. 1.0d0) go to 110
            s = 1.0d0/dabs(z(k))
            Call dscal(n,s,z,1)
  110    Continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 Continue
      s = 1.0d0/dasum(n,z,1)
      Call dscal(n,s,z,1)
!
      ynorm = 1.0d0
!
!     solve l*v = y
!
      Do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         If (k .Lt. n) Call daxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         If (dabs(z(k)) .Le. 1.0d0) go to 130
            s = 1.0d0/dabs(z(k))
            Call dscal(n,s,z,1)
            ynorm = s*ynorm
  130    Continue
  140 Continue
      s = 1.0d0/dasum(n,z,1)
      Call dscal(n,s,z,1)
      ynorm = s*ynorm
!
!     solve  u*z = v
!
      Do 160 kb = 1, n
         k = n + 1 - kb
         If (dabs(z(k)) .Le. dabs(a(k,k))) go to 150
            s = dabs(a(k,k))/dabs(z(k))
            Call dscal(n,s,z,1)
            ynorm = s*ynorm
  150    Continue
         If (a(k,k) .Ne. 0.0d0) z(k) = z(k)/a(k,k)
         If (a(k,k) .Eq. 0.0d0) z(k) = 1.0d0
         t = -z(k)
         Call daxpy(k-1,t,a(1,k),1,z(1),1)
  160 Continue
!     make znorm = 1.0
      s = 1.0d0/dasum(n,z,1)
      Call dscal(n,s,z,1)
      ynorm = s*ynorm
!
      If (anorm .Ne. 0.0d0) rcond = ynorm/anorm
      If (anorm .Eq. 0.0d0) rcond = 0.0d0
      Return
      End Subroutine dgeco

!***************************************************************************************************
      Subroutine dgefa(a,lda,n,ipvt,info)
      Integer lda,n,ipvt(1),info
      Real(Kind=rk) :: a(lda,*)
!
!     dgefa factors a double precision matrix by gaussian elimination.
!
!     dgefa is usually called by dgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
!
!     on entry
!
!        a       double precision(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that dgesl or dgedi will divide by zero
!                     if called.  use  rcond  in dgeco for a reliable
!                     indication of singularity.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,dscal,idamax
!
!     internal variables
!
      Real(Kind=rk) :: t
!Kuester      integer idamax,j,k,kp1,l,nm1
      Integer j,k,kp1,l,nm1
!
!
!     gaussian elimination with partial pivoting
!
      info = 0
      nm1 = n - 1
      If (nm1 .Lt. 1) go to 70
      Do 60 k = 1, nm1
         kp1 = k + 1
!
!        find l = pivot index
!
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
!
!        zero pivot implies this column already triangularized
!
         If (a(l,k) .Eq. 0.0d0) go to 40
!
!           interchange if necessary
!
            If (l .Eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       Continue
!
!           compute multipliers
!
            t = -1.0d0/a(k,k)
            Call dscal(n-k,t,a(k+1,k),1)
!
!           row elimination with column indexing
!
            Do 30 j = kp1, n
               t = a(l,j)
               If (l .Eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          Continue
               Call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       Continue
         go to 50
   40    Continue
            info = k
   50    Continue
   60 Continue
   70 Continue
      ipvt(n) = n
      If (a(n,n) .Eq. 0.0d0) info = n
      Return
      End Subroutine dgefa

!***************************************************************************************************

End Module
