module eispack

  implicit none

  Integer, Parameter, Private :: rk=8

! brought to f90 module
! declarations of pythag removed
! declarations of epslon removed
! kind=rk for all reals and all constants
! generic intrinsics max,min,sqrt,sign,abs
! continuation lines collapsed
! declared k1,...,k8 in rsm

contains

!*********************************************************************************
!*********************************************************************************
      subroutine rg(nm,n,a,wr,wi,matz,z,iv1,fv1,ierr)
!
      integer n,nm,is1,is2,ierr,matz
      real(kind=rk)  :: a(nm,n),wr(n),wi(n),z(nm,n),fv1(n)
      integer iv1(n)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a real general matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        a  contains the real general matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        wr  and  wi  contain the real and imaginary parts,
!        respectively, of the eigenvalues.  complex conjugate
!        pairs of eigenvalues appear consecutively with the
!        eigenvalue having the positive imaginary part first.
!
!        z  contains the real and imaginary parts of the eigenvectors
!        if matz is not zero.  if the j-th eigenvalue is real, the
!        j-th column of  z  contains its eigenvector.  if the j-th
!        eigenvalue is complex with positive imaginary part, the
!        j-th and (j+1)-th columns of  z  contain the real and
!        imaginary parts of its eigenvector.  the conjugate of this
!        vector is the eigenvector for the conjugate eigenvalue.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for hqr
!           and hqr2.  the normal completion code is zero.
!
!        iv1  and  fv1  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
!
!      .... To turn on balancing .....
   10 call  balanc(nm,n,a,is1,is2,fv1)
!      write(    *      ,fmt_str)'Matrix balancing done !'

!      .... To turn of balancing .....
!  10 is1 = 1
!     is2 = n
!
!      write(    *      ,fmt_str_AdNo)'Reduceing matrix to upper Hessenberg form .... '
      call  elmhes(nm,n,is1,is2,a,iv1)
!      write(    *      ,fmt_str)' done !'
!
      if (matz .ne. 0) go to 20
!     .......... find eigenvalues only ..........
!      write(    *      ,fmt_str_AdNo)'Calculation of eigenvalues with qr-method .... '
      call  hqr(nm,n,is1,is2,a,wr,wi,ierr)
!      write(    *      ,fmt_str)' done !'
!
      go to 50
!     .......... find both eigenvalues and eigenvectors ..........
   20 call  eltran(nm,n,is1,is2,a,iv1,z)
!      write(    *       ,fmt_str_AdNo)'Calculation of eigenvalues and eigenvectors with qr-method .... '
      call  hqr2(nm,n,is1,is2,a,wr,wi,z,ierr)
      if (ierr .ne. 0) go to 50
!     .... if balancing is turned of comment out next line ....
      call  balbak(nm,n,is1,is2,fv1,n,z)
!      write(    *      ,fmt_str)' done !'
!
   50 return
      end subroutine
!*********************************************************************************
!*********************************************************************************




!*********************************************************************************
      subroutine cdiv(ar,ai,br,bi,cr,ci)
      real(kind=rk)  :: ar,ai,br,bi,cr,ci
!
!     complex division, (cr,ci) = (ar,ai)/(br,bi)
!
      real(kind=rk)  :: s,ars,ais,brs,bis
      s = abs(br) + abs(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2 + bis**2
      cr = (ars*brs + ais*bis)/s
      ci = (ais*brs - ars*bis)/s
      return
      end subroutine
!*********************************************************************************
      subroutine csroot(xr,xi,yr,yi)
      real(kind=rk)  :: xr,xi,yr,yi
!
!     (yr,yi) = complex sqrt(xr,xi)
!     branch chosen so that yr .ge. 0.0 and sign(yi) .eq. sign(xi)
!
      real(kind=rk)  :: s,tr,ti
      tr = xr
      ti = xi
      s = sqrt(0.5e0_rk*(pythag(tr,ti) + abs(tr)))
      if (tr .ge. 0.0e0_rk) yr = s
      if (ti .lt. 0.0e0_rk) s = -s
      if (tr .le. 0.0e0_rk) yi = s
      if (tr .lt. 0.0e0_rk) yr = 0.5e0_rk*(ti/yi)
      if (tr .gt. 0.0e0_rk) yi = 0.5e0_rk*(ti/yr)
      return
      end subroutine
!*********************************************************************************
      real(kind=rk)   function epslon (x)
      real(kind=rk)  :: x
!
!     estimate unit roundoff in quantities of size x.
!
      real(kind=rk)  :: a,b,c,eps
!
!     this program should function properly on all systems
!     satisfying the following two assumptions,
!        1.  the base used in representing floating point
!            numbers is not a power of three.
!        2.  the quantity  a  in statement 10 is represented to
!            the accuracy used in floating point variables
!            that are stored in memory.
!     the statement number 10 and the go to 10 are intended to
!     force optimizing compilers to generate code satisfying
!     assumption 2.
!     under these assumptions, it should be true that,
!            a  is not exactly equal to four-thirds,
!            b  has a zero for its last bit or digit,
!            c  is not exactly equal to one,
!            eps  measures the separation of 1.0 from
!                 the next larger floating point number.
!     the developers of eispack would appreciate being informed
!     about any systems where these assumptions do not hold.
!
!     this version dated 4/6/83.
!
      a = 4.0e0_rk/3.0e0_rk
   10 b = a - 1.0e0_rk
      c = b + b + b
      eps = abs(c-1.0e0_rk)
      if (eps .eq. 0.0e0_rk) go to 10
      epslon = eps*abs(x)
      return
      end function epslon
!*********************************************************************************
      real(kind=rk)  function pythag(a,b)
      real(kind=rk)  :: a,b
!
!     finds sqrt(a**2+b**2) without overflow or destructive underflow
!
      real(kind=rk)  :: p,r,s,t,u
      p = max(abs(a),abs(b))
      if (p .eq. 0.0e0_rk) go to 20
      r = (min(abs(a),abs(b))/p)**2
   10 continue
         t = 4.0e0_rk + r
         if (t .eq. 4.0e0_rk) go to 20
         s = r/t
         u = 1.0e0_rk + 2.0e0_rk*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end function pythag
!*********************************************************************************
      subroutine bakvec(nm,n,t,e,m,z,ierr)
!
      integer i,j,m,n,nm,ierr
      real(kind=rk)  :: t(nm,3),e(n),z(nm,m)
!
!     this subroutine forms the eigenvectors of a nonsymmetric
!     tridiagonal matrix by back transforming those of the
!     corresponding symmetric matrix determined by  figi.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        t contains the nonsymmetric matrix.  its subdiagonal is
!          stored in the last n-1 positions of the first column,
!          its diagonal in the n positions of the second column,
!          and its superdiagonal in the first n-1 positions of
!          the third column.  t(1,1) and t(n,3) are arbitrary.
!
!        e contains the subdiagonal elements of the symmetric
!          matrix in its last n-1 positions.  e(1) is arbitrary.
!
!        m is the number of eigenvectors to be back transformed.
!
!        z contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        t is unaltered.
!
!        e is destroyed.
!
!        z contains the transformed eigenvectors
!          in its first m columns.
!
!        ierr is set to
!          zero       for normal return,
!          2*n+i      if e(i) is zero with t(i,1) or t(i-1,3) non-zero.
!                     in this case, the symmetric matrix is not similar
!                     to the original matrix, and the eigenvectors
!                     cannot be found by this program.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      if (m .eq. 0) go to 1001
      e(1) = 1.0e0_rk
      if (n .eq. 1) go to 1001
!
      do 100 i = 2, n
         if (e(i) .ne. 0.0e0_rk) go to 80
         if (t(i,1) .ne. 0.0e0_rk .or. t(i-1,3) .ne. 0.0e0_rk) go to 1000
         e(i) = 1.0e0_rk
         go to 100
   80    e(i) = e(i-1) * e(i) / t(i-1,3)
  100 continue
!
      do 120 j = 1, m
!
         do 120 i = 2, n
         z(i,j) = z(i,j) * e(i)
  120 continue
!
      go to 1001
!     .......... set error -- eigenvectors cannot be
!                found by this program ..........
 1000 ierr = 2 * n + i
 1001 return
      end subroutine
!*********************************************************************************
      subroutine balanc(nm,n,a,low,igh,scale)
!
      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      real(kind=rk)  :: a(nm,n),scale(n)
      real(kind=rk)  :: c,f,g,r,s,b2,radix
      logical noconv
!
!     this subroutine is a translation of the algol procedure balance,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!
!     this subroutine balances a real matrix and isolates
!     eigenvalues whenever possible.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the input matrix to be balanced.
!
!     on output
!
!        a contains the balanced matrix.
!
!        low and igh are two integers such that a(i,j)
!          is equal to zero if
!           (1) i is greater than j and
!           (2) j=1,...,low-1 or i=igh+1,...,n.
!
!        scale contains information determining the
!           permutations and scaling factors used.
!
!     suppose that the principal submatrix in rows low through igh
!     has been balanced, that p(j) denotes the index interchanged
!     with j during the permutation step, and that the elements
!     of the diagonal matrix used are denoted by d(i,j).  then
!        scale(j) = p(j),    for j = 1,...,low-1
!                 = d(j,j),      j = low,...,igh
!                 = p(j)         j = igh+1,...,n.
!     the order in which the interchanges are made is n to igh+1,
!     then 1 to low-1.
!
!     note that 1 is returned for igh if igh is zero formally.
!
!     the algol procedure exc contained in balance appears in
!     balanc  in line.  (note that the algol roles of identifiers
!     k,l have been reversed.)
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      radix = 16.0e0_rk
!
      b2 = radix * radix
      k = 1
      l = n
      go to 100
!     .......... in-line procedure for row and
!                column exchange ..........
   20 scale(m) = j
      if (j .eq. m) go to 50
!
      do 30 i = 1, l
         f = a(i,j)
         a(i,j) = a(i,m)
         a(i,m) = f
   30 continue
!
      do 40 i = k, n
         f = a(j,i)
         a(j,i) = a(m,i)
         a(m,i) = f
   40 continue
!
   50 go to (80,130), iexc
!     .......... search for rows isolating an eigenvalue
!                and push them down ..........
   80 if (l .eq. 1) go to 280
      l = l - 1
!     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l
         j = l + 1 - jj
!
         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (a(j,i) .ne. 0.0e0_rk) go to 120
  110    continue
!
         m = l
         iexc = 1
         go to 20
  120 continue
!
      go to 140
!     .......... search for columns isolating an eigenvalue
!                and push them left ..........
  130 k = k + 1
!
  140 do 170 j = k, l
!
         do 150 i = k, l
            if (i .eq. j) go to 150
            if (a(i,j) .ne. 0.0e0_rk) go to 170
  150    continue
!
         m = k
         iexc = 2
         go to 20
  170 continue
!     .......... now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 scale(i) = 1.0e0_rk
!     .......... iterative loop for norm reduction ..........
  190 noconv = .false.
!
      do 270 i = k, l
         c = 0.0e0_rk
         r = 0.0e0_rk
!
         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + abs(a(j,i))
            r = r + abs(a(i,j))
  200    continue
!     .......... guard against zero c or r due to underflow ..........
         if (c .eq. 0.0e0_rk .or. r .eq. 0.0e0_rk) go to 270
         g = r / radix
         f = 1.0e0_rk
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
!     .......... now balance ..........
  240    if ((c + r) / f .ge. 0.95e0_rk * s) go to 270
         g = 1.0e0_rk / f
         scale(i) = scale(i) * f
         noconv = .true.
!
         do 250 j = k, n
  250    a(i,j) = a(i,j) * g
!
         do 260 j = 1, l
  260    a(j,i) = a(j,i) * f
!
  270 continue
!
      if (noconv) go to 190
!
  280 low = k
      igh = l
      return
      end subroutine
!*********************************************************************************
      subroutine balbak(nm,n,low,igh,scale,m,z)
!
      integer i,j,k,m,n,ii,nm,igh,low
      real(kind=rk)  :: scale(n),z(nm,m)
      real(kind=rk)  :: s
!
!     this subroutine is a translation of the algol procedure balbak,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!
!     this subroutine forms the eigenvectors of a real general
!     matrix by back transforming those of the corresponding
!     balanced matrix determined by  balanc.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by  balanc.
!
!        scale contains information determining the permutations
!          and scaling factors used by  balanc.
!
!        m is the number of columns of z to be back transformed.
!
!        z contains the real and imaginary parts of the eigen-
!          vectors to be back transformed in its first m columns.
!
!     on output
!
!        z contains the real and imaginary parts of the
!          transformed eigenvectors in its first m columns.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
!
      do 110 i = low, igh
         s = scale(i)
!     .......... left hand eigenvectors are back transformed
!                if the foregoing statement is replaced by
!                s=1.0e0_rk/scale(i). ..........
         do 100 j = 1, m
  100    z(i,j) = z(i,j) * s
!
  110 continue
!     ......... for i=low-1 step -1 until 1,
!               igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
!
         do 130 j = 1, m
            s = z(i,j)
            z(i,j) = z(k,j)
            z(k,j) = s
  130    continue
!
  140 continue
!
  200 return
      end subroutine
!*********************************************************************************
      subroutine bandr(nm,n,mb,a,d,e,e2,matz,z)
!
      integer j,k,l,n,r,i1,i2,j1,j2,kr,mb,mr,m1,nm,n2,r1,ugl,maxl,maxr
      real(kind=rk)  :: a(nm,mb),d(n),e(n),e2(n),z(nm,n)
      real(kind=rk)  :: g,u,b1,b2,c2,f1,f2,s2,dmin,dminrt
      logical matz
!
!     this subroutine is a translation of the algol procedure bandrd,
!     num. math. 12, 231-241(1968) by schwarz.
!     handbook for auto. comp., vol.ii-linear algebra, 273-283(1971).
!
!     this subroutine reduces a real symmetric band matrix
!     to a symmetric tridiagonal matrix using and optionally
!     accumulating orthogonal similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        mb is the (half) band width of the matrix, defined as the
!          number of adjacent diagonals, including the principal
!          diagonal, required to specify the non-zero portion of the
!          lower triangle of the matrix.
!
!        a contains the lower triangle of the symmetric band input
!          matrix stored as an n by mb array.  its lowest subdiagonal
!          is stored in the last n+1-mb positions of the first column,
!          its next subdiagonal in the last n+2-mb positions of the
!          second column, further subdiagonals similarly, and finally
!          its principal diagonal in the n positions of the last column.
!          contents of storages not part of the matrix are arbitrary.
!
!        matz should be set to .true. if the transformation matrix is
!          to be accumulated, and to .false. otherwise.
!
!     on output
!
!        a has been destroyed, except for its last two columns which
!          contain a copy of the tridiagonal matrix.
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
!        z contains the orthogonal transformation matrix produced in
!          the reduction if matz has been set to .true.  otherwise, z
!          is not referenced.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      dmin = 2.0e0_rk**(-64)
      dminrt = 2.0e0_rk**(-32)
!     .......... initialize diagonal scaling matrix ..........
      do 30 j = 1, n
   30 d(j) = 1.0e0_rk
!
      if (.not. matz) go to 60
!
      do 50 j = 1, n
!
         do 40 k = 1, n
   40    z(j,k) = 0.0e0_rk
!
         z(j,j) = 1.0e0_rk
   50 continue
!
   60 m1 = mb - 1
      if (m1 - 1) 900, 800, 70
   70 n2 = n - 2
!
      do 700 k = 1, n2
         maxr = min0(m1,n-k)
!     .......... for r=maxr step -1 until 2 do -- ..........
         do 600 r1 = 2, maxr
            r = maxr + 2 - r1
            kr = k + r
            mr = mb - r
            g = a(kr,mr)
            a(kr-1,1) = a(kr-1,mr+1)
            ugl = k
!
            do 500 j = kr, n, m1
               j1 = j - 1
               j2 = j1 - 1
               if (g .eq. 0.0e0_rk) go to 600
               b1 = a(j1,1) / g
               b2 = b1 * d(j1) / d(j)
               s2 = 1.0e0_rk / (1.0e0_rk + b1 * b2)
               if (s2 .ge. 0.5e0_rk ) go to 450
               b1 = g / a(j1,1)
               b2 = b1 * d(j) / d(j1)
               c2 = 1.0e0_rk - s2
               d(j1) = c2 * d(j1)
               d(j) = c2 * d(j)
               f1 = 2.0e0_rk * a(j,m1)
               f2 = b1 * a(j1,mb)
               a(j,m1) = -b2 * (b1 * a(j,m1) - a(j,mb)) - f2 + a(j,m1)
               a(j1,mb) = b2 * (b2 * a(j,mb) + f1) + a(j1,mb)
               a(j,mb) = b1 * (f2 - f1) + a(j,mb)
!
               do 200 l = ugl, j2
                  i2 = mb - j + l
                  u = a(j1,i2+1) + b2 * a(j,i2)
                  a(j,i2) = -b1 * a(j1,i2+1) + a(j,i2)
                  a(j1,i2+1) = u
  200          continue
!
               ugl = j
               a(j1,1) = a(j1,1) + b2 * g
               if (j .eq. n) go to 350
               maxl = min0(m1,n-j1)
!
               do 300 l = 2, maxl
                  i1 = j1 + l
                  i2 = mb - l
                  u = a(i1,i2) + b2 * a(i1,i2+1)
                  a(i1,i2+1) = -b1 * a(i1,i2) + a(i1,i2+1)
                  a(i1,i2) = u
  300          continue
!
               i1 = j + m1
               if (i1 .gt. n) go to 350
               g = b2 * a(i1,1)
  350          if (.not. matz) go to 500
!
               do 400 l = 1, n
                  u = z(l,j1) + b2 * z(l,j)
                  z(l,j) = -b1 * z(l,j1) + z(l,j)
                  z(l,j1) = u
  400          continue
!
               go to 500
!
  450          u = d(j1)
               d(j1) = s2 * d(j)
               d(j) = s2 * u
               f1 = 2.0e0_rk * a(j,m1)
               f2 = b1 * a(j,mb)
               u = b1 * (f2 - f1) + a(j1,mb)
               a(j,m1) = b2 * (b1 * a(j,m1) - a(j1,mb)) + f2 - a(j,m1)
               a(j1,mb) = b2 * (b2 * a(j1,mb) + f1) + a(j,mb)
               a(j,mb) = u
!
               do 460 l = ugl, j2
                  i2 = mb - j + l
                  u = b2 * a(j1,i2+1) + a(j,i2)
                  a(j,i2) = -a(j1,i2+1) + b1 * a(j,i2)
                  a(j1,i2+1) = u
  460          continue
!
               ugl = j
               a(j1,1) = b2 * a(j1,1) + g
               if (j .eq. n) go to 480
               maxl = min0(m1,n-j1)
!
               do 470 l = 2, maxl
                  i1 = j1 + l
                  i2 = mb - l
                  u = b2 * a(i1,i2) + a(i1,i2+1)
                  a(i1,i2+1) = -a(i1,i2) + b1 * a(i1,i2+1)
                  a(i1,i2) = u
  470          continue
!
               i1 = j + m1
               if (i1 .gt. n) go to 480
               g = a(i1,1)
               a(i1,1) = b1 * a(i1,1)
  480          if (.not. matz) go to 500
!
               do 490 l = 1, n
                  u = b2 * z(l,j1) + z(l,j)
                  z(l,j) = -z(l,j1) + b1 * z(l,j)
                  z(l,j1) = u
  490          continue
!
  500       continue
!
  600    continue
!
         if (mod(k,64) .ne. 0) go to 700
!     .......... rescale to avoid underflow or overflow ..........
         do 650 j = k, n
            if (d(j) .ge. dmin) go to 650
            maxl = max0(1,mb+1-j)
!
            do 610 l = maxl, m1
  610       a(j,l) = dminrt * a(j,l)
!
            if (j .eq. n) go to 630
            maxl = min0(m1,n-j)
!
            do 620 l = 1, maxl
               i1 = j + l
               i2 = mb - l
               a(i1,i2) = dminrt * a(i1,i2)
  620       continue
!
  630       if (.not. matz) go to 645
!
            do 640 l = 1, n
  640       z(l,j) = dminrt * z(l,j)
!
  645       a(j,mb) = dmin * a(j,mb)
            d(j) = d(j) / dmin
  650    continue
!
  700 continue
!     .......... form square root of scaling matrix ..........
  800 do 810 j = 2, n
  810 e(j) = sqrt(d(j))
!
      if (.not. matz) go to 840
!
      do 830 j = 1, n
!
         do 820 k = 2, n
  820    z(j,k) = e(k) * z(j,k)
!
  830 continue
!
  840 u = 1.0e0_rk
!
      do 850 j = 2, n
         a(j,m1) = u * e(j) * a(j,m1)
         u = e(j)
         e2(j) = a(j,m1) ** 2
         a(j,mb) = d(j) * a(j,mb)
         d(j) = a(j,mb)
         e(j) = a(j,m1)
  850 continue
!
      d(1) = a(1,mb)
      e(1) = 0.0e0_rk
      e2(1) = 0.0e0_rk
      go to 1001
!
  900 do 950 j = 1, n
         d(j) = a(j,mb)
         e(j) = 0.0e0_rk
         e2(j) = 0.0e0_rk
  950 continue
!
 1001 return
      end subroutine
!*********************************************************************************
      subroutine bandv(nm,n,mbw,a,e21,m,w,z,ierr,nv,rv,rv6)
!
      integer i,j,k,m,n,r,ii,ij,jj,kj,mb,m1,nm,nv,ij1,its,kj1,mbw,m21,ierr,maxj,maxk,group
      real(kind=rk)  :: a(nm,mbw),w(m),z(nm,m),rv(nv),rv6(n)
      real(kind=rk)  :: u,v,uk,xu,x0,x1,e21,eps2,eps3,eps4,norm,order
!
!     this subroutine finds those eigenvectors of a real symmetric
!     band matrix corresponding to specified eigenvalues, using inverse
!     iteration.  the subroutine may also be used to solve systems
!     of linear equations with a symmetric or non-symmetric band
!     coefficient matrix.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        mbw is the number of columns of the array a used to store the
!          band matrix.  if the matrix is symmetric, mbw is its (half)
!          band width, denoted mb and defined as the number of adjacent
!          diagonals, including the principal diagonal, required to
!          specify the non-zero portion of the lower triangle of the
!          matrix.  if the subroutine is being used to solve systems
!          of linear equations and the coefficient matrix is not
!          symmetric, it must however have the same number of adjacent
!          diagonals above the main diagonal as below, and in this
!          case, mbw=2*mb-1.
!
!        a contains the lower triangle of the symmetric band input
!          matrix stored as an n by mb array.  its lowest subdiagonal
!          is stored in the last n+1-mb positions of the first column,
!          its next subdiagonal in the last n+2-mb positions of the
!          second column, further subdiagonals similarly, and finally
!          its principal diagonal in the n positions of column mb.
!          if the subroutine is being used to solve systems of linear
!          equations and the coefficient matrix is not symmetric, a is
!          n by 2*mb-1 instead with lower triangle as above and with
!          its first superdiagonal stored in the first n-1 positions of
!          column mb+1, its second superdiagonal in the first n-2
!          positions of column mb+2, further superdiagonals similarly,
!          and finally its highest superdiagonal in the first n+1-mb
!          positions of the last column.
!          contents of storages not part of the matrix are arbitrary.
!
!        e21 specifies the ordering of the eigenvalues and contains
!            0.0e0_rk if the eigenvalues are in ascending order, or
!            2.0e0_rk if the eigenvalues are in descending order.
!          if the subroutine is being used to solve systems of linear
!          equations, e21 should be set to 1.0e0_rk if the coefficient
!          matrix is symmetric and to -1.0e0_rk if not.
!
!        m is the number of specified eigenvalues or the number of
!          systems of linear equations.
!
!        w contains the m eigenvalues in ascending or descending order.
!          if the subroutine is being used to solve systems of linear
!          equations (a-w(r)*i)*x(r)=b(r), where i is the identity
!          matrix, w(r) should be set accordingly, for r=1,2,...,m.
!
!        z contains the constant matrix columns (b(r),r=1,2,...,m), if
!          the subroutine is used to solve systems of linear equations.
!
!        nv must be set to the dimension of the array parameter rv
!          as declared in the calling program dimension statement.
!
!     on output
!
!        a and w are unaltered.
!
!        z contains the associated set of orthogonal eigenvectors.
!          any vector which fails to converge is set to zero.  if the
!          subroutine is used to solve systems of linear equations,
!          z contains the solution matrix columns (x(r),r=1,2,...,m).
!
!        ierr is set to
!          zero       for normal return,
!          -r         if the eigenvector corresponding to the r-th
!                     eigenvalue fails to converge, or if the r-th
!                     system of linear equations is nearly singular.
!
!        rv and rv6 are temporary storage arrays.  note that rv is
!          of dimension at least n*(2*mb-1).  if the subroutine
!          is being used to solve systems of linear equations, the
!          determinant (up to sign) of a-w(m)*i is available, upon
!          return, as the product of the first n elements of rv.
!
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      if (m .eq. 0) go to 1001
      mb = mbw
      if (e21 .lt. 0.0e0_rk) mb = (mbw + 1) / 2
      m1 = mb - 1
      m21 = m1 + mb
      order = 1.0e0_rk - abs(e21)
!     .......... find vectors by inverse iteration ..........
      do 920 r = 1, m
         its = 1
         x1 = w(r)
         if (r .ne. 1) go to 100
!     .......... compute norm of matrix ..........
         norm = 0.0e0_rk
!
         do 60 j = 1, mb
            jj = mb + 1 - j
            kj = jj + m1
            ij = 1
            v = 0.0e0_rk
!
            do 40 i = jj, n
               v = v + abs(a(i,j))
               if (e21 .ge. 0.0e0_rk) go to 40
               v = v + abs(a(ij,kj))
               ij = ij + 1
   40       continue
!
            norm = max(norm,v)
   60    continue
!
         if (e21 .lt. 0.0e0_rk) norm = 0.5e0_rk * norm
!     .......... eps2 is the criterion for grouping,
!                eps3 replaces zero pivots and equal
!                roots are modified by eps3,
!                eps4 is taken very small to avoid overflow ..........
         if (norm .eq. 0.0e0_rk) norm = 1.0e0_rk
         eps2 = 1.0e-3_rk * norm * abs(order)
         eps3 = epslon(norm)
         uk = n
         uk = sqrt(uk)
         eps4 = uk * eps3
   80    group = 0
         go to 120
!     .......... look for close or coincident roots ..........
  100    if (abs(x1-x0) .ge. eps2) go to 80
         group = group + 1
         if (order * (x1 - x0) .le. 0.0e0_rk) x1 = x0 + order * eps3
!     .......... expand matrix, subtract eigenvalue,
!                and initialize vector ..........
  120    do 200 i = 1, n
            ij = i + min0(0,i-m1) * n
            kj = ij + mb * n
            ij1 = kj + m1 * n
            if (m1 .eq. 0) go to 180
!
            do 150 j = 1, m1
               if (ij .gt. m1) go to 125
               if (ij .gt. 0) go to 130
               rv(ij1) = 0.0e0_rk
               ij1 = ij1 + n
               go to 130
  125          rv(ij) = a(i,j)
  130          ij = ij + n
               ii = i + j
               if (ii .gt. n) go to 150
               jj = mb - j
               if (e21 .ge. 0.0e0_rk) go to 140
               ii = i
               jj = mb + j
  140          rv(kj) = a(ii,jj)
               kj = kj + n
  150       continue
!
  180       rv(ij) = a(i,mb) - x1
            rv6(i) = eps4
            if (order .eq. 0.0e0_rk) rv6(i) = z(i,r)
  200    continue
!
         if (m1 .eq. 0) go to 600
!     .......... elimination with interchanges ..........
         do 580 i = 1, n
            ii = i + 1
            maxk = min0(i+m1-1,n)
            maxj = min0(n-i,m21-2) * n
!
            do 360 k = i, maxk
               kj1 = k
               j = kj1 + n
               jj = j + maxj
!
               do 340 kj = j, jj, n
                  rv(kj1) = rv(kj)
                  kj1 = kj
  340          continue
!
               rv(kj1) = 0.0e0_rk
  360       continue
!
            if (i .eq. n) go to 580
            u = 0.0e0_rk
            maxk = min0(i+m1,n)
            maxj = min0(n-ii,m21-2) * n
!
            do 450 j = i, maxk
               if (abs(rv(j)) .lt. abs(u)) go to 450
               u = rv(j)
               k = j
  450       continue
!
            j = i + n
            jj = j + maxj
            if (k .eq. i) go to 520
            kj = k
!
            do 500 ij = i, jj, n
               v = rv(ij)
               rv(ij) = rv(kj)
               rv(kj) = v
               kj = kj + n
  500       continue
!
            if (order .ne. 0.0e0_rk) go to 520
            v = rv6(i)
            rv6(i) = rv6(k)
            rv6(k) = v
  520       if (u .eq. 0.0e0_rk) go to 580
!
            do 560 k = ii, maxk
               v = rv(k) / u
               kj = k
!
               do 540 ij = j, jj, n
                  kj = kj + n
                  rv(kj) = rv(kj) - v * rv(ij)
  540          continue
!
               if (order .eq. 0.0e0_rk) rv6(k) = rv6(k) - v * rv6(i)
  560       continue
!
  580    continue
!     .......... back substitution
!                for i=n step -1 until 1 do -- ..........
  600    do 630 ii = 1, n
            i = n + 1 - ii
            maxj = min0(ii,m21)
            if (maxj .eq. 1) go to 620
            ij1 = i
            j = ij1 + n
            jj = j + (maxj - 2) * n
!
            do 610 ij = j, jj, n
               ij1 = ij1 + 1
               rv6(i) = rv6(i) - rv(ij) * rv6(ij1)
  610       continue
!
  620       v = rv(i)
            if (abs(v) .ge. eps3) go to 625
!     .......... set error -- nearly singular linear system ..........
            if (order .eq. 0.0e0_rk) ierr = -r
            v = sign(eps3,v)
  625       rv6(i) = rv6(i) / v
  630    continue
!
         xu = 1.0e0_rk
         if (order .eq. 0.0e0_rk) go to 870
!     .......... orthogonalize with respect to previous
!                members of group ..........
         if (group .eq. 0) go to 700
!
         do 680 jj = 1, group
            j = r - group - 1 + jj
            xu = 0.0e0_rk
!
            do 640 i = 1, n
  640       xu = xu + rv6(i) * z(i,j)
!
            do 660 i = 1, n
  660       rv6(i) = rv6(i) - xu * z(i,j)
!
  680    continue
!
  700    norm = 0.0e0_rk
!
         do 720 i = 1, n
  720    norm = norm + abs(rv6(i))
!
         if (norm .ge. 0.1e0_rk) go to 840
!     .......... in-line procedure for choosing
!                a new starting vector ..........
         if (its .ge. n) go to 830
         its = its + 1
         xu = eps4 / (uk + 1.0e0_rk)
         rv6(1) = eps4
!
         do 760 i = 2, n
  760    rv6(i) = xu
!
         rv6(its) = rv6(its) - eps4 * uk
         go to 600
!     .......... set error -- non-converged eigenvector ..........
  830    ierr = -r
         xu = 0.0e0_rk
         go to 870
!     .......... normalize so that sum of squares is
!                1 and expand to full order ..........
  840    u = 0.0e0_rk
!
         do 860 i = 1, n
  860    u = pythag(u,rv6(i))
!
         xu = 1.0e0_rk / u
!
  870    do 900 i = 1, n
  900    z(i,r) = rv6(i) * xu
!
         x0 = x1
  920 continue
!
 1001 return
      end subroutine
!*********************************************************************************
      subroutine bisect(n,eps1,d,e,e2,lb,ub,mm,m,w,ind,ierr,rv4,rv5)
!
      integer i,j,k,l,m,n,p,q,r,s,ii,mm,m1,m2,tag,ierr,isturm
      real(kind=rk)  :: d(n),e(n),e2(n),w(mm),rv4(n),rv5(n)
      real(kind=rk)  :: u,v,lb,t1,t2,ub,xu,x0,x1,eps1,tst1,tst2
      integer ind(mm)
!
!     this subroutine is a translation of the bisection technique
!     in the algol procedure tristurm by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
!
!     this subroutine finds those eigenvalues of a tridiagonal
!     symmetric matrix which lie in a specified interval,
!     using bisection.
!
!     on input
!
!        n is the order of the matrix.
!
!        eps1 is an absolute error tolerance for the computed
!          eigenvalues.  if the input eps1 is non-positive,
!          it is reset for each submatrix to a default value,
!          namely, minus the product of the relative machine
!          precision and the 1-norm of the submatrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2(1) is arbitrary.
!
!        lb and ub define the interval to be searched for eigenvalues.
!          if lb is not less than ub, no eigenvalues will be found.
!
!        mm should be set to an upper bound for the number of
!          eigenvalues in the interval.  warning. if more than
!          mm eigenvalues are determined to lie in the interval,
!          an error return is made with no eigenvalues found.
!
!     on output
!
!        eps1 is unaltered unless it has been reset to its
!          (last) default value.
!
!        d and e are unaltered.
!
!        elements of e2, corresponding to elements of e regarded
!          as negligible, have been replaced by zero causing the
!          matrix to split into a direct sum of submatrices.
!          e2(1) is also set to zero.
!
!        m is the number of eigenvalues determined to lie in (lb,ub).
!
!        w contains the m eigenvalues in ascending order.
!
!        ind contains in its first m positions the submatrix indices
!          associated with the corresponding eigenvalues in w --
!          1 for eigenvalues belonging to the first submatrix from
!          the top, 2 for those belonging to the second submatrix, etc..
!
!        ierr is set to
!          zero       for normal return,
!          3*n+1      if m exceeds mm.
!
!        rv4 and rv5 are temporary storage arrays.
!
!     the algol procedure sturmcnt contained in tristurm
!     appears in bisect in-line.
!
!     note that subroutine tql1 or imtql1 is generally faster than
!     bisect, if more than n/4 eigenvalues are to be found.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      tag = 0
      t1 = lb
      t2 = ub
!     .......... look for small sub-diagonal entries ..........
      do 40 i = 1, n
         if (i .eq. 1) go to 20
         tst1 = abs(d(i)) + abs(d(i-1))
         tst2 = tst1 + abs(e(i))
         if (tst2 .gt. tst1) go to 40
   20    e2(i) = 0.0e0_rk
   40 continue
!     .......... determine the number of eigenvalues
!                in the interval ..........
      p = 1
      q = n
      x1 = ub
      isturm = 1
      go to 320
   60 m = s
      x1 = lb
      isturm = 2
      go to 320
   80 m = m - s
      if (m .gt. mm) go to 980
      q = 0
      r = 0
!     .......... establish and process next submatrix, refining
!                interval by the gerschgorin bounds ..........
  100 if (r .eq. m) go to 1001
      tag = tag + 1
      p = q + 1
      xu = d(p)
      x0 = d(p)
      u = 0.0e0_rk
!
      do 120 q = p, n
         x1 = u
         u = 0.0e0_rk
         v = 0.0e0_rk
         if (q .eq. n) go to 110
         u = abs(e(q+1))
         v = e2(q+1)
  110    xu = min(d(q)-(x1+u),xu)
         x0 = max(d(q)+(x1+u),x0)
         if (v .eq. 0.0e0_rk) go to 140
  120 continue
!
  140 x1 = epslon(max(abs(xu),abs(x0)))
      if (eps1 .le. 0.0e0_rk) eps1 = -x1
      if (p .ne. q) go to 180
!     .......... check for isolated root within interval ..........
      if (t1 .gt. d(p) .or. d(p) .ge. t2) go to 940
      m1 = p
      m2 = p
      rv5(p) = d(p)
      go to 900
  180 x1 = x1 * (q - p + 1)
      lb = max(t1,xu-x1)
      ub = min(t2,x0+x1)
      x1 = lb
      isturm = 3
      go to 320
  200 m1 = s + 1
      x1 = ub
      isturm = 4
      go to 320
  220 m2 = s
      if (m1 .gt. m2) go to 940
!     .......... find roots by bisection ..........
      x0 = ub
      isturm = 5
!
      do 240 i = m1, m2
         rv5(i) = ub
         rv4(i) = lb
  240 continue
!     .......... loop for k-th eigenvalue
!                for k=m2 step -1 until m1 do --
!                (-do- not used to legalize -computed go to-) ..........
      k = m2
  250    xu = lb
!     .......... for i=k step -1 until m1 do -- ..........
         do 260 ii = m1, k
            i = m1 + k - ii
            if (xu .ge. rv4(i)) go to 260
            xu = rv4(i)
            go to 280
  260    continue
!
  280    if (x0 .gt. rv5(k)) x0 = rv5(k)
!     .......... next bisection step ..........
  300    x1 = (xu + x0) * 0.5e0_rk
         if ((x0 - xu) .le. abs(eps1)) go to 420
         tst1 = 2.0e0_rk * (abs(xu) + abs(x0))
         tst2 = tst1 + (x0 - xu)
         if (tst2 .eq. tst1) go to 420
!     .......... in-line procedure for sturm sequence ..........
  320    s = p - 1
         u = 1.0e0_rk
!
         do 340 i = p, q
            if (u .ne. 0.0e0_rk) go to 325
            v = abs(e(i)) / epslon(1.0e0_rk)
            if (e2(i) .eq. 0.0e0_rk) v = 0.0e0_rk
            go to 330
  325       v = e2(i) / u
  330       u = d(i) - x1 - v
            if (u .lt. 0.0e0_rk) s = s + 1
  340    continue
!
         go to (60,80,200,220,360), isturm
!     .......... refine intervals ..........
  360    if (s .ge. k) go to 400
         xu = x1
         if (s .ge. m1) go to 380
         rv4(m1) = x1
         go to 300
  380    rv4(s+1) = x1
         if (rv5(s) .gt. x1) rv5(s) = x1
         go to 300
  400    x0 = x1
         go to 300
!     .......... k-th eigenvalue found ..........
  420    rv5(k) = x1
      k = k - 1
      if (k .ge. m1) go to 250
!     .......... order eigenvalues tagged with their
!                submatrix associations ..........
  900 s = r
      r = r + m2 - m1 + 1
      j = 1
      k = m1
!
      do 920 l = 1, r
         if (j .gt. s) go to 910
         if (k .gt. m2) go to 940
         if (rv5(k) .ge. w(l)) go to 915
!
         do 905 ii = j, s
            i = l + s - ii
            w(i+1) = w(i)
            ind(i+1) = ind(i)
  905    continue
!
  910    w(l) = rv5(k)
         ind(l) = tag
         k = k + 1
         go to 920
  915    j = j + 1
  920 continue
!
  940 if (q .lt. n) go to 100
      go to 1001
!     .......... set error -- underestimate of number of
!                eigenvalues in interval ..........
  980 ierr = 3 * n + 1
 1001 lb = t1
      ub = t2
      return
      end subroutine
!*********************************************************************************
      subroutine bqr(nm,n,mb,a,t,r,ierr,nv,rv)
!
      integer i,j,k,l,m,n,ii,ik,jk,jm,kj,kk,km,ll,mb,mk,mn,mz,m1,m2,m3,m4,ni,nm,nv,its,kj1,m21,m31,ierr,imult
      real(kind=rk)  :: a(nm,mb),rv(nv)
      real(kind=rk)  :: f,g,q,r,s,t,tst1,tst2,scale
!
!     this subroutine is a translation of the algol procedure bqr,
!     num. math. 16, 85-92(1970) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol ii-linear algebra, 266-272(1971).
!
!     this subroutine finds the eigenvalue of smallest (usually)
!     magnitude of a real symmetric band matrix using the
!     qr algorithm with shifts of origin.  consecutive calls
!     can be made to find further eigenvalues.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        mb is the (half) band width of the matrix, defined as the
!          number of adjacent diagonals, including the principal
!          diagonal, required to specify the non-zero portion of the
!          lower triangle of the matrix.
!
!        a contains the lower triangle of the symmetric band input
!          matrix stored as an n by mb array.  its lowest subdiagonal
!          is stored in the last n+1-mb positions of the first column,
!          its next subdiagonal in the last n+2-mb positions of the
!          second column, further subdiagonals similarly, and finally
!          its principal diagonal in the n positions of the last column.
!          contents of storages not part of the matrix are arbitrary.
!          on a subsequent call, its output contents from the previous
!          call should be passed.
!
!        t specifies the shift (of eigenvalues) applied to the diagonal
!          of a in forming the input matrix. what is actually determined
!          is the eigenvalue of a+ti (i is the identity matrix) nearest
!          to t.  on a subsequent call, the output value of t from the
!          previous call should be passed if the next nearest eigenvalue
!          is sought.
!
!        r should be specified as zero on the first call, and as its
!          output value from the previous call on a subsequent call.
!          it is used to determine when the last row and column of
!          the transformed band matrix can be regarded as negligible.
!
!        nv must be set to the dimension of the array parameter rv
!          as declared in the calling program dimension statement.
!
!     on output
!
!        a contains the transformed band matrix.  the matrix a+ti
!          derived from the output parameters is similar to the
!          input a+ti to within rounding errors.  its last row and
!          column are null (if ierr is zero).
!
!        t contains the computed eigenvalue of a+ti (if ierr is zero).
!
!        r contains the maximum of its input value and the norm of the
!          last column of the input matrix a.
!
!        ierr is set to
!          zero       for normal return,
!          n          if the eigenvalue has not been
!                     determined after 30 iterations.
!
!        rv is a temporary storage array of dimension at least
!          (2*mb**2+4*mb-3).  the first (3*mb-2) locations correspond
!          to the algol array b, the next (2*mb-1) locations correspond
!          to the algol array h, and the final (2*mb**2-mb) locations
!          correspond to the mb by (2*mb-1) algol array u.
!
!     note. for a subsequent call, n should be replaced by n-1, but
!     mb should not be altered even when it exceeds the current n.
!
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      m1 = min0(mb,n)
      m = m1 - 1
      m2 = m + m
      m21 = m2 + 1
      m3 = m21 + m
      m31 = m3 + 1
      m4 = m31 + m2
      mn = m + n
      mz = mb - m1
      its = 0
!     .......... test for convergence ..........
   40 g = a(n,mb)
      if (m .eq. 0) go to 360
      f = 0.0e0_rk
!
      do 50 k = 1, m
         mk = k + mz
         f = f + abs(a(n,mk))
   50 continue
!
      if (its .eq. 0 .and. f .gt. r) r = f
      tst1 = r
      tst2 = tst1 + f
      if (tst2 .le. tst1) go to 360
      if (its .eq. 30) go to 1000
      its = its + 1
!     .......... form shift from bottom 2 by 2 minor ..........
      if (f .gt. 0.25e0_rk * r .and. its .lt. 5) go to 90
      f = a(n,mb-1)
      if (f .eq. 0.0e0_rk) go to 70
      q = (a(n-1,mb) - g) / (2.0e0_rk * f)
      s = pythag(q,1.0e0_rk)
      g = g - f / (q + sign(s,q))
   70 t = t + g
!
      do 80 i = 1, n
   80 a(i,mb) = a(i,mb) - g
!
   90 do 100 k = m31, m4
  100 rv(k) = 0.0e0_rk
!
      do 350 ii = 1, mn
         i = ii - m
         ni = n - ii
         if (ni .lt. 0) go to 230
!     .......... form column of shifted matrix a-g*i ..........
         l = max0(1,2-i)
!
         do 110 k = 1, m3
  110    rv(k) = 0.0e0_rk
!
         do 120 k = l, m1
            km = k + m
            mk = k + mz
            rv(km) = a(ii,mk)
  120    continue
!
         ll = min0(m,ni)
         if (ll .eq. 0) go to 135
!
         do 130 k = 1, ll
            km = k + m21
            ik = ii + k
            mk = mb - k
            rv(km) = a(ik,mk)
  130    continue
!     .......... pre-multiply with householder reflections ..........
  135    ll = m2
         imult = 0
!     .......... multiplication procedure ..........
  140    kj = m4 - m1
!
         do 170 j = 1, ll
            kj = kj + m1
            jm = j + m3
            if (rv(jm) .eq. 0.0e0_rk) go to 170
            f = 0.0e0_rk
!
            do 150 k = 1, m1
               kj = kj + 1
               jk = j + k - 1
               f = f + rv(kj) * rv(jk)
  150       continue
!
            f = f / rv(jm)
            kj = kj - m1
!
            do 160 k = 1, m1
               kj = kj + 1
               jk = j + k - 1
               rv(jk) = rv(jk) - rv(kj) * f
  160       continue
!
            kj = kj - m1
  170    continue
!
         if (imult .ne. 0) go to 280
!     .......... householder reflection ..........
         f = rv(m21)
         s = 0.0e0_rk
         rv(m4) = 0.0e0_rk
         scale = 0.0e0_rk
!
         do 180 k = m21, m3
  180    scale = scale + abs(rv(k))
!
         if (scale .eq. 0.0e0_rk) go to 210
!
         do 190 k = m21, m3
  190    s = s + (rv(k)/scale)**2
!
         s = scale * scale * s
         g = -sign(sqrt(s),f)
         rv(m21) = g
         rv(m4) = s - f * g
         kj = m4 + m2 * m1 + 1
         rv(kj) = f - g
!
         do 200 k = 2, m1
            kj = kj + 1
            km = k + m2
            rv(kj) = rv(km)
  200    continue
!     .......... save column of triangular factor r ..........
  210    do 220 k = l, m1
            km = k + m
            mk = k + mz
            a(ii,mk) = rv(km)
  220    continue
!
  230    l = max0(1,m1+1-i)
         if (i .le. 0) go to 300
!     .......... perform additional steps ..........
         do 240 k = 1, m21
  240    rv(k) = 0.0e0_rk
!
         ll = min0(m1,ni+m1)
!     .......... get row of triangular factor r ..........
         do 250 kk = 1, ll
            k = kk - 1
            km = k + m1
            ik = i + k
            mk = mb - k
            rv(km) = a(ik,mk)
  250    continue
!     .......... post-multiply with householder reflections ..........
         ll = m1
         imult = 1
         go to 140
!     .......... store column of new a matrix ..........
  280    do 290 k = l, m1
            mk = k + mz
            a(i,mk) = rv(k)
  290    continue
!     .......... update householder reflections ..........
  300    if (l .gt. 1) l = l - 1
         kj1 = m4 + l * m1
!
         do 320 j = l, m2
            jm = j + m3
            rv(jm) = rv(jm+1)
!
            do 320 k = 1, m1
               kj1 = kj1 + 1
               kj = kj1 - m1
               rv(kj) = rv(kj1)
  320    continue
!
  350 continue
!
      go to 40
!     .......... convergence ..........
  360 t = t + g
!
      do 380 i = 1, n
  380 a(i,mb) = a(i,mb) - g
!
      do 400 k = 1, m1
         mk = k + mz
         a(n,mk) = 0.0e0_rk
  400 continue
!
      go to 1001
!     .......... set error -- no convergence to
!                eigenvalue after 30 iterations ..........
 1000 ierr = n
 1001 return
      end subroutine
!*********************************************************************************
      subroutine cbabk2(nm,n,low,igh,scale,m,zr,zi)
!
      integer i,j,k,m,n,ii,nm,igh,low
      real(kind=rk)  :: scale(n),zr(nm,m),zi(nm,m)
      real(kind=rk)  :: s
!
!     this subroutine is a translation of the algol procedure
!     cbabk2, which is a complex version of balbak,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!
!     this subroutine forms the eigenvectors of a complex general
!     matrix by back transforming those of the corresponding
!     balanced matrix determined by  cbal.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by  cbal.
!
!        scale contains information determining the permutations
!          and scaling factors used by  cbal.
!
!        m is the number of eigenvectors to be back transformed.
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the eigenvectors to be
!          back transformed in their first m columns.
!
!     on output
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the transformed eigenvectors
!          in their first m columns.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
!
      do 110 i = low, igh
         s = scale(i)
!     .......... left hand eigenvectors are back transformed
!                if the foregoing statement is replaced by
!                s=1.0e0_rk/scale(i). ..........
         do 100 j = 1, m
            zr(i,j) = zr(i,j) * s
            zi(i,j) = zi(i,j) * s
  100    continue
!
  110 continue
!     .......... for i=low-1 step -1 until 1,
!                igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
!
         do 130 j = 1, m
            s = zr(i,j)
            zr(i,j) = zr(k,j)
            zr(k,j) = s
            s = zi(i,j)
            zi(i,j) = zi(k,j)
            zi(k,j) = s
  130    continue
!
  140 continue
!
  200 return
      end subroutine
!*********************************************************************************
      subroutine cbal(nm,n,ar,ai,low,igh,scale)
!
      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      real(kind=rk)  :: ar(nm,n),ai(nm,n),scale(n)
      real(kind=rk)  :: c,f,g,r,s,b2,radix
      logical noconv
!
!     this subroutine is a translation of the algol procedure
!     cbalance, which is a complex version of balance,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!
!     this subroutine balances a complex matrix and isolates
!     eigenvalues whenever possible.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the complex matrix to be balanced.
!
!     on output
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the balanced matrix.
!
!        low and igh are two integers such that ar(i,j) and ai(i,j)
!          are equal to zero if
!           (1) i is greater than j and
!           (2) j=1,...,low-1 or i=igh+1,...,n.
!
!        scale contains information determining the
!           permutations and scaling factors used.
!
!     suppose that the principal submatrix in rows low through igh
!     has been balanced, that p(j) denotes the index interchanged
!     with j during the permutation step, and that the elements
!     of the diagonal matrix used are denoted by d(i,j).  then
!        scale(j) = p(j),    for j = 1,...,low-1
!                 = d(j,j)       j = low,...,igh
!                 = p(j)         j = igh+1,...,n.
!     the order in which the interchanges are made is n to igh+1,
!     then 1 to low-1.
!
!     note that 1 is returned for igh if igh is zero formally.
!
!     the algol procedure exc contained in cbalance appears in
!     cbal  in line.  (note that the algol roles of identifiers
!     k,l have been reversed.)
!
!     arithmetic is real throughout.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      radix = 16.0e0_rk
!
      b2 = radix * radix
      k = 1
      l = n
      go to 100
!     .......... in-line procedure for row and
!                column exchange ..........
   20 scale(m) = j
      if (j .eq. m) go to 50
!
      do 30 i = 1, l
         f = ar(i,j)
         ar(i,j) = ar(i,m)
         ar(i,m) = f
         f = ai(i,j)
         ai(i,j) = ai(i,m)
         ai(i,m) = f
   30 continue
!
      do 40 i = k, n
         f = ar(j,i)
         ar(j,i) = ar(m,i)
         ar(m,i) = f
         f = ai(j,i)
         ai(j,i) = ai(m,i)
         ai(m,i) = f
   40 continue
!
   50 go to (80,130), iexc
!     .......... search for rows isolating an eigenvalue
!                and push them down ..........
   80 if (l .eq. 1) go to 280
      l = l - 1
!     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l
         j = l + 1 - jj
!
         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (ar(j,i) .ne. 0.0e0_rk .or. ai(j,i) .ne. 0.0e0_rk) go to 120
  110    continue
!
         m = l
         iexc = 1
         go to 20
  120 continue
!
      go to 140
!     .......... search for columns isolating an eigenvalue
!                and push them left ..........
  130 k = k + 1
!
  140 do 170 j = k, l
!
         do 150 i = k, l
            if (i .eq. j) go to 150
            if (ar(i,j) .ne. 0.0e0_rk .or. ai(i,j) .ne. 0.0e0_rk) go to 170
  150    continue
!
         m = k
         iexc = 2
         go to 20
  170 continue
!     .......... now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 scale(i) = 1.0e0_rk
!     .......... iterative loop for norm reduction ..........
  190 noconv = .false.
!
      do 270 i = k, l
         c = 0.0e0_rk
         r = 0.0e0_rk
!
         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + abs(ar(j,i)) + abs(ai(j,i))
            r = r + abs(ar(i,j)) + abs(ai(i,j))
  200    continue
!     .......... guard against zero c or r due to underflow ..........
         if (c .eq. 0.0e0_rk .or. r .eq. 0.0e0_rk) go to 270
         g = r / radix
         f = 1.0e0_rk
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
!     .......... now balance ..........
  240    if ((c + r) / f .ge. 0.95e0_rk * s) go to 270
         g = 1.0e0_rk / f
         scale(i) = scale(i) * f
         noconv = .true.
!
         do 250 j = k, n
            ar(i,j) = ar(i,j) * g
            ai(i,j) = ai(i,j) * g
  250    continue
!
         do 260 j = 1, l
            ar(j,i) = ar(j,i) * f
            ai(j,i) = ai(j,i) * f
  260    continue
!
  270 continue
!
      if (noconv) go to 190
!
  280 low = k
      igh = l
      return
      end subroutine
!*********************************************************************************
      subroutine cg(nm,n,ar,ai,wr,wi,matz,zr,zi,fv1,fv2,fv3,ierr)
!
      integer n,nm,is1,is2,ierr,matz
      real(kind=rk)  :: ar(nm,n),ai(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),fv1(n),fv2(n),fv3(n)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a complex general matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a=(ar,ai).
!
!        ar  and  ai  contain the real and imaginary parts,
!        respectively, of the complex general matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        wr  and  wi  contain the real and imaginary parts,
!        respectively, of the eigenvalues.
!
!        zr  and  zi  contain the real and imaginary parts,
!        respectively, of the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for comqr
!           and comqr2.  the normal completion code is zero.
!
!        fv1, fv2, and  fv3  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
!
   10 call  cbal(nm,n,ar,ai,is1,is2,fv1)
      call  corth(nm,n,is1,is2,ar,ai,fv2,fv3)
      if (matz .ne. 0) go to 20
!     .......... find eigenvalues only ..........
      call  comqr(nm,n,is1,is2,ar,ai,wr,wi,ierr)
      go to 50
!     .......... find both eigenvalues and eigenvectors ..........
   20 call  comqr2(nm,n,is1,is2,fv2,fv3,ar,ai,wr,wi,zr,zi,ierr)
      if (ierr .ne. 0) go to 50
      call  cbabk2(nm,n,is1,is2,fv1,n,zr,zi)
   50 return
      end subroutine
!*********************************************************************************
      subroutine ch(nm,n,ar,ai,w,matz,zr,zi,fv1,fv2,fm1,ierr)
!
      integer i,j,n,nm,ierr,matz
      real(kind=rk)  :: ar(nm,n),ai(nm,n),w(n),zr(nm,n),zi(nm,n),fv1(n),fv2(n),fm1(2,n)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a complex hermitian matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a=(ar,ai).
!
!        ar  and  ai  contain the real and imaginary parts,
!        respectively, of the complex hermitian matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        zr  and  zi  contain the real and imaginary parts,
!        respectively, of the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and tql2.  the normal completion code is zero.
!
!        fv1, fv2, and  fm1  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
!
   10 call  htridi(nm,n,ar,ai,w,fv1,fv2,fm1)
      if (matz .ne. 0) go to 20
!     .......... find eigenvalues only ..........
      call  tqlrat(n,w,fv2,ierr)
      go to 50
!     .......... find both eigenvalues and eigenvectors ..........
   20 do 40 i = 1, n
!
         do 30 j = 1, n
            zr(j,i) = 0.0e0_rk
   30    continue
!
         zr(i,i) = 1.0e0_rk
   40 continue
!
      call  tql2(nm,n,w,fv1,zr,ierr)
      if (ierr .ne. 0) go to 50
      call  htribk(nm,n,ar,ai,fm1,n,zr,zi)
   50 return
      end subroutine
!*********************************************************************************
      subroutine cinvit(nm,n,ar,ai,wr,wi,select,mm,m,zr,zi,ierr,rm1,rm2,rv1,rv2)
!
      integer i,j,k,m,n,s,ii,mm,mp,nm,uk,ip1,its,km1,ierr
      real(kind=rk)  :: ar(nm,n),ai(nm,n),wr(n),wi(n),zr(nm,mm),zi(nm,mm),rm1(n,n),rm2(n,n),rv1(n),rv2(n)
      real(kind=rk)  :: x,y,eps3,norm,normv,growto,ilambd,rlambd,ukroot
      logical select(n)
!
!     this subroutine is a translation of the algol procedure cx invit
!     by peters and wilkinson.
!     handbook for auto. comp. vol.ii-linear algebra, 418-439(1971).
!
!     this subroutine finds those eigenvectors of a complex upper
!     hessenberg matrix corresponding to specified eigenvalues,
!     using inverse iteration.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the hessenberg matrix.
!
!        wr and wi contain the real and imaginary parts, respectively,
!          of the eigenvalues of the matrix.  the eigenvalues must be
!          stored in a manner identical to that of subroutine  comlr,
!          which recognizes possible splitting of the matrix.
!
!        select specifies the eigenvectors to be found.  the
!          eigenvector corresponding to the j-th eigenvalue is
!          specified by setting select(j) to .true..
!
!        mm should be set to an upper bound for the number of
!          eigenvectors to be found.
!
!     on output
!
!        ar, ai, wi, and select are unaltered.
!
!        wr may have been altered since close eigenvalues are perturbed
!          slightly in searching for independent eigenvectors.
!
!        m is the number of eigenvectors actually found.
!
!        zr and zi contain the real and imaginary parts, respectively,
!          of the eigenvectors.  the eigenvectors are normalized
!          so that the component of largest magnitude is 1.
!          any vector which fails the acceptance test is set to zero.
!
!        ierr is set to
!          zero       for normal return,
!          -(2*n+1)   if more than mm eigenvectors have been specified,
!          -k         if the iteration corresponding to the k-th
!                     value fails,
!          -(n+k)     if both error situations occur.
!
!        rm1, rm2, rv1, and rv2 are temporary storage arrays.
!
!     the algol procedure guessvec appears in cinvit in line.
!
!     calls cdiv for complex division.
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      uk = 0
      s = 1
!
      do 980 k = 1, n
         if (.not. select(k)) go to 980
         if (s .gt. mm) go to 1000
         if (uk .ge. k) go to 200
!     .......... check for possible splitting ..........
         do 120 uk = k, n
            if (uk .eq. n) go to 140
            if (ar(uk+1,uk) .eq. 0.0e0_rk .and. ai(uk+1,uk) .eq. 0.0e0_rk) go to 140
  120    continue
!     .......... compute infinity norm of leading uk by uk
!                (hessenberg) matrix ..........
  140    norm = 0.0e0_rk
         mp = 1
!
         do 180 i = 1, uk
            x = 0.0e0_rk
!
            do 160 j = mp, uk
  160       x = x + pythag(ar(i,j),ai(i,j))
!
            if (x .gt. norm) norm = x
            mp = i
  180    continue
!     .......... eps3 replaces zero pivot in decomposition
!                and close roots are modified by eps3 ..........
         if (norm .eq. 0.0e0_rk) norm = 1.0e0_rk
         eps3 = epslon(norm)
!     .......... growto is the criterion for growth ..........
         ukroot = uk
         ukroot = sqrt(ukroot)
         growto = 0.1e0_rk / ukroot
  200    rlambd = wr(k)
         ilambd = wi(k)
         if (k .eq. 1) go to 280
         km1 = k - 1
         go to 240
!     .......... perturb eigenvalue if it is close
!                to any previous eigenvalue ..........
  220    rlambd = rlambd + eps3
!     .......... for i=k-1 step -1 until 1 do -- ..........
  240    do 260 ii = 1, km1
            i = k - ii
            if (select(i) .and. abs(wr(i)-rlambd) .lt. eps3 .and. abs(wi(i)-ilambd) .lt. eps3) go to 220
  260    continue
!
         wr(k) = rlambd
!     .......... form upper hessenberg (ar,ai)-(rlambd,ilambd)*i
!                and initial complex vector ..........
  280    mp = 1
!
         do 320 i = 1, uk
!
            do 300 j = mp, uk
               rm1(i,j) = ar(i,j)
               rm2(i,j) = ai(i,j)
  300       continue
!
            rm1(i,i) = rm1(i,i) - rlambd
            rm2(i,i) = rm2(i,i) - ilambd
            mp = i
            rv1(i) = eps3
  320    continue
!     .......... triangular decomposition with interchanges,
!                replacing zero pivots by eps3 ..........
         if (uk .eq. 1) go to 420
!
         do 400 i = 2, uk
            mp = i - 1
            if (pythag(rm1(i,mp),rm2(i,mp)) .le. pythag(rm1(mp,mp),rm2(mp,mp))) go to 360
!
            do 340 j = mp, uk
               y = rm1(i,j)
               rm1(i,j) = rm1(mp,j)
               rm1(mp,j) = y
               y = rm2(i,j)
               rm2(i,j) = rm2(mp,j)
               rm2(mp,j) = y
  340       continue
!
  360       if (rm1(mp,mp) .eq. 0.0e0_rk .and. rm2(mp,mp) .eq. 0.0e0_rk) rm1(mp,mp) = eps3
            call cdiv(rm1(i,mp),rm2(i,mp),rm1(mp,mp),rm2(mp,mp),x,y)
            if (x .eq. 0.0e0_rk .and. y .eq. 0.0e0_rk) go to 400
!
            do 380 j = i, uk
               rm1(i,j) = rm1(i,j) - x * rm1(mp,j) + y * rm2(mp,j)
               rm2(i,j) = rm2(i,j) - x * rm2(mp,j) - y * rm1(mp,j)
  380       continue
!
  400    continue
!
  420    if (rm1(uk,uk) .eq. 0.0e0_rk .and. rm2(uk,uk) .eq. 0.0e0_rk) rm1(uk,uk) = eps3
         its = 0
!     .......... back substitution
!                for i=uk step -1 until 1 do -- ..........
  660    do 720 ii = 1, uk
            i = uk + 1 - ii
            x = rv1(i)
            y = 0.0e0_rk
            if (i .eq. uk) go to 700
            ip1 = i + 1
!
            do 680 j = ip1, uk
               x = x - rm1(i,j) * rv1(j) + rm2(i,j) * rv2(j)
               y = y - rm1(i,j) * rv2(j) - rm2(i,j) * rv1(j)
  680       continue
!
  700       call cdiv(x,y,rm1(i,i),rm2(i,i),rv1(i),rv2(i))
  720    continue
!     .......... acceptance test for eigenvector
!                and normalization ..........
         its = its + 1
         norm = 0.0e0_rk
         normv = 0.0e0_rk
!
         do 780 i = 1, uk
            x = pythag(rv1(i),rv2(i))
            if (normv .ge. x) go to 760
            normv = x
            j = i
  760       norm = norm + x
  780    continue
!
         if (norm .lt. growto) go to 840
!     .......... accept vector ..........
         x = rv1(j)
         y = rv2(j)
!
         do 820 i = 1, uk
            call cdiv(rv1(i),rv2(i),x,y,zr(i,s),zi(i,s))
  820    continue
!
         if (uk .eq. n) go to 940
         j = uk + 1
         go to 900
!     .......... in-line procedure for choosing
!                a new starting vector ..........
  840    if (its .ge. uk) go to 880
         x = ukroot
         y = eps3 / (x + 1.0e0_rk)
         rv1(1) = eps3
!
         do 860 i = 2, uk
  860    rv1(i) = y
!
         j = uk - its + 1
         rv1(j) = rv1(j) - eps3 * x
         go to 660
!     .......... set error -- unaccepted eigenvector ..........
  880    j = 1
         ierr = -k
!     .......... set remaining vector components to zero ..........
  900    do 920 i = j, n
            zr(i,s) = 0.0e0_rk
            zi(i,s) = 0.0e0_rk
  920    continue
!
  940    s = s + 1
  980 continue
!
      go to 1001
!     .......... set error -- underestimate of eigenvector
!                space required ..........
 1000 if (ierr .ne. 0) ierr = ierr - n
      if (ierr .eq. 0) ierr = -(2 * n + 1)
 1001 m = s - 1
      return
      end subroutine
!*********************************************************************************
      subroutine combak(nm,low,igh,ar,ai,int,m,zr,zi)
!
      integer i,j,m,la,mm,mp,nm,igh,kp1,low,mp1
      real(kind=rk)  :: ar(nm,igh),ai(nm,igh),zr(nm,m),zi(nm,m)
      real(kind=rk)  :: xr,xi
      integer int(igh)
!
!     this subroutine is a translation of the algol procedure combak,
!     num. math. 12, 349-368(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     this subroutine forms the eigenvectors of a complex general
!     matrix by back transforming those of the corresponding
!     upper hessenberg matrix determined by  comhes.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        low and igh are integers determined by the balancing
!          subroutine  cbal.  if  cbal  has not been used,
!          set low=1 and igh equal to the order of the matrix.
!
!        ar and ai contain the multipliers which were used in the
!          reduction by  comhes  in their lower triangles
!          below the subdiagonal.
!
!        int contains information on the rows and columns
!          interchanged in the reduction by  comhes.
!          only elements low through igh are used.
!
!        m is the number of eigenvectors to be back transformed.
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the eigenvectors to be
!          back transformed in their first m columns.
!
!     on output
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the transformed eigenvectors
!          in their first m columns.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (m .eq. 0) go to 200
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
!     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = kp1, la
         mp = low + igh - mm
         mp1 = mp + 1
!
         do 110 i = mp1, igh
            xr = ar(i,mp-1)
            xi = ai(i,mp-1)
            if (xr .eq. 0.0e0_rk .and. xi .eq. 0.0e0_rk) go to 110
!
            do 100 j = 1, m
               zr(i,j) = zr(i,j) + xr * zr(mp,j) - xi * zi(mp,j)
               zi(i,j) = zi(i,j) + xr * zi(mp,j) + xi * zr(mp,j)
  100       continue
!
  110    continue
!
         i = int(mp)
         if (i .eq. mp) go to 140
!
         do 130 j = 1, m
            xr = zr(i,j)
            zr(i,j) = zr(mp,j)
            zr(mp,j) = xr
            xi = zi(i,j)
            zi(i,j) = zi(mp,j)
            zi(mp,j) = xi
  130    continue
!
  140 continue
!
  200 return
      end subroutine
!*********************************************************************************
      subroutine comhes(nm,n,low,igh,ar,ai,int)
!
      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
      real(kind=rk)  :: ar(nm,n),ai(nm,n)
      real(kind=rk)  :: xr,xi,yr,yi
      integer int(igh)
!
!     this subroutine is a translation of the algol procedure comhes,
!     num. math. 12, 349-368(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     given a complex general matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     low through igh to upper hessenberg form by
!     stabilized elementary similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  cbal.  if  cbal  has not been used,
!          set low=1, igh=n.
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the complex input matrix.
!
!     on output
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the hessenberg matrix.  the
!          multipliers which were used in the reduction
!          are stored in the remaining triangles under the
!          hessenberg matrix.
!
!        int contains information on the rows and columns
!          interchanged in the reduction.
!          only elements low through igh are used.
!
!     calls cdiv for complex division.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
!
      do 180 m = kp1, la
         mm1 = m - 1
         xr = 0.0e0_rk
         xi = 0.0e0_rk
         i = m
!
         do 100 j = m, igh
            if (abs(ar(j,mm1)) + abs(ai(j,mm1)) .le. abs(xr) + abs(xi)) go to 100
            xr = ar(j,mm1)
            xi = ai(j,mm1)
            i = j
  100    continue
!
         int(m) = i
         if (i .eq. m) go to 130
!     .......... interchange rows and columns of ar and ai ..........
         do 110 j = mm1, n
            yr = ar(i,j)
            ar(i,j) = ar(m,j)
            ar(m,j) = yr
            yi = ai(i,j)
            ai(i,j) = ai(m,j)
            ai(m,j) = yi
  110    continue
!
         do 120 j = 1, igh
            yr = ar(j,i)
            ar(j,i) = ar(j,m)
            ar(j,m) = yr
            yi = ai(j,i)
            ai(j,i) = ai(j,m)
            ai(j,m) = yi
  120    continue
!     .......... end interchange ..........
  130    if (xr .eq. 0.0e0_rk .and. xi .eq. 0.0e0_rk) go to 180
         mp1 = m + 1
!
         do 160 i = mp1, igh
            yr = ar(i,mm1)
            yi = ai(i,mm1)
            if (yr .eq. 0.0e0_rk .and. yi .eq. 0.0e0_rk) go to 160
            call cdiv(yr,yi,xr,xi,yr,yi)
            ar(i,mm1) = yr
            ai(i,mm1) = yi
!
            do 140 j = m, n
               ar(i,j) = ar(i,j) - yr * ar(m,j) + yi * ai(m,j)
               ai(i,j) = ai(i,j) - yr * ai(m,j) - yi * ar(m,j)
  140       continue
!
            do 150 j = 1, igh
               ar(j,m) = ar(j,m) + yr * ar(j,i) - yi * ai(j,i)
               ai(j,m) = ai(j,m) + yr * ai(j,i) + yi * ar(j,i)
  150       continue
!
  160    continue
!
  180 continue
!
  200 return
      end subroutine
!*********************************************************************************
      subroutine comlr(nm,n,low,igh,hr,hi,wr,wi,ierr)
!
      integer i,j,l,m,n,en,ll,mm,nm,igh,im1,itn,its,low,mp1,enm1,ierr
      real(kind=rk)  :: hr(nm,n),hi(nm,n),wr(n),wi(n)
      real(kind=rk)  :: si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,tst1,tst2
!
!     this subroutine is a translation of the algol procedure comlr,
!     num. math. 12, 369-376(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
!
!     this subroutine finds the eigenvalues of a complex
!     upper hessenberg matrix by the modified lr method.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  cbal.  if  cbal  has not been used,
!          set low=1, igh=n.
!
!        hr and hi contain the real and imaginary parts,
!          respectively, of the complex upper hessenberg matrix.
!          their lower triangles below the subdiagonal contain the
!          multipliers which were used in the reduction by  comhes,
!          if performed.
!
!     on output
!
!        the upper hessenberg portions of hr and hi have been
!          destroyed.  therefore, they must be saved before
!          calling  comlr  if subsequent calculation of
!          eigenvectors is to be performed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  if an error
!          exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
!     calls cdiv for complex division.
!     calls csroot for complex square root.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
!     .......... store roots isolated by cbal ..........
      do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
!
      en = igh
      tr = 0.0e0_rk
      ti = 0.0e0_rk
      itn = 30*n
!     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 1001
      its = 0
      enm1 = en - 1
!     .......... look for single small sub-diagonal element
!                for l=en step -1 until low d0 -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1)) + abs(hr(l,l)) + abs(hi(l,l))
         tst2 = tst1 + abs(hr(l,l-1)) + abs(hi(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
!     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1) - hi(enm1,en) * hi(en,enm1)
      xi = hr(enm1,en) * hi(en,enm1) + hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0e0_rk .and. xi .eq. 0.0e0_rk) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0e0_rk
      yi = (hi(enm1,enm1) - si) / 2.0e0_rk
      call csroot(yr**2-yi**2+xr,2.0e0_rk*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0e0_rk) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
!     .......... form exceptional shift ..........
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
      si = abs(hi(en,enm1)) + abs(hi(enm1,en-2))
!
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
!
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
!     .......... look for two consecutive small
!                sub-diagonal elements ..........
      xr = abs(hr(enm1,enm1)) + abs(hi(enm1,enm1))
      yr = abs(hr(en,enm1)) + abs(hi(en,enm1))
      zzr = abs(hr(en,en)) + abs(hi(en,en))
!     .......... for m=en-1 step -1 until l do -- ..........
      do 380 mm = l, enm1
         m = enm1 + l - mm
         if (m .eq. l) go to 420
         yi = yr
         yr = abs(hr(m,m-1)) + abs(hi(m,m-1))
         xi = zzr
         zzr = xr
         xr = abs(hr(m-1,m-1)) + abs(hi(m-1,m-1))
         tst1 = zzr / yi * (zzr + xr + xi)
         tst2 = tst1 + yr
         if (tst2 .eq. tst1) go to 420
  380 continue
!     .......... triangular decomposition h=l*r ..........
  420 mp1 = m + 1
!
      do 520 i = mp1, en
         im1 = i - 1
         xr = hr(im1,im1)
         xi = hi(im1,im1)
         yr = hr(i,im1)
         yi = hi(i,im1)
         if (abs(xr) + abs(xi) .ge. abs(yr) + abs(yi)) go to 460
!     .......... interchange rows of hr and hi ..........
         do 440 j = im1, en
            zzr = hr(im1,j)
            hr(im1,j) = hr(i,j)
            hr(i,j) = zzr
            zzi = hi(im1,j)
            hi(im1,j) = hi(i,j)
            hi(i,j) = zzi
  440    continue
!
         call cdiv(xr,xi,yr,yi,zzr,zzi)
         wr(i) = 1.0e0_rk
         go to 480
  460    call cdiv(yr,yi,xr,xi,zzr,zzi)
         wr(i) = -1.0e0_rk
  480    hr(i,im1) = zzr
         hi(i,im1) = zzi
!
         do 500 j = i, en
            hr(i,j) = hr(i,j) - zzr * hr(im1,j) + zzi * hi(im1,j)
            hi(i,j) = hi(i,j) - zzr * hi(im1,j) - zzi * hr(im1,j)
  500    continue
!
  520 continue
!     .......... composition r*l=h ..........
      do 640 j = mp1, en
         xr = hr(j,j-1)
         xi = hi(j,j-1)
         hr(j,j-1) = 0.0e0_rk
         hi(j,j-1) = 0.0e0_rk
!     .......... interchange columns of hr and hi,
!                if necessary ..........
         if (wr(j) .le. 0.0e0_rk) go to 580
!
         do 540 i = l, j
            zzr = hr(i,j-1)
            hr(i,j-1) = hr(i,j)
            hr(i,j) = zzr
            zzi = hi(i,j-1)
            hi(i,j-1) = hi(i,j)
            hi(i,j) = zzi
  540    continue
!
  580    do 600 i = l, j
            hr(i,j-1) = hr(i,j-1) + xr * hr(i,j) - xi * hi(i,j)
            hi(i,j-1) = hi(i,j-1) + xr * hi(i,j) + xi * hr(i,j)
  600    continue
!
  640 continue
!
      go to 240
!     .......... a root found ..........
  660 wr(en) = hr(en,en) + tr
      wi(en) = hi(en,en) + ti
      en = enm1
      go to 220
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end subroutine
!*********************************************************************************
      subroutine comlr2(nm,n,low,igh,int,hr,hi,wr,wi,zr,zi,ierr)
!
      integer i,j,k,l,m,n,en,ii,jj,ll,mm,nm,nn,igh,im1,ip1,itn,its,low,mp1,enm1,iend,ierr
      real(kind=rk)  :: hr(nm,n),hi(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n)
      real(kind=rk)  :: si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2
      integer int(igh)
!
!     this subroutine is a translation of the algol procedure comlr2,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a complex upper hessenberg matrix by the modified lr
!     method.  the eigenvectors of a complex general matrix
!     can also be found if  comhes  has been used to reduce
!     this general matrix to hessenberg form.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  cbal.  if  cbal  has not been used,
!          set low=1, igh=n.
!
!        int contains information on the rows and columns interchanged
!          in the reduction by  comhes, if performed.  only elements
!          low through igh are used.  if the eigenvectors of the hessen-
!          berg matrix are desired, set int(j)=j for these elements.
!
!        hr and hi contain the real and imaginary parts,
!          respectively, of the complex upper hessenberg matrix.
!          their lower triangles below the subdiagonal contain the
!          multipliers which were used in the reduction by  comhes,
!          if performed.  if the eigenvectors of the hessenberg
!          matrix are desired, these elements must be set to zero.
!
!     on output
!
!        the upper hessenberg portions of hr and hi have been
!          destroyed, but the location hr(1,1) contains the norm
!          of the triangularized matrix.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  if an error
!          exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the eigenvectors.  the eigenvectors
!          are unnormalized.  if an error exit is made, none of
!          the eigenvectors has been found.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
!
!     calls cdiv for complex division.
!     calls csroot for complex square root.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
!     .......... initialize eigenvector matrix ..........
      do 100 i = 1, n
!
         do 100 j = 1, n
            zr(i,j) = 0.0e0_rk
            zi(i,j) = 0.0e0_rk
            if (i .eq. j) zr(i,j) = 1.0e0_rk
  100 continue
!     .......... form the matrix of accumulated transformations
!                from the information left by comhes ..........
      iend = igh - low - 1
      if (iend .le. 0) go to 180
!     .......... for i=igh-1 step -1 until low+1 do -- ..........
      do 160 ii = 1, iend
         i = igh - ii
         ip1 = i + 1
!
         do 120 k = ip1, igh
            zr(k,i) = hr(k,i-1)
            zi(k,i) = hi(k,i-1)
  120    continue
!
         j = int(i)
         if (i .eq. j) go to 160
!
         do 140 k = i, igh
            zr(i,k) = zr(j,k)
            zi(i,k) = zi(j,k)
            zr(j,k) = 0.0e0_rk
            zi(j,k) = 0.0e0_rk
  140    continue
!
         zr(j,i) = 1.0e0_rk
  160 continue
!     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
!
      en = igh
      tr = 0.0e0_rk
      ti = 0.0e0_rk
      itn = 30*n
!     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 680
      its = 0
      enm1 = en - 1
!     .......... look for single small sub-diagonal element
!                for l=en step -1 until low do -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1)) + abs(hr(l,l)) + abs(hi(l,l))
         tst2 = tst1 + abs(hr(l,l-1)) + abs(hi(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
!     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1) - hi(enm1,en) * hi(en,enm1)
      xi = hr(enm1,en) * hi(en,enm1) + hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0e0_rk .and. xi .eq. 0.0e0_rk) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0e0_rk
      yi = (hi(enm1,enm1) - si) / 2.0e0_rk
      call csroot(yr**2-yi**2+xr,2.0e0_rk*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0e0_rk) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
!     .......... form exceptional shift ..........
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
      si = abs(hi(en,enm1)) + abs(hi(enm1,en-2))
!
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
!
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
!     .......... look for two consecutive small
!                sub-diagonal elements ..........
      xr = abs(hr(enm1,enm1)) + abs(hi(enm1,enm1))
      yr = abs(hr(en,enm1)) + abs(hi(en,enm1))
      zzr = abs(hr(en,en)) + abs(hi(en,en))
!     .......... for m=en-1 step -1 until l do -- ..........
      do 380 mm = l, enm1
         m = enm1 + l - mm
         if (m .eq. l) go to 420
         yi = yr
         yr = abs(hr(m,m-1)) + abs(hi(m,m-1))
         xi = zzr
         zzr = xr
         xr = abs(hr(m-1,m-1)) + abs(hi(m-1,m-1))
         tst1 = zzr / yi * (zzr + xr + xi)
         tst2 = tst1 + yr
         if (tst2 .eq. tst1) go to 420
  380 continue
!     .......... triangular decomposition h=l*r ..........
  420 mp1 = m + 1
!
      do 520 i = mp1, en
         im1 = i - 1
         xr = hr(im1,im1)
         xi = hi(im1,im1)
         yr = hr(i,im1)
         yi = hi(i,im1)
         if (abs(xr) + abs(xi) .ge. abs(yr) + abs(yi)) go to 460
!     .......... interchange rows of hr and hi ..........
         do 440 j = im1, n
            zzr = hr(im1,j)
            hr(im1,j) = hr(i,j)
            hr(i,j) = zzr
            zzi = hi(im1,j)
            hi(im1,j) = hi(i,j)
            hi(i,j) = zzi
  440    continue
!
         call cdiv(xr,xi,yr,yi,zzr,zzi)
         wr(i) = 1.0e0_rk
         go to 480
  460    call cdiv(yr,yi,xr,xi,zzr,zzi)
         wr(i) = -1.0e0_rk
  480    hr(i,im1) = zzr
         hi(i,im1) = zzi
!
         do 500 j = i, n
            hr(i,j) = hr(i,j) - zzr * hr(im1,j) + zzi * hi(im1,j)
            hi(i,j) = hi(i,j) - zzr * hi(im1,j) - zzi * hr(im1,j)
  500    continue
!
  520 continue
!     .......... composition r*l=h ..........
      do 640 j = mp1, en
         xr = hr(j,j-1)
         xi = hi(j,j-1)
         hr(j,j-1) = 0.0e0_rk
         hi(j,j-1) = 0.0e0_rk
!     .......... interchange columns of hr, hi, zr, and zi,
!                if necessary ..........
         if (wr(j) .le. 0.0e0_rk) go to 580
!
         do 540 i = 1, j
            zzr = hr(i,j-1)
            hr(i,j-1) = hr(i,j)
            hr(i,j) = zzr
            zzi = hi(i,j-1)
            hi(i,j-1) = hi(i,j)
            hi(i,j) = zzi
  540    continue
!
         do 560 i = low, igh
            zzr = zr(i,j-1)
            zr(i,j-1) = zr(i,j)
            zr(i,j) = zzr
            zzi = zi(i,j-1)
            zi(i,j-1) = zi(i,j)
            zi(i,j) = zzi
  560    continue
!
  580    do 600 i = 1, j
            hr(i,j-1) = hr(i,j-1) + xr * hr(i,j) - xi * hi(i,j)
            hi(i,j-1) = hi(i,j-1) + xr * hi(i,j) + xi * hr(i,j)
  600    continue
!     .......... accumulate transformations ..........
         do 620 i = low, igh
            zr(i,j-1) = zr(i,j-1) + xr * zr(i,j) - xi * zi(i,j)
            zi(i,j-1) = zi(i,j-1) + xr * zi(i,j) + xi * zr(i,j)
  620    continue
!
  640 continue
!
      go to 240
!     .......... a root found ..........
  660 hr(en,en) = hr(en,en) + tr
      wr(en) = hr(en,en)
      hi(en,en) = hi(en,en) + ti
      wi(en) = hi(en,en)
      en = enm1
      go to 220
!     .......... all roots found.  backsubstitute to find
!                vectors of upper triangular form ..........
  680 norm = 0.0e0_rk
!
      do 720 i = 1, n
!
         do 720 j = i, n
            tr = abs(hr(i,j)) + abs(hi(i,j))
            if (tr .gt. norm) norm = tr
  720 continue
!
      hr(1,1) = norm
      if (n .eq. 1 .or. norm .eq. 0.0e0_rk) go to 1001
!     .......... for en=n step -1 until 2 do -- ..........
      do 800 nn = 2, n
         en = n + 2 - nn
         xr = wr(en)
         xi = wi(en)
         hr(en,en) = 1.0e0_rk
         hi(en,en) = 0.0e0_rk
         enm1 = en - 1
!     .......... for i=en-1 step -1 until 1 do -- ..........
         do 780 ii = 1, enm1
            i = en - ii
            zzr = 0.0e0_rk
            zzi = 0.0e0_rk
            ip1 = i + 1
!
            do 740 j = ip1, en
               zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
               zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
  740       continue
!
            yr = xr - wr(i)
            yi = xi - wi(i)
            if (yr .ne. 0.0e0_rk .or. yi .ne. 0.0e0_rk) go to 765
               tst1 = norm
               yr = tst1
  760          yr = 0.01e0_rk * yr
               tst2 = norm + yr
               if (tst2 .gt. tst1) go to 760
  765       continue
            call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
!     .......... overflow control ..........
            tr = abs(hr(i,en)) + abs(hi(i,en))
            if (tr .eq. 0.0e0_rk) go to 780
            tst1 = tr
            tst2 = tst1 + 1.0e0_rk/tst1
            if (tst2 .gt. tst1) go to 780
            do 770 j = i, en
               hr(j,en) = hr(j,en)/tr
               hi(j,en) = hi(j,en)/tr
  770       continue
!
  780    continue
!
  800 continue
!     .......... end backsubstitution ..........
      enm1 = n - 1
!     .......... vectors of isolated roots ..........
      do  840 i = 1, enm1
         if (i .ge. low .and. i .le. igh) go to 840
         ip1 = i + 1
!
         do 820 j = ip1, n
            zr(i,j) = hr(i,j)
            zi(i,j) = hi(i,j)
  820    continue
!
  840 continue
!     .......... multiply by transformation matrix to give
!                vectors of original full matrix.
!                for j=n step -1 until low+1 do -- ..........
      do 880 jj = low, enm1
         j = n + low - jj
         m = min0(j,igh)
!
         do 880 i = low, igh
            zzr = 0.0e0_rk
            zzi = 0.0e0_rk
!
            do 860 k = low, m
               zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
               zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
  860       continue
!
            zr(i,j) = zzr
            zi(i,j) = zzi
  880 continue
!
      go to 1001
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end subroutine
!*********************************************************************************
      subroutine comqr(nm,n,low,igh,hr,hi,wr,wi,ierr)
!
      integer i,j,l,n,en,ll,nm,igh,itn,its,low,lp1,enm1,ierr
      real(kind=rk)  :: hr(nm,n),hi(nm,n),wr(n),wi(n)
      real(kind=rk)  :: si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2
!
!     this subroutine is a translation of a unitary analogue of the
!     algol procedure  comlr, num. math. 12, 369-376(1968) by martin
!     and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
!     the unitary analogue substitutes the qr algorithm of francis
!     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
!
!     this subroutine finds the eigenvalues of a complex
!     upper hessenberg matrix by the qr method.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  cbal.  if  cbal  has not been used,
!          set low=1, igh=n.
!
!        hr and hi contain the real and imaginary parts,
!          respectively, of the complex upper hessenberg matrix.
!          their lower triangles below the subdiagonal contain
!          information about the unitary transformations used in
!          the reduction by  corth, if performed.
!
!     on output
!
!        the upper hessenberg portions of hr and hi have been
!          destroyed.  therefore, they must be saved before
!          calling  comqr  if subsequent calculation of
!          eigenvectors is to be performed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  if an error
!          exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
!     calls cdiv for complex division.
!     calls csroot for complex square root.
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      if (low .eq. igh) go to 180
!     .......... create real subdiagonal elements ..........
      l = low + 1
!
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0e0_rk) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0e0_rk
!
         do 155 j = i, igh
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
!
         do 160 j = low, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
!
  170 continue
!     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
!
      en = igh
      tr = 0.0e0_rk
      ti = 0.0e0_rk
      itn = 30*n
!     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 1001
      its = 0
      enm1 = en - 1
!     .......... look for single small sub-diagonal element
!                for l=en step -1 until low d0 -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1)) + abs(hr(l,l)) + abs(hi(l,l))
         tst2 = tst1 + abs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
!     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0e0_rk .and. xi .eq. 0.0e0_rk) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0e0_rk
      yi = (hi(enm1,enm1) - si) / 2.0e0_rk
      call csroot(yr**2-yi**2+xr,2.0e0_rk*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0e0_rk) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
!     .......... form exceptional shift ..........
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
      si = 0.0e0_rk
!
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
!
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
!     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
!
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0e0_rk
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0e0_rk
         hi(i,i-1) = sr / norm
!
         do 490 j = i, en
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
!
  500 continue
!
      si = hi(en,en)
      if (si .eq. 0.0e0_rk) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0e0_rk
!     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
!
         do 580 i = l, j
            yr = hr(i,j-1)
            yi = 0.0e0_rk
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
!
  600 continue
!
      if (si .eq. 0.0e0_rk) go to 240
!
      do 630 i = l, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
!
      go to 240
!     .......... a root found ..........
  660 wr(en) = hr(en,en) + tr
      wi(en) = hi(en,en) + ti
      en = enm1
      go to 220
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end subroutine
!*********************************************************************************
      subroutine comqr2(nm,n,low,igh,ortr,orti,hr,hi,wr,wi,zr,zi,ierr)
!
      integer i,j,k,l,m,n,en,ii,jj,ll,nm,nn,igh,ip1, itn,its,low,lp1,enm1,iend,ierr
      real(kind=rk)  :: hr(nm,n),hi(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n), ortr(igh),orti(igh)
      real(kind=rk)  :: si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2
!
!     this subroutine is a translation of a unitary analogue of the
!     algol procedure  comlr2, num. math. 16, 181-204(1970) by peters
!     and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!     the unitary analogue substitutes the qr algorithm of francis
!     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a complex upper hessenberg matrix by the qr
!     method.  the eigenvectors of a complex general matrix
!     can also be found if  corth  has been used to reduce
!     this general matrix to hessenberg form.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  cbal.  if  cbal  has not been used,
!          set low=1, igh=n.
!
!        ortr and orti contain information about the unitary trans-
!          formations used in the reduction by  corth, if performed.
!          only elements low through igh are used.  if the eigenvectors
!          of the hessenberg matrix are desired, set ortr(j) and
!          orti(j) to 0.0e0_rk for these elements.
!
!        hr and hi contain the real and imaginary parts,
!          respectively, of the complex upper hessenberg matrix.
!          their lower triangles below the subdiagonal contain further
!          information about the transformations which were used in the
!          reduction by  corth, if performed.  if the eigenvectors of
!          the hessenberg matrix are desired, these elements may be
!          arbitrary.
!
!     on output
!
!        ortr, orti, and the upper hessenberg portions of hr and hi
!          have been destroyed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  if an error
!          exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the eigenvectors.  the eigenvectors
!          are unnormalized.  if an error exit is made, none of
!          the eigenvectors has been found.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
!     calls cdiv for complex division.
!     calls csroot for complex square root.
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
!     .......... initialize eigenvector matrix ..........
      do 101 j = 1, n
!
         do 100 i = 1, n
            zr(i,j) = 0.0e0_rk
            zi(i,j) = 0.0e0_rk
  100    continue
         zr(j,j) = 1.0e0_rk
  101 continue
!     .......... form the matrix of accumulated transformations
!                from the information left by corth ..........
      iend = igh - low - 1
      if (iend) 180, 150, 105
!     .......... for i=igh-1 step -1 until low+1 do -- ..........
  105 do 140 ii = 1, iend
         i = igh - ii
         if (ortr(i) .eq. 0.0e0_rk .and. orti(i) .eq. 0.0e0_rk) go to 140
         if (hr(i,i-1) .eq. 0.0e0_rk .and. hi(i,i-1) .eq. 0.0e0_rk) go to 140
!     .......... norm below is negative of h formed in corth ..........
         norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)
         ip1 = i + 1
!
         do 110 k = ip1, igh
            ortr(k) = hr(k,i-1)
            orti(k) = hi(k,i-1)
  110    continue
!
         do 130 j = i, igh
            sr = 0.0e0_rk
            si = 0.0e0_rk
!
            do 115 k = i, igh
               sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
               si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
  115       continue
!
            sr = sr / norm
            si = si / norm
!
            do 120 k = i, igh
               zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
               zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
  120       continue
!
  130    continue
!
  140 continue
!     .......... create real subdiagonal elements ..........
  150 l = low + 1
!
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0e0_rk) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0e0_rk
!
         do 155 j = i, n
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
!
         do 160 j = 1, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
!
         do 165 j = low, igh
            si = yr * zi(j,i) + yi * zr(j,i)
            zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
            zi(j,i) = si
  165    continue
!
  170 continue
!     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
!
      en = igh
      tr = 0.0e0_rk
      ti = 0.0e0_rk
      itn = 30*n
!     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 680
      its = 0
      enm1 = en - 1
!     .......... look for single small sub-diagonal element
!                for l=en step -1 until low do -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1)) + abs(hr(l,l)) + abs(hi(l,l))
         tst2 = tst1 + abs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
!     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0e0_rk .and. xi .eq. 0.0e0_rk) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0e0_rk
      yi = (hi(enm1,enm1) - si) / 2.0e0_rk
      call csroot(yr**2-yi**2+xr,2.0e0_rk*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0e0_rk) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
!     .......... form exceptional shift ..........
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
      si = 0.0e0_rk
!
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
!
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
!     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
!
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0e0_rk
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0e0_rk
         hi(i,i-1) = sr / norm
!
         do 490 j = i, n
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
!
  500 continue
!
      si = hi(en,en)
      if (si .eq. 0.0e0_rk) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0e0_rk
      if (en .eq. n) go to 540
      ip1 = en + 1
!
      do 520 j = ip1, n
         yr = hr(en,j)
         yi = hi(en,j)
         hr(en,j) = sr * yr + si * yi
         hi(en,j) = sr * yi - si * yr
  520 continue
!     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
!
         do 580 i = 1, j
            yr = hr(i,j-1)
            yi = 0.0e0_rk
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
!
         do 590 i = low, igh
            yr = zr(i,j-1)
            yi = zi(i,j-1)
            zzr = zr(i,j)
            zzi = zi(i,j)
            zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
            zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  590    continue
!
  600 continue
!
      if (si .eq. 0.0e0_rk) go to 240
!
      do 630 i = 1, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
!
      do 640 i = low, igh
         yr = zr(i,en)
         yi = zi(i,en)
         zr(i,en) = sr * yr - si * yi
         zi(i,en) = sr * yi + si * yr
  640 continue
!
      go to 240
!     .......... a root found ..........
  660 hr(en,en) = hr(en,en) + tr
      wr(en) = hr(en,en)
      hi(en,en) = hi(en,en) + ti
      wi(en) = hi(en,en)
      en = enm1
      go to 220
!     .......... all roots found.  backsubstitute to find
!                vectors of upper triangular form ..........
  680 norm = 0.0e0_rk
!
      do 720 i = 1, n
!
         do 720 j = i, n
            tr = abs(hr(i,j)) + abs(hi(i,j))
            if (tr .gt. norm) norm = tr
  720 continue
!
      if (n .eq. 1 .or. norm .eq. 0.0e0_rk) go to 1001
!     .......... for en=n step -1 until 2 do -- ..........
      do 800 nn = 2, n
         en = n + 2 - nn
         xr = wr(en)
         xi = wi(en)
         hr(en,en) = 1.0e0_rk
         hi(en,en) = 0.0e0_rk
         enm1 = en - 1
!     .......... for i=en-1 step -1 until 1 do -- ..........
         do 780 ii = 1, enm1
            i = en - ii
            zzr = 0.0e0_rk
            zzi = 0.0e0_rk
            ip1 = i + 1
!
            do 740 j = ip1, en
               zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
               zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
  740       continue
!
            yr = xr - wr(i)
            yi = xi - wi(i)
            if (yr .ne. 0.0e0_rk .or. yi .ne. 0.0e0_rk) go to 765
               tst1 = norm
               yr = tst1
  760          yr = 0.01e0_rk * yr
               tst2 = norm + yr
               if (tst2 .gt. tst1) go to 760
  765       continue
            call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
!     .......... overflow control ..........
            tr = abs(hr(i,en)) + abs(hi(i,en))
            if (tr .eq. 0.0e0_rk) go to 780
            tst1 = tr
            tst2 = tst1 + 1.0e0_rk/tst1
            if (tst2 .gt. tst1) go to 780
            do 770 j = i, en
               hr(j,en) = hr(j,en)/tr
               hi(j,en) = hi(j,en)/tr
  770       continue
!
  780    continue
!
  800 continue
!     .......... end backsubstitution ..........
      enm1 = n - 1
!     .......... vectors of isolated roots ..........
      do  840 i = 1, enm1
         if (i .ge. low .and. i .le. igh) go to 840
         ip1 = i + 1
!
         do 820 j = ip1, n
            zr(i,j) = hr(i,j)
            zi(i,j) = hi(i,j)
  820    continue
!
  840 continue
!     .......... multiply by transformation matrix to give
!                vectors of original full matrix.
!                for j=n step -1 until low+1 do -- ..........
      do 880 jj = low, enm1
         j = n + low - jj
         m = min0(j,igh)
!
         do 880 i = low, igh
            zzr = 0.0e0_rk
            zzi = 0.0e0_rk
!
            do 860 k = low, m
               zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
               zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
  860       continue
!
            zr(i,j) = zzr
            zi(i,j) = zzi
  880 continue
!
      go to 1001
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end subroutine
!*********************************************************************************
      subroutine cortb(nm,low,igh,ar,ai,ortr,orti,m,zr,zi)
!
      integer i,j,m,la,mm,mp,nm,igh,kp1,low,mp1
      real(kind=rk)  :: ar(nm,igh),ai(nm,igh),ortr(igh),orti(igh),zr(nm,m),zi(nm,m)
      real(kind=rk)  :: h,gi,gr
!
!     this subroutine is a translation of a complex analogue of
!     the algol procedure ortbak, num. math. 12, 349-368(1968)
!     by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     this subroutine forms the eigenvectors of a complex general
!     matrix by back transforming those of the corresponding
!     upper hessenberg matrix determined by  corth.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        low and igh are integers determined by the balancing
!          subroutine  cbal.  if  cbal  has not been used,
!          set low=1 and igh equal to the order of the matrix.
!
!        ar and ai contain information about the unitary
!          transformations used in the reduction by  corth
!          in their strict lower triangles.
!
!        ortr and orti contain further information about the
!          transformations used in the reduction by  corth.
!          only elements low through igh are used.
!
!        m is the number of columns of zr and zi to be back transformed.
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the eigenvectors to be
!          back transformed in their first m columns.
!
!     on output
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the transformed eigenvectors
!          in their first m columns.
!
!        ortr and orti have been altered.
!
!     note that cortb preserves vector euclidean norms.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (m .eq. 0) go to 200
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
!     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = kp1, la
         mp = low + igh - mm
         if (ar(mp,mp-1) .eq. 0.0e0_rk .and. ai(mp,mp-1) .eq. 0.0e0_rk) go to 140
!     .......... h below is negative of h formed in corth ..........
         h = ar(mp,mp-1) * ortr(mp) + ai(mp,mp-1) * orti(mp)
         mp1 = mp + 1
!
         do 100 i = mp1, igh
            ortr(i) = ar(i,mp-1)
            orti(i) = ai(i,mp-1)
  100    continue
!
         do 130 j = 1, m
            gr = 0.0e0_rk
            gi = 0.0e0_rk
!
            do 110 i = mp, igh
               gr = gr + ortr(i) * zr(i,j) + orti(i) * zi(i,j)
               gi = gi + ortr(i) * zi(i,j) - orti(i) * zr(i,j)
  110       continue
!
            gr = gr / h
            gi = gi / h
!
            do 120 i = mp, igh
               zr(i,j) = zr(i,j) + gr * ortr(i) - gi * orti(i)
               zi(i,j) = zi(i,j) + gr * orti(i) + gi * ortr(i)
  120       continue
!
  130    continue
!
  140 continue
!
  200 return
      end subroutine
!*********************************************************************************
      subroutine corth(nm,n,low,igh,ar,ai,ortr,orti)
!
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      real(kind=rk)  :: ar(nm,n),ai(nm,n),ortr(igh),orti(igh)
      real(kind=rk)  :: f,g,h,fi,fr,scale
!
!     this subroutine is a translation of a complex analogue of
!     the algol procedure orthes, num. math. 12, 349-368(1968)
!     by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     given a complex general matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     low through igh to upper hessenberg form by
!     unitary similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  cbal.  if  cbal  has not been used,
!          set low=1, igh=n.
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the complex input matrix.
!
!     on output
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the hessenberg matrix.  information
!          about the unitary transformations used in the reduction
!          is stored in the remaining triangles under the
!          hessenberg matrix.
!
!        ortr and orti contain further information about the
!          transformations.  only elements low through igh are used.
!
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
!
      do 180 m = kp1, la
         h = 0.0e0_rk
         ortr(m) = 0.0e0_rk
         orti(m) = 0.0e0_rk
         scale = 0.0e0_rk
!     .......... scale column (algol tol then not needed) ..........
         do 90 i = m, igh
   90    scale = scale + abs(ar(i,m-1)) + abs(ai(i,m-1))
!
         if (scale .eq. 0.0e0_rk) go to 180
         mp = m + igh
!     .......... for i=igh step -1 until m do -- ..........
         do 100 ii = m, igh
            i = mp - ii
            ortr(i) = ar(i,m-1) / scale
            orti(i) = ai(i,m-1) / scale
            h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
  100    continue
!
         g = sqrt(h)
         f = pythag(ortr(m),orti(m))
         if (f .eq. 0.0e0_rk) go to 103
         h = h + f * g
         g = g / f
         ortr(m) = (1.0e0_rk + g) * ortr(m)
         orti(m) = (1.0e0_rk + g) * orti(m)
         go to 105
!
  103    ortr(m) = g
         ar(m,m-1) = scale
!     .......... form (i-(u*ut)/h) * a ..........
  105    do 130 j = m, n
            fr = 0.0e0_rk
            fi = 0.0e0_rk
!     .......... for i=igh step -1 until m do -- ..........
            do 110 ii = m, igh
               i = mp - ii
               fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
               fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
  110       continue
!
            fr = fr / h
            fi = fi / h
!
            do 120 i = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)
               ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)
  120       continue
!
  130    continue
!     .......... form (i-(u*ut)/h)*a*(i-(u*ut)/h) ..........
         do 160 i = 1, igh
            fr = 0.0e0_rk
            fi = 0.0e0_rk
!     .......... for j=igh step -1 until m do -- ..........
            do 140 jj = m, igh
               j = mp - jj
               fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
               fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
  140       continue
!
            fr = fr / h
            fi = fi / h
!
            do 150 j = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)
               ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)
  150       continue
!
  160    continue
!
         ortr(m) = scale * ortr(m)
         orti(m) = scale * orti(m)
         ar(m,m-1) = -g * ar(m,m-1)
         ai(m,m-1) = -g * ai(m,m-1)
  180 continue
!
  200 return
      end subroutine
!*********************************************************************************
      subroutine elmbak(nm,low,igh,a,int,m,z)
!
      integer i,j,m,la,mm,mp,nm,igh,kp1,low,mp1
      real(kind=rk)  :: a(nm,igh),z(nm,m)
      real(kind=rk)  :: x
      integer int(igh)
!
!     this subroutine is a translation of the algol procedure elmbak,
!     num. math. 12, 349-368(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     this subroutine forms the eigenvectors of a real general
!     matrix by back transforming those of the corresponding
!     upper hessenberg matrix determined by  elmhes.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1 and igh equal to the order of the matrix.
!
!        a contains the multipliers which were used in the
!          reduction by  elmhes  in its lower triangle
!          below the subdiagonal.
!
!        int contains information on the rows and columns
!          interchanged in the reduction by  elmhes.
!          only elements low through igh are used.
!
!        m is the number of columns of z to be back transformed.
!
!        z contains the real and imaginary parts of the eigen-
!          vectors to be back transformed in its first m columns.
!
!     on output
!
!        z contains the real and imaginary parts of the
!          transformed eigenvectors in its first m columns.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (m .eq. 0) go to 200
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
!     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = kp1, la
         mp = low + igh - mm
         mp1 = mp + 1
!
         do 110 i = mp1, igh
            x = a(i,mp-1)
            if (x .eq. 0.0e0_rk) go to 110
!
            do 100 j = 1, m
  100       z(i,j) = z(i,j) + x * z(mp,j)
!
  110    continue
!
         i = int(mp)
         if (i .eq. mp) go to 140
!
         do 130 j = 1, m
            x = z(i,j)
            z(i,j) = z(mp,j)
            z(mp,j) = x
  130    continue
!
  140 continue
!
  200 return
      end subroutine
!*********************************************************************************
      subroutine elmhes(nm,n,low,igh,a,int)
!
      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
      real(kind=rk)  :: a(nm,n)
      real(kind=rk)  :: x,y
      integer int(igh)
!
!     this subroutine is a translation of the algol procedure elmhes,
!     num. math. 12, 349-368(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     given a real general matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     low through igh to upper hessenberg form by
!     stabilized elementary similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        a contains the input matrix.
!
!     on output
!
!        a contains the hessenberg matrix.  the multipliers
!          which were used in the reduction are stored in the
!          remaining triangle under the hessenberg matrix.
!
!        int contains information on the rows and columns
!          interchanged in the reduction.
!          only elements low through igh are used.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
!
      do 180 m = kp1, la
         mm1 = m - 1
         x = 0.0e0_rk
         i = m
!
         do 100 j = m, igh
            if (abs(a(j,mm1)) .le. abs(x)) go to 100
            x = a(j,mm1)
            i = j
  100    continue
!
         int(m) = i
         if (i .eq. m) go to 130
!     .......... interchange rows and columns of a ..........
         do 110 j = mm1, n
            y = a(i,j)
            a(i,j) = a(m,j)
            a(m,j) = y
  110    continue
!
         do 120 j = 1, igh
            y = a(j,i)
            a(j,i) = a(j,m)
            a(j,m) = y
  120    continue
!     .......... end interchange ..........
  130    if (x .eq. 0.0e0_rk) go to 180
         mp1 = m + 1
!
         do 160 i = mp1, igh
            y = a(i,mm1)
            if (y .eq. 0.0e0_rk) go to 160
            y = y / x
            a(i,mm1) = y
!
            do 140 j = m, n
  140       a(i,j) = a(i,j) - y * a(m,j)
!
            do 150 j = 1, igh
  150       a(j,m) = a(j,m) + y * a(j,i)
!
  160    continue
!
  180 continue
!
  200 return
      end subroutine
!*********************************************************************************
      subroutine eltran(nm,n,low,igh,a,int,z)
!
      integer i,j,n,kl,mm,mp,nm,igh,low,mp1
      real(kind=rk)  :: a(nm,igh),z(nm,n)
      integer int(igh)
!
!     this subroutine is a translation of the algol procedure elmtrans,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine accumulates the stabilized elementary
!     similarity transformations used in the reduction of a
!     real general matrix to upper hessenberg form by  elmhes.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        a contains the multipliers which were used in the
!          reduction by  elmhes  in its lower triangle
!          below the subdiagonal.
!
!        int contains information on the rows and columns
!          interchanged in the reduction by  elmhes.
!          only elements low through igh are used.
!
!     on output
!
!        z contains the transformation matrix produced in the
!          reduction by  elmhes.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
!     .......... initialize z to identity matrix ..........
      do 80 j = 1, n
!
         do 60 i = 1, n
   60    z(i,j) = 0.0e0_rk
!
         z(j,j) = 1.0e0_rk
   80 continue
!
      kl = igh - low - 1
      if (kl .lt. 1) go to 200
!     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = 1, kl
         mp = igh - mm
         mp1 = mp + 1
!
         do 100 i = mp1, igh
  100    z(i,mp) = a(i,mp-1)
!
         i = int(mp)
         if (i .eq. mp) go to 140
!
         do 130 j = mp, igh
            z(mp,j) = z(i,j)
            z(i,j) = 0.0e0_rk
  130    continue
!
         z(i,mp) = 1.0e0_rk
  140 continue
!
  200 return
      end subroutine
!*********************************************************************************
      subroutine figi(nm,n,t,d,e,e2,ierr)
!
      integer i,n,nm,ierr
      real(kind=rk)  :: t(nm,3),d(n),e(n),e2(n)
!
!     given a nonsymmetric tridiagonal matrix such that the products
!     of corresponding pairs of off-diagonal elements are all
!     non-negative, this subroutine reduces it to a symmetric
!     tridiagonal matrix with the same eigenvalues.  if, further,
!     a zero product only occurs when both factors are zero,
!     the reduced matrix is similar to the original matrix.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        t contains the input matrix.  its subdiagonal is
!          stored in the last n-1 positions of the first column,
!          its diagonal in the n positions of the second column,
!          and its superdiagonal in the first n-1 positions of
!          the third column.  t(1,1) and t(n,3) are arbitrary.
!
!     on output
!
!        t is unaltered.
!
!        d contains the diagonal elements of the symmetric matrix.
!
!        e contains the subdiagonal elements of the symmetric
!          matrix in its last n-1 positions.  e(1) is not set.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
!        ierr is set to
!          zero       for normal return,
!          n+i        if t(i,1)*t(i-1,3) is negative,
!          -(3*n+i)   if t(i,1)*t(i-1,3) is zero with one factor
!                     non-zero.  in this case, the eigenvectors of
!                     the symmetric matrix are not simply related
!                     to those of  t  and should not be sought.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
!
      do 100 i = 1, n
         if (i .eq. 1) go to 90
         e2(i) = t(i,1) * t(i-1,3)
         if (e2(i)) 1000, 60, 80
   60    if (t(i,1) .eq. 0.0e0_rk .and. t(i-1,3) .eq. 0.0e0_rk) go to 80
!     .......... set error -- product of some pair of off-diagonal
!                elements is zero with one member non-zero ..........
         ierr = -(3 * n + i)
   80    e(i) = sqrt(e2(i))
   90    d(i) = t(i,2)
  100 continue
!
      go to 1001
!     .......... set error -- product of some pair of off-diagonal
!                elements is negative ..........
 1000 ierr = n + i
 1001 return
      end subroutine
!*********************************************************************************
      subroutine figi2(nm,n,t,d,e,z,ierr)
!
      integer i,j,n,nm,ierr
      real(kind=rk)  :: t(nm,3),d(n),e(n),z(nm,n)
      real(kind=rk)  :: h
!
!     given a nonsymmetric tridiagonal matrix such that the products
!     of corresponding pairs of off-diagonal elements are all
!     non-negative, and zero only when both factors are zero, this
!     subroutine reduces it to a symmetric tridiagonal matrix
!     using and accumulating diagonal similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        t contains the input matrix.  its subdiagonal is
!          stored in the last n-1 positions of the first column,
!          its diagonal in the n positions of the second column,
!          and its superdiagonal in the first n-1 positions of
!          the third column.  t(1,1) and t(n,3) are arbitrary.
!
!     on output
!
!        t is unaltered.
!
!        d contains the diagonal elements of the symmetric matrix.
!
!        e contains the subdiagonal elements of the symmetric
!          matrix in its last n-1 positions.  e(1) is not set.
!
!        z contains the transformation matrix produced in
!          the reduction.
!
!        ierr is set to
!          zero       for normal return,
!          n+i        if t(i,1)*t(i-1,3) is negative,
!          2*n+i      if t(i,1)*t(i-1,3) is zero with
!                     one factor non-zero.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
!
      do 100 i = 1, n
!
         do 50 j = 1, n
   50    z(i,j) = 0.0e0_rk
!
         if (i .eq. 1) go to 70
         h = t(i,1) * t(i-1,3)
         if (h) 900, 60, 80
   60    if (t(i,1) .ne. 0.0e0_rk .or. t(i-1,3) .ne. 0.0e0_rk) go to 1000
         e(i) = 0.0e0_rk
   70    z(i,i) = 1.0e0_rk
         go to 90
   80    e(i) = sqrt(h)
         z(i,i) = z(i-1,i-1) * e(i) / t(i-1,3)
   90    d(i) = t(i,2)
  100 continue
!
      go to 1001
!     .......... set error -- product of some pair of off-diagonal
!                elements is negative ..........
  900 ierr = n + i
      go to 1001
!     .......... set error -- product of some pair of off-diagonal
!                elements is zero with one member non-zero ..........
 1000 ierr = 2 * n + i
 1001 return
      end subroutine
!*********************************************************************************
      subroutine hqr(nm,n,low,igh,h,wr,wi,ierr)
!
      integer i,j,k,l,m,n,en,ll,mm,na,nm,igh,itn,its,low,mp2,enm2,ierr
      real(kind=rk)  :: h(nm,n),wr(n),wi(n)
      real(kind=rk)  :: p,q,r,s,t,w,x,y,zz,norm,tst1,tst2
      logical notlas
!
!     this subroutine is a translation of the algol procedure hqr,
!     num. math. 14, 219-231(1970) by martin, peters, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).
!
!     this subroutine finds the eigenvalues of a real
!     upper hessenberg matrix by the qr method.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        h contains the upper hessenberg matrix.  information about
!          the transformations used in the reduction to hessenberg
!          form by  elmhes  or  orthes, if performed, is stored
!          in the remaining triangle under the hessenberg matrix.
!
!     on output
!
!        h has been destroyed.  therefore, it must be saved
!          before calling  hqr  if subsequent calculation and
!          back transformation of eigenvectors is to be performed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  the eigenvalues
!          are unordered except that complex conjugate pairs
!          of values appear consecutively with the eigenvalue
!          having the positive imaginary part first.  if an
!          error exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      norm = 0.0e0_rk
      k = 1
!     .......... store roots isolated by balanc
!                and compute matrix norm ..........
      do 50 i = 1, n
!
         do 40 j = k, n
   40    norm = norm + abs(h(i,j))
!
         k = i
         if (i .ge. low .and. i .le. igh) go to 50
         wr(i) = h(i,i)
         wi(i) = 0.0e0_rk
   50 continue
!
      en = igh
      t = 0.0e0_rk
      itn = 30*n
!     .......... search for next eigenvalues ..........
   60 if (en .lt. low) go to 1001
      its = 0
      na = en - 1
      enm2 = na - 1
!     .......... look for single small sub-diagonal element
!                for l=en step -1 until low do -- ..........
   70 do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 100
         s = abs(h(l-1,l-1)) + abs(h(l,l))
         if (s .eq. 0.0e0_rk) s = norm
         tst1 = s
         tst2 = tst1 + abs(h(l,l-1))
         if (tst2 .eq. tst1) go to 100
   80 continue
!     .......... form shift ..........
  100 x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (itn .eq. 0) go to 1000
      if (its .ne. 10 .and. its .ne. 20) go to 130
!     .......... form exceptional shift ..........
      t = t + x
!
      do 120 i = low, en
  120 h(i,i) = h(i,i) - x
!
      s = abs(h(en,na)) + abs(h(na,enm2))
      x = 0.75e0_rk * s
      y = x
      w = -0.4375e0_rk * s * s
  130 its = its + 1
      itn = itn - 1
!     .......... look for two consecutive small
!                sub-diagonal elements.
!                for m=en-2 step -1 until l do -- ..........
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = abs(p) + abs(q) + abs(r)
         p = p / s
         q = q / s
         r = r / s
         if (m .eq. l) go to 150
         tst1 = abs(p)*(abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1)))
         tst2 = tst1 + abs(h(m,m-1))*(abs(q) + abs(r))
         if (tst2 .eq. tst1) go to 150
  140 continue
!
  150 mp2 = m + 2
!
      do 160 i = mp2, en
         h(i,i-2) = 0.0e0_rk
         if (i .eq. mp2) go to 160
         h(i,i-3) = 0.0e0_rk
  160 continue
!     .......... double qr step involving rows l to en and
!                columns m to en ..........
      do 260 k = m, na
         notlas = k .ne. na
         if (k .eq. m) go to 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0.0e0_rk
         if (notlas) r = h(k+2,k-1)
         x = abs(p) + abs(q) + abs(r)
         if (x .eq. 0.0e0_rk) go to 260
         p = p / x
         q = q / x
         r = r / x
  170    s = sign(sqrt(p*p+q*q+r*r),p)
         if (k .eq. m) go to 180
         h(k,k-1) = -s * x
         go to 190
  180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
         if (notlas) go to 225
!     .......... row modification ..........
         do 200 j = k, n
            p = h(k,j) + q * h(k+1,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
  200    continue
!
         j = min0(en,k+3)
!     .......... column modification ..........
         do 210 i = 1, j
            p = x * h(i,k) + y * h(i,k+1)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
  210    continue
         go to 255
  225    continue
!     .......... row modification ..........
         do 230 j = k, n
            p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
            h(k+2,j) = h(k+2,j) - p * zz
  230    continue
!
         j = min0(en,k+3)
!     .......... column modification ..........
         do 240 i = 1, j
            p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
            h(i,k+2) = h(i,k+2) - p * r
  240    continue
  255    continue
!
  260 continue
!
      go to 70
!     .......... one root found ..........
  270 wr(en) = x + t
      wi(en) = 0.0e0_rk
      en = na
      go to 60
!     .......... two roots found ..........
  280 p = (y - x) / 2.0e0_rk
      q = p * p + w
      zz = sqrt(abs(q))
      x = x + t
      if (q .lt. 0.0e0_rk) go to 320
!     .......... real pair ..........
      zz = p + sign(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0.0e0_rk) wr(en) = x - w / zz
      wi(na) = 0.0e0_rk
      wi(en) = 0.0e0_rk
      go to 330
!     .......... complex pair ..........
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      go to 60
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end subroutine
!*********************************************************************************
      subroutine hqr2(nm,n,low,igh,h,wr,wi,z,ierr)
!
      integer i,j,k,l,m,n,en,ii,jj,ll,mm,na,nm,nn, igh,itn,its,low,mp2,enm2,ierr
      real(kind=rk)  :: h(nm,n),wr(n),wi(n),z(nm,n)
      real(kind=rk)  :: p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,tst1,tst2
      logical notlas
!
!     this subroutine is a translation of the algol procedure hqr2,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a real upper hessenberg matrix by the qr method.  the
!     eigenvectors of a real general matrix can also be found
!     if  elmhes  and  eltran  or  orthes  and  ortran  have
!     been used to reduce this general matrix to hessenberg form
!     and to accumulate the similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        h contains the upper hessenberg matrix.
!
!        z contains the transformation matrix produced by  eltran
!          after the reduction by  elmhes, or by  ortran  after the
!          reduction by  orthes, if performed.  if the eigenvectors
!          of the hessenberg matrix are desired, z must contain the
!          identity matrix.
!
!     on output
!
!        h has been destroyed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  the eigenvalues
!          are unordered except that complex conjugate pairs
!          of values appear consecutively with the eigenvalue
!          having the positive imaginary part first.  if an
!          error exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        z contains the real and imaginary parts of the eigenvectors.
!          if the i-th eigenvalue is real, the i-th column of z
!          contains its eigenvector.  if the i-th eigenvalue is complex
!          with positive imaginary part, the i-th and (i+1)-th
!          columns of z contain the real and imaginary parts of its
!          eigenvector.  the eigenvectors are unnormalized.  if an
!          error exit is made, none of the eigenvectors has been found.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
!     calls cdiv for complex division.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      norm = 0.0e0_rk
      k = 1
!     .......... store roots isolated by balanc
!                and compute matrix norm ..........
      do 50 i = 1, n
!
         do 40 j = k, n
   40    norm = norm + abs(h(i,j))
!
         k = i
         if (i .ge. low .and. i .le. igh) go to 50
         wr(i) = h(i,i)
         wi(i) = 0.0e0_rk
   50 continue
!
      en = igh
      t = 0.0e0_rk
      itn = 30*n
!     .......... search for next eigenvalues ..........
   60 if (en .lt. low) go to 340
      its = 0
      na = en - 1
      enm2 = na - 1
!     .......... look for single small sub-diagonal element
!                for l=en step -1 until low do -- ..........
   70 do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 100
         s = abs(h(l-1,l-1)) + abs(h(l,l))
         if (s .eq. 0.0e0_rk) s = norm
         tst1 = s
         tst2 = tst1 + abs(h(l,l-1))
         if (tst2 .eq. tst1) go to 100
   80 continue
!     .......... form shift ..........
  100 x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (itn .eq. 0) go to 1000
      if (its .ne. 10 .and. its .ne. 20) go to 130
!     .......... form exceptional shift ..........
      t = t + x
!
      do 120 i = low, en
  120 h(i,i) = h(i,i) - x
!
      s = abs(h(en,na)) + abs(h(na,enm2))
      x = 0.75e0_rk * s
      y = x
      w = -0.4375e0_rk * s * s
  130 its = its + 1
      itn = itn - 1
!     .......... look for two consecutive small
!                sub-diagonal elements.
!                for m=en-2 step -1 until l do -- ..........
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = abs(p) + abs(q) + abs(r)
         p = p / s
         q = q / s
         r = r / s
         if (m .eq. l) go to 150
         tst1 = abs(p)*(abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1)))
         tst2 = tst1 + abs(h(m,m-1))*(abs(q) + abs(r))
         if (tst2 .eq. tst1) go to 150
  140 continue
!
  150 mp2 = m + 2
!
      do 160 i = mp2, en
         h(i,i-2) = 0.0e0_rk
         if (i .eq. mp2) go to 160
         h(i,i-3) = 0.0e0_rk
  160 continue
!     .......... double qr step involving rows l to en and
!                columns m to en ..........
      do 260 k = m, na
         notlas = k .ne. na
         if (k .eq. m) go to 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0.0e0_rk
         if (notlas) r = h(k+2,k-1)
         x = abs(p) + abs(q) + abs(r)
         if (x .eq. 0.0e0_rk) go to 260
         p = p / x
         q = q / x
         r = r / x
  170    s = sign(sqrt(p*p+q*q+r*r),p)
         if (k .eq. m) go to 180
         h(k,k-1) = -s * x
         go to 190
  180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
         if (notlas) go to 225
!     .......... row modification ..........
         do 200 j = k, n
            p = h(k,j) + q * h(k+1,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
  200    continue
!
         j = min0(en,k+3)
!     .......... column modification ..........
         do 210 i = 1, j
            p = x * h(i,k) + y * h(i,k+1)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
  210    continue
!     .......... accumulate transformations ..........
         do 220 i = low, igh
            p = x * z(i,k) + y * z(i,k+1)
            z(i,k) = z(i,k) - p
            z(i,k+1) = z(i,k+1) - p * q
  220    continue
         go to 255
  225    continue
!     .......... row modification ..........
         do 230 j = k, n
            p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
            h(k+2,j) = h(k+2,j) - p * zz
  230    continue
!
         j = min0(en,k+3)
!     .......... column modification ..........
         do 240 i = 1, j
            p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
            h(i,k+2) = h(i,k+2) - p * r
  240    continue
!     .......... accumulate transformations ..........
         do 250 i = low, igh
            p = x * z(i,k) + y * z(i,k+1) + zz * z(i,k+2)
            z(i,k) = z(i,k) - p
            z(i,k+1) = z(i,k+1) - p * q
            z(i,k+2) = z(i,k+2) - p * r
  250    continue
  255    continue
!
  260 continue
!
      go to 70
!     .......... one root found ..........
  270 h(en,en) = x + t
      wr(en) = h(en,en)
      wi(en) = 0.0e0_rk
      en = na
      go to 60
!     .......... two roots found ..........
  280 p = (y - x) / 2.0e0_rk
      q = p * p + w
      zz = sqrt(abs(q))
      h(en,en) = x + t
      x = h(en,en)
      h(na,na) = y + t
      if (q .lt. 0.0e0_rk) go to 320
!     .......... real pair ..........
      zz = p + sign(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0.0e0_rk) wr(en) = x - w / zz
      wi(na) = 0.0e0_rk
      wi(en) = 0.0e0_rk
      x = h(en,na)
      s = abs(x) + abs(zz)
      p = x / s
      q = zz / s
      r = sqrt(p*p+q*q)
      p = p / r
      q = q / r
!     .......... row modification ..........
      do 290 j = na, n
         zz = h(na,j)
         h(na,j) = q * zz + p * h(en,j)
         h(en,j) = q * h(en,j) - p * zz
  290 continue
!     .......... column modification ..........
      do 300 i = 1, en
         zz = h(i,na)
         h(i,na) = q * zz + p * h(i,en)
         h(i,en) = q * h(i,en) - p * zz
  300 continue
!     .......... accumulate transformations ..........
      do 310 i = low, igh
         zz = z(i,na)
         z(i,na) = q * zz + p * z(i,en)
         z(i,en) = q * z(i,en) - p * zz
  310 continue
!
      go to 330
!     .......... complex pair ..........
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      go to 60
!     .......... all roots found.  backsubstitute to find
!                vectors of upper triangular form ..........
  340 if (norm .eq. 0.0e0_rk) go to 1001
!     .......... for en=n step -1 until 1 do -- ..........
      do 800 nn = 1, n
         en = n + 1 - nn
         p = wr(en)
         q = wi(en)
         na = en - 1
         if (q) 710, 600, 800
!     .......... real vector ..........
  600    m = en
         h(en,en) = 1.0e0_rk
         if (na .eq. 0) go to 800
!     .......... for i=en-1 step -1 until 1 do -- ..........
         do 700 ii = 1, na
            i = en - ii
            w = h(i,i) - p
            r = 0.0e0_rk
!
            do 610 j = m, en
  610       r = r + h(i,j) * h(j,en)
!
            if (wi(i) .ge. 0.0e0_rk) go to 630
            zz = w
            s = r
            go to 700
  630       m = i
            if (wi(i) .ne. 0.0e0_rk) go to 640
            t = w
            if (t .ne. 0.0e0_rk) go to 635
               tst1 = norm
               t = tst1
  632          t = 0.01e0_rk * t
               tst2 = norm + t
               if (tst2 .gt. tst1) go to 632
  635       h(i,en) = -r / t
            go to 680
!     .......... solve real equations ..........
  640       x = h(i,i+1)
            y = h(i+1,i)
            q = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i)
            t = (x * s - zz * r) / q
            h(i,en) = t
            if (abs(x) .le. abs(zz)) go to 650
            h(i+1,en) = (-r - w * t) / x
            go to 680
  650       h(i+1,en) = (-s - y * t) / zz
!
!     .......... overflow control ..........
  680       t = abs(h(i,en))
            if (t .eq. 0.0e0_rk) go to 700
            tst1 = t
            tst2 = tst1 + 1.0e0_rk/tst1
            if (tst2 .gt. tst1) go to 700
            do 690 j = i, en
               h(j,en) = h(j,en)/t
  690       continue
!
  700    continue
!     .......... end real vector ..........
         go to 800
!     .......... complex vector ..........
  710    m = na
!     .......... last vector component chosen imaginary so that
!                eigenvector matrix is triangular ..........
         if (abs(h(en,na)) .le. abs(h(na,en))) go to 720
         h(na,na) = q / h(en,na)
         h(na,en) = -(h(en,en) - p) / h(en,na)
         go to 730
  720    call cdiv(0.0e0_rk,-h(na,en),h(na,na)-p,q,h(na,na),h(na,en))
  730    h(en,na) = 0.0e0_rk
         h(en,en) = 1.0e0_rk
         enm2 = na - 1
         if (enm2 .eq. 0) go to 800
!     .......... for i=en-2 step -1 until 1 do -- ..........
         do 795 ii = 1, enm2
            i = na - ii
            w = h(i,i) - p
            ra = 0.0e0_rk
            sa = 0.0e0_rk
!
            do 760 j = m, en
               ra = ra + h(i,j) * h(j,na)
               sa = sa + h(i,j) * h(j,en)
  760       continue
!
            if (wi(i) .ge. 0.0e0_rk) go to 770
            zz = w
            r = ra
            s = sa
            go to 795
  770       m = i
            if (wi(i) .ne. 0.0e0_rk) go to 780
            call cdiv(-ra,-sa,w,q,h(i,na),h(i,en))
            go to 790
!     .......... solve complex equations ..........
  780       x = h(i,i+1)
            y = h(i+1,i)
            vr = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) - q * q
            vi = (wr(i) - p) * 2.0e0_rk * q
            if (vr .ne. 0.0e0_rk .or. vi .ne. 0.0e0_rk) go to 784
               tst1 = norm * (abs(w) + abs(q) + abs(x) + abs(y) + abs(zz))
               vr = tst1
  783          vr = 0.01e0_rk * vr
               tst2 = tst1 + vr
               if (tst2 .gt. tst1) go to 783
  784       call cdiv(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi, h(i,na),h(i,en))
            if (abs(x) .le. abs(zz) + abs(q)) go to 785
            h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
            h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
            go to 790
  785       call cdiv(-r-y*h(i,na),-s-y*h(i,en),zz,q, h(i+1,na),h(i+1,en))
!
!     .......... overflow control ..........
  790       t = max(abs(h(i,na)), abs(h(i,en)))
            if (t .eq. 0.0e0_rk) go to 795
            tst1 = t
            tst2 = tst1 + 1.0e0_rk/tst1
            if (tst2 .gt. tst1) go to 795
            do 792 j = i, en
               h(j,na) = h(j,na)/t
               h(j,en) = h(j,en)/t
  792       continue
!
  795    continue
!     .......... end complex vector ..........
  800 continue
!     .......... end back substitution.
!                vectors of isolated roots ..........
      do 840 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 840
!
         do 820 j = i, n
  820    z(i,j) = h(i,j)
!
  840 continue
!     .......... multiply by transformation matrix to give
!                vectors of original full matrix.
!                for j=n step -1 until low do -- ..........
      do 880 jj = low, n
         j = n + low - jj
         m = min0(j,igh)
!
         do 880 i = low, igh
            zz = 0.0e0_rk
!
            do 860 k = low, m
  860       zz = zz + z(i,k) * h(k,j)
!
            z(i,j) = zz
  880 continue
!
      go to 1001
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end subroutine
!*********************************************************************************
      subroutine htrib3(nm,n,a,tau,m,zr,zi)
!
      integer i,j,k,l,m,n,nm
      real(kind=rk)  :: a(nm,n),tau(2,n),zr(nm,m),zi(nm,m)
      real(kind=rk)  :: h,s,si
!
!     this subroutine is a translation of a complex analogue of
!     the algol procedure trbak3, num. math. 11, 181-195(1968)
!     by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine forms the eigenvectors of a complex hermitian
!     matrix by back transforming those of the corresponding
!     real symmetric tridiagonal matrix determined by  htrid3.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains information about the unitary transformations
!          used in the reduction by  htrid3.
!
!        tau contains further information about the transformations.
!
!        m is the number of eigenvectors to be back transformed.
!
!        zr contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the transformed eigenvectors
!          in their first m columns.
!
!     note that the last component of each returned vector
!     is real and that vector euclidean norms are preserved.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (m .eq. 0) go to 200
!     .......... transform the eigenvectors of the real symmetric
!                tridiagonal matrix to those of the hermitian
!                tridiagonal matrix. ..........
      do 50 k = 1, n
!
         do 50 j = 1, m
            zi(k,j) = -zr(k,j) * tau(2,k)
            zr(k,j) = zr(k,j) * tau(1,k)
   50 continue
!
      if (n .eq. 1) go to 200
!     .......... recover and apply the householder matrices ..........
      do 140 i = 2, n
         l = i - 1
         h = a(i,i)
         if (h .eq. 0.0e0_rk) go to 140
!
         do 130 j = 1, m
            s = 0.0e0_rk
            si = 0.0e0_rk
!
            do 110 k = 1, l
               s = s + a(i,k) * zr(k,j) - a(k,i) * zi(k,j)
               si = si + a(i,k) * zi(k,j) + a(k,i) * zr(k,j)
  110       continue
!     .......... double divisions avoid possible underflow ..........
            s = (s / h) / h
            si = (si / h) / h
!
            do 120 k = 1, l
               zr(k,j) = zr(k,j) - s * a(i,k) - si * a(k,i)
               zi(k,j) = zi(k,j) - si * a(i,k) + s * a(k,i)
  120       continue
!
  130    continue
!
  140 continue
!
  200 return
      end subroutine
!*********************************************************************************
      subroutine htribk(nm,n,ar,ai,tau,m,zr,zi)
!
      integer i,j,k,l,m,n,nm
      real(kind=rk)  :: ar(nm,n),ai(nm,n),tau(2,n),zr(nm,m),zi(nm,m)
      real(kind=rk)  :: h,s,si
!
!     this subroutine is a translation of a complex analogue of
!     the algol procedure trbak1, num. math. 11, 181-195(1968)
!     by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine forms the eigenvectors of a complex hermitian
!     matrix by back transforming those of the corresponding
!     real symmetric tridiagonal matrix determined by  htridi.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        ar and ai contain information about the unitary trans-
!          formations used in the reduction by  htridi  in their
!          full lower triangles except for the diagonal of ar.
!
!        tau contains further information about the transformations.
!
!        m is the number of eigenvectors to be back transformed.
!
!        zr contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the transformed eigenvectors
!          in their first m columns.
!
!     note that the last component of each returned vector
!     is real and that vector euclidean norms are preserved.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (m .eq. 0) go to 200
!     .......... transform the eigenvectors of the real symmetric
!                tridiagonal matrix to those of the hermitian
!                tridiagonal matrix. ..........
      do 50 k = 1, n
!
         do 50 j = 1, m
            zi(k,j) = -zr(k,j) * tau(2,k)
            zr(k,j) = zr(k,j) * tau(1,k)
   50 continue
!
      if (n .eq. 1) go to 200
!     .......... recover and apply the householder matrices ..........
      do 140 i = 2, n
         l = i - 1
         h = ai(i,i)
         if (h .eq. 0.0e0_rk) go to 140
!
         do 130 j = 1, m
            s = 0.0e0_rk
            si = 0.0e0_rk
!
            do 110 k = 1, l
               s = s + ar(i,k) * zr(k,j) - ai(i,k) * zi(k,j)
               si = si + ar(i,k) * zi(k,j) + ai(i,k) * zr(k,j)
  110       continue
!     .......... double divisions avoid possible underflow ..........
            s = (s / h) / h
            si = (si / h) / h
!
            do 120 k = 1, l
               zr(k,j) = zr(k,j) - s * ar(i,k) - si * ai(i,k)
               zi(k,j) = zi(k,j) - si * ar(i,k) + s * ai(i,k)
  120       continue
!
  130    continue
!
  140 continue
!
  200 return
      end subroutine
!*********************************************************************************
      subroutine htrid3(nm,n,a,d,e,e2,tau)
!
      integer i,j,k,l,n,ii,nm,jm1,jp1
      real(kind=rk)  :: a(nm,n),d(n),e(n),e2(n),tau(2,n)
      real(kind=rk)  :: f,g,h,fi,gi,hh,si,scale
!
!     this subroutine is a translation of a complex analogue of
!     the algol procedure tred3, num. math. 11, 181-195(1968)
!     by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a complex hermitian matrix, stored as
!     a single square array, to a real symmetric tridiagonal matrix
!     using unitary similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the lower triangle of the complex hermitian input
!          matrix.  the real parts of the matrix elements are stored
!          in the full lower triangle of a, and the imaginary parts
!          are stored in the transposed positions of the strict upper
!          triangle of a.  no storage is required for the zero
!          imaginary parts of the diagonal elements.
!
!     on output
!
!        a contains information about the unitary transformations
!          used in the reduction.
!
!        d contains the diagonal elements of the the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
!        tau contains further information about the transformations.
!
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      tau(1,n) = 1.0e0_rk
      tau(2,n) = 0.0e0_rk
!     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0e0_rk
         scale = 0.0e0_rk
         if (l .lt. 1) go to 130
!     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + abs(a(i,k)) + abs(a(k,i))
!
         if (scale .ne. 0.0e0_rk) go to 140
         tau(1,l) = 1.0e0_rk
         tau(2,l) = 0.0e0_rk
  130    e(i) = 0.0e0_rk
         e2(i) = 0.0e0_rk
         go to 290
!
  140    do 150 k = 1, l
            a(i,k) = a(i,k) / scale
            a(k,i) = a(k,i) / scale
            h = h + a(i,k) * a(i,k) + a(k,i) * a(k,i)
  150    continue
!
         e2(i) = scale * scale * h
         g = sqrt(h)
         e(i) = scale * g
         f = pythag(a(i,l),a(l,i))
!     .......... form next diagonal element of matrix t ..........
         if (f .eq. 0.0e0_rk) go to 160
         tau(1,l) = (a(l,i) * tau(2,i) - a(i,l) * tau(1,i)) / f
         si = (a(i,l) * tau(2,i) + a(l,i) * tau(1,i)) / f
         h = h + f * g
         g = 1.0e0_rk + g / f
         a(i,l) = g * a(i,l)
         a(l,i) = g * a(l,i)
         if (l .eq. 1) go to 270
         go to 170
  160    tau(1,l) = -tau(1,i)
         si = tau(2,i)
         a(i,l) = g
  170    f = 0.0e0_rk
!
         do 240 j = 1, l
            g = 0.0e0_rk
            gi = 0.0e0_rk
            if (j .eq. 1) go to 190
            jm1 = j - 1
!     .......... form element of a*u ..........
            do 180 k = 1, jm1
               g = g + a(j,k) * a(i,k) + a(k,j) * a(k,i)
               gi = gi - a(j,k) * a(k,i) + a(k,j) * a(i,k)
  180       continue
!
  190       g = g + a(j,j) * a(i,j)
            gi = gi - a(j,j) * a(j,i)
            jp1 = j + 1
            if (l .lt. jp1) go to 220
!
            do 200 k = jp1, l
               g = g + a(k,j) * a(i,k) - a(j,k) * a(k,i)
               gi = gi - a(k,j) * a(k,i) - a(j,k) * a(i,k)
  200       continue
!     .......... form element of p ..........
  220       e(j) = g / h
            tau(2,j) = gi / h
            f = f + e(j) * a(i,j) - tau(2,j) * a(j,i)
  240    continue
!
         hh = f / (h + h)
!     .......... form reduced a ..........
         do 260 j = 1, l
            f = a(i,j)
            g = e(j) - hh * f
            e(j) = g
            fi = -a(j,i)
            gi = tau(2,j) - hh * fi
            tau(2,j) = -gi
            a(j,j) = a(j,j) - 2.0e0_rk * (f * g + fi * gi)
            if (j .eq. 1) go to 260
            jm1 = j - 1
!
            do 250 k = 1, jm1
               a(j,k) = a(j,k) - f * e(k) - g * a(i,k) + fi * tau(2,k) + gi * a(k,i)
               a(k,j) = a(k,j) - f * tau(2,k) - g * a(k,i) - fi * e(k) - gi * a(i,k)
  250       continue
!
  260    continue
!
  270    do 280 k = 1, l
            a(i,k) = scale * a(i,k)
            a(k,i) = scale * a(k,i)
  280    continue
!
         tau(2,l) = -si
  290    d(i) = a(i,i)
         a(i,i) = scale * sqrt(h)
  300 continue
!
      return
      end subroutine
!*********************************************************************************
      subroutine htridi(nm,n,ar,ai,d,e,e2,tau)
!
      integer i,j,k,l,n,ii,nm,jp1
      real(kind=rk)  :: ar(nm,n),ai(nm,n),d(n),e(n),e2(n),tau(2,n)
      real(kind=rk)  :: f,g,h,fi,gi,hh,si,scale
!
!     this subroutine is a translation of a complex analogue of
!     the algol procedure tred1, num. math. 11, 181-195(1968)
!     by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a complex hermitian matrix
!     to a real symmetric tridiagonal matrix using
!     unitary similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the complex hermitian input matrix.
!          only the lower triangle of the matrix need be supplied.
!
!     on output
!
!        ar and ai contain information about the unitary trans-
!          formations used in the reduction in their full lower
!          triangles.  their strict upper triangles and the
!          diagonal of ar are unaltered.
!
!        d contains the diagonal elements of the the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
!        tau contains further information about the transformations.
!
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      tau(1,n) = 1.0e0_rk
      tau(2,n) = 0.0e0_rk
!
      do 100 i = 1, n
  100 d(i) = ar(i,i)
!     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0e0_rk
         scale = 0.0e0_rk
         if (l .lt. 1) go to 130
!     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + abs(ar(i,k)) + abs(ai(i,k))
!
         if (scale .ne. 0.0e0_rk) go to 140
         tau(1,l) = 1.0e0_rk
         tau(2,l) = 0.0e0_rk
  130    e(i) = 0.0e0_rk
         e2(i) = 0.0e0_rk
         go to 290
!
  140    do 150 k = 1, l
            ar(i,k) = ar(i,k) / scale
            ai(i,k) = ai(i,k) / scale
            h = h + ar(i,k) * ar(i,k) + ai(i,k) * ai(i,k)
  150    continue
!
         e2(i) = scale * scale * h
         g = sqrt(h)
         e(i) = scale * g
         f = pythag(ar(i,l),ai(i,l))
!     .......... form next diagonal element of matrix t ..........
         if (f .eq. 0.0e0_rk) go to 160
         tau(1,l) = (ai(i,l) * tau(2,i) - ar(i,l) * tau(1,i)) / f
         si = (ar(i,l) * tau(2,i) + ai(i,l) * tau(1,i)) / f
         h = h + f * g
         g = 1.0e0_rk + g / f
         ar(i,l) = g * ar(i,l)
         ai(i,l) = g * ai(i,l)
         if (l .eq. 1) go to 270
         go to 170
  160    tau(1,l) = -tau(1,i)
         si = tau(2,i)
         ar(i,l) = g
  170    f = 0.0e0_rk
!
         do 240 j = 1, l
            g = 0.0e0_rk
            gi = 0.0e0_rk
!     .......... form element of a*u ..........
            do 180 k = 1, j
               g = g + ar(j,k) * ar(i,k) + ai(j,k) * ai(i,k)
               gi = gi - ar(j,k) * ai(i,k) + ai(j,k) * ar(i,k)
  180       continue
!
            jp1 = j + 1
            if (l .lt. jp1) go to 220
!
            do 200 k = jp1, l
               g = g + ar(k,j) * ar(i,k) - ai(k,j) * ai(i,k)
               gi = gi - ar(k,j) * ai(i,k) - ai(k,j) * ar(i,k)
  200       continue
!     .......... form element of p ..........
  220       e(j) = g / h
            tau(2,j) = gi / h
            f = f + e(j) * ar(i,j) - tau(2,j) * ai(i,j)
  240    continue
!
         hh = f / (h + h)
!     .......... form reduced a ..........
         do 260 j = 1, l
            f = ar(i,j)
            g = e(j) - hh * f
            e(j) = g
            fi = -ai(i,j)
            gi = tau(2,j) - hh * fi
            tau(2,j) = -gi
!
            do 260 k = 1, j
               ar(j,k) = ar(j,k) - f * e(k) - g * ar(i,k) + fi * tau(2,k) + gi * ai(i,k)
               ai(j,k) = ai(j,k) - f * tau(2,k) - g * ai(i,k) - fi * e(k) - gi * ar(i,k)
  260    continue
!
  270    do 280 k = 1, l
            ar(i,k) = scale * ar(i,k)
            ai(i,k) = scale * ai(i,k)
  280    continue
!
         tau(2,l) = -si
  290    hh = d(i)
         d(i) = ar(i,i)
         ar(i,i) = hh
         ai(i,i) = scale * sqrt(h)
  300 continue
!
      return
      end subroutine
!*********************************************************************************
      subroutine imtql1(n,d,e,ierr)
!
      integer i,j,l,m,n,ii,mml,ierr
      real(kind=rk)  :: d(n),e(n)
      real(kind=rk)  :: b,c,f,g,p,r,s,tst1,tst2
!
!     this subroutine is a translation of the algol procedure imtql1,
!     num. math. 12, 377-383(1968) by martin and wilkinson,
!     as modified in num. math. 15, 450(1970) by dubrulle.
!     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
!
!     this subroutine finds the eigenvalues of a symmetric
!     tridiagonal matrix by the implicit ql method.
!
!     on input
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct and
!          ordered for indices 1,2,...ierr-1, but may not be
!          the smallest eigenvalues.
!
!        e has been destroyed.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      if (n .eq. 1) go to 1001
!
      do 100 i = 2, n
  100 e(i-1) = e(i)
!
      e(n) = 0.0e0_rk
!
      do 290 l = 1, n
         j = 0
!     .......... look for small sub-diagonal element ..........
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            tst1 = abs(d(m)) + abs(d(m+1))
            tst2 = tst1 + abs(e(m))
            if (tst2 .eq. tst1) go to 120
  110    continue
!
  120    p = d(l)
         if (m .eq. l) go to 215
         if (j .eq. 30) go to 1000
         j = j + 1
!     .......... form shift ..........
         g = (d(l+1) - p) / (2.0e0_rk * e(l))
         r = pythag(g,1.0e0_rk)
         g = d(m) - p + e(l) / (g + sign(r,g))
         s = 1.0e0_rk
         c = 1.0e0_rk
         p = 0.0e0_rk
         mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            r = pythag(f,g)
            e(i+1) = r
            if (r .eq. 0.0e0_rk) go to 210
            s = f / r
            c = g / r
            g = d(i+1) - p
            r = (d(i) - g) * s + 2.0e0_rk * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
  200    continue
!
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0e0_rk
         go to 105
!     .......... recover from underflow ..........
  210    d(i+1) = d(i+1) - p
         e(m) = 0.0e0_rk
         go to 105
!     .......... order eigenvalues ..........
  215    if (l .eq. 1) go to 250
!     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
!
  250    i = 1
  270    d(i) = p
  290 continue
!
      go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end subroutine
!*********************************************************************************
      subroutine imtql2(nm,n,d,e,z,ierr)
!
      integer i,j,k,l,m,n,ii,nm,mml,ierr
      real(kind=rk)  :: d(n),e(n),z(nm,n)
      real(kind=rk)  :: b,c,f,g,p,r,s,tst1,tst2
!
!     this subroutine is a translation of the algol procedure imtql2,
!     num. math. 12, 377-383(1968) by martin and wilkinson,
!     as modified in num. math. 15, 450(1970) by dubrulle.
!     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a symmetric tridiagonal matrix by the implicit ql method.
!     the eigenvectors of a full symmetric matrix can also
!     be found if  tred2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        z contains the transformation matrix produced in the
!          reduction by  tred2, if performed.  if the eigenvectors
!          of the tridiagonal matrix are desired, z must contain
!          the identity matrix.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,ierr-1.
!
!        e has been destroyed.
!
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      if (n .eq. 1) go to 1001
!
      do 100 i = 2, n
  100 e(i-1) = e(i)
!
      e(n) = 0.0e0_rk
!
      do 240 l = 1, n
         j = 0
!     .......... look for small sub-diagonal element ..........
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            tst1 = abs(d(m)) + abs(d(m+1))
            tst2 = tst1 + abs(e(m))
            if (tst2 .eq. tst1) go to 120
  110    continue
!
  120    p = d(l)
         if (m .eq. l) go to 240
         if (j .eq. 30) go to 1000
         j = j + 1
!     .......... form shift ..........
         g = (d(l+1) - p) / (2.0e0_rk * e(l))
         r = pythag(g,1.0e0_rk)
         g = d(m) - p + e(l) / (g + sign(r,g))
         s = 1.0e0_rk
         c = 1.0e0_rk
         p = 0.0e0_rk
         mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            r = pythag(f,g)
            e(i+1) = r
            if (r .eq. 0.0e0_rk) go to 210
            s = f / r
            c = g / r
            g = d(i+1) - p
            r = (d(i) - g) * s + 2.0e0_rk * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
!     .......... form vector ..........
            do 180 k = 1, n
               f = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * f
               z(k,i) = c * z(k,i) - s * f
  180       continue
!
  200    continue
!
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0e0_rk
         go to 105
!     .......... recover from underflow ..........
  210    d(i+1) = d(i+1) - p
         e(m) = 0.0e0_rk
         go to 105
  240 continue
!     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
!
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
!
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
!
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
!
  300 continue
!
      go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end subroutine
!*********************************************************************************
      subroutine imtqlv(n,d,e,e2,w,ind,ierr,rv1)
!
      integer i,j,k,l,m,n,ii,mml,tag,ierr
      real(kind=rk)  :: d(n),e(n),e2(n),w(n),rv1(n)
      real(kind=rk)  :: b,c,f,g,p,r,s,tst1,tst2
      integer ind(n)
!
!     this subroutine is a variant of  imtql1  which is a translation of
!     algol procedure imtql1, num. math. 12, 377-383(1968) by martin and
!     wilkinson, as modified in num. math. 15, 450(1970) by dubrulle.
!     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
!
!     this subroutine finds the eigenvalues of a symmetric tridiagonal
!     matrix by the implicit ql method and associates with them
!     their corresponding submatrix indices.
!
!     on input
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2(1) is arbitrary.
!
!     on output
!
!        d and e are unaltered.
!
!        elements of e2, corresponding to elements of e regarded
!          as negligible, have been replaced by zero causing the
!          matrix to split into a direct sum of submatrices.
!          e2(1) is also set to zero.
!
!        w contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct and
!          ordered for indices 1,2,...ierr-1, but may not be
!          the smallest eigenvalues.
!
!        ind contains the submatrix indices associated with the
!          corresponding eigenvalues in w -- 1 for eigenvalues
!          belonging to the first submatrix from the top,
!          2 for those belonging to the second submatrix, etc..
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!        rv1 is a temporary storage array.
!
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      k = 0
      tag = 0
!
      do 100 i = 1, n
         w(i) = d(i)
         if (i .ne. 1) rv1(i-1) = e(i)
  100 continue
!
      e2(1) = 0.0e0_rk
      rv1(n) = 0.0e0_rk
!
      do 290 l = 1, n
         j = 0
!     .......... look for small sub-diagonal element ..........
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            tst1 = abs(w(m)) + abs(w(m+1))
            tst2 = tst1 + abs(rv1(m))
            if (tst2 .eq. tst1) go to 120
!     .......... guard against underflowed element of e2 ..........
            if (e2(m+1) .eq. 0.0e0_rk) go to 125
  110    continue
!
  120    if (m .le. k) go to 130
         if (m .ne. n) e2(m+1) = 0.0e0_rk
  125    k = m
         tag = tag + 1
  130    p = w(l)
         if (m .eq. l) go to 215
         if (j .eq. 30) go to 1000
         j = j + 1
!     .......... form shift ..........
         g = (w(l+1) - p) / (2.0e0_rk * rv1(l))
         r = pythag(g,1.0e0_rk)
         g = w(m) - p + rv1(l) / (g + sign(r,g))
         s = 1.0e0_rk
         c = 1.0e0_rk
         p = 0.0e0_rk
         mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            f = s * rv1(i)
            b = c * rv1(i)
            r = pythag(f,g)
            rv1(i+1) = r
            if (r .eq. 0.0e0_rk) go to 210
            s = f / r
            c = g / r
            g = w(i+1) - p
            r = (w(i) - g) * s + 2.0e0_rk * c * b
            p = s * r
            w(i+1) = g + p
            g = c * r - b
  200    continue
!
         w(l) = w(l) - p
         rv1(l) = g
         rv1(m) = 0.0e0_rk
         go to 105
!     .......... recover from underflow ..........
  210    w(i+1) = w(i+1) - p
         rv1(m) = 0.0e0_rk
         go to 105
!     .......... order eigenvalues ..........
  215    if (l .eq. 1) go to 250
!     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. w(i-1)) go to 270
            w(i) = w(i-1)
            ind(i) = ind(i-1)
  230    continue
!
  250    i = 1
  270    w(i) = p
         ind(i) = tag
  290 continue
!
      go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end subroutine
!*********************************************************************************
      subroutine invit(nm,n,a,wr,wi,select,mm,m,z,ierr,rm1,rv1,rv2)
!
      integer i,j,k,l,m,n,s,ii,ip,mm,mp,nm,ns,n1,uk,ip1,its,km1,ierr
      real(kind=rk)  :: a(nm,n),wr(n),wi(n),z(nm,mm),rm1(n,n), rv1(n),rv2(n)
      real(kind=rk)  :: t,w,x,y,eps3,norm,normv,growto,ilambd,rlambd,ukroot
      logical select(n)
!
!     this subroutine is a translation of the algol procedure invit
!     by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
!
!     this subroutine finds those eigenvectors of a real upper
!     hessenberg matrix corresponding to specified eigenvalues,
!     using inverse iteration.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the hessenberg matrix.
!
!        wr and wi contain the real and imaginary parts, respectively,
!          of the eigenvalues of the matrix.  the eigenvalues must be
!          stored in a manner identical to that of subroutine  hqr,
!          which recognizes possible splitting of the matrix.
!
!        select specifies the eigenvectors to be found. the
!          eigenvector corresponding to the j-th eigenvalue is
!          specified by setting select(j) to .true..
!
!        mm should be set to an upper bound for the number of
!          columns required to store the eigenvectors to be found.
!          note that two columns are required to store the
!          eigenvector corresponding to a complex eigenvalue.
!
!     on output
!
!        a and wi are unaltered.
!
!        wr may have been altered since close eigenvalues are perturbed
!          slightly in searching for independent eigenvectors.
!
!        select may have been altered.  if the elements corresponding
!          to a pair of conjugate complex eigenvalues were each
!          initially set to .true., the program resets the second of
!          the two elements to .false..
!
!        m is the number of columns actually used to store
!          the eigenvectors.
!
!        z contains the real and imaginary parts of the eigenvectors.
!          if the next selected eigenvalue is real, the next column
!          of z contains its eigenvector.  if the eigenvalue is
!          complex, the next two columns of z contain the real and
!          imaginary parts of its eigenvector.  the eigenvectors are
!          normalized so that the component of largest magnitude is 1.
!          any vector which fails the acceptance test is set to zero.
!
!        ierr is set to
!          zero       for normal return,
!          -(2*n+1)   if more than mm columns of z are necessary
!                     to store the eigenvectors corresponding to
!                     the specified eigenvalues.
!          -k         if the iteration corresponding to the k-th
!                     value fails,
!          -(n+k)     if both error situations occur.
!
!        rm1, rv1, and rv2 are temporary storage arrays.  note that rm1
!          is square of dimension n by n and, augmented by two columns
!          of z, is the transpose of the corresponding algol b array.
!
!     the algol procedure guessvec appears in invit in line.
!
!     calls cdiv for complex division.
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      uk = 0
      s = 1
!     .......... ip = 0, real eigenvalue
!                     1, first of conjugate complex pair
!                    -1, second of conjugate complex pair ..........
      ip = 0
      n1 = n - 1
!
      do 980 k = 1, n
         if (wi(k) .eq. 0.0e0_rk .or. ip .lt. 0) go to 100
         ip = 1
         if (select(k) .and. select(k+1)) select(k+1) = .false.
  100    if (.not. select(k)) go to 960
         if (wi(k) .ne. 0.0e0_rk) s = s + 1
         if (s .gt. mm) go to 1000
         if (uk .ge. k) go to 200
!     .......... check for possible splitting ..........
         do 120 uk = k, n
            if (uk .eq. n) go to 140
            if (a(uk+1,uk) .eq. 0.0e0_rk) go to 140
  120    continue
!     .......... compute infinity norm of leading uk by uk
!                (hessenberg) matrix ..........
  140    norm = 0.0e0_rk
         mp = 1
!
         do 180 i = 1, uk
            x = 0.0e0_rk
!
            do 160 j = mp, uk
  160       x = x + abs(a(i,j))
!
            if (x .gt. norm) norm = x
            mp = i
  180    continue
!     .......... eps3 replaces zero pivot in decomposition
!                and close roots are modified by eps3 ..........
         if (norm .eq. 0.0e0_rk) norm = 1.0e0_rk
         eps3 = epslon(norm)
!     .......... growto is the criterion for the growth ..........
         ukroot = uk
         ukroot = sqrt(ukroot)
         growto = 0.1e0_rk / ukroot
  200    rlambd = wr(k)
         ilambd = wi(k)
         if (k .eq. 1) go to 280
         km1 = k - 1
         go to 240
!     .......... perturb eigenvalue if it is close
!                to any previous eigenvalue ..........
  220    rlambd = rlambd + eps3
!     .......... for i=k-1 step -1 until 1 do -- ..........
  240    do 260 ii = 1, km1
            i = k - ii
            if (select(i) .and. abs(wr(i)-rlambd) .lt. eps3 .and. abs(wi(i)-ilambd) .lt. eps3) go to 220
  260    continue
!
         wr(k) = rlambd
!     .......... perturb conjugate eigenvalue to match ..........
         ip1 = k + ip
         wr(ip1) = rlambd
!     .......... form upper hessenberg a-rlambd*i (transposed)
!                and initial real vector ..........
  280    mp = 1
!
         do 320 i = 1, uk
!
            do 300 j = mp, uk
  300       rm1(j,i) = a(i,j)
!
            rm1(i,i) = rm1(i,i) - rlambd
            mp = i
            rv1(i) = eps3
  320    continue
!
         its = 0
         if (ilambd .ne. 0.0e0_rk) go to 520
!     .......... real eigenvalue.
!                triangular decomposition with interchanges,
!                replacing zero pivots by eps3 ..........
         if (uk .eq. 1) go to 420
!
         do 400 i = 2, uk
            mp = i - 1
            if (abs(rm1(mp,i)) .le. abs(rm1(mp,mp))) go to 360
!
            do 340 j = mp, uk
               y = rm1(j,i)
               rm1(j,i) = rm1(j,mp)
               rm1(j,mp) = y
  340       continue
!
  360       if (rm1(mp,mp) .eq. 0.0e0_rk) rm1(mp,mp) = eps3
            x = rm1(mp,i) / rm1(mp,mp)
            if (x .eq. 0.0e0_rk) go to 400
!
            do 380 j = i, uk
  380       rm1(j,i) = rm1(j,i) - x * rm1(j,mp)
!
  400    continue
!
  420    if (rm1(uk,uk) .eq. 0.0e0_rk) rm1(uk,uk) = eps3
!     .......... back substitution for real vector
!                for i=uk step -1 until 1 do -- ..........
  440    do 500 ii = 1, uk
            i = uk + 1 - ii
            y = rv1(i)
            if (i .eq. uk) go to 480
            ip1 = i + 1
!
            do 460 j = ip1, uk
  460       y = y - rm1(j,i) * rv1(j)
!
  480       rv1(i) = y / rm1(i,i)
  500    continue
!
         go to 740
!     .......... complex eigenvalue.
!                triangular decomposition with interchanges,
!                replacing zero pivots by eps3.  store imaginary
!                parts in upper triangle starting at (1,3) ..........
  520    ns = n - s
         z(1,s-1) = -ilambd
         z(1,s) = 0.0e0_rk
         if (n .eq. 2) go to 550
         rm1(1,3) = -ilambd
         z(1,s-1) = 0.0e0_rk
         if (n .eq. 3) go to 550
!
         do 540 i = 4, n
  540    rm1(1,i) = 0.0e0_rk
!
  550    do 640 i = 2, uk
            mp = i - 1
            w = rm1(mp,i)
            if (i .lt. n) t = rm1(mp,i+1)
            if (i .eq. n) t = z(mp,s-1)
            x = rm1(mp,mp) * rm1(mp,mp) + t * t
            if (w * w .le. x) go to 580
            x = rm1(mp,mp) / w
            y = t / w
            rm1(mp,mp) = w
            if (i .lt. n) rm1(mp,i+1) = 0.0e0_rk
            if (i .eq. n) z(mp,s-1) = 0.0e0_rk
!
            do 560 j = i, uk
               w = rm1(j,i)
               rm1(j,i) = rm1(j,mp) - x * w
               rm1(j,mp) = w
               if (j .lt. n1) go to 555
               l = j - ns
               z(i,l) = z(mp,l) - y * w
               z(mp,l) = 0.0e0_rk
               go to 560
  555          rm1(i,j+2) = rm1(mp,j+2) - y * w
               rm1(mp,j+2) = 0.0e0_rk
  560       continue
!
            rm1(i,i) = rm1(i,i) - y * ilambd
            if (i .lt. n1) go to 570
            l = i - ns
            z(mp,l) = -ilambd
            z(i,l) = z(i,l) + x * ilambd
            go to 640
  570       rm1(mp,i+2) = -ilambd
            rm1(i,i+2) = rm1(i,i+2) + x * ilambd
            go to 640
  580       if (x .ne. 0.0e0_rk) go to 600
            rm1(mp,mp) = eps3
            if (i .lt. n) rm1(mp,i+1) = 0.0e0_rk
            if (i .eq. n) z(mp,s-1) = 0.0e0_rk
            t = 0.0e0_rk
            x = eps3 * eps3
  600       w = w / x
            x = rm1(mp,mp) * w
            y = -t * w
!
            do 620 j = i, uk
               if (j .lt. n1) go to 610
               l = j - ns
               t = z(mp,l)
               z(i,l) = -x * t - y * rm1(j,mp)
               go to 615
  610          t = rm1(mp,j+2)
               rm1(i,j+2) = -x * t - y * rm1(j,mp)
  615          rm1(j,i) = rm1(j,i) - x * rm1(j,mp) + y * t
  620       continue
!
            if (i .lt. n1) go to 630
            l = i - ns
            z(i,l) = z(i,l) - ilambd
            go to 640
  630       rm1(i,i+2) = rm1(i,i+2) - ilambd
  640    continue
!
         if (uk .lt. n1) go to 650
         l = uk - ns
         t = z(uk,l)
         go to 655
  650    t = rm1(uk,uk+2)
  655    if (rm1(uk,uk) .eq. 0.0e0_rk .and. t .eq. 0.0e0_rk) rm1(uk,uk) = eps3
!     .......... back substitution for complex vector
!                for i=uk step -1 until 1 do -- ..........
  660    do 720 ii = 1, uk
            i = uk + 1 - ii
            x = rv1(i)
            y = 0.0e0_rk
            if (i .eq. uk) go to 700
            ip1 = i + 1
!
            do 680 j = ip1, uk
               if (j .lt. n1) go to 670
               l = j - ns
               t = z(i,l)
               go to 675
  670          t = rm1(i,j+2)
  675          x = x - rm1(j,i) * rv1(j) + t * rv2(j)
               y = y - rm1(j,i) * rv2(j) - t * rv1(j)
  680       continue
!
  700       if (i .lt. n1) go to 710
            l = i - ns
            t = z(i,l)
            go to 715
  710       t = rm1(i,i+2)
  715       call cdiv(x,y,rm1(i,i),t,rv1(i),rv2(i))
  720    continue
!     .......... acceptance test for real or complex
!                eigenvector and normalization ..........
  740    its = its + 1
         norm = 0.0e0_rk
         normv = 0.0e0_rk
!
         do 780 i = 1, uk
            if (ilambd .eq. 0.0e0_rk) x = abs(rv1(i))
            if (ilambd .ne. 0.0e0_rk) x = pythag(rv1(i),rv2(i))
            if (normv .ge. x) go to 760
            normv = x
            j = i
  760       norm = norm + x
  780    continue
!
         if (norm .lt. growto) go to 840
!     .......... accept vector ..........
         x = rv1(j)
         if (ilambd .eq. 0.0e0_rk) x = 1.0e0_rk / x
         if (ilambd .ne. 0.0e0_rk) y = rv2(j)
!
         do 820 i = 1, uk
            if (ilambd .ne. 0.0e0_rk) go to 800
            z(i,s) = rv1(i) * x
            go to 820
  800       call cdiv(rv1(i),rv2(i),x,y,z(i,s-1),z(i,s))
  820    continue
!
         if (uk .eq. n) go to 940
         j = uk + 1
         go to 900
!     .......... in-line procedure for choosing
!                a new starting vector ..........
  840    if (its .ge. uk) go to 880
         x = ukroot
         y = eps3 / (x + 1.0e0_rk)
         rv1(1) = eps3
!
         do 860 i = 2, uk
  860    rv1(i) = y
!
         j = uk - its + 1
         rv1(j) = rv1(j) - eps3 * x
         if (ilambd .eq. 0.0e0_rk) go to 440
         go to 660
!     .......... set error -- unaccepted eigenvector ..........
  880    j = 1
         ierr = -k
!     .......... set remaining vector components to zero ..........
  900    do 920 i = j, n
            z(i,s) = 0.0e0_rk
            if (ilambd .ne. 0.0e0_rk) z(i,s-1) = 0.0e0_rk
  920    continue
!
  940    s = s + 1
  960    if (ip .eq. (-1)) ip = 0
         if (ip .eq. 1) ip = -1
  980 continue
!
      go to 1001
!     .......... set error -- underestimate of eigenvector
!                space required ..........
 1000 if (ierr .ne. 0) ierr = ierr - n
      if (ierr .eq. 0) ierr = -(2 * n + 1)
 1001 m = s - 1 - iabs(ip)
      return
      end subroutine
!*********************************************************************************
      subroutine minfit(nm,m,n,a,w,ip,b,ierr,rv1)
!
      integer i,j,k,l,m,n,ii,ip,i1,kk,k1,ll,l1,m1,nm,its,ierr
      real(kind=rk)  :: a(nm,n),w(n),b(nm,ip),rv1(n)
      real(kind=rk)  :: c,f,g,h,s,x,y,z,tst1,tst2,scale
!
!     this subroutine is a translation of the algol procedure minfit,
!     num. math. 14, 403-420(1970) by golub and reinsch.
!     handbook for auto. comp., vol ii-linear algebra, 134-151(1971).
!
!     this subroutine determines, towards the solution of the linear
!                                                        t
!     system ax=b, the singular value decomposition a=usv  of a real
!                                         t
!     m by n rectangular matrix, forming u b rather than u.  householder
!     bidiagonalization and a variant of the qr algorithm are used.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.  note that nm must be at least
!          as large as the maximum of m and n.
!
!        m is the number of rows of a and b.
!
!        n is the number of columns of a and the order of v.
!
!        a contains the rectangular coefficient matrix of the system.
!
!        ip is the number of columns of b.  ip can be zero.
!
!        b contains the constant column matrix of the system
!          if ip is not zero.  otherwise b is not referenced.
!
!     on output
!
!        a has been overwritten by the matrix v (orthogonal) of the
!          decomposition in its first n rows and columns.  if an
!          error exit is made, the columns of v corresponding to
!          indices of correct singular values should be correct.
!
!        w contains the n (non-negative) singular values of a (the
!          diagonal elements of s).  they are unordered.  if an
!          error exit is made, the singular values should be correct
!          for indices ierr+1,ierr+2,...,n.
!
!                                   t
!        b has been overwritten by u b.  if an error exit is made,
!                       t
!          the rows of u b corresponding to indices of correct
!          singular values should be correct.
!
!        ierr is set to
!          zero       for normal return,
!          k          if the k-th singular value has not been
!                     determined after 30 iterations.
!
!        rv1 is a temporary storage array.
!
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
!     .......... householder reduction to bidiagonal form ..........
      g = 0.0e0_rk
      scale = 0.0e0_rk
      x = 0.0e0_rk
!
      do 300 i = 1, n
         l = i + 1
         rv1(i) = scale * g
         g = 0.0e0_rk
         s = 0.0e0_rk
         scale = 0.0e0_rk
         if (i .gt. m) go to 210
!
         do 120 k = i, m
  120    scale = scale + abs(a(k,i))
!
         if (scale .eq. 0.0e0_rk) go to 210
!
         do 130 k = i, m
            a(k,i) = a(k,i) / scale
            s = s + a(k,i)**2
  130    continue
!
         f = a(i,i)
         g = -sign(sqrt(s),f)
         h = f * g - s
         a(i,i) = f - g
         if (i .eq. n) go to 160
!
         do 150 j = l, n
            s = 0.0e0_rk
!
            do 140 k = i, m
  140       s = s + a(k,i) * a(k,j)
!
            f = s / h
!
            do 150 k = i, m
               a(k,j) = a(k,j) + f * a(k,i)
  150    continue
!
  160    if (ip .eq. 0) go to 190
!
         do 180 j = 1, ip
            s = 0.0e0_rk
!
            do 170 k = i, m
  170       s = s + a(k,i) * b(k,j)
!
            f = s / h
!
            do 180 k = i, m
               b(k,j) = b(k,j) + f * a(k,i)
  180    continue
!
  190    do 200 k = i, m
  200    a(k,i) = scale * a(k,i)
!
  210    w(i) = scale * g
         g = 0.0e0_rk
         s = 0.0e0_rk
         scale = 0.0e0_rk
         if (i .gt. m .or. i .eq. n) go to 290
!
         do 220 k = l, n
  220    scale = scale + abs(a(i,k))
!
         if (scale .eq. 0.0e0_rk) go to 290
!
         do 230 k = l, n
            a(i,k) = a(i,k) / scale
            s = s + a(i,k)**2
  230    continue
!
         f = a(i,l)
         g = -sign(sqrt(s),f)
         h = f * g - s
         a(i,l) = f - g
!
         do 240 k = l, n
  240    rv1(k) = a(i,k) / h
!
         if (i .eq. m) go to 270
!
         do 260 j = l, m
            s = 0.0e0_rk
!
            do 250 k = l, n
  250       s = s + a(j,k) * a(i,k)
!
            do 260 k = l, n
               a(j,k) = a(j,k) + s * rv1(k)
  260    continue
!
  270    do 280 k = l, n
  280    a(i,k) = scale * a(i,k)
!
  290    x = max(x,abs(w(i))+abs(rv1(i)))
  300 continue
!     .......... accumulation of right-hand transformations.
!                for i=n step -1 until 1 do -- ..........
      do 400 ii = 1, n
         i = n + 1 - ii
         if (i .eq. n) go to 390
         if (g .eq. 0.0e0_rk) go to 360
!
         do 320 j = l, n
!     .......... double division avoids possible underflow ..........
  320    a(j,i) = (a(i,j) / a(i,l)) / g
!
         do 350 j = l, n
            s = 0.0e0_rk
!
            do 340 k = l, n
  340       s = s + a(i,k) * a(k,j)
!
            do 350 k = l, n
               a(k,j) = a(k,j) + s * a(k,i)
  350    continue
!
  360    do 380 j = l, n
            a(i,j) = 0.0e0_rk
            a(j,i) = 0.0e0_rk
  380    continue
!
  390    a(i,i) = 1.0e0_rk
         g = rv1(i)
         l = i
  400 continue
!
      if (m .ge. n .or. ip .eq. 0) go to 510
      m1 = m + 1
!
      do 500 i = m1, n
!
         do 500 j = 1, ip
            b(i,j) = 0.0e0_rk
  500 continue
!     .......... diagonalization of the bidiagonal form ..........
  510 tst1 = x
!     .......... for k=n step -1 until 1 do -- ..........
      do 700 kk = 1, n
         k1 = n - kk
         k = k1 + 1
         its = 0
!     .......... test for splitting.
!                for l=k step -1 until 1 do -- ..........
  520    do 530 ll = 1, k
            l1 = k - ll
            l = l1 + 1
            tst2 = tst1 + abs(rv1(l))
            if (tst2 .eq. tst1) go to 565
!     .......... rv1(1) is always zero, so there is no exit
!                through the bottom of the loop ..........
            tst2 = tst1 + abs(w(l1))
            if (tst2 .eq. tst1) go to 540
  530    continue
!     .......... cancellation of rv1(l) if l greater than 1 ..........
  540    c = 0.0e0_rk
         s = 1.0e0_rk
!
         do 560 i = l, k
            f = s * rv1(i)
            rv1(i) = c * rv1(i)
            tst2 = tst1 + abs(f)
            if (tst2 .eq. tst1) go to 565
            g = w(i)
            h = pythag(f,g)
            w(i) = h
            c = g / h
            s = -f / h
            if (ip .eq. 0) go to 560
!
            do 550 j = 1, ip
               y = b(l1,j)
               z = b(i,j)
               b(l1,j) = y * c + z * s
               b(i,j) = -y * s + z * c
  550       continue
!
  560    continue
!     .......... test for convergence ..........
  565    z = w(k)
         if (l .eq. k) go to 650
!     .......... shift from bottom 2 by 2 minor ..........
         if (its .eq. 30) go to 1000
         its = its + 1
         x = w(l)
         y = w(k1)
         g = rv1(k1)
         h = rv1(k)
         f = 0.5e0_rk * (((g + z) / h) * ((g - z) / y) + y / h - h / y)
         g = pythag(f,1.0e0_rk)
         f = x - (z / x) * z + (h / x) * (y / (f + sign(g,f)) - h)
!     .......... next qr transformation ..........
         c = 1.0e0_rk
         s = 1.0e0_rk
!
         do 600 i1 = l, k1
            i = i1 + 1
            g = rv1(i)
            y = w(i)
            h = s * g
            g = c * g
            z = pythag(f,h)
            rv1(i1) = z
            c = f / z
            s = h / z
            f = x * c + g * s
            g = -x * s + g * c
            h = y * s
            y = y * c
!
            do 570 j = 1, n
               x = a(j,i1)
               z = a(j,i)
               a(j,i1) = x * c + z * s
               a(j,i) = -x * s + z * c
  570       continue
!
            z = pythag(f,h)
            w(i1) = z
!     .......... rotation can be arbitrary if z is zero ..........
            if (z .eq. 0.0e0_rk) go to 580
            c = f / z
            s = h / z
  580       f = c * g + s * y
            x = -s * g + c * y
            if (ip .eq. 0) go to 600
!
            do 590 j = 1, ip
               y = b(i1,j)
               z = b(i,j)
               b(i1,j) = y * c + z * s
               b(i,j) = -y * s + z * c
  590       continue
!
  600    continue
!
         rv1(l) = 0.0e0_rk
         rv1(k) = f
         w(k) = x
         go to 520
!     .......... convergence ..........
  650    if (z .ge. 0.0e0_rk) go to 700
!     .......... w(k) is made non-negative ..........
         w(k) = -z
!
         do 690 j = 1, n
  690    a(j,k) = -a(j,k)
!
  700 continue
!
      go to 1001
!     .......... set error -- no convergence to a
!                singular value after 30 iterations ..........
 1000 ierr = k
 1001 return
      end subroutine
!*********************************************************************************
      subroutine ortbak(nm,low,igh,a,ort,m,z)
!
      integer i,j,m,la,mm,mp,nm,igh,kp1,low,mp1
      real(kind=rk)  :: a(nm,igh),ort(igh),z(nm,m)
      real(kind=rk)  :: g
!
!     this subroutine is a translation of the algol procedure ortbak,
!     num. math. 12, 349-368(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     this subroutine forms the eigenvectors of a real general
!     matrix by back transforming those of the corresponding
!     upper hessenberg matrix determined by  orthes.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1 and igh equal to the order of the matrix.
!
!        a contains information about the orthogonal trans-
!          formations used in the reduction by  orthes
!          in its strict lower triangle.
!
!        ort contains further information about the trans-
!          formations used in the reduction by  orthes.
!          only elements low through igh are used.
!
!        m is the number of columns of z to be back transformed.
!
!        z contains the real and imaginary parts of the eigen-
!          vectors to be back transformed in its first m columns.
!
!     on output
!
!        z contains the real and imaginary parts of the
!          transformed eigenvectors in its first m columns.
!
!        ort has been altered.
!
!     note that ortbak preserves vector euclidean norms.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (m .eq. 0) go to 200
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
!     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = kp1, la
         mp = low + igh - mm
         if (a(mp,mp-1) .eq. 0.0e0_rk) go to 140
         mp1 = mp + 1
!
         do 100 i = mp1, igh
  100    ort(i) = a(i,mp-1)
!
         do 130 j = 1, m
            g = 0.0e0_rk
!
            do 110 i = mp, igh
  110       g = g + ort(i) * z(i,j)
!     .......... divisor below is negative of h formed in orthes.
!                double division avoids possible underflow ..........
            g = (g / ort(mp)) / a(mp,mp-1)
!
            do 120 i = mp, igh
  120       z(i,j) = z(i,j) + g * ort(i)
!
  130    continue
!
  140 continue
!
  200 return
      end subroutine
!*********************************************************************************
      subroutine orthes(nm,n,low,igh,a,ort)
!
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      real(kind=rk)  :: a(nm,n),ort(igh)
      real(kind=rk)  :: f,g,h,scale
!
!     this subroutine is a translation of the algol procedure orthes,
!     num. math. 12, 349-368(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     given a real general matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     low through igh to upper hessenberg form by
!     orthogonal similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        a contains the input matrix.
!
!     on output
!
!        a contains the hessenberg matrix.  information about
!          the orthogonal transformations used in the reduction
!          is stored in the remaining triangle under the
!          hessenberg matrix.
!
!        ort contains further information about the transformations.
!          only elements low through igh are used.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
!
      do 180 m = kp1, la
         h = 0.0e0_rk
         ort(m) = 0.0e0_rk
         scale = 0.0e0_rk
!     .......... scale column (algol tol then not needed) ..........
         do 90 i = m, igh
   90    scale = scale + abs(a(i,m-1))
!
         if (scale .eq. 0.0e0_rk) go to 180
         mp = m + igh
!     .......... for i=igh step -1 until m do -- ..........
         do 100 ii = m, igh
            i = mp - ii
            ort(i) = a(i,m-1) / scale
            h = h + ort(i) * ort(i)
  100    continue
!
         g = -sign(sqrt(h),ort(m))
         h = h - ort(m) * g
         ort(m) = ort(m) - g
!     .......... form (i-(u*ut)/h) * a ..........
         do 130 j = m, n
            f = 0.0e0_rk
!     .......... for i=igh step -1 until m do -- ..........
            do 110 ii = m, igh
               i = mp - ii
               f = f + ort(i) * a(i,j)
  110       continue
!
            f = f / h
!
            do 120 i = m, igh
  120       a(i,j) = a(i,j) - f * ort(i)
!
  130    continue
!     .......... form (i-(u*ut)/h)*a*(i-(u*ut)/h) ..........
         do 160 i = 1, igh
            f = 0.0e0_rk
!     .......... for j=igh step -1 until m do -- ..........
            do 140 jj = m, igh
               j = mp - jj
               f = f + ort(j) * a(i,j)
  140       continue
!
            f = f / h
!
            do 150 j = m, igh
  150       a(i,j) = a(i,j) - f * ort(j)
!
  160    continue
!
         ort(m) = scale * ort(m)
         a(m,m-1) = scale * g
  180 continue
!
  200 return
      end subroutine
!*********************************************************************************
      subroutine ortran(nm,n,low,igh,a,ort,z)
!
      integer i,j,n,kl,mm,mp,nm,igh,low,mp1
      real(kind=rk)  :: a(nm,igh),ort(igh),z(nm,n)
      real(kind=rk)  :: g
!
!     this subroutine is a translation of the algol procedure ortrans,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine accumulates the orthogonal similarity
!     transformations used in the reduction of a real general
!     matrix to upper hessenberg form by  orthes.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        a contains information about the orthogonal trans-
!          formations used in the reduction by  orthes
!          in its strict lower triangle.
!
!        ort contains further information about the trans-
!          formations used in the reduction by  orthes.
!          only elements low through igh are used.
!
!     on output
!
!        z contains the transformation matrix produced in the
!          reduction by  orthes.
!
!        ort has been altered.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
!     .......... initialize z to identity matrix ..........
      do 80 j = 1, n
!
         do 60 i = 1, n
   60    z(i,j) = 0.0e0_rk
!
         z(j,j) = 1.0e0_rk
   80 continue
!
      kl = igh - low - 1
      if (kl .lt. 1) go to 200
!     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = 1, kl
         mp = igh - mm
         if (a(mp,mp-1) .eq. 0.0e0_rk) go to 140
         mp1 = mp + 1
!
         do 100 i = mp1, igh
  100    ort(i) = a(i,mp-1)
!
         do 130 j = mp, igh
            g = 0.0e0_rk
!
            do 110 i = mp, igh
  110       g = g + ort(i) * z(i,j)
!     .......... divisor below is negative of h formed in orthes.
!                double division avoids possible underflow ..........
            g = (g / ort(mp)) / a(mp,mp-1)
!
            do 120 i = mp, igh
  120       z(i,j) = z(i,j) + g * ort(i)
!
  130    continue
!
  140 continue
!
  200 return
      end subroutine
!*********************************************************************************
      subroutine qzhes(nm,n,a,b,matz,z)
!
      integer i,j,k,l,n,lb,l1,nm,nk1,nm1,nm2
      real(kind=rk)  :: a(nm,n),b(nm,n),z(nm,n)
      real(kind=rk)  :: r,s,t,u1,u2,v1,v2,rho
      logical matz
!
!     this subroutine is the first step of the qz algorithm
!     for solving generalized matrix eigenvalue problems,
!     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
!
!     this subroutine accepts a pair of real general matrices and
!     reduces one of them to upper hessenberg form and the other
!     to upper triangular form using orthogonal transformations.
!     it is usually followed by  qzit,  qzval  and, possibly,  qzvec.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrices.
!
!        a contains a real general matrix.
!
!        b contains a real general matrix.
!
!        matz should be set to .true. if the right hand transformations
!          are to be accumulated for later use in computing
!          eigenvectors, and to .false. otherwise.
!
!     on output
!
!        a has been reduced to upper hessenberg form.  the elements
!          below the first subdiagonal have been set to zero.
!
!        b has been reduced to upper triangular form.  the elements
!          below the main diagonal have been set to zero.
!
!        z contains the product of the right hand transformations if
!          matz has been set to .true.  otherwise, z is not referenced.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
!     .......... initialize z ..........
      if (.not. matz) go to 10
!
      do 3 j = 1, n
!
         do 2 i = 1, n
            z(i,j) = 0.0e0_rk
    2    continue
!
         z(j,j) = 1.0e0_rk
    3 continue
!     .......... reduce b to upper triangular form ..........
   10 if (n .le. 1) go to 170
      nm1 = n - 1
!
      do 100 l = 1, nm1
         l1 = l + 1
         s = 0.0e0_rk
!
         do 20 i = l1, n
            s = s + abs(b(i,l))
   20    continue
!
         if (s .eq. 0.0e0_rk) go to 100
         s = s + abs(b(l,l))
         r = 0.0e0_rk
!
         do 25 i = l, n
            b(i,l) = b(i,l) / s
            r = r + b(i,l)**2
   25    continue
!
         r = sign(sqrt(r),b(l,l))
         b(l,l) = b(l,l) + r
         rho = r * b(l,l)
!
         do 50 j = l1, n
            t = 0.0e0_rk
!
            do 30 i = l, n
               t = t + b(i,l) * b(i,j)
   30       continue
!
            t = -t / rho
!
            do 40 i = l, n
               b(i,j) = b(i,j) + t * b(i,l)
   40       continue
!
   50    continue
!
         do 80 j = 1, n
            t = 0.0e0_rk
!
            do 60 i = l, n
               t = t + b(i,l) * a(i,j)
   60       continue
!
            t = -t / rho
!
            do 70 i = l, n
               a(i,j) = a(i,j) + t * b(i,l)
   70       continue
!
   80    continue
!
         b(l,l) = -s * r
!
         do 90 i = l1, n
            b(i,l) = 0.0e0_rk
   90    continue
!
  100 continue
!     .......... reduce a to upper hessenberg form, while
!                keeping b triangular ..........
      if (n .eq. 2) go to 170
      nm2 = n - 2
!
      do 160 k = 1, nm2
         nk1 = nm1 - k
!     .......... for l=n-1 step -1 until k+1 do -- ..........
         do 150 lb = 1, nk1
            l = n - lb
            l1 = l + 1
!     .......... zero a(l+1,k) ..........
            s = abs(a(l,k)) + abs(a(l1,k))
            if (s .eq. 0.0e0_rk) go to 150
            u1 = a(l,k) / s
            u2 = a(l1,k) / s
            r = sign(sqrt(u1*u1+u2*u2),u1)
            v1 =  -(u1 + r) / r
            v2 = -u2 / r
            u2 = v2 / v1
!
            do 110 j = k, n
               t = a(l,j) + u2 * a(l1,j)
               a(l,j) = a(l,j) + t * v1
               a(l1,j) = a(l1,j) + t * v2
  110       continue
!
            a(l1,k) = 0.0e0_rk
!
            do 120 j = l, n
               t = b(l,j) + u2 * b(l1,j)
               b(l,j) = b(l,j) + t * v1
               b(l1,j) = b(l1,j) + t * v2
  120       continue
!     .......... zero b(l+1,l) ..........
            s = abs(b(l1,l1)) + abs(b(l1,l))
            if (s .eq. 0.0e0_rk) go to 150
            u1 = b(l1,l1) / s
            u2 = b(l1,l) / s
            r = sign(sqrt(u1*u1+u2*u2),u1)
            v1 =  -(u1 + r) / r
            v2 = -u2 / r
            u2 = v2 / v1
!
            do 130 i = 1, l1
               t = b(i,l1) + u2 * b(i,l)
               b(i,l1) = b(i,l1) + t * v1
               b(i,l) = b(i,l) + t * v2
  130       continue
!
            b(l1,l) = 0.0e0_rk
!
            do 140 i = 1, n
               t = a(i,l1) + u2 * a(i,l)
               a(i,l1) = a(i,l1) + t * v1
               a(i,l) = a(i,l) + t * v2
  140       continue
!
            if (.not. matz) go to 150
!
            do 145 i = 1, n
               t = z(i,l1) + u2 * z(i,l)
               z(i,l1) = z(i,l1) + t * v1
               z(i,l) = z(i,l) + t * v2
  145       continue
!
  150    continue
!
  160 continue
!
  170 return
      end subroutine
!*********************************************************************************
      subroutine qzit(nm,n,a,b,eps1,matz,z,ierr)
!
      integer i,j,k,l,n,en,k1,k2,ld,ll,l1,na,nm,ish,itn,its,km1,lm1, enm2,ierr,lor1,enorn
      real(kind=rk)  :: a(nm,n),b(nm,n),z(nm,n)
      real(kind=rk)  :: r,s,t,a1,a2,a3,ep,sh,u1,u2,u3,v1,v2,v3,ani,a11
      real(kind=rk)  :: a12,a21,a22,a33,a34,a43,a44,bni,b11,b12,b22,b33,b34
      real(kind=rk)  :: b44,epsa,epsb,eps1,anorm,bnorm
      logical matz,notlas
!
!     this subroutine is the second step of the qz algorithm
!     for solving generalized matrix eigenvalue problems,
!     siam j. numer. anal. 10, 241-256(1973) by moler and stewart,
!     as modified in technical note nasa tn d-7305(1973) by ward.
!
!     this subroutine accepts a pair of real matrices, one of them
!     in upper hessenberg form and the other in upper triangular form.
!     it reduces the hessenberg matrix to quasi-triangular form using
!     orthogonal transformations while maintaining the triangular form
!     of the other matrix.  it is usually preceded by  qzhes  and
!     followed by  qzval  and, possibly,  qzvec.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrices.
!
!        a contains a real upper hessenberg matrix.
!
!        b contains a real upper triangular matrix.
!
!        eps1 is a tolerance used to determine negligible elements.
!          eps1 = 0.0 (or negative) may be input, in which case an
!          element will be neglected only if it is less than roundoff
!          error times the norm of its matrix.  if the input eps1 is
!          positive, then an element will be considered negligible
!          if it is less than eps1 times the norm of its matrix.  a
!          positive value of eps1 may result in faster execution,
!          but less accurate results.
!
!        matz should be set to .true. if the right hand transformations
!          are to be accumulated for later use in computing
!          eigenvectors, and to .false. otherwise.
!
!        z contains, if matz has been set to .true., the
!          transformation matrix produced in the reduction
!          by  qzhes, if performed, or else the identity matrix.
!          if matz has been set to .false., z is not referenced.
!
!     on output
!
!        a has been reduced to quasi-triangular form.  the elements
!          below the first subdiagonal are still zero and no two
!          consecutive subdiagonal elements are nonzero.
!
!        b is still in upper triangular form, although its elements
!          have been altered.  the location b(n,1) is used to store
!          eps1 times the norm of b for later use by  qzval  and  qzvec.
!
!        z contains the product of the right hand transformations
!          (for both steps) if matz has been set to .true..
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
!     .......... compute epsa,epsb ..........
      anorm = 0.0e0_rk
      bnorm = 0.0e0_rk
!
      do 30 i = 1, n
         ani = 0.0e0_rk
         if (i .ne. 1) ani = abs(a(i,i-1))
         bni = 0.0e0_rk
!
         do 20 j = i, n
            ani = ani + abs(a(i,j))
            bni = bni + abs(b(i,j))
   20    continue
!
         if (ani .gt. anorm) anorm = ani
         if (bni .gt. bnorm) bnorm = bni
   30 continue
!
      if (anorm .eq. 0.0e0_rk) anorm = 1.0e0_rk
      if (bnorm .eq. 0.0e0_rk) bnorm = 1.0e0_rk
      ep = eps1
      if (ep .gt. 0.0e0_rk) go to 50
!     .......... use roundoff level if eps1 is zero ..........
      ep = epslon(1.0e0_rk)
   50 epsa = ep * anorm
      epsb = ep * bnorm
!     .......... reduce a to quasi-triangular form, while
!                keeping b triangular ..........
      lor1 = 1
      enorn = n
      en = n
      itn = 30*n
!     .......... begin qz step ..........
   60 if (en .le. 2) go to 1001
      if (.not. matz) enorn = en
      its = 0
      na = en - 1
      enm2 = na - 1
   70 ish = 2
!     .......... check for convergence or reducibility.
!                for l=en step -1 until 1 do -- ..........
      do 80 ll = 1, en
         lm1 = en - ll
         l = lm1 + 1
         if (l .eq. 1) go to 95
         if (abs(a(l,lm1)) .le. epsa) go to 90
   80 continue
!
   90 a(l,lm1) = 0.0e0_rk
      if (l .lt. na) go to 95
!     .......... 1-by-1 or 2-by-2 block isolated ..........
      en = lm1
      go to 60
!     .......... check for small top of b ..........
   95 ld = l
  100 l1 = l + 1
      b11 = b(l,l)
      if (abs(b11) .gt. epsb) go to 120
      b(l,l) = 0.0e0_rk
      s = abs(a(l,l)) + abs(a(l1,l))
      u1 = a(l,l) / s
      u2 = a(l1,l) / s
      r = sign(sqrt(u1*u1+u2*u2),u1)
      v1 = -(u1 + r) / r
      v2 = -u2 / r
      u2 = v2 / v1
!
      do 110 j = l, enorn
         t = a(l,j) + u2 * a(l1,j)
         a(l,j) = a(l,j) + t * v1
         a(l1,j) = a(l1,j) + t * v2
         t = b(l,j) + u2 * b(l1,j)
         b(l,j) = b(l,j) + t * v1
         b(l1,j) = b(l1,j) + t * v2
  110 continue
!
      if (l .ne. 1) a(l,lm1) = -a(l,lm1)
      lm1 = l
      l = l1
      go to 90
  120 a11 = a(l,l) / b11
      a21 = a(l1,l) / b11
      if (ish .eq. 1) go to 140
!     .......... iteration strategy ..........
      if (itn .eq. 0) go to 1000
      if (its .eq. 10) go to 155
!     .......... determine type of shift ..........
      b22 = b(l1,l1)
      if (abs(b22) .lt. epsb) b22 = epsb
      b33 = b(na,na)
      if (abs(b33) .lt. epsb) b33 = epsb
      b44 = b(en,en)
      if (abs(b44) .lt. epsb) b44 = epsb
      a33 = a(na,na) / b33
      a34 = a(na,en) / b44
      a43 = a(en,na) / b33
      a44 = a(en,en) / b44
      b34 = b(na,en) / b44
      t = 0.5e0_rk * (a43 * b34 - a33 - a44)
      r = t * t + a34 * a43 - a33 * a44
      if (r .lt. 0.0e0_rk) go to 150
!     .......... determine single shift zeroth column of a ..........
      ish = 1
      r = sqrt(r)
      sh = -t + r
      s = -t - r
      if (abs(s-a44) .lt. abs(sh-a44)) sh = s
!     .......... look for two consecutive small
!                sub-diagonal elements of a.
!                for l=en-2 step -1 until ld do -- ..........
      do 130 ll = ld, enm2
         l = enm2 + ld - ll
         if (l .eq. ld) go to 140
         lm1 = l - 1
         l1 = l + 1
         t = a(l,l)
         if (abs(b(l,l)) .gt. epsb) t = t - sh * b(l,l)
         if (abs(a(l,lm1)) .le. abs(t/a(l1,l)) * epsa) go to 100
  130 continue
!
  140 a1 = a11 - sh
      a2 = a21
      if (l .ne. ld) a(l,lm1) = -a(l,lm1)
      go to 160
!     .......... determine double shift zeroth column of a ..........
  150 a12 = a(l,l1) / b22
      a22 = a(l1,l1) / b22
      b12 = b(l,l1) / b22
      a1 = ((a33 - a11) * (a44 - a11) - a34 * a43 + a43 * b34 * a11) / a21 + a12 - a11 * b12
      a2 = (a22 - a11) - a21 * b12 - (a33 - a11) - (a44 - a11)+ a43 * b34
      a3 = a(l1+1,l1) / b22
      go to 160
!     .......... ad hoc shift ..........
  155 a1 = 0.0e0_rk
      a2 = 1.0e0_rk
      a3 = 1.1605e0_rk
  160 its = its + 1
      itn = itn - 1
      if (.not. matz) lor1 = ld
!     .......... main loop ..........
      do 260 k = l, na
         notlas = k .ne. na .and. ish .eq. 2
         k1 = k + 1
         k2 = k + 2
         km1 = max0(k-1,l)
         ll = min0(en,k1+ish)
         if (notlas) go to 190
!     .......... zero a(k+1,k-1) ..........
         if (k .eq. l) go to 170
         a1 = a(k,km1)
         a2 = a(k1,km1)
  170    s = abs(a1) + abs(a2)
         if (s .eq. 0.0e0_rk) go to 70
         u1 = a1 / s
         u2 = a2 / s
         r = sign(sqrt(u1*u1+u2*u2),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         u2 = v2 / v1
!
         do 180 j = km1, enorn
            t = a(k,j) + u2 * a(k1,j)
            a(k,j) = a(k,j) + t * v1
            a(k1,j) = a(k1,j) + t * v2
            t = b(k,j) + u2 * b(k1,j)
            b(k,j) = b(k,j) + t * v1
            b(k1,j) = b(k1,j) + t * v2
  180    continue
!
         if (k .ne. l) a(k1,km1) = 0.0e0_rk
         go to 240
!     .......... zero a(k+1,k-1) and a(k+2,k-1) ..........
  190    if (k .eq. l) go to 200
         a1 = a(k,km1)
         a2 = a(k1,km1)
         a3 = a(k2,km1)
  200    s = abs(a1) + abs(a2) + abs(a3)
         if (s .eq. 0.0e0_rk) go to 260
         u1 = a1 / s
         u2 = a2 / s
         u3 = a3 / s
         r = sign(sqrt(u1*u1+u2*u2+u3*u3),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         v3 = -u3 / r
         u2 = v2 / v1
         u3 = v3 / v1
!
         do 210 j = km1, enorn
            t = a(k,j) + u2 * a(k1,j) + u3 * a(k2,j)
            a(k,j) = a(k,j) + t * v1
            a(k1,j) = a(k1,j) + t * v2
            a(k2,j) = a(k2,j) + t * v3
            t = b(k,j) + u2 * b(k1,j) + u3 * b(k2,j)
            b(k,j) = b(k,j) + t * v1
            b(k1,j) = b(k1,j) + t * v2
            b(k2,j) = b(k2,j) + t * v3
  210    continue
!
         if (k .eq. l) go to 220
         a(k1,km1) = 0.0e0_rk
         a(k2,km1) = 0.0e0_rk
!     .......... zero b(k+2,k+1) and b(k+2,k) ..........
  220    s = abs(b(k2,k2)) + abs(b(k2,k1)) + abs(b(k2,k))
         if (s .eq. 0.0e0_rk) go to 240
         u1 = b(k2,k2) / s
         u2 = b(k2,k1) / s
         u3 = b(k2,k) / s
         r = sign(sqrt(u1*u1+u2*u2+u3*u3),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         v3 = -u3 / r
         u2 = v2 / v1
         u3 = v3 / v1
!
         do 230 i = lor1, ll
            t = a(i,k2) + u2 * a(i,k1) + u3 * a(i,k)
            a(i,k2) = a(i,k2) + t * v1
            a(i,k1) = a(i,k1) + t * v2
            a(i,k) = a(i,k) + t * v3
            t = b(i,k2) + u2 * b(i,k1) + u3 * b(i,k)
            b(i,k2) = b(i,k2) + t * v1
            b(i,k1) = b(i,k1) + t * v2
            b(i,k) = b(i,k) + t * v3
  230    continue
!
         b(k2,k) = 0.0e0_rk
         b(k2,k1) = 0.0e0_rk
         if (.not. matz) go to 240
!
         do 235 i = 1, n
            t = z(i,k2) + u2 * z(i,k1) + u3 * z(i,k)
            z(i,k2) = z(i,k2) + t * v1
            z(i,k1) = z(i,k1) + t * v2
            z(i,k) = z(i,k) + t * v3
  235    continue
!     .......... zero b(k+1,k) ..........
  240    s = abs(b(k1,k1)) + abs(b(k1,k))
         if (s .eq. 0.0e0_rk) go to 260
         u1 = b(k1,k1) / s
         u2 = b(k1,k) / s
         r = sign(sqrt(u1*u1+u2*u2),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         u2 = v2 / v1
!
         do 250 i = lor1, ll
            t = a(i,k1) + u2 * a(i,k)
            a(i,k1) = a(i,k1) + t * v1
            a(i,k) = a(i,k) + t * v2
            t = b(i,k1) + u2 * b(i,k)
            b(i,k1) = b(i,k1) + t * v1
            b(i,k) = b(i,k) + t * v2
  250    continue
!
         b(k1,k) = 0.0e0_rk
         if (.not. matz) go to 260
!
         do 255 i = 1, n
            t = z(i,k1) + u2 * z(i,k)
            z(i,k1) = z(i,k1) + t * v1
            z(i,k) = z(i,k) + t * v2
  255    continue
!
  260 continue
!     .......... end qz step ..........
      go to 70
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 ierr = en
!     .......... save epsb for use by qzval and qzvec ..........
 1001 if (n .gt. 1) b(n,1) = epsb
      return
      end subroutine
!*********************************************************************************
      subroutine qzval(nm,n,a,b,alfr,alfi,beta,matz,z)
!
      integer i,j,n,en,na,nm,nn,isw
      real(kind=rk)  :: a(nm,n),b(nm,n),alfr(n),alfi(n),beta(n),z(nm,n)
      real(kind=rk)  :: c,d,e,r,s,t,an,a1,a2,bn,cq,cz,di,dr,ei,ti,tr,u1
      real(kind=rk)  :: u2,v1,v2,a1i,a11,a12,a2i,a21,a22,b11,b12,b22,sqi,sqr
      real(kind=rk)  :: ssi,ssr,szi,szr,a11i,a11r,a12i,a12r,a22i,a22r,epsb
      logical matz
!
!     this subroutine is the third step of the qz algorithm
!     for solving generalized matrix eigenvalue problems,
!     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
!
!     this subroutine accepts a pair of real matrices, one of them
!     in quasi-triangular form and the other in upper triangular form.
!     it reduces the quasi-triangular matrix further, so that any
!     remaining 2-by-2 blocks correspond to pairs of complex
!     eigenvalues, and returns quantities whose ratios give the
!     generalized eigenvalues.  it is usually preceded by  qzhes
!     and  qzit  and may be followed by  qzvec.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrices.
!
!        a contains a real upper quasi-triangular matrix.
!
!        b contains a real upper triangular matrix.  in addition,
!          location b(n,1) contains the tolerance quantity (epsb)
!          computed and saved in  qzit.
!
!        matz should be set to .true. if the right hand transformations
!          are to be accumulated for later use in computing
!          eigenvectors, and to .false. otherwise.
!
!        z contains, if matz has been set to .true., the
!          transformation matrix produced in the reductions by qzhes
!          and qzit, if performed, or else the identity matrix.
!          if matz has been set to .false., z is not referenced.
!
!     on output
!
!        a has been reduced further to a quasi-triangular matrix
!          in which all nonzero subdiagonal elements correspond to
!          pairs of complex eigenvalues.
!
!        b is still in upper triangular form, although its elements
!          have been altered.  b(n,1) is unaltered.
!
!        alfr and alfi contain the real and imaginary parts of the
!          diagonal elements of the triangular matrix that would be
!          obtained if a were reduced completely to triangular form
!          by unitary transformations.  non-zero values of alfi occur
!          in pairs, the first member positive and the second negative.
!
!        beta contains the diagonal elements of the corresponding b,
!          normalized to be real and non-negative.  the generalized
!          eigenvalues are then the ratios ((alfr+i*alfi)/beta).
!
!        z contains the product of the right hand transformations
!          (for all three steps) if matz has been set to .true.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      epsb = b(n,1)
      isw = 1
!     .......... find eigenvalues of quasi-triangular matrices.
!                for en=n step -1 until 1 do -- ..........
      do 510 nn = 1, n
         en = n + 1 - nn
         na = en - 1
         if (isw .eq. 2) go to 505
         if (en .eq. 1) go to 410
         if (a(en,na) .ne. 0.0e0_rk) go to 420
!     .......... 1-by-1 block, one real root ..........
  410    alfr(en) = a(en,en)
         if (b(en,en) .lt. 0.0e0_rk) alfr(en) = -alfr(en)
         beta(en) = abs(b(en,en))
         alfi(en) = 0.0e0_rk
         go to 510
!     .......... 2-by-2 block ..........
  420    if (abs(b(na,na)) .le. epsb) go to 455
         if (abs(b(en,en)) .gt. epsb) go to 430
         a1 = a(en,en)
         a2 = a(en,na)
         bn = 0.0e0_rk
         go to 435
  430    an = abs(a(na,na)) + abs(a(na,en)) + abs(a(en,na)) + abs(a(en,en))
         bn = abs(b(na,na)) + abs(b(na,en)) + abs(b(en,en))
         a11 = a(na,na) / an
         a12 = a(na,en) / an
         a21 = a(en,na) / an
         a22 = a(en,en) / an
         b11 = b(na,na) / bn
         b12 = b(na,en) / bn
         b22 = b(en,en) / bn
         e = a11 / b11
         ei = a22 / b22
         s = a21 / (b11 * b22)
         t = (a22 - e * b22) / b22
         if (abs(e) .le. abs(ei)) go to 431
         e = ei
         t = (a11 - e * b11) / b11
  431    c = 0.5e0_rk * (t - s * b12)
         d = c * c + s * (a12 - e * b12)
         if (d .lt. 0.0e0_rk) go to 480
!     .......... two real roots.
!                zero both a(en,na) and b(en,na) ..........
         e = e + (c + sign(sqrt(d),c))
         a11 = a11 - e * b11
         a12 = a12 - e * b12
         a22 = a22 - e * b22
         if (abs(a11) + abs(a12) .lt. abs(a21) + abs(a22)) go to 432
         a1 = a12
         a2 = a11
         go to 435
  432    a1 = a22
         a2 = a21
!     .......... choose and apply real z ..........
  435    s = abs(a1) + abs(a2)
         u1 = a1 / s
         u2 = a2 / s
         r = sign(sqrt(u1*u1+u2*u2),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         u2 = v2 / v1
!
         do 440 i = 1, en
            t = a(i,en) + u2 * a(i,na)
            a(i,en) = a(i,en) + t * v1
            a(i,na) = a(i,na) + t * v2
            t = b(i,en) + u2 * b(i,na)
            b(i,en) = b(i,en) + t * v1
            b(i,na) = b(i,na) + t * v2
  440    continue
!
         if (.not. matz) go to 450
!
         do 445 i = 1, n
            t = z(i,en) + u2 * z(i,na)
            z(i,en) = z(i,en) + t * v1
            z(i,na) = z(i,na) + t * v2
  445    continue
!
  450    if (bn .eq. 0.0e0_rk) go to 475
         if (an .lt. abs(e) * bn) go to 455
         a1 = b(na,na)
         a2 = b(en,na)
         go to 460
  455    a1 = a(na,na)
         a2 = a(en,na)
!     .......... choose and apply real q ..........
  460    s = abs(a1) + abs(a2)
         if (s .eq. 0.0e0_rk) go to 475
         u1 = a1 / s
         u2 = a2 / s
         r = sign(sqrt(u1*u1+u2*u2),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         u2 = v2 / v1
!
         do 470 j = na, n
            t = a(na,j) + u2 * a(en,j)
            a(na,j) = a(na,j) + t * v1
            a(en,j) = a(en,j) + t * v2
            t = b(na,j) + u2 * b(en,j)
            b(na,j) = b(na,j) + t * v1
            b(en,j) = b(en,j) + t * v2
  470    continue
!
  475    a(en,na) = 0.0e0_rk
         b(en,na) = 0.0e0_rk
         alfr(na) = a(na,na)
         alfr(en) = a(en,en)
         if (b(na,na) .lt. 0.0e0_rk) alfr(na) = -alfr(na)
         if (b(en,en) .lt. 0.0e0_rk) alfr(en) = -alfr(en)
         beta(na) = abs(b(na,na))
         beta(en) = abs(b(en,en))
         alfi(en) = 0.0e0_rk
         alfi(na) = 0.0e0_rk
         go to 505
!     .......... two complex roots ..........
  480    e = e + c
         ei = sqrt(-d)
         a11r = a11 - e * b11
         a11i = ei * b11
         a12r = a12 - e * b12
         a12i = ei * b12
         a22r = a22 - e * b22
         a22i = ei * b22
         if (abs(a11r) + abs(a11i) + abs(a12r) + abs(a12i) .lt.abs(a21) + abs(a22r) + abs(a22i)) go to 482
         a1 = a12r
         a1i = a12i
         a2 = -a11r
         a2i = -a11i
         go to 485
  482    a1 = a22r
         a1i = a22i
         a2 = -a21
         a2i = 0.0e0_rk
!     .......... choose complex z ..........
  485    cz = sqrt(a1*a1+a1i*a1i)
         if (cz .eq. 0.0e0_rk) go to 487
         szr = (a1 * a2 + a1i * a2i) / cz
         szi = (a1 * a2i - a1i * a2) / cz
         r = sqrt(cz*cz+szr*szr+szi*szi)
         cz = cz / r
         szr = szr / r
         szi = szi / r
         go to 490
  487    szr = 1.0e0_rk
         szi = 0.0e0_rk
  490    if (an .lt. (abs(e) + ei) * bn) go to 492
         a1 = cz * b11 + szr * b12
         a1i = szi * b12
         a2 = szr * b22
         a2i = szi * b22
         go to 495
  492    a1 = cz * a11 + szr * a12
         a1i = szi * a12
         a2 = cz * a21 + szr * a22
         a2i = szi * a22
!     .......... choose complex q ..........
  495    cq = sqrt(a1*a1+a1i*a1i)
         if (cq .eq. 0.0e0_rk) go to 497
         sqr = (a1 * a2 + a1i * a2i) / cq
         sqi = (a1 * a2i - a1i * a2) / cq
         r = sqrt(cq*cq+sqr*sqr+sqi*sqi)
         cq = cq / r
         sqr = sqr / r
         sqi = sqi / r
         go to 500
  497    sqr = 1.0e0_rk
         sqi = 0.0e0_rk
!     .......... compute diagonal elements that would result
!                if transformations were applied ..........
  500    ssr = sqr * szr + sqi * szi
         ssi = sqr * szi - sqi * szr
         i = 1
         tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21 + ssr * a22
         ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22
         dr = cq * cz * b11 + cq * szr * b12 + ssr * b22
         di = cq * szi * b12 + ssi * b22
         go to 503
  502    i = 2
         tr = ssr * a11 - sqr * cz * a12 - cq * szr * a21 + cq * cz * a22
         ti = -ssi * a11 - sqi * cz * a12 + cq * szi * a21
         dr = ssr * b11 - sqr * cz * b12 + cq * cz * b22
         di = -ssi * b11 - sqi * cz * b12
  503    t = ti * dr - tr * di
         j = na
         if (t .lt. 0.0e0_rk) j = en
         r = sqrt(dr*dr+di*di)
         beta(j) = bn * r
         alfr(j) = an * (tr * dr + ti * di) / r
         alfi(j) = an * t / r
         if (i .eq. 1) go to 502
  505    isw = 3 - isw
  510 continue
      b(n,1) = epsb
!
      return
      end subroutine
!*********************************************************************************
      subroutine qzvec(nm,n,a,b,alfr,alfi,beta,z)
!
      integer i,j,k,m,n,en,ii,jj,na,nm,nn,isw,enm2
      real(kind=rk)  :: a(nm,n),b(nm,n),alfr(n),alfi(n),beta(n),z(nm,n)
      real(kind=rk)  :: d,q,r,s,t,w,x,y,di,dr,ra,rr,sa,ti,tr,t1,t2,w1,x1,zz,z1,alfm,almi,almr,betm,epsb
!
!     this subroutine is the optional fourth step of the qz algorithm
!     for solving generalized matrix eigenvalue problems,
!     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
!
!     this subroutine accepts a pair of real matrices, one of them in
!     quasi-triangular form (in which each 2-by-2 block corresponds to
!     a pair of complex eigenvalues) and the other in upper triangular
!     form.  it computes the eigenvectors of the triangular problem and
!     transforms the results back to the original coordinate system.
!     it is usually preceded by  qzhes,  qzit, and  qzval.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrices.
!
!        a contains a real upper quasi-triangular matrix.
!
!        b contains a real upper triangular matrix.  in addition,
!          location b(n,1) contains the tolerance quantity (epsb)
!          computed and saved in  qzit.
!
!        alfr, alfi, and beta  are vectors with components whose
!          ratios ((alfr+i*alfi)/beta) are the generalized
!          eigenvalues.  they are usually obtained from  qzval.
!
!        z contains the transformation matrix produced in the
!          reductions by  qzhes,  qzit, and  qzval, if performed.
!          if the eigenvectors of the triangular problem are
!          desired, z must contain the identity matrix.
!
!     on output
!
!        a is unaltered.  its subdiagonal elements provide information
!           about the storage of the complex eigenvectors.
!
!        b has been destroyed.
!
!        alfr, alfi, and beta are unaltered.
!
!        z contains the real and imaginary parts of the eigenvectors.
!          if alfi(i) .eq. 0.0, the i-th eigenvalue is real and
!            the i-th column of z contains its eigenvector.
!          if alfi(i) .ne. 0.0, the i-th eigenvalue is complex.
!            if alfi(i) .gt. 0.0, the eigenvalue is the first of
!              a complex pair and the i-th and (i+1)-th columns
!              of z contain its eigenvector.
!            if alfi(i) .lt. 0.0, the eigenvalue is the second of
!              a complex pair and the (i-1)-th and i-th columns
!              of z contain the conjugate of its eigenvector.
!          each eigenvector is normalized so that the modulus
!          of its largest component is 1.0 .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      epsb = b(n,1)
      isw = 1
!     .......... for en=n step -1 until 1 do -- ..........
      do 800 nn = 1, n
         en = n + 1 - nn
         na = en - 1
         if (isw .eq. 2) go to 795
         if (alfi(en) .ne. 0.0e0_rk) go to 710
!     .......... real vector ..........
         m = en
         b(en,en) = 1.0e0_rk
         if (na .eq. 0) go to 800
         alfm = alfr(m)
         betm = beta(m)
!     .......... for i=en-1 step -1 until 1 do -- ..........
         do 700 ii = 1, na
            i = en - ii
            w = betm * a(i,i) - alfm * b(i,i)
            r = 0.0e0_rk
!
            do 610 j = m, en
  610       r = r + (betm * a(i,j) - alfm * b(i,j)) * b(j,en)
!
            if (i .eq. 1 .or. isw .eq. 2) go to 630
            if (betm * a(i,i-1) .eq. 0.0e0_rk) go to 630
            zz = w
            s = r
            go to 690
  630       m = i
            if (isw .eq. 2) go to 640
!     .......... real 1-by-1 block ..........
            t = w
            if (w .eq. 0.0e0_rk) t = epsb
            b(i,en) = -r / t
            go to 700
!     .......... real 2-by-2 block ..........
  640       x = betm * a(i,i+1) - alfm * b(i,i+1)
            y = betm * a(i+1,i)
            q = w * zz - x * y
            t = (x * s - zz * r) / q
            b(i,en) = t
            if (abs(x) .le. abs(zz)) go to 650
            b(i+1,en) = (-r - w * t) / x
            go to 690
  650       b(i+1,en) = (-s - y * t) / zz
  690       isw = 3 - isw
  700    continue
!     .......... end real vector ..........
         go to 800
!     .......... complex vector ..........
  710    m = na
         almr = alfr(m)
         almi = alfi(m)
         betm = beta(m)
!     .......... last vector component chosen imaginary so that
!                eigenvector matrix is triangular ..........
         y = betm * a(en,na)
         b(na,na) = -almi * b(en,en) / y
         b(na,en) = (almr * b(en,en) - betm * a(en,en)) / y
         b(en,na) = 0.0e0_rk
         b(en,en) = 1.0e0_rk
         enm2 = na - 1
         if (enm2 .eq. 0) go to 795
!     .......... for i=en-2 step -1 until 1 do -- ..........
         do 790 ii = 1, enm2
            i = na - ii
            w = betm * a(i,i) - almr * b(i,i)
            w1 = -almi * b(i,i)
            ra = 0.0e0_rk
            sa = 0.0e0_rk
!
            do 760 j = m, en
               x = betm * a(i,j) - almr * b(i,j)
               x1 = -almi * b(i,j)
               ra = ra + x * b(j,na) - x1 * b(j,en)
               sa = sa + x * b(j,en) + x1 * b(j,na)
  760       continue
!
            if (i .eq. 1 .or. isw .eq. 2) go to 770
            if (betm * a(i,i-1) .eq. 0.0e0_rk) go to 770
            zz = w
            z1 = w1
            r = ra
            s = sa
            isw = 2
            go to 790
  770       m = i
            if (isw .eq. 2) go to 780
!     .......... complex 1-by-1 block ..........
            tr = -ra
            ti = -sa
  773       dr = w
            di = w1
!     .......... complex divide (t1,t2) = (tr,ti) / (dr,di) ..........
  775       if (abs(di) .gt. abs(dr)) go to 777
            rr = di / dr
            d = dr + di * rr
            t1 = (tr + ti * rr) / d
            t2 = (ti - tr * rr) / d
            go to (787,782), isw
  777       rr = dr / di
            d = dr * rr + di
            t1 = (tr * rr + ti) / d
            t2 = (ti * rr - tr) / d
            go to (787,782), isw
!     .......... complex 2-by-2 block ..........
  780       x = betm * a(i,i+1) - almr * b(i,i+1)
            x1 = -almi * b(i,i+1)
            y = betm * a(i+1,i)
            tr = y * ra - w * r + w1 * s
            ti = y * sa - w * s - w1 * r
            dr = w * zz - w1 * z1 - x * y
            di = w * z1 + w1 * zz - x1 * y
            if (dr .eq. 0.0e0_rk .and. di .eq. 0.0e0_rk) dr = epsb
            go to 775
  782       b(i+1,na) = t1
            b(i+1,en) = t2
            isw = 1
            if (abs(y) .gt. abs(w) + abs(w1)) go to 785
            tr = -ra - x * b(i+1,na) + x1 * b(i+1,en)
            ti = -sa - x * b(i+1,en) - x1 * b(i+1,na)
            go to 773
  785       t1 = (-r - zz * b(i+1,na) + z1 * b(i+1,en)) / y
            t2 = (-s - zz * b(i+1,en) - z1 * b(i+1,na)) / y
  787       b(i,na) = t1
            b(i,en) = t2
  790    continue
!     .......... end complex vector ..........
  795    isw = 3 - isw
  800 continue
!     .......... end back substitution.
!                transform to original coordinate system.
!                for j=n step -1 until 1 do -- ..........
      do 880 jj = 1, n
         j = n + 1 - jj
!
         do 880 i = 1, n
            zz = 0.0e0_rk
!
            do 860 k = 1, j
  860       zz = zz + z(i,k) * b(k,j)
!
            z(i,j) = zz
  880 continue
!     .......... normalize so that modulus of largest
!                component of each vector is 1.
!                (isw is 1 initially from before) ..........
      do 950 j = 1, n
         d = 0.0e0_rk
         if (isw .eq. 2) go to 920
         if (alfi(j) .ne. 0.0e0_rk) go to 945
!
         do 890 i = 1, n
            if (abs(z(i,j)) .gt. d) d = abs(z(i,j))
  890    continue
!
         do 900 i = 1, n
  900    z(i,j) = z(i,j) / d
!
         go to 950
!
  920    do 930 i = 1, n
            r = abs(z(i,j-1)) + abs(z(i,j))
            if (r .ne. 0.0e0_rk) r = r * sqrt((z(i,j-1)/r)**2 +(z(i,j)/r)**2)
            if (r .gt. d) d = r
  930    continue
!
         do 940 i = 1, n
            z(i,j-1) = z(i,j-1) / d
            z(i,j) = z(i,j) / d
  940    continue
!
  945    isw = 3 - isw
  950 continue
!
      return
      end subroutine
!*********************************************************************************
      subroutine ratqr(n,eps1,d,e,e2,m,w,ind,bd,type,idef,ierr)
!
      integer i,j,k,m,n,ii,jj,k1,idef,ierr,jdef
      real(kind=rk)  :: d(n),e(n),e2(n),w(n),bd(n)
      real(kind=rk)  :: f,p,q,r,s,ep,qp,err,tot,eps1,delta
      integer ind(n)
      logical type
!
!     this subroutine is a translation of the algol procedure ratqr,
!     num. math. 11, 264-272(1968) by reinsch and bauer.
!     handbook for auto. comp., vol.ii-linear algebra, 257-265(1971).
!
!     this subroutine finds the algebraically smallest or largest
!     eigenvalues of a symmetric tridiagonal matrix by the
!     rational qr method with newton corrections.
!
!     on input
!
!        n is the order of the matrix.
!
!        eps1 is a theoretical absolute error tolerance for the
!          computed eigenvalues.  if the input eps1 is non-positive,
!          or indeed smaller than its default value, it is reset
!          at each iteration to the respective default value,
!          namely, the product of the relative machine precision
!          and the magnitude of the current eigenvalue iterate.
!          the theoretical absolute error in the k-th eigenvalue
!          is usually not greater than k times eps1.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2(1) is arbitrary.
!
!        m is the number of eigenvalues to be found.
!
!        idef should be set to 1 if the input matrix is known to be
!          positive definite, to -1 if the input matrix is known to
!          be negative definite, and to 0 otherwise.
!
!        type should be set to .true. if the smallest eigenvalues
!          are to be found, and to .false. if the largest eigenvalues
!          are to be found.
!
!     on output
!
!        eps1 is unaltered unless it has been reset to its
!          (last) default value.
!
!        d and e are unaltered (unless w overwrites d).
!
!        elements of e2, corresponding to elements of e regarded
!          as negligible, have been replaced by zero causing the
!          matrix to split into a direct sum of submatrices.
!          e2(1) is set to 0.0e0_rk if the smallest eigenvalues have been
!          found, and to 2.0e0_rk if the largest eigenvalues have been
!          found.  e2 is otherwise unaltered (unless overwritten by bd).
!
!        w contains the m algebraically smallest eigenvalues in
!          ascending order, or the m largest eigenvalues in
!          descending order.  if an error exit is made because of
!          an incorrect specification of idef, no eigenvalues
!          are found.  if the newton iterates for a particular
!          eigenvalue are not monotone, the best estimate obtained
!          is returned and ierr is set.  w may coincide with d.
!
!        ind contains in its first m positions the submatrix indices
!          associated with the corresponding eigenvalues in w --
!          1 for eigenvalues belonging to the first submatrix from
!          the top, 2 for those belonging to the second submatrix, etc..
!
!        bd contains refined bounds for the theoretical errors of the
!          corresponding eigenvalues in w.  these bounds are usually
!          within the tolerance specified by eps1.  bd may coincide
!          with e2.
!
!        ierr is set to
!          zero       for normal return,
!          6*n+1      if  idef  is set to 1 and  type  to .true.
!                     when the matrix is not positive definite, or
!                     if  idef  is set to -1 and  type  to .false.
!                     when the matrix is not negative definite,
!          5*n+k      if successive iterates to the k-th eigenvalue
!                     are not monotone increasing, where k refers
!                     to the last such occurrence.
!
!     note that subroutine tridib is generally faster and more
!     accurate than ratqr if the eigenvalues are clustered.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      jdef = idef
!     .......... copy d array into w ..........
      do 20 i = 1, n
   20 w(i) = d(i)
!
      if (type) go to 40
      j = 1
      go to 400
   40 err = 0.0e0_rk
      s = 0.0e0_rk
!     .......... look for small sub-diagonal entries and define
!                initial shift from lower gerschgorin bound.
!                copy e2 array into bd ..........
      tot = w(1)
      q = 0.0e0_rk
      j = 0
!
      do 100 i = 1, n
         p = q
         if (i .eq. 1) go to 60
         if (p .gt. epslon(abs(d(i)) + abs(d(i-1)))) go to 80
   60    e2(i) = 0.0e0_rk
   80    bd(i) = e2(i)
!     .......... count also if element of e2 has underflowed ..........
         if (e2(i) .eq. 0.0e0_rk) j = j + 1
         ind(i) = j
         q = 0.0e0_rk
         if (i .ne. n) q = abs(e(i+1))
         tot = min(w(i)-p-q,tot)
  100 continue
!
      if (jdef .eq. 1 .and. tot .lt. 0.0e0_rk) go to 140
!
      do 110 i = 1, n
  110 w(i) = w(i) - tot
!
      go to 160
  140 tot = 0.0e0_rk
!
  160 do 360 k = 1, m
!     .......... next qr transformation ..........
  180    tot = tot + s
         delta = w(n) - s
         i = n
         f = abs(epslon(tot))
         if (eps1 .lt. f) eps1 = f
         if (delta .gt. eps1) go to 190
         if (delta .lt. (-eps1)) go to 1000
         go to 300
!     .......... replace small sub-diagonal squares by zero
!                to reduce the incidence of underflows ..........
  190    if (k .eq. n) go to 210
         k1 = k + 1
         do 200 j = k1, n
            if (bd(j) .le. (epslon(w(j)+w(j-1))) ** 2) bd(j) = 0.0e0_rk
  200    continue
!
  210    f = bd(n) / delta
         qp = delta + f
         p = 1.0e0_rk
         if (k .eq. n) go to 260
         k1 = n - k
!     .......... for i=n-1 step -1 until k do -- ..........
         do 240 ii = 1, k1
            i = n - ii
            q = w(i) - s - f
            r = q / qp
            p = p * r + 1.0e0_rk
            ep = f * r
            w(i+1) = qp + ep
            delta = q - ep
            if (delta .gt. eps1) go to 220
            if (delta .lt. (-eps1)) go to 1000
            go to 300
  220       f = bd(i) / q
            qp = delta + f
            bd(i+1) = qp * ep
  240    continue
!
  260    w(k) = qp
         s = qp / p
         if (tot + s .gt. tot) go to 180
!     .......... set error -- irregular end of iteration.
!                deflate minimum diagonal element ..........
         ierr = 5 * n + k
         s = 0.0e0_rk
         delta = qp
!
         do 280 j = k, n
            if (w(j) .gt. delta) go to 280
            i = j
            delta = w(j)
  280    continue
!     .......... convergence ..........
  300    if (i .lt. n) bd(i+1) = bd(i) * f / qp
         ii = ind(i)
         if (i .eq. k) go to 340
         k1 = i - k
!     .......... for j=i-1 step -1 until k do -- ..........
         do 320 jj = 1, k1
            j = i - jj
            w(j+1) = w(j) - s
            bd(j+1) = bd(j)
            ind(j+1) = ind(j)
  320    continue
!
  340    w(k) = tot
         err = err + abs(delta)
         bd(k) = err
         ind(k) = ii
  360 continue
!
      if (type) go to 1001
      f = bd(1)
      e2(1) = 2.0e0_rk
      bd(1) = f
      j = 2
!     .......... negate elements of w for largest values ..........
  400 do 500 i = 1, n
  500 w(i) = -w(i)
!
      jdef = -jdef
      go to (40,1001), j
!     .......... set error -- idef specified incorrectly ..........
 1000 ierr = 6 * n + 1
 1001 return
      end subroutine
!*********************************************************************************
      subroutine rebak(nm,n,b,dl,m,z)
!
      integer i,j,k,m,n,i1,ii,nm
      real(kind=rk)  :: b(nm,n),dl(n),z(nm,m)
      real(kind=rk)  :: x
!
!     this subroutine is a translation of the algol procedure rebaka,
!     num. math. 11, 99-110(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
!
!     this subroutine forms the eigenvectors of a generalized
!     symmetric eigensystem by back transforming those of the
!     derived symmetric matrix determined by  reduc.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix system.
!
!        b contains information about the similarity transformation
!          (cholesky decomposition) used in the reduction by  reduc
!          in its strict lower triangle.
!
!        dl contains further information about the transformation.
!
!        m is the number of eigenvectors to be back transformed.
!
!        z contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        z contains the transformed eigenvectors
!          in its first m columns.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (m .eq. 0) go to 200
!
      do 100 j = 1, m
!     .......... for i=n step -1 until 1 do -- ..........
         do 100 ii = 1, n
            i = n + 1 - ii
            i1 = i + 1
            x = z(i,j)
            if (i .eq. n) go to 80
!
            do 60 k = i1, n
   60       x = x - b(k,i) * z(k,j)
!
   80       z(i,j) = x / dl(i)
  100 continue
!
  200 return
      end subroutine
!*********************************************************************************
      subroutine rebakb(nm,n,b,dl,m,z)
!
      integer i,j,k,m,n,i1,ii,nm
      real(kind=rk)  :: b(nm,n),dl(n),z(nm,m)
      real(kind=rk)  :: x
!
!     this subroutine is a translation of the algol procedure rebakb,
!     num. math. 11, 99-110(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
!
!     this subroutine forms the eigenvectors of a generalized
!     symmetric eigensystem by back transforming those of the
!     derived symmetric matrix determined by  reduc2.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix system.
!
!        b contains information about the similarity transformation
!          (cholesky decomposition) used in the reduction by  reduc2
!          in its strict lower triangle.
!
!        dl contains further information about the transformation.
!
!        m is the number of eigenvectors to be back transformed.
!
!        z contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        z contains the transformed eigenvectors
!          in its first m columns.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (m .eq. 0) go to 200
!
      do 100 j = 1, m
!     .......... for i=n step -1 until 1 do -- ..........
         do 100 ii = 1, n
            i1 = n - ii
            i = i1 + 1
            x = dl(i) * z(i,j)
            if (i .eq. 1) go to 80
!
            do 60 k = 1, i1
   60       x = x + b(i,k) * z(k,j)
!
   80       z(i,j) = x
  100 continue
!
  200 return
      end subroutine
!*********************************************************************************
      subroutine reduc(nm,n,a,b,dl,ierr)
!
      integer i,j,k,n,i1,j1,nm,nn,ierr
      real(kind=rk)  :: a(nm,n),b(nm,n),dl(n)
      real(kind=rk)  :: x,y
!
!     this subroutine is a translation of the algol procedure reduc1,
!     num. math. 11, 99-110(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
!
!     this subroutine reduces the generalized symmetric eigenproblem
!     ax=(lambda)bx, where b is positive definite, to the standard
!     symmetric eigenproblem using the cholesky factorization of b.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrices a and b.  if the cholesky
!          factor l of b is already available, n should be prefixed
!          with a minus sign.
!
!        a and b contain the real symmetric input matrices.  only the
!          full upper triangles of the matrices need be supplied.  if
!          n is negative, the strict lower triangle of b contains,
!          instead, the strict lower triangle of its cholesky factor l.
!
!        dl contains, if n is negative, the diagonal elements of l.
!
!     on output
!
!        a contains in its full lower triangle the full lower triangle
!          of the symmetric matrix derived from the reduction to the
!          standard form.  the strict upper triangle of a is unaltered.
!
!        b contains in its strict lower triangle the strict lower
!          triangle of its cholesky factor l.  the full upper
!          triangle of b is unaltered.
!
!        dl contains the diagonal elements of l.
!
!        ierr is set to
!          zero       for normal return,
!          7*n+1      if b is not positive definite.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      nn = iabs(n)
      if (n .lt. 0) go to 100
!     .......... form l in the arrays b and dl ..........
      do 80 i = 1, n
         i1 = i - 1
!
         do 80 j = i, n
            x = b(i,j)
            if (i .eq. 1) go to 40
!
            do 20 k = 1, i1
   20       x = x - b(i,k) * b(j,k)
!
   40       if (j .ne. i) go to 60
            if (x .le. 0.0e0_rk) go to 1000
            y = sqrt(x)
            dl(i) = y
            go to 80
   60       b(j,i) = x / y
   80 continue
!     .......... form the transpose of the upper triangle of inv(l)*a
!                in the lower triangle of the array a ..........
  100 do 200 i = 1, nn
         i1 = i - 1
         y = dl(i)
!
         do 200 j = i, nn
            x = a(i,j)
            if (i .eq. 1) go to 180
!
            do 160 k = 1, i1
  160       x = x - b(i,k) * a(j,k)
!
  180       a(j,i) = x / y
  200 continue
!     .......... pre-multiply by inv(l) and overwrite ..........
      do 300 j = 1, nn
         j1 = j - 1
!
         do 300 i = j, nn
            x = a(i,j)
            if (i .eq. j) go to 240
            i1 = i - 1
!
            do 220 k = j, i1
  220       x = x - a(k,j) * b(i,k)
!
  240       if (j .eq. 1) go to 280
!
            do 260 k = 1, j1
  260       x = x - a(j,k) * b(i,k)
!
  280       a(i,j) = x / dl(i)
  300 continue
!
      go to 1001
!     .......... set error -- b is not positive definite ..........
 1000 ierr = 7 * n + 1
 1001 return
      end subroutine
!*********************************************************************************
      subroutine reduc2(nm,n,a,b,dl,ierr)
!
      integer i,j,k,n,i1,j1,nm,nn,ierr
      real(kind=rk)  :: a(nm,n),b(nm,n),dl(n)
      real(kind=rk)  :: x,y
!
!     this subroutine is a translation of the algol procedure reduc2,
!     num. math. 11, 99-110(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
!
!     this subroutine reduces the generalized symmetric eigenproblems
!     abx=(lambda)x or bay=(lambda)y, where b is positive definite,
!     to the standard symmetric eigenproblem using the cholesky
!     factorization of b.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrices a and b.  if the cholesky
!          factor l of b is already available, n should be prefixed
!          with a minus sign.
!
!        a and b contain the real symmetric input matrices.  only the
!          full upper triangles of the matrices need be supplied.  if
!          n is negative, the strict lower triangle of b contains,
!          instead, the strict lower triangle of its cholesky factor l.
!
!        dl contains, if n is negative, the diagonal elements of l.
!
!     on output
!
!        a contains in its full lower triangle the full lower triangle
!          of the symmetric matrix derived from the reduction to the
!          standard form.  the strict upper triangle of a is unaltered.
!
!        b contains in its strict lower triangle the strict lower
!          triangle of its cholesky factor l.  the full upper
!          triangle of b is unaltered.
!
!        dl contains the diagonal elements of l.
!
!        ierr is set to
!          zero       for normal return,
!          7*n+1      if b is not positive definite.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      nn = iabs(n)
      if (n .lt. 0) go to 100
!     .......... form l in the arrays b and dl ..........
      do 80 i = 1, n
         i1 = i - 1
!
         do 80 j = i, n
            x = b(i,j)
            if (i .eq. 1) go to 40
!
            do 20 k = 1, i1
   20       x = x - b(i,k) * b(j,k)
!
   40       if (j .ne. i) go to 60
            if (x .le. 0.0e0_rk) go to 1000
            y = sqrt(x)
            dl(i) = y
            go to 80
   60       b(j,i) = x / y
   80 continue
!     .......... form the lower triangle of a*l
!                in the lower triangle of the array a ..........
  100 do 200 i = 1, nn
         i1 = i + 1
!
         do 200 j = 1, i
            x = a(j,i) * dl(j)
            if (j .eq. i) go to 140
            j1 = j + 1
!
            do 120 k = j1, i
  120       x = x + a(k,i) * b(k,j)
!
  140       if (i .eq. nn) go to 180
!
            do 160 k = i1, nn
  160       x = x + a(i,k) * b(k,j)
!
  180       a(i,j) = x
  200 continue
!     .......... pre-multiply by transpose(l) and overwrite ..........
      do 300 i = 1, nn
         i1 = i + 1
         y = dl(i)
!
         do 300 j = 1, i
            x = y * a(i,j)
            if (i .eq. nn) go to 280
!
            do 260 k = i1, nn
  260       x = x + a(k,j) * b(k,i)
!
  280       a(i,j) = x
  300 continue
!
      go to 1001
!     .......... set error -- b is not positive definite ..........
 1000 ierr = 7 * n + 1
 1001 return
      end subroutine
!*********************************************************************************
      subroutine rgg(nm,n,a,b,alfr,alfi,beta,matz,z,ierr)
!
      integer n,nm,ierr,matz
      real(kind=rk)  :: a(nm,n),b(nm,n),alfr(n),alfi(n),beta(n),z(nm,n)
      logical tf
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     for the real general generalized eigenproblem  ax = (lambda)bx.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrices  a  and  b.
!
!        a  contains a real general matrix.
!
!        b  contains a real general matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        alfr  and  alfi  contain the real and imaginary parts,
!        respectively, of the numerators of the eigenvalues.
!
!        beta  contains the denominators of the eigenvalues,
!        which are thus given by the ratios  (alfr+i*alfi)/beta.
!        complex conjugate pairs of eigenvalues appear consecutively
!        with the eigenvalue having the positive imaginary part first.
!
!        z  contains the real and imaginary parts of the eigenvectors
!        if matz is not zero.  if the j-th eigenvalue is real, the
!        j-th column of  z  contains its eigenvector.  if the j-th
!        eigenvalue is complex with positive imaginary part, the
!        j-th and (j+1)-th columns of  z  contain the real and
!        imaginary parts of its eigenvector.  the conjugate of this
!        vector is the eigenvector for the conjugate eigenvalue.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for qzit.
!           the normal completion code is zero.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
!
   10 if (matz .ne. 0) go to 20
!     .......... find eigenvalues only ..........
      tf = .false.
      call  qzhes(nm,n,a,b,tf,z)
      call  qzit(nm,n,a,b,0.0e0_rk,tf,z,ierr)
      call  qzval(nm,n,a,b,alfr,alfi,beta,tf,z)
      go to 50
!     .......... find both eigenvalues and eigenvectors ..........
   20 tf = .true.
      call  qzhes(nm,n,a,b,tf,z)
      call  qzit(nm,n,a,b,0.0e0_rk,tf,z,ierr)
      call  qzval(nm,n,a,b,alfr,alfi,beta,tf,z)
      if (ierr .ne. 0) go to 50
      call  qzvec(nm,n,a,b,alfr,alfi,beta,z)
   50 return
      end subroutine
!*********************************************************************************
      subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
!
      integer n,nm,ierr,matz
      real(kind=rk)  :: a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a real symmetric matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        a  contains the real symmetric matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and tql2.  the normal completion code is zero.
!
!        fv1  and  fv2  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
!
   10 if (matz .ne. 0) go to 20
!     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fv1,fv2)
      call  tqlrat(n,w,fv2,ierr)
      go to 50
!     .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
   50 return
      end subroutine
!*********************************************************************************
      subroutine rsb(nm,n,mb,a,w,matz,z,fv1,fv2,ierr)
!
      integer n,mb,nm,ierr,matz
      real(kind=rk)  :: a(nm,mb),w(n),z(nm,n),fv1(n),fv2(n)
      logical tf
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a real symmetric band matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        mb  is the half band width of the matrix, defined as the
!        number of adjacent diagonals, including the principal
!        diagonal, required to specify the non-zero portion of the
!        lower triangle of the matrix.
!
!        a  contains the lower triangle of the real symmetric
!        band matrix.  its lowest subdiagonal is stored in the
!        last  n+1-mb  positions of the first column, its next
!        subdiagonal in the last  n+2-mb  positions of the
!        second column, further subdiagonals similarly, and
!        finally its principal diagonal in the  n  positions
!        of the last column.  contents of storages not part
!        of the matrix are arbitrary.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and tql2.  the normal completion code is zero.
!
!        fv1  and  fv2  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (n .le. nm) go to 5
      ierr = 10 * n
      go to 50
    5 if (mb .gt. 0) go to 10
      ierr = 12 * n
      go to 50
   10 if (mb .le. n) go to 15
      ierr = 12 * n
      go to 50
!
   15 if (matz .ne. 0) go to 20
!     .......... find eigenvalues only ..........
      tf = .false.
      call  bandr(nm,n,mb,a,w,fv1,fv2,tf,z)
      call  tqlrat(n,w,fv2,ierr)
      go to 50
!     .......... find both eigenvalues and eigenvectors ..........
   20 tf = .true.
      call  bandr(nm,n,mb,a,w,fv1,fv1,tf,z)
      call  tql2(nm,n,w,fv1,z,ierr)
   50 return
      end subroutine
!*********************************************************************************
      subroutine rsg(nm,n,a,b,w,matz,z,fv1,fv2,ierr)
!
      integer n,nm,ierr,matz
      real(kind=rk)  :: a(nm,n),b(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     for the real symmetric generalized eigenproblem  ax = (lambda)bx.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrices  a  and  b.
!
!        a  contains a real symmetric matrix.
!
!        b  contains a positive definite real symmetric matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and tql2.  the normal completion code is zero.
!
!        fv1  and  fv2  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
!
   10 call  reduc(nm,n,a,b,fv2,ierr)
      if (ierr .ne. 0) go to 50
      if (matz .ne. 0) go to 20
!     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fv1,fv2)
      call  tqlrat(n,w,fv2,ierr)
      go to 50
!     .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
      if (ierr .ne. 0) go to 50
      call  rebak(nm,n,b,fv2,n,z)
   50 return
      end subroutine
!*********************************************************************************
      subroutine rsgab(nm,n,a,b,w,matz,z,fv1,fv2,ierr)
!
      integer n,nm,ierr,matz
      real(kind=rk)  :: a(nm,n),b(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     for the real symmetric generalized eigenproblem  abx = (lambda)x.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrices  a  and  b.
!
!        a  contains a real symmetric matrix.
!
!        b  contains a positive definite real symmetric matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and tql2.  the normal completion code is zero.
!
!        fv1  and  fv2  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
!
   10 call  reduc2(nm,n,a,b,fv2,ierr)
      if (ierr .ne. 0) go to 50
      if (matz .ne. 0) go to 20
!     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fv1,fv2)
      call  tqlrat(n,w,fv2,ierr)
      go to 50
!     .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
      if (ierr .ne. 0) go to 50
      call  rebak(nm,n,b,fv2,n,z)
   50 return
      end subroutine
!*********************************************************************************
      subroutine rsgba(nm,n,a,b,w,matz,z,fv1,fv2,ierr)
!
      integer n,nm,ierr,matz
      real(kind=rk)  :: a(nm,n),b(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     for the real symmetric generalized eigenproblem  bax = (lambda)x.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrices  a  and  b.
!
!        a  contains a real symmetric matrix.
!
!        b  contains a positive definite real symmetric matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and tql2.  the normal completion code is zero.
!
!        fv1  and  fv2  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
!
   10 call  reduc2(nm,n,a,b,fv2,ierr)
      if (ierr .ne. 0) go to 50
      if (matz .ne. 0) go to 20
!     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fv1,fv2)
      call  tqlrat(n,w,fv2,ierr)
      go to 50
!     .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
      if (ierr .ne. 0) go to 50
      call  rebakb(nm,n,b,fv2,n,z)
   50 return
      end subroutine
!*********************************************************************************
      subroutine rsm(nm,n,a,w,m,z,fwork,iwork,ierr)
!
      integer n,nm,m,iwork(n),ierr
      real(kind=rk)  :: a(nm,n),w(n),z(nm,m),fwork(1)
      integer   :: k1,k2,k3,k4,k5,k6,k7,k8
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find all of the eigenvalues and some of the eigenvectors
!     of a real symmetric matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        a  contains the real symmetric matrix.
!
!        m  the eigenvectors corresponding to the first m eigenvalues
!           are to be computed.
!           if m = 0 then no eigenvectors are computed.
!           if m = n then all of the eigenvectors are computed.
!
!     on output
!
!        w  contains all n eigenvalues in ascending order.
!
!        z  contains the orthonormal eigenvectors associated with
!           the first m eigenvalues.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat,
!           imtqlv and tinvit.  the normal completion code is zero.
!
!        fwork  is a temporary storage array of dimension 8*n.
!
!        iwork  is an integer temporary storage array of dimension n.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 10 * n
      if (n .gt. nm .or. m .gt. nm) go to 50
      k1 = 1
      k2 = k1 + n
      k3 = k2 + n
      k4 = k3 + n
      k5 = k4 + n
      k6 = k5 + n
      k7 = k6 + n
      k8 = k7 + n
      if (m .gt. 0) go to 10
!     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fwork(k1),fwork(k2))
      call  tqlrat(n,w,fwork(k2),ierr)
      go to 50
!     .......... find all eigenvalues and m eigenvectors ..........
   10 call  tred1(nm,n,a,fwork(k1),fwork(k2),fwork(k3))
      call  imtqlv(n,fwork(k1),fwork(k2),fwork(k3),w,iwork, ierr,fwork(k4))
      call  tinvit(nm,n,fwork(k1),fwork(k2),fwork(k3),m,w,iwork,z,ierr, fwork(k4),fwork(k5),fwork(k6),fwork(k7),fwork(k8))
      call  trbak1(nm,n,a,fwork(k2),m,z)
   50 return
      end subroutine
!*********************************************************************************
      subroutine rsp(nm,n,nv,a,w,matz,z,fv1,fv2,ierr)
!
      integer i,j,n,nm,nv,ierr,matz
      real(kind=rk)  :: a(nv),w(n),z(nm,n),fv1(n),fv2(n)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a real symmetric packed matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        nv  is an integer variable set equal to the
!        dimension of the array  a  as specified for
!        a  in the calling program.  nv  must not be
!        less than  n*(n+1)/2.
!
!        a  contains the lower triangle of the real symmetric
!        packed matrix stored row-wise.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and tql2.  the normal completion code is zero.
!
!        fv1  and  fv2  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (n .le. nm) go to 5
      ierr = 10 * n
      go to 50
    5 if (nv .ge. (n * (n + 1)) / 2) go to 10
      ierr = 20 * n
      go to 50
!
   10 call  tred3(n,nv,a,w,fv1,fv2)
      if (matz .ne. 0) go to 20
!     .......... find eigenvalues only ..........
      call  tqlrat(n,w,fv2,ierr)
      go to 50
!     .......... find both eigenvalues and eigenvectors ..........
   20 do 40 i = 1, n
!
         do 30 j = 1, n
            z(j,i) = 0.0e0_rk
   30    continue
!
         z(i,i) = 1.0e0_rk
   40 continue
!
      call  tql2(nm,n,w,fv1,z,ierr)
      if (ierr .ne. 0) go to 50
      call  trbak3(nm,n,nv,a,n,z)
   50 return
      end subroutine
!*********************************************************************************
      subroutine rst(nm,n,w,e,matz,z,ierr)
!
      integer i,j,n,nm,ierr,matz
      real(kind=rk)  :: w(n),e(n),z(nm,n)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a real symmetric tridiagonal matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix.
!
!        w  contains the diagonal elements of the real
!        symmetric tridiagonal matrix.
!
!        e  contains the subdiagonal elements of the matrix in
!        its last n-1 positions.  e(1) is arbitrary.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for imtql1
!           and imtql2.  the normal completion code is zero.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
!
   10 if (matz .ne. 0) go to 20
!     .......... find eigenvalues only ..........
      call  imtql1(n,w,e,ierr)
      go to 50
!     .......... find both eigenvalues and eigenvectors ..........
   20 do 40 i = 1, n
!
         do 30 j = 1, n
            z(j,i) = 0.0e0_rk
   30    continue
!
         z(i,i) = 1.0e0_rk
   40 continue
!
      call  imtql2(nm,n,w,e,z,ierr)
   50 return
      end subroutine
!*********************************************************************************
      subroutine rt(nm,n,a,w,matz,z,fv1,ierr)
!
      integer n,nm,ierr,matz
      real(kind=rk)  :: a(nm,3),w(n),z(nm,n),fv1(n)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a special real tridiagonal matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        a  contains the special real tridiagonal matrix in its
!        first three columns.  the subdiagonal elements are stored
!        in the last  n-1  positions of the first column, the
!        diagonal elements in the second column, and the superdiagonal
!        elements in the first  n-1  positions of the third column.
!        elements  a(1,1)  and  a(n,3)  are arbitrary.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for imtql1
!           and imtql2.  the normal completion code is zero.
!
!        fv1  is a temporary storage array.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
!
   10 if (matz .ne. 0) go to 20
!     .......... find eigenvalues only ..........
      call  figi(nm,n,a,w,fv1,fv1,ierr)
      if (ierr .gt. 0) go to 50
      call  imtql1(n,w,fv1,ierr)
      go to 50
!     .......... find both eigenvalues and eigenvectors ..........
   20 call  figi2(nm,n,a,w,fv1,z,ierr)
      if (ierr .ne. 0) go to 50
      call  imtql2(nm,n,w,fv1,z,ierr)
   50 return
      end subroutine
!*********************************************************************************
      subroutine svd(nm,m,n,a,w,matu,u,matv,v,ierr,rv1)
!
      integer i,j,k,l,m,n,ii,i1,kk,k1,ll,l1,mn,nm,its,ierr
      real(kind=rk)  :: a(nm,n),w(n),u(nm,n),v(nm,n),rv1(n)
      real(kind=rk)  :: c,f,g,h,s,x,y,z,tst1,tst2,scale
      logical matu,matv
!
!     this subroutine is a translation of the algol procedure svd,
!     num. math. 14, 403-420(1970) by golub and reinsch.
!     handbook for auto. comp., vol ii-linear algebra, 134-151(1971).
!
!     this subroutine determines the singular value decomposition
!          t
!     a=usv  of a real m by n rectangular matrix.  householder
!     bidiagonalization and a variant of the qr algorithm are used.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.  note that nm must be at least
!          as large as the maximum of m and n.
!
!        m is the number of rows of a (and u).
!
!        n is the number of columns of a (and u) and the order of v.
!
!        a contains the rectangular input matrix to be decomposed.
!
!        matu should be set to .true. if the u matrix in the
!          decomposition is desired, and to .false. otherwise.
!
!        matv should be set to .true. if the v matrix in the
!          decomposition is desired, and to .false. otherwise.
!
!     on output
!
!        a is unaltered (unless overwritten by u or v).
!
!        w contains the n (non-negative) singular values of a (the
!          diagonal elements of s).  they are unordered.  if an
!          error exit is made, the singular values should be correct
!          for indices ierr+1,ierr+2,...,n.
!
!        u contains the matrix u (orthogonal column vectors) of the
!          decomposition if matu has been set to .true.  otherwise
!          u is used as a temporary array.  u may coincide with a.
!          if an error exit is made, the columns of u corresponding
!          to indices of correct singular values should be correct.
!
!        v contains the matrix v (orthogonal) of the decomposition if
!          matv has been set to .true.  otherwise v is not referenced.
!          v may also coincide with a if u is not needed.  if an error
!          exit is made, the columns of v corresponding to indices of
!          correct singular values should be correct.
!
!        ierr is set to
!          zero       for normal return,
!          k          if the k-th singular value has not been
!                     determined after 30 iterations.
!
!        rv1 is a temporary storage array.
!
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
!
      do 100 i = 1, m
!
         do 100 j = 1, n
            u(i,j) = a(i,j)
  100 continue
!     .......... householder reduction to bidiagonal form ..........
      g = 0.0e0_rk
      scale = 0.0e0_rk
      x = 0.0e0_rk
!
      do 300 i = 1, n
         l = i + 1
         rv1(i) = scale * g
         g = 0.0e0_rk
         s = 0.0e0_rk
         scale = 0.0e0_rk
         if (i .gt. m) go to 210
!
         do 120 k = i, m
  120    scale = scale + abs(u(k,i))
!
         if (scale .eq. 0.0e0_rk) go to 210
!
         do 130 k = i, m
            u(k,i) = u(k,i) / scale
            s = s + u(k,i)**2
  130    continue
!
         f = u(i,i)
         g = -sign(sqrt(s),f)
         h = f * g - s
         u(i,i) = f - g
         if (i .eq. n) go to 190
!
         do 150 j = l, n
            s = 0.0e0_rk
!
            do 140 k = i, m
  140       s = s + u(k,i) * u(k,j)
!
            f = s / h
!
            do 150 k = i, m
               u(k,j) = u(k,j) + f * u(k,i)
  150    continue
!
  190    do 200 k = i, m
  200    u(k,i) = scale * u(k,i)
!
  210    w(i) = scale * g
         g = 0.0e0_rk
         s = 0.0e0_rk
         scale = 0.0e0_rk
         if (i .gt. m .or. i .eq. n) go to 290
!
         do 220 k = l, n
  220    scale = scale + abs(u(i,k))
!
         if (scale .eq. 0.0e0_rk) go to 290
!
         do 230 k = l, n
            u(i,k) = u(i,k) / scale
            s = s + u(i,k)**2
  230    continue
!
         f = u(i,l)
         g = -sign(sqrt(s),f)
         h = f * g - s
         u(i,l) = f - g
!
         do 240 k = l, n
  240    rv1(k) = u(i,k) / h
!
         if (i .eq. m) go to 270
!
         do 260 j = l, m
            s = 0.0e0_rk
!
            do 250 k = l, n
  250       s = s + u(j,k) * u(i,k)
!
            do 260 k = l, n
               u(j,k) = u(j,k) + s * rv1(k)
  260    continue
!
  270    do 280 k = l, n
  280    u(i,k) = scale * u(i,k)
!
  290    x = max(x,abs(w(i))+abs(rv1(i)))
  300 continue
!     .......... accumulation of right-hand transformations ..........
      if (.not. matv) go to 410
!     .......... for i=n step -1 until 1 do -- ..........
      do 400 ii = 1, n
         i = n + 1 - ii
         if (i .eq. n) go to 390
         if (g .eq. 0.0e0_rk) go to 360
!
         do 320 j = l, n
!     .......... double division avoids possible underflow ..........
  320    v(j,i) = (u(i,j) / u(i,l)) / g
!
         do 350 j = l, n
            s = 0.0e0_rk
!
            do 340 k = l, n
  340       s = s + u(i,k) * v(k,j)
!
            do 350 k = l, n
               v(k,j) = v(k,j) + s * v(k,i)
  350    continue
!
  360    do 380 j = l, n
            v(i,j) = 0.0e0_rk
            v(j,i) = 0.0e0_rk
  380    continue
!
  390    v(i,i) = 1.0e0_rk
         g = rv1(i)
         l = i
  400 continue
!     .......... accumulation of left-hand transformations ..........
  410 if (.not. matu) go to 510
!     ..........for i=min(m,n) step -1 until 1 do -- ..........
      mn = n
      if (m .lt. n) mn = m
!
      do 500 ii = 1, mn
         i = mn + 1 - ii
         l = i + 1
         g = w(i)
         if (i .eq. n) go to 430
!
         do 420 j = l, n
  420    u(i,j) = 0.0e0_rk
!
  430    if (g .eq. 0.0e0_rk) go to 475
         if (i .eq. mn) go to 460
!
         do 450 j = l, n
            s = 0.0e0_rk
!
            do 440 k = l, m
  440       s = s + u(k,i) * u(k,j)
!     .......... double division avoids possible underflow ..........
            f = (s / u(i,i)) / g
!
            do 450 k = i, m
               u(k,j) = u(k,j) + f * u(k,i)
  450    continue
!
  460    do 470 j = i, m
  470    u(j,i) = u(j,i) / g
!
         go to 490
!
  475    do 480 j = i, m
  480    u(j,i) = 0.0e0_rk
!
  490    u(i,i) = u(i,i) + 1.0e0_rk
  500 continue
!     .......... diagonalization of the bidiagonal form ..........
  510 tst1 = x
!     .......... for k=n step -1 until 1 do -- ..........
      do 700 kk = 1, n
         k1 = n - kk
         k = k1 + 1
         its = 0
!     .......... test for splitting.
!                for l=k step -1 until 1 do -- ..........
  520    do 530 ll = 1, k
            l1 = k - ll
            l = l1 + 1
            tst2 = tst1 + abs(rv1(l))
            if (tst2 .eq. tst1) go to 565
!     .......... rv1(1) is always zero, so there is no exit
!                through the bottom of the loop ..........
            tst2 = tst1 + abs(w(l1))
            if (tst2 .eq. tst1) go to 540
  530    continue
!     .......... cancellation of rv1(l) if l greater than 1 ..........
  540    c = 0.0e0_rk
         s = 1.0e0_rk
!
         do 560 i = l, k
            f = s * rv1(i)
            rv1(i) = c * rv1(i)
            tst2 = tst1 + abs(f)
            if (tst2 .eq. tst1) go to 565
            g = w(i)
            h = pythag(f,g)
            w(i) = h
            c = g / h
            s = -f / h
            if (.not. matu) go to 560
!
            do 550 j = 1, m
               y = u(j,l1)
               z = u(j,i)
               u(j,l1) = y * c + z * s
               u(j,i) = -y * s + z * c
  550       continue
!
  560    continue
!     .......... test for convergence ..........
  565    z = w(k)
         if (l .eq. k) go to 650
!     .......... shift from bottom 2 by 2 minor ..........
         if (its .eq. 30) go to 1000
         its = its + 1
         x = w(l)
         y = w(k1)
         g = rv1(k1)
         h = rv1(k)
         f = 0.5e0_rk * (((g + z) / h) * ((g - z) / y) + y / h - h / y)
         g = pythag(f,1.0e0_rk)
         f = x - (z / x) * z + (h / x) * (y / (f + sign(g,f)) - h)
!     .......... next qr transformation ..........
         c = 1.0e0_rk
         s = 1.0e0_rk
!
         do 600 i1 = l, k1
            i = i1 + 1
            g = rv1(i)
            y = w(i)
            h = s * g
            g = c * g
            z = pythag(f,h)
            rv1(i1) = z
            c = f / z
            s = h / z
            f = x * c + g * s
            g = -x * s + g * c
            h = y * s
            y = y * c
            if (.not. matv) go to 575
!
            do 570 j = 1, n
               x = v(j,i1)
               z = v(j,i)
               v(j,i1) = x * c + z * s
               v(j,i) = -x * s + z * c
  570       continue
!
  575       z = pythag(f,h)
            w(i1) = z
!     .......... rotation can be arbitrary if z is zero ..........
            if (z .eq. 0.0e0_rk) go to 580
            c = f / z
            s = h / z
  580       f = c * g + s * y
            x = -s * g + c * y
            if (.not. matu) go to 600
!
            do 590 j = 1, m
               y = u(j,i1)
               z = u(j,i)
               u(j,i1) = y * c + z * s
               u(j,i) = -y * s + z * c
  590       continue
!
  600    continue
!
         rv1(l) = 0.0e0_rk
         rv1(k) = f
         w(k) = x
         go to 520
!     .......... convergence ..........
  650    if (z .ge. 0.0e0_rk) go to 700
!     .......... w(k) is made non-negative ..........
         w(k) = -z
         if (.not. matv) go to 700
!
         do 690 j = 1, n
  690    v(j,k) = -v(j,k)
!
  700 continue
!
      go to 1001
!     .......... set error -- no convergence to a
!                singular value after 30 iterations ..........
 1000 ierr = k
 1001 return
      end subroutine
!*********************************************************************************
      subroutine tinvit(nm,n,d,e,e2,m,w,ind,z, ierr,rv1,rv2,rv3,rv4,rv6)
!
      integer i,j,m,n,p,q,r,s,ii,ip,jj,nm,its,tag,ierr,group
      real(kind=rk)  :: d(n),e(n),e2(n),w(m),z(nm,m),rv1(n),rv2(n),rv3(n),rv4(n),rv6(n)
      real(kind=rk)  :: u,v,uk,xu,x0,x1,eps2,eps3,eps4,norm,order
      integer ind(m)
!
!     this subroutine is a translation of the inverse iteration tech-
!     nique in the algol procedure tristurm by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
!
!     this subroutine finds those eigenvectors of a tridiagonal
!     symmetric matrix corresponding to specified eigenvalues,
!     using inverse iteration.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        e2 contains the squares of the corresponding elements of e,
!          with zeros corresponding to negligible elements of e.
!          e(i) is considered negligible if it is not larger than
!          the product of the relative machine precision and the sum
!          of the magnitudes of d(i) and d(i-1).  e2(1) must contain
!          0.0e0_rk if the eigenvalues are in ascending order, or 2.0e0_rk
!          if the eigenvalues are in descending order.  if  bisect,
!          tridib, or  imtqlv  has been used to find the eigenvalues,
!          their output e2 array is exactly what is expected here.
!
!        m is the number of specified eigenvalues.
!
!        w contains the m eigenvalues in ascending or descending order.
!
!        ind contains in its first m positions the submatrix indices
!          associated with the corresponding eigenvalues in w --
!          1 for eigenvalues belonging to the first submatrix from
!          the top, 2 for those belonging to the second submatrix, etc.
!
!     on output
!
!        all input arrays are unaltered.
!
!        z contains the associated set of orthonormal eigenvectors.
!          any vector which fails to converge is set to zero.
!
!        ierr is set to
!          zero       for normal return,
!          -r         if the eigenvector corresponding to the r-th
!                     eigenvalue fails to converge in 5 iterations.
!
!        rv1, rv2, rv3, rv4, and rv6 are temporary storage arrays.
!
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      if (m .eq. 0) go to 1001
      tag = 0
      order = 1.0e0_rk - e2(1)
      q = 0
!     .......... establish and process next submatrix ..........
  100 p = q + 1
!
      do 120 q = p, n
         if (q .eq. n) go to 140
         if (e2(q+1) .eq. 0.0e0_rk) go to 140
  120 continue
!     .......... find vectors by inverse iteration ..........
  140 tag = tag + 1
      s = 0
!
      do 920 r = 1, m
         if (ind(r) .ne. tag) go to 920
         its = 1
         x1 = w(r)
         if (s .ne. 0) go to 510
!     .......... check for isolated root ..........
         xu = 1.0e0_rk
         if (p .ne. q) go to 490
         rv6(p) = 1.0e0_rk
         go to 870
  490    norm = abs(d(p))
         ip = p + 1
!
         do 500 i = ip, q
  500    norm = max(norm, abs(d(i))+abs(e(i)))
!     .......... eps2 is the criterion for grouping,
!                eps3 replaces zero pivots and equal
!                roots are modified by eps3,
!                eps4 is taken very small to avoid overflow ..........
         eps2 = 1.0e-3_rk * norm
         eps3 = epslon(norm)
         uk = q - p + 1
         eps4 = uk * eps3
         uk = eps4 / sqrt(uk)
         s = p
  505    group = 0
         go to 520
!     .......... look for close or coincident roots ..........
  510    if (abs(x1-x0) .ge. eps2) go to 505
         group = group + 1
         if (order * (x1 - x0) .le. 0.0e0_rk) x1 = x0 + order * eps3
!     .......... elimination with interchanges and
!                initialization of vector ..........
  520    v = 0.0e0_rk
!
         do 580 i = p, q
            rv6(i) = uk
            if (i .eq. p) go to 560
            if (abs(e(i)) .lt. abs(u)) go to 540
!     .......... warning -- a divide check may occur here if
!                e2 array has not been specified correctly ..........
            xu = u / e(i)
            rv4(i) = xu
            rv1(i-1) = e(i)
            rv2(i-1) = d(i) - x1
            rv3(i-1) = 0.0e0_rk
            if (i .ne. q) rv3(i-1) = e(i+1)
            u = v - xu * rv2(i-1)
            v = -xu * rv3(i-1)
            go to 580
  540       xu = e(i) / u
            rv4(i) = xu
            rv1(i-1) = u
            rv2(i-1) = v
            rv3(i-1) = 0.0e0_rk
  560       u = d(i) - x1 - xu * v
            if (i .ne. q) v = e(i+1)
  580    continue
!
         if (u .eq. 0.0e0_rk) u = eps3
         rv1(q) = u
         rv2(q) = 0.0e0_rk
         rv3(q) = 0.0e0_rk
!     .......... back substitution
!                for i=q step -1 until p do -- ..........
  600    do 620 ii = p, q
            i = p + q - ii
            rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i)
            v = u
            u = rv6(i)
  620    continue
!     .......... orthogonalize with respect to previous
!                members of group ..........
         if (group .eq. 0) go to 700
         j = r
!
         do 680 jj = 1, group
  630       j = j - 1
            if (ind(j) .ne. tag) go to 630
            xu = 0.0e0_rk
!
            do 640 i = p, q
  640       xu = xu + rv6(i) * z(i,j)
!
            do 660 i = p, q
  660       rv6(i) = rv6(i) - xu * z(i,j)
!
  680    continue
!
  700    norm = 0.0e0_rk
!
         do 720 i = p, q
  720    norm = norm + abs(rv6(i))
!
         if (norm .ge. 1.0e0_rk) go to 840
!     .......... forward substitution ..........
         if (its .eq. 5) go to 830
         if (norm .ne. 0.0e0_rk) go to 740
         rv6(s) = eps4
         s = s + 1
         if (s .gt. q) s = p
         go to 780
  740    xu = eps4 / norm
!
         do 760 i = p, q
  760    rv6(i) = rv6(i) * xu
!     .......... elimination operations on next vector
!                iterate ..........
  780    do 820 i = ip, q
            u = rv6(i)
!     .......... if rv1(i-1) .eq. e(i), a row interchange
!                was performed earlier in the
!                triangularization process ..........
            if (rv1(i-1) .ne. e(i)) go to 800
            u = rv6(i-1)
            rv6(i-1) = rv6(i)
  800       rv6(i) = u - rv4(i) * rv6(i-1)
  820    continue
!
         its = its + 1
         go to 600
!     .......... set error -- non-converged eigenvector ..........
  830    ierr = -r
         xu = 0.0e0_rk
         go to 870
!     .......... normalize so that sum of squares is
!                1 and expand to full order ..........
  840    u = 0.0e0_rk
!
         do 860 i = p, q
  860    u = pythag(u,rv6(i))
!
         xu = 1.0e0_rk / u
!
  870    do 880 i = 1, n
  880    z(i,r) = 0.0e0_rk
!
         do 900 i = p, q
  900    z(i,r) = rv6(i) * xu
!
         x0 = x1
  920 continue
!
      if (q .lt. n) go to 100
 1001 return
      end subroutine
!*********************************************************************************
      subroutine tql1(n,d,e,ierr)
!
      integer i,j,l,m,n,ii,l1,l2,mml,ierr
      real(kind=rk)  :: d(n),e(n)
      real(kind=rk)  :: c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2
!
!     this subroutine is a translation of the algol procedure tql1,
!     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!     wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
!     this subroutine finds the eigenvalues of a symmetric
!     tridiagonal matrix by the ql method.
!
!     on input
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct and
!          ordered for indices 1,2,...ierr-1, but may not be
!          the smallest eigenvalues.
!
!        e has been destroyed.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      if (n .eq. 1) go to 1001
!
      do 100 i = 2, n
  100 e(i-1) = e(i)
!
      f = 0.0e0_rk
      tst1 = 0.0e0_rk
      e(n) = 0.0e0_rk
!
      do 290 l = 1, n
         j = 0
         h = abs(d(l)) + abs(e(l))
         if (tst1 .lt. h) tst1 = h
!     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + abs(e(m))
            if (tst2 .eq. tst1) go to 120
!     .......... e(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
  110    continue
!
  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
!     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0e0_rk * e(l))
         r = pythag(p,1.0e0_rk)
         d(l) = e(l) / (p + sign(r,p))
         d(l1) = e(l) * (p + sign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
!
         do 140 i = l2, n
  140    d(i) = d(i) - h
!
  145    f = f + h
!     .......... ql transformation ..........
         p = d(m)
         c = 1.0e0_rk
         c2 = c
         el1 = e(l1)
         s = 0.0e0_rk
         mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
  200    continue
!
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + abs(e(l))
         if (tst2 .gt. tst1) go to 130
  210    p = d(l) + f
!     .......... order eigenvalues ..........
         if (l .eq. 1) go to 250
!     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
!
  250    i = 1
  270    d(i) = p
  290 continue
!
      go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end subroutine
!*********************************************************************************
      subroutine tql2(nm,n,d,e,z,ierr)
!
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      real(kind=rk)  :: d(n),e(n),z(nm,n)
      real(kind=rk)  :: c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2
!
!     this subroutine is a translation of the algol procedure tql2,
!     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!     wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a symmetric tridiagonal matrix by the ql method.
!     the eigenvectors of a full symmetric matrix can also
!     be found if  tred2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        z contains the transformation matrix produced in the
!          reduction by  tred2, if performed.  if the eigenvectors
!          of the tridiagonal matrix are desired, z must contain
!          the identity matrix.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,ierr-1.
!
!        e has been destroyed.
!
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      if (n .eq. 1) go to 1001
!
      do 100 i = 2, n
  100 e(i-1) = e(i)
!
      f = 0.0e0_rk
      tst1 = 0.0e0_rk
      e(n) = 0.0e0_rk
!
      do 240 l = 1, n
         j = 0
         h = abs(d(l)) + abs(e(l))
         if (tst1 .lt. h) tst1 = h
!     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + abs(e(m))
            if (tst2 .eq. tst1) go to 120
!     .......... e(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
  110    continue
!
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
!     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0e0_rk * e(l))
         r = pythag(p,1.0e0_rk)
         d(l) = e(l) / (p + sign(r,p))
         d(l1) = e(l) * (p + sign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
!
         do 140 i = l2, n
  140    d(i) = d(i) - h
!
  145    f = f + h
!     .......... ql transformation ..........
         p = d(m)
         c = 1.0e0_rk
         c2 = c
         el1 = e(l1)
         s = 0.0e0_rk
         mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
!     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
!
  200    continue
!
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + abs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
!     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
!
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
!
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
!
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
!
  300 continue
!
      go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end subroutine
!*********************************************************************************
      subroutine tqlrat(n,d,e2,ierr)
!
      integer i,j,l,m,n,ii,l1,mml,ierr
      real(kind=rk)  :: d(n),e2(n)
      real(kind=rk)  :: b,c,f,g,h,p,r,s,t
!
!     this subroutine is a translation of the algol procedure tqlrat,
!     algorithm 464, comm. acm 16, 689(1973) by reinsch.
!
!     this subroutine finds the eigenvalues of a symmetric
!     tridiagonal matrix by the rational ql method.
!
!     on input
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e2 contains the squares of the subdiagonal elements of the
!          input matrix in its last n-1 positions.  e2(1) is arbitrary.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct and
!          ordered for indices 1,2,...ierr-1, but may not be
!          the smallest eigenvalues.
!
!        e2 has been destroyed.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      if (n .eq. 1) go to 1001
!
      do 100 i = 2, n
  100 e2(i-1) = e2(i)
!
      f = 0.0e0_rk
      t = 0.0e0_rk
      e2(n) = 0.0e0_rk
!
      do 290 l = 1, n
         j = 0
         h = abs(d(l)) + sqrt(e2(l))
         if (t .gt. h) go to 105
         t = h
         b = epslon(t)
         c = b * b
!     .......... look for small squared sub-diagonal element ..........
  105    do 110 m = l, n
            if (e2(m) .le. c) go to 120
!     .......... e2(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
  110    continue
!
  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
!     .......... form shift ..........
         l1 = l + 1
         s = sqrt(e2(l))
         g = d(l)
         p = (d(l1) - g) / (2.0e0_rk * s)
         r = pythag(p,1.0e0_rk)
         d(l) = s / (p + sign(r,p))
         h = g - d(l)
!
         do 140 i = l1, n
  140    d(i) = d(i) - h
!
         f = f + h
!     .......... rational ql transformation ..........
         g = d(m)
         if (g .eq. 0.0e0_rk) g = b
         h = g
         s = 0.0e0_rk
         mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            p = g * h
            r = p + e2(i)
            e2(i+1) = s * r
            s = e2(i) / r
            d(i+1) = h + s * (h + d(i))
            g = d(i) - e2(i) / g
            if (g .eq. 0.0e0_rk) g = b
            h = g * p / r
  200    continue
!
         e2(l) = s * g
         d(l) = h
!     .......... guard against underflow in convergence test ..........
         if (h .eq. 0.0e0_rk) go to 210
         if (abs(e2(l)) .le. abs(c/h)) go to 210
         e2(l) = h * e2(l)
         if (e2(l) .ne. 0.0e0_rk) go to 130
  210    p = d(l) + f
!     .......... order eigenvalues ..........
         if (l .eq. 1) go to 250
!     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
!
  250    i = 1
  270    d(i) = p
  290 continue
!
      go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end subroutine
!*********************************************************************************
      subroutine trbak1(nm,n,a,e,m,z)
!
      integer i,j,k,l,m,n,nm
      real(kind=rk)  :: a(nm,n),e(n),z(nm,m)
      real(kind=rk)  :: s
!
!     this subroutine is a translation of the algol procedure trbak1,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine forms the eigenvectors of a real symmetric
!     matrix by back transforming those of the corresponding
!     symmetric tridiagonal matrix determined by  tred1.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains information about the orthogonal trans-
!          formations used in the reduction by  tred1
!          in its strict lower triangle.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is arbitrary.
!
!        m is the number of eigenvectors to be back transformed.
!
!        z contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        z contains the transformed eigenvectors
!          in its first m columns.
!
!     note that trbak1 preserves vector euclidean norms.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (m .eq. 0) go to 200
      if (n .eq. 1) go to 200
!
      do 140 i = 2, n
         l = i - 1
         if (e(i) .eq. 0.0e0_rk) go to 140
!
         do 130 j = 1, m
            s = 0.0e0_rk
!
            do 110 k = 1, l
  110       s = s + a(i,k) * z(k,j)
!     .......... divisor below is negative of h formed in tred1.
!                double division avoids possible underflow ..........
            s = (s / a(i,l)) / e(i)
!
            do 120 k = 1, l
  120       z(k,j) = z(k,j) + s * a(i,k)
!
  130    continue
!
  140 continue
!
  200 return
      end subroutine
!*********************************************************************************
      subroutine trbak3(nm,n,nv,a,m,z)
!
      integer i,j,k,l,m,n,ik,iz,nm,nv
      real(kind=rk)  :: a(nv),z(nm,m)
      real(kind=rk)  :: h,s
!
!     this subroutine is a translation of the algol procedure trbak3,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine forms the eigenvectors of a real symmetric
!     matrix by back transforming those of the corresponding
!     symmetric tridiagonal matrix determined by  tred3.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        nv must be set to the dimension of the array parameter a
!          as declared in the calling program dimension statement.
!
!        a contains information about the orthogonal transformations
!          used in the reduction by  tred3  in its first
!          n*(n+1)/2 positions.
!
!        m is the number of eigenvectors to be back transformed.
!
!        z contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        z contains the transformed eigenvectors
!          in its first m columns.
!
!     note that trbak3 preserves vector euclidean norms.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (m .eq. 0) go to 200
      if (n .eq. 1) go to 200
!
      do 140 i = 2, n
         l = i - 1
         iz = (i * l) / 2
         ik = iz + i
         h = a(ik)
         if (h .eq. 0.0e0_rk) go to 140
!
         do 130 j = 1, m
            s = 0.0e0_rk
            ik = iz
!
            do 110 k = 1, l
               ik = ik + 1
               s = s + a(ik) * z(k,j)
  110       continue
!     .......... double division avoids possible underflow ..........
            s = (s / h) / h
            ik = iz
!
            do 120 k = 1, l
               ik = ik + 1
               z(k,j) = z(k,j) - s * a(ik)
  120       continue
!
  130    continue
!
  140 continue
!
  200 return
      end subroutine
!*********************************************************************************
      subroutine tred1(nm,n,a,d,e,e2)
!
      integer i,j,k,l,n,ii,nm,jp1
      real(kind=rk)  :: a(nm,n),d(n),e(n),e2(n)
      real(kind=rk)  :: f,g,h,scale
!
!     this subroutine is a translation of the algol procedure tred1,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix
!     to a symmetric tridiagonal matrix using
!     orthogonal similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the real symmetric input matrix.  only the
!          lower triangle of the matrix need be supplied.
!
!     on output
!
!        a contains information about the orthogonal trans-
!          formations used in the reduction in its strict lower
!          triangle.  the full upper triangle of a is unaltered.
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      do 100 i = 1, n
         d(i) = a(n,i)
         a(n,i) = a(i,i)
  100 continue
!     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0e0_rk
         scale = 0.0e0_rk
         if (l .lt. 1) go to 130
!     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + abs(d(k))
!
         if (scale .ne. 0.0e0_rk) go to 140
!
         do 125 j = 1, l
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = 0.0e0_rk
  125    continue
!
  130    e(i) = 0.0e0_rk
         e2(i) = 0.0e0_rk
         go to 300
!
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
!
         e2(i) = scale * scale * h
         f = d(l)
         g = -sign(sqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         if (l .eq. 1) go to 285
!     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0e0_rk
!
         do 240 j = 1, l
            f = d(j)
            g = e(j) + a(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
!
            do 200 k = jp1, l
               g = g + a(k,j) * d(k)
               e(k) = e(k) + a(k,j) * f
  200       continue
!
  220       e(j) = g
  240    continue
!     .......... form p ..........
         f = 0.0e0_rk
!
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
!
         h = f / (h + h)
!     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - h * d(j)
!     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
!
            do 260 k = j, l
  260       a(k,j) = a(k,j) - f * e(k) - g * d(k)
!
  280    continue
!
  285    do 290 j = 1, l
            f = d(j)
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = f * scale
  290    continue
!
  300 continue
!
      return
      end subroutine
!*********************************************************************************
      subroutine tred2(nm,n,a,d,e,z)
!
      integer i,j,k,l,n,ii,nm,jp1
      real(kind=rk)  :: a(nm,n),d(n),e(n),z(nm,n)
      real(kind=rk)  :: f,g,h,hh,scale
!
!     this subroutine is a translation of the algol procedure tred2,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix to a
!     symmetric tridiagonal matrix using and accumulating
!     orthogonal similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the real symmetric input matrix.  only the
!          lower triangle of the matrix need be supplied.
!
!     on output
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        z contains the orthogonal transformation matrix
!          produced in the reduction.
!
!        a and z may coincide.  if distinct, a is unaltered.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      do 100 i = 1, n
!
         do 80 j = i, n
   80    z(j,i) = a(j,i)
!
         d(i) = a(n,i)
  100 continue
!
      if (n .eq. 1) go to 510
!     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0e0_rk
         scale = 0.0e0_rk
         if (l .lt. 2) go to 130
!     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + abs(d(k))
!
         if (scale .ne. 0.0e0_rk) go to 140
  130    e(i) = d(l)
!
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0e0_rk
            z(j,i) = 0.0e0_rk
  135    continue
!
         go to 290
!
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
!
         f = d(l)
         g = -sign(sqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
!     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0e0_rk
!
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
!
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
!
  220       e(j) = g
  240    continue
!     .......... form p ..........
         f = 0.0e0_rk
!
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
!
         hh = f / (h + h)
!     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
!     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
!
            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
!
            d(j) = z(l,j)
            z(i,j) = 0.0e0_rk
  280    continue
!
  290    d(i) = h
  300 continue
!     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0e0_rk
         h = d(i)
         if (h .eq. 0.0e0_rk) go to 380
!
         do 330 k = 1, l
  330    d(k) = z(k,i) / h
!
         do 360 j = 1, l
            g = 0.0e0_rk
!
            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)
!
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue
!
  380    do 400 k = 1, l
  400    z(k,i) = 0.0e0_rk
!
  500 continue
!
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0e0_rk
  520 continue
!
      z(n,n) = 1.0e0_rk
      e(1) = 0.0e0_rk
      return
      end subroutine
!*********************************************************************************
      subroutine tred3(n,nv,a,d,e,e2)
!
      integer i,j,k,l,n,ii,iz,jk,nv,jm1
      real(kind=rk)  :: a(nv),d(n),e(n),e2(n)
      real(kind=rk)  :: f,g,h,hh,scale
!
!     this subroutine is a translation of the algol procedure tred3,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix, stored as
!     a one-dimensional array, to a symmetric tridiagonal matrix
!     using orthogonal similarity transformations.
!
!     on input
!
!        n is the order of the matrix.
!
!        nv must be set to the dimension of the array parameter a
!          as declared in the calling program dimension statement.
!
!        a contains the lower triangle of the real symmetric
!          input matrix, stored row-wise as a one-dimensional
!          array, in its first n*(n+1)/2 positions.
!
!     on output
!
!        a contains information about the orthogonal
!          transformations used in the reduction.
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
!     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         iz = (i * l) / 2
         h = 0.0e0_rk
         scale = 0.0e0_rk
         if (l .lt. 1) go to 130
!     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
            iz = iz + 1
            d(k) = a(iz)
            scale = scale + abs(d(k))
  120    continue
!
         if (scale .ne. 0.0e0_rk) go to 140
  130    e(i) = 0.0e0_rk
         e2(i) = 0.0e0_rk
         go to 290
!
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
!
         e2(i) = scale * scale * h
         f = d(l)
         g = -sign(sqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         a(iz) = scale * d(l)
         if (l .eq. 1) go to 290
         jk = 1
!
         do 240 j = 1, l
            f = d(j)
            g = 0.0e0_rk
            jm1 = j - 1
            if (jm1 .lt. 1) go to 220
!
            do 200 k = 1, jm1
               g = g + a(jk) * d(k)
               e(k) = e(k) + a(jk) * f
               jk = jk + 1
  200       continue
!
  220       e(j) = g + a(jk) * f
            jk = jk + 1
  240    continue
!     .......... form p ..........
         f = 0.0e0_rk
!
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
!
         hh = f / (h + h)
!     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
!
         jk = 1
!     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
!
            do 260 k = 1, j
               a(jk) = a(jk) - f * e(k) - g * d(k)
               jk = jk + 1
  260       continue
!
  280    continue
!
  290    d(i) = a(iz+1)
         a(iz+1) = scale * sqrt(h)
  300 continue
!
      return
      end subroutine
!*********************************************************************************
      subroutine tridib(n,eps1,d,e,e2,lb,ub,m11,m,w,ind,ierr,rv4,rv5)
!
      integer i,j,k,l,m,n,p,q,r,s,ii,m1,m2,m11,m22,tag,ierr,isturm
      real(kind=rk)  :: d(n),e(n),e2(n),w(m),rv4(n),rv5(n)
      real(kind=rk)  :: u,v,lb,t1,t2,ub,xu,x0,x1,eps1,tst1,tst2
      integer ind(m)
!
!     this subroutine is a translation of the algol procedure bisect,
!     num. math. 9, 386-393(1967) by barth, martin, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 249-256(1971).
!
!     this subroutine finds those eigenvalues of a tridiagonal
!     symmetric matrix between specified boundary indices,
!     using bisection.
!
!     on input
!
!        n is the order of the matrix.
!
!        eps1 is an absolute error tolerance for the computed
!          eigenvalues.  if the input eps1 is non-positive,
!          it is reset for each submatrix to a default value,
!          namely, minus the product of the relative machine
!          precision and the 1-norm of the submatrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2(1) is arbitrary.
!
!        m11 specifies the lower boundary index for the desired
!          eigenvalues.
!
!        m specifies the number of eigenvalues desired.  the upper
!          boundary index m22 is then obtained as m22=m11+m-1.
!
!     on output
!
!        eps1 is unaltered unless it has been reset to its
!          (last) default value.
!
!        d and e are unaltered.
!
!        elements of e2, corresponding to elements of e regarded
!          as negligible, have been replaced by zero causing the
!          matrix to split into a direct sum of submatrices.
!          e2(1) is also set to zero.
!
!        lb and ub define an interval containing exactly the desired
!          eigenvalues.
!
!        w contains, in its first m positions, the eigenvalues
!          between indices m11 and m22 in ascending order.
!
!        ind contains in its first m positions the submatrix indices
!          associated with the corresponding eigenvalues in w --
!          1 for eigenvalues belonging to the first submatrix from
!          the top, 2 for those belonging to the second submatrix, etc..
!
!        ierr is set to
!          zero       for normal return,
!          3*n+1      if multiple eigenvalues at index m11 make
!                     unique selection impossible,
!          3*n+2      if multiple eigenvalues at index m22 make
!                     unique selection impossible.
!
!        rv4 and rv5 are temporary storage arrays.
!
!     note that subroutine tql1, imtql1, or tqlrat is generally faster
!     than tridib, if more than n/4 eigenvalues are to be found.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      tag = 0
      xu = d(1)
      x0 = d(1)
      u = 0.0e0_rk
!     .......... look for small sub-diagonal entries and determine an
!                interval containing all the eigenvalues ..........
      do 40 i = 1, n
         x1 = u
         u = 0.0e0_rk
         if (i .ne. n) u = abs(e(i+1))
         xu = min(d(i)-(x1+u),xu)
         x0 = max(d(i)+(x1+u),x0)
         if (i .eq. 1) go to 20
         tst1 = abs(d(i)) + abs(d(i-1))
         tst2 = tst1 + abs(e(i))
         if (tst2 .gt. tst1) go to 40
   20    e2(i) = 0.0e0_rk
   40 continue
!
      x1 = n
      x1 = x1 * epslon(max(abs(xu),abs(x0)))
      xu = xu - x1
      t1 = xu
      x0 = x0 + x1
      t2 = x0
!     .......... determine an interval containing exactly
!                the desired eigenvalues ..........
      p = 1
      q = n
      m1 = m11 - 1
      if (m1 .eq. 0) go to 75
      isturm = 1
   50 v = x1
      x1 = xu + (x0 - xu) * 0.5e0_rk
      if (x1 .eq. v) go to 980
      go to 320
   60 if (s - m1) 65, 73, 70
   65 xu = x1
      go to 50
   70 x0 = x1
      go to 50
   73 xu = x1
      t1 = x1
   75 m22 = m1 + m
      if (m22 .eq. n) go to 90
      x0 = t2
      isturm = 2
      go to 50
   80 if (s - m22) 65, 85, 70
   85 t2 = x1
   90 q = 0
      r = 0
!     .......... establish and process next submatrix, refining
!                interval by the gerschgorin bounds ..........
  100 if (r .eq. m) go to 1001
      tag = tag + 1
      p = q + 1
      xu = d(p)
      x0 = d(p)
      u = 0.0e0_rk
!
      do 120 q = p, n
         x1 = u
         u = 0.0e0_rk
         v = 0.0e0_rk
         if (q .eq. n) go to 110
         u = abs(e(q+1))
         v = e2(q+1)
  110    xu = min(d(q)-(x1+u),xu)
         x0 = max(d(q)+(x1+u),x0)
         if (v .eq. 0.0e0_rk) go to 140
  120 continue
!
  140 x1 = epslon(max(abs(xu),abs(x0)))
      if (eps1 .le. 0.0e0_rk) eps1 = -x1
      if (p .ne. q) go to 180
!     .......... check for isolated root within interval ..........
      if (t1 .gt. d(p) .or. d(p) .ge. t2) go to 940
      m1 = p
      m2 = p
      rv5(p) = d(p)
      go to 900
  180 x1 = x1 * (q - p + 1)
      lb = max(t1,xu-x1)
      ub = min(t2,x0+x1)
      x1 = lb
      isturm = 3
      go to 320
  200 m1 = s + 1
      x1 = ub
      isturm = 4
      go to 320
  220 m2 = s
      if (m1 .gt. m2) go to 940
!     .......... find roots by bisection ..........
      x0 = ub
      isturm = 5
!
      do 240 i = m1, m2
         rv5(i) = ub
         rv4(i) = lb
  240 continue
!     .......... loop for k-th eigenvalue
!                for k=m2 step -1 until m1 do --
!                (-do- not used to legalize -computed go to-) ..........
      k = m2
  250    xu = lb
!     .......... for i=k step -1 until m1 do -- ..........
         do 260 ii = m1, k
            i = m1 + k - ii
            if (xu .ge. rv4(i)) go to 260
            xu = rv4(i)
            go to 280
  260    continue
!
  280    if (x0 .gt. rv5(k)) x0 = rv5(k)
!     .......... next bisection step ..........
  300    x1 = (xu + x0) * 0.5e0_rk
         if ((x0 - xu) .le. abs(eps1)) go to 420
         tst1 = 2.0e0_rk * (abs(xu) + abs(x0))
         tst2 = tst1 + (x0 - xu)
         if (tst2 .eq. tst1) go to 420
!     .......... in-line procedure for sturm sequence ..........
  320    s = p - 1
         u = 1.0e0_rk
!
         do 340 i = p, q
            if (u .ne. 0.0e0_rk) go to 325
            v = abs(e(i)) / epslon(1.0e0_rk)
            if (e2(i) .eq. 0.0e0_rk) v = 0.0e0_rk
            go to 330
  325       v = e2(i) / u
  330       u = d(i) - x1 - v
            if (u .lt. 0.0e0_rk) s = s + 1
  340    continue
!
         go to (60,80,200,220,360), isturm
!     .......... refine intervals ..........
  360    if (s .ge. k) go to 400
         xu = x1
         if (s .ge. m1) go to 380
         rv4(m1) = x1
         go to 300
  380    rv4(s+1) = x1
         if (rv5(s) .gt. x1) rv5(s) = x1
         go to 300
  400    x0 = x1
         go to 300
!     .......... k-th eigenvalue found ..........
  420    rv5(k) = x1
      k = k - 1
      if (k .ge. m1) go to 250
!     .......... order eigenvalues tagged with their
!                submatrix associations ..........
  900 s = r
      r = r + m2 - m1 + 1
      j = 1
      k = m1
!
      do 920 l = 1, r
         if (j .gt. s) go to 910
         if (k .gt. m2) go to 940
         if (rv5(k) .ge. w(l)) go to 915
!
         do 905 ii = j, s
            i = l + s - ii
            w(i+1) = w(i)
            ind(i+1) = ind(i)
  905    continue
!
  910    w(l) = rv5(k)
         ind(l) = tag
         k = k + 1
         go to 920
  915    j = j + 1
  920 continue
!
  940 if (q .lt. n) go to 100
      go to 1001
!     .......... set error -- interval cannot be found containing
!                exactly the desired eigenvalues ..........
  980 ierr = 3 * n + isturm
 1001 lb = t1
      ub = t2
      return
      end subroutine
!*********************************************************************************
      subroutine tsturm(nm,n,eps1,d,e,e2,lb,ub,mm,m,w,z, ierr,rv1,rv2,rv3,rv4,rv5,rv6)
!
      integer i,j,k,m,n,p,q,r,s,ii,ip,jj,mm,m1,m2,nm,its,ierr,group,isturm
      real(kind=rk)  :: d(n),e(n),e2(n),w(mm),z(nm,mm), rv1(n),rv2(n),rv3(n),rv4(n),rv5(n),rv6(n)
      real(kind=rk)  :: u,v,lb,t1,t2,ub,uk,xu,x0,x1,eps1,eps2,eps3,eps4,norm,tst1,tst2
!
!     this subroutine is a translation of the algol procedure tristurm
!     by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
!
!     this subroutine finds those eigenvalues of a tridiagonal
!     symmetric matrix which lie in a specified interval and their
!     associated eigenvectors, using bisection and inverse iteration.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        eps1 is an absolute error tolerance for the computed
!          eigenvalues.  it should be chosen commensurate with
!          relative perturbations in the matrix elements of the
!          order of the relative machine precision.  if the
!          input eps1 is non-positive, it is reset for each
!          submatrix to a default value, namely, minus the
!          product of the relative machine precision and the
!          1-norm of the submatrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2(1) is arbitrary.
!
!        lb and ub define the interval to be searched for eigenvalues.
!          if lb is not less than ub, no eigenvalues will be found.
!
!        mm should be set to an upper bound for the number of
!          eigenvalues in the interval.  warning. if more than
!          mm eigenvalues are determined to lie in the interval,
!          an error return is made with no values or vectors found.
!
!     on output
!
!        eps1 is unaltered unless it has been reset to its
!          (last) default value.
!
!        d and e are unaltered.
!
!        elements of e2, corresponding to elements of e regarded
!          as negligible, have been replaced by zero causing the
!          matrix to split into a direct sum of submatrices.
!          e2(1) is also set to zero.
!
!        m is the number of eigenvalues determined to lie in (lb,ub).
!
!        w contains the m eigenvalues in ascending order if the matrix
!          does not split.  if the matrix splits, the eigenvalues are
!          in ascending order for each submatrix.  if a vector error
!          exit is made, w contains those values already found.
!
!        z contains the associated set of orthonormal eigenvectors.
!          if an error exit is made, z contains those vectors
!          already found.
!
!        ierr is set to
!          zero       for normal return,
!          3*n+1      if m exceeds mm.
!          4*n+r      if the eigenvector corresponding to the r-th
!                     eigenvalue fails to converge in 5 iterations.
!
!        rv1, rv2, rv3, rv4, rv5, and rv6 are temporary storage arrays.
!
!     the algol procedure sturmcnt contained in tristurm
!     appears in tsturm in-line.
!
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      t1 = lb
      t2 = ub
!     .......... look for small sub-diagonal entries ..........
      do 40 i = 1, n
         if (i .eq. 1) go to 20
         tst1 = abs(d(i)) + abs(d(i-1))
         tst2 = tst1 + abs(e(i))
         if (tst2 .gt. tst1) go to 40
   20    e2(i) = 0.0e0_rk
   40 continue
!     .......... determine the number of eigenvalues
!                in the interval ..........
      p = 1
      q = n
      x1 = ub
      isturm = 1
      go to 320
   60 m = s
      x1 = lb
      isturm = 2
      go to 320
   80 m = m - s
      if (m .gt. mm) go to 980
      q = 0
      r = 0
!     .......... establish and process next submatrix, refining
!                interval by the gerschgorin bounds ..........
  100 if (r .eq. m) go to 1001
      p = q + 1
      xu = d(p)
      x0 = d(p)
      u = 0.0e0_rk
!
      do 120 q = p, n
         x1 = u
         u = 0.0e0_rk
         v = 0.0e0_rk
         if (q .eq. n) go to 110
         u = abs(e(q+1))
         v = e2(q+1)
  110    xu = min(d(q)-(x1+u),xu)
         x0 = max(d(q)+(x1+u),x0)
         if (v .eq. 0.0e0_rk) go to 140
  120 continue
!
  140 x1 = epslon(max(abs(xu),abs(x0)))
      if (eps1 .le. 0.0e0_rk) eps1 = -x1
      if (p .ne. q) go to 180
!     .......... check for isolated root within interval ..........
      if (t1 .gt. d(p) .or. d(p) .ge. t2) go to 940
      r = r + 1
!
      do 160 i = 1, n
  160 z(i,r) = 0.0e0_rk
!
      w(r) = d(p)
      z(p,r) = 1.0e0_rk
      go to 940
  180 u = q-p+1
      x1 = u * x1
      lb = max(t1,xu-x1)
      ub = min(t2,x0+x1)
      x1 = lb
      isturm = 3
      go to 320
  200 m1 = s + 1
      x1 = ub
      isturm = 4
      go to 320
  220 m2 = s
      if (m1 .gt. m2) go to 940
!     .......... find roots by bisection ..........
      x0 = ub
      isturm = 5
!
      do 240 i = m1, m2
         rv5(i) = ub
         rv4(i) = lb
  240 continue
!     .......... loop for k-th eigenvalue
!                for k=m2 step -1 until m1 do --
!                (-do- not used to legalize -computed go to-) ..........
      k = m2
  250    xu = lb
!     .......... for i=k step -1 until m1 do -- ..........
         do 260 ii = m1, k
            i = m1 + k - ii
            if (xu .ge. rv4(i)) go to 260
            xu = rv4(i)
            go to 280
  260    continue
!
  280    if (x0 .gt. rv5(k)) x0 = rv5(k)
!     .......... next bisection step ..........
  300    x1 = (xu + x0) * 0.5e0_rk
         if ((x0 - xu) .le. abs(eps1)) go to 420
         tst1 = 2.0e0_rk * (abs(xu) + abs(x0))
         tst2 = tst1 + (x0 - xu)
         if (tst2 .eq. tst1) go to 420
!     .......... in-line procedure for sturm sequence ..........
  320    s = p - 1
         u = 1.0e0_rk
!
         do 340 i = p, q
            if (u .ne. 0.0e0_rk) go to 325
            v = abs(e(i)) / epslon(1.0e0_rk)
            if (e2(i) .eq. 0.0e0_rk) v = 0.0e0_rk
            go to 330
  325       v = e2(i) / u
  330       u = d(i) - x1 - v
            if (u .lt. 0.0e0_rk) s = s + 1
  340    continue
!
         go to (60,80,200,220,360), isturm
!     .......... refine intervals ..........
  360    if (s .ge. k) go to 400
         xu = x1
         if (s .ge. m1) go to 380
         rv4(m1) = x1
         go to 300
  380    rv4(s+1) = x1
         if (rv5(s) .gt. x1) rv5(s) = x1
         go to 300
  400    x0 = x1
         go to 300
!     .......... k-th eigenvalue found ..........
  420    rv5(k) = x1
      k = k - 1
      if (k .ge. m1) go to 250
!     .......... find vectors by inverse iteration ..........
      norm = abs(d(p))
      ip = p + 1
!
      do 500 i = ip, q
  500 norm = max(norm, abs(d(i)) + abs(e(i)))
!     .......... eps2 is the criterion for grouping,
!                eps3 replaces zero pivots and equal
!                roots are modified by eps3,
!                eps4 is taken very small to avoid overflow ..........
      eps2 = 1.0e-3_rk * norm
      eps3 = epslon(norm)
      uk = q - p + 1
      eps4 = uk * eps3
      uk = eps4 / sqrt(uk)
      group = 0
      s = p
!
      do 920 k = m1, m2
         r = r + 1
         its = 1
         w(r) = rv5(k)
         x1 = rv5(k)
!     .......... look for close or coincident roots ..........
         if (k .eq. m1) go to 520
         if (x1 - x0 .ge. eps2) group = -1
         group = group + 1
         if (x1 .le. x0) x1 = x0 + eps3
!     .......... elimination with interchanges and
!                initialization of vector ..........
  520    v = 0.0e0_rk
!
         do 580 i = p, q
            rv6(i) = uk
            if (i .eq. p) go to 560
            if (abs(e(i)) .lt. abs(u)) go to 540
            xu = u / e(i)
            rv4(i) = xu
            rv1(i-1) = e(i)
            rv2(i-1) = d(i) - x1
            rv3(i-1) = 0.0e0_rk
            if (i .ne. q) rv3(i-1) = e(i+1)
            u = v - xu * rv2(i-1)
            v = -xu * rv3(i-1)
            go to 580
  540       xu = e(i) / u
            rv4(i) = xu
            rv1(i-1) = u
            rv2(i-1) = v
            rv3(i-1) = 0.0e0_rk
  560       u = d(i) - x1 - xu * v
            if (i .ne. q) v = e(i+1)
  580    continue
!
         if (u .eq. 0.0e0_rk) u = eps3
         rv1(q) = u
         rv2(q) = 0.0e0_rk
         rv3(q) = 0.0e0_rk
!     .......... back substitution
!                for i=q step -1 until p do -- ..........
  600    do 620 ii = p, q
            i = p + q - ii
            rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i)
            v = u
            u = rv6(i)
  620    continue
!     .......... orthogonalize with respect to previous
!                members of group ..........
         if (group .eq. 0) go to 700
!
         do 680 jj = 1, group
            j = r - group - 1 + jj
            xu = 0.0e0_rk
!
            do 640 i = p, q
  640       xu = xu + rv6(i) * z(i,j)
!
            do 660 i = p, q
  660       rv6(i) = rv6(i) - xu * z(i,j)
!
  680    continue
!
  700    norm = 0.0e0_rk
!
         do 720 i = p, q
  720    norm = norm + abs(rv6(i))
!
         if (norm .ge. 1.0e0_rk) go to 840
!     .......... forward substitution ..........
         if (its .eq. 5) go to 960
         if (norm .ne. 0.0e0_rk) go to 740
         rv6(s) = eps4
         s = s + 1
         if (s .gt. q) s = p
         go to 780
  740    xu = eps4 / norm
!
         do 760 i = p, q
  760    rv6(i) = rv6(i) * xu
!     .......... elimination operations on next vector
!                iterate ..........
  780    do 820 i = ip, q
            u = rv6(i)
!     .......... if rv1(i-1) .eq. e(i), a row interchange
!                was performed earlier in the
!                triangularization process ..........
            if (rv1(i-1) .ne. e(i)) go to 800
            u = rv6(i-1)
            rv6(i-1) = rv6(i)
  800       rv6(i) = u - rv4(i) * rv6(i-1)
  820    continue
!
         its = its + 1
         go to 600
!     .......... normalize so that sum of squares is
!                1 and expand to full order ..........
  840    u = 0.0e0_rk
!
         do 860 i = p, q
  860    u = pythag(u,rv6(i))
!
         xu = 1.0e0_rk / u
!
         do 880 i = 1, n
  880    z(i,r) = 0.0e0_rk
!
         do 900 i = p, q
  900    z(i,r) = rv6(i) * xu
!
         x0 = x1
  920 continue
!
  940 if (q .lt. n) go to 100
      go to 1001
!     .......... set error -- non-converged eigenvector ..........
  960 ierr = 4 * n + r
      go to 1001
!     .......... set error -- underestimate of number of
!                eigenvalues in interval ..........
  980 ierr = 3 * n + 1
 1001 lb = t1
      ub = t2
      return
      end subroutine
!*********************************************************************************
end module eispack
