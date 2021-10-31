Module blas_1

  USE global_std
  
  Implicit None

Contains

!*******************************************************************************************
      Subroutine  dcopy(n,dx,incx,dy,incy)
!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      Double Precision dx(*),dy(*)
      Integer i,incx,incy,ix,iy,m,mp1,n
!
      If(n.Le.0)Return
      If(incx.Eq.1.And.incy.Eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      If(incx.Lt.0)ix = (-n+1)*incx + 1
      If(incy.Lt.0)iy = (-n+1)*incy + 1
      Do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 Continue
      Return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = Mod(n,7)
      If( m .Eq. 0 ) go to 40
      Do 30 i = 1,m
        dy(i) = dx(i)
   30 Continue
      If( n .Lt. 7 ) Return
   40 mp1 = m + 1
      Do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 Continue
      Return
      End Subroutine dcopy
!*******************************************************************************************
      Subroutine drotg(da,db,c,s)
!
!     construct givens plane rotation.
!     jack dongarra, linpack, 3/11/78.
!
      Real(Kind=rk) :: da,db,c,s,roe,scale,r,z
!
      roe = db
      If( dabs(da) .Gt. dabs(db) ) roe = da
      scale = dabs(da) + dabs(db)
      If( scale .Ne. 0.0d0 ) go to 10
         c = 1.0d0
         s = 0.0d0
         r = 0.0d0
         z = 0.0d0
         go to 20
   10 r = scale*dsqrt((da/scale)**2 + (db/scale)**2)
      r = dsign(1.0d0,roe)*r
      c = da/r
      s = db/r
      z = 1.0d0
      If( dabs(da) .Gt. dabs(db) ) z = s
      If( dabs(db) .Ge. dabs(da) .And. c .Ne. 0.0d0 ) z = 1.0d0/c
   20 da = r
      db = z
      Return
      End Subroutine drotg
!*******************************************************************************************
      Subroutine  drot (n,dx,incx,dy,incy,c,s)
!
!     applies a plane rotation.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      Real(Kind=rk) :: dx(*),dy(*),dtemp,c,s
      Integer i,incx,incy,ix,iy,n
!
      If(n.Le.0)Return
      If(incx.Eq.1.And.incy.Eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      If(incx.Lt.0)ix = (-n+1)*incx + 1
      If(incy.Lt.0)iy = (-n+1)*incy + 1
      Do 10 i = 1,n
        dtemp = c*dx(ix) + s*dy(iy)
        dy(iy) = c*dy(iy) - s*dx(ix)
        dx(ix) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 Continue
      Return
!
!       code for both increments equal to 1
!
   20 Do 30 i = 1,n
        dtemp = c*dx(i) + s*dy(i)
        dy(i) = c*dy(i) - s*dx(i)
        dx(i) = dtemp
   30 Continue
      Return
      End Subroutine drot
!*******************************************************************************************
      Real(Kind=rk) Function dnrm2 ( n, x, incx )
!     .. scalar arguments ..
      Integer                           incx, n
!     .. array arguments ..
      Real(Kind=rk) ::                  x( * )
!     ..
!
!  dnrm2 returns the euclidean norm of a vector via the function
!  name, so that
!
!     dnrm2 := sqrt( x'*x )
!
!
!
!  -- This version written on 25-October-1982.
!     Modified on 14-October-1993 to inline the call to DLASSQ.
!     Sven Hammarling, Nag Ltd.
!
!
!     .. Parameters ..
      Real(Kind=rk) ::      one         , zero
      Parameter           ( one = 1.0d+0, zero = 0.0d+0 )
!     .. local scalars ..
      Integer               ix
      Real(Kind=rk) ::      absxi, norm, scale, ssq
!     .. intrinsic functions ..
      Intrinsic             abs, sqrt
!     ..
!     .. executable statements ..
      If( n.Lt.1 .Or. incx.Lt.1 )Then
         norm  = zero
      Else If( n.Eq.1 )Then
         norm  = Abs( x( 1 ) )
      Else
         scale = zero
         ssq   = one
!        the following loop is equivalent to this call to the lapack
!        auxiliary routine:
!        call dlassq( n, x, incx, scale, ssq )
!
         Do 10, ix = 1, 1 + ( n - 1 )*incx, incx
            If( x( ix ).Ne.zero )Then
               absxi = Abs( x( ix ) )
               If( scale.Lt.absxi )Then
                  ssq   = one   + ssq*( scale/absxi )**2
                  scale = absxi
               Else
                  ssq   = ssq   +     ( absxi/scale )**2
               End If
            End If
   10    Continue
         norm  = scale * Sqrt( ssq )
      End If
!
      dnrm2 = norm
      Return
!
!     end of dnrm2.
!
      End Function dnrm2
!*******************************************************************************************
      Integer Function idamax(n,dx,incx)
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      Real(Kind=rk) :: dx(*),dmax
      Integer i,incx,ix,n
!
      idamax = 0
      If( n.Lt.1 .Or. incx.Le.0 ) Return
      idamax = 1
      If(n.Eq.1)Return
      If(incx.Eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      Do 10 i = 2,n
         If(dabs(dx(ix)).Le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 Continue
      Return
!
!        code for increment equal to 1
!
   20 dmax = dabs(dx(1))
      Do 30 i = 2,n
         If(dabs(dx(i)).Le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 Continue
      Return
      End Function idamax
!*******************************************************************************************
      Real(Kind=rk) Function dasum(n,dx,incx)
!
!     takes the sum of the absolute values.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      Real(Kind=rk) :: dx(*),dtemp
      Integer i,incx,m,mp1,n,nincx
!
      dasum = 0.0d0
      dtemp = 0.0d0
      If( n.Le.0 .Or. incx.Le.0 )Return
      If(incx.Eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      Do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
   10 Continue
      dasum = dtemp
      Return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 m = Mod(n,6)
      If( m .Eq. 0 ) go to 40
      Do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 Continue
      If( n .Lt. 6 ) go to 60
   40 mp1 = m + 1
      Do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))+ dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 Continue
   60 dasum = dtemp
      Return
      End Function dasum
!*******************************************************************************************
      Subroutine daxpy(n,da,dx,incx,dy,incy)
!
!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      Real(Kind=rk) :: dx(*),dy(*),da
      Integer i,incx,incy,ix,iy,m,mp1,n
!
      If(n.Le.0)Return
      If (da .Eq. 0.0d0) Return
      If(incx.Eq.1.And.incy.Eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      If(incx.Lt.0)ix = (-n+1)*incx + 1
      If(incy.Lt.0)iy = (-n+1)*incy + 1
      Do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 Continue
      Return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = Mod(n,4)
      If( m .Eq. 0 ) go to 40
      Do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 Continue
      If( n .Lt. 4 ) Return
   40 mp1 = m + 1
      Do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 Continue
      Return
      End Subroutine daxpy
!*******************************************************************************************
      Real(Kind=rk) Function ddot(n,dx,incx,dy,incy)
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      Real(Kind=rk) :: dx(*),dy(*),dtemp
      Integer i,incx,incy,ix,iy,m,mp1,n
!
      ddot = 0.0d0
      dtemp = 0.0d0
      If(n.Le.0)Return
      If(incx.Eq.1.And.incy.Eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      If(incx.Lt.0)ix = (-n+1)*incx + 1
      If(incy.Lt.0)iy = (-n+1)*incy + 1
      Do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 Continue
      ddot = dtemp
      Return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = Mod(n,5)
      If( m .Eq. 0 ) go to 40
      Do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 Continue
      If( n .Lt. 5 ) go to 60
   40 mp1 = m + 1
      Do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 Continue
   60 ddot = dtemp
      Return
      End Function ddot
!*******************************************************************************************
      Subroutine  dscal(n,da,dx,incx)
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      Real(Kind=rk) :: da,dx(*)
      Integer i,incx,m,mp1,n,nincx
!
      If( n.Le.0 .Or. incx.Le.0 )Return
      If(incx.Eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      Do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 Continue
      Return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 m = Mod(n,5)
      If( m .Eq. 0 ) go to 40
      Do 30 i = 1,m
        dx(i) = da*dx(i)
   30 Continue
      If( n .Lt. 5 ) Return
   40 mp1 = m + 1
      Do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 Continue
      Return
      End Subroutine dscal
!*******************************************************************************************
      Subroutine  dswap (n,dx,incx,dy,incy)
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      Real(Kind=rk) :: dx(*),dy(*),dtemp
      Integer i,incx,incy,ix,iy,m,mp1,n
!
      If(n.Le.0)Return
      If(incx.Eq.1.And.incy.Eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      If(incx.Lt.0)ix = (-n+1)*incx + 1
      If(incy.Lt.0)iy = (-n+1)*incy + 1
      Do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 Continue
      Return
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
   20 m = Mod(n,3)
      If( m .Eq. 0 ) go to 40
      Do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 Continue
      If( n .Lt. 3 ) Return
   40 mp1 = m + 1
      Do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 Continue
      Return
      End Subroutine dswap

!*******************************************************************************************

End Module blas_1
