      subroutine internal_df_dx(N,M,oldx,oldy,myindex,sharp,df_dx)
Cf2py intent(in) N
Cf2py intent(in) M
Cf2py intent(in) oldx
Cf2py intent(in) oldy
Cf2py intent(in) myindex
Cf2py intent(in) sharp
Cf2py intent(out) df_dx
      
      integer N,M,myindex
      real*8 oldx(N),oldy(M)
      logical sharp
      real*8 df_dx
      
      real*8 xm1,fm1,x,f
      real*8 xp1,fp1
      
      xm1 = oldx(myindex-1)
      fm1 = oldy(myindex-1)
      x   = oldx(myindex)
      f   = oldy(myindex)
      xp1 = oldx(myindex+1)
      fp1 = oldy(myindex+1)
      
      
      if (sharp.eqV..true.) then
        df_dx = min( (fp1-f)/(xp1-x), (f-fm1)/(x-xm1)  )
      else
        df_dx=1.0/(xp1-xm1)*(fp1*(x-xm1)/(xp1-x)-fm1*(xp1-x)/(x-xm1))+
     &       f*(xp1-2.0**x+xm1) / ( (x-xm1)*(xp1-x))
      end if
      
      RETURN
      END
      
     
      subroutine local_interpolation(N,M,newx,oldx,oldy,newy)     
Cf2py intent(in) N
Cf2py intent(in) M
Cf2py intent(in) newx
Cf2py intent(in) oldx
Cf2py intent(in) oldy
Cf2py intent(out) newy

      integer N,M,myindex,i
      real*8 oldx(M),oldy(M)
      real*8 oldxn(M+1),oldyn(M+1)
      real*8 newx(N),newy(N)
      real*8 numerator,denominator,sharpness,dfdx0,dfdx1
      real*8 x0,f0,x1,f1,x2,f2,P1,P2,P3,P4
      logical sharp
     
      df_dx = 0.0
C     extend axis to be able to interpolate last point
      oldxn(+1) = 2*oldx(M)-oldx(M-1)
      oldyn(M+1) = (oldy(M)-oldy(M-1))/(oldx(M)-oldx(M-1))*
     &             (oldx(M)-lastx) + oldy(M)
      do 10, i=1,M
        oldxn(i) = oldx(i)
        oldyn(i) = oldy(i)
   10 continue 
      
      myindex = 0
      do 20, i=1,N
   90   if (oldxn(myindex).LE.newx(i)) then
          myindex = myindex + 1
          goto 90
        endif
      x0 = oldxn(myindex-1)
      f0 = oldyn(myindex-1)
      x1 = oldxn(myindex)
      f1 = oldyn(myindex)
      x2 = oldxn(myindex+1)
      f2 = oldyn(myindex+1)
      numerator = ((x1-x0)*(f1-f2))
      denominator = ((f1-f0)*(x1-x2))
      
      if (denominator.eq.0.0) then
        sharpness = 0
      else
        sharpness = numerator/denominator
      endif
      
      sharp = (0.2.le.sharpness).AND.(sharpness.le.0.5) 
      
      call internal_df_dx(M+1,N,oldxn,oldyn,myindex-1,sharp,dfdx0)
      call internal_df_dx(M+1,N,oldxn,oldyn,myindex,sharp,dfdx1)
      
      P1 =  (x-x1)**2.0 * (2.0*x-3.0*x0+x1) / (x1-x0)**3.0
      P2 = -(x-x0)**2.0 * (2.0*x-3.0*x1+x0) / (x1-x0)**3.0
      P3 =  (x-x0)    * (x-x1)**2.0     / (x1-x0)**2.0
      P4 =  (x-x0)**2.0 * (x-x1)        / (x1-x0)**2.0
      
      newy(i) = f0*P1 + f1*P2 + dfdx0*P3 + dfdx1*P4
   20 continue
      RETURN
      END
      
      

      PROGRAM dummy
      END