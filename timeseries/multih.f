C FILE: MULTIH.F
      subroutine sfou(kk,t,f,ll,fr0,frs,nh,mode,th)
c   Purpose: Computes periodogram by fitting multiharmonic Fourier series
c   Method: Orthogonal projections on base polynomials in complex space
c   Computing cost: time scales with nh*ll, memory scales with 8*mk (words)
c   Please quote in publications the reference:
c         Schwarzenberg-Czerny, A., 1996, ApJ 460, L107.
c   (C) A.Schwarzenberg-Czerny, e-mail: alex@camk.edu.pl
c
c   Input parameters:
c   kk - number of observations (kk<mk)
c   t(kk) - times of observations	(double precision)
c   f(kk) - values of observations
c   ll - number of frequencies in periodogram
c   fr0 - first frequency (units of 1/time unit)
c   frs - frequency step (must satisfy frs<0.5/(t_end-t_start) )
c   nh - number of harmonics in use (nh<mn/2)
c   mode - output mode (recommended mode=1)
c   Note: To prevent undersampling caused by excessively large frs
c   remember that resolution in the periodogram grows by factor nh
c   with respect to that in power spectrum.
c
c   Output parameter:
c   th(ll) - values of the periodogram, depend on mode:
c      - AoV Fisher-Snedecor F(df1,df2) statistic for mode=1
c      - total power fitted in all harmonics for mode=2
c      - 'total amplitude=sqrt(total power)' for mode=3
c        Note: for mode=3 amplitude requires some renormalization
c        depending on ll (not implemented)
c   Message singular frequency means that phase coverage
c   is insufficient at this frequency for Fourier series fit

      implicit none
      integer k,kk,mode,n,nn,nh,l,ll,mk,mn
      complex z,zn,cf,p,pc,san,sad,al,a,sn,sc,c
      real*8 f,th,fr0,frs,df1,df2,sav,s2,tol,ave
      double precision t,om,pi2,fr0d,frsd
      logical sing
      parameter (mk=500000,mn=100,tol=1e-6)
      dimension t(*),f(*),cf(mk),p(mk),z(mk),zn(mk),
     $   a(0:mn),c(0:mn),th(*)

Cf2py intent(in) kk
Cf2py intent(in) t
Cf2py intent(in) f
Cf2py intent(in) ll
Cf2py intent(in) fr0
Cf2py intent(in) frs
Cf2py intent(in) nh
Cf2py intent(in) mode
Cf2py intent(in) th
Cf2py intent(out) th

      nn=nh*2
      if (nn.gt.mn.or.kk.gt.mk) then
        write(*,*) ' SFOU: Wrong dimensions'
        stop
      endif
      pi2=2d0*atan2(0d0,-1d0)
      fr0d=fr0
      frsd=frs
      df1=kk-nn-1
      df2=nn
      ave=0.
      do 1 k=1,kk
        ave=ave+f(k)
 1    continue
      ave=ave/kk
      s2=0.
      do 2 k=1,kk
        sav=f(k)-ave
        s2=s2+sav*sav
 2    continue

      do 10 l=1,ll
        sing=.false.
        th(l)=0.
        al=0.
        om=pi2*(fr0d+(l-1)*frsd)
        do 11 k=1,kk
          sav=mod(om*t(k),pi2)
          z(k)=exp(cmplx(0.,sav))
          zn(k)=conjg(z(k))
          p(k)=zn(k)
          sav=nh*sav
          cf(k)=(f(k)-ave)*exp(cmplx(0.,sav))
 11     continue

        do 12 n=0,nn
          san=0.
          sad=0.
          sn=0.
          sc=0.

          do 13 k=1,kk
            pc=p(k)
            pc=z(k)*pc-al*zn(k)*conjg(pc)
            p(k)=pc
            pc=conjg(pc)
            sav=pc*conjg(pc)
            sn=sn+sav
            sc=sc+cf(k)*pc
            zn(k)=zn(k)*z(k)
            sad=sad+zn(k)*pc*pc
            san=san+z(k)*sav
 13       continue
          if (abs(sad)+abs(sn).le.tol) then
            if (.not.sing) then
              write(*,*) ' SFOU: Singular frequency',om/pi2
            endif
            al=0.
            c(n)=0.
            sing=.true.
          else
            al=san/sad
            c(n)=sc/sqrt(sn)
          endif
          a(n)=al
          th(l)=th(l)+c(n)*conjg(c(n))
 12     continue
        if (mode.eq.2) then
          th(l)=th(l)
        elseif (mode.eq.3) then
          th(l)=sqrt(th(l))
        elseif (mode.eq.4) then
c To convert F-values into probabilities use Numerical Recipes f77 routines
c betai,betacf,gammln. This is rather slow and you may prefere doing it for
c the max values only.
c         th(l)=betai(half*df2,half*df1,df2/(df2+df1*th(l)))
        else
          th(l)=df1*th(l)/(df2*max(s2-th(l),1e-32))
        endif
 10   continue

      end        

      complex function dotp(nz,z,w)
      integer nz,iz
      complex z,w,sum
      dimension z(*),w(*)
      sum=0.
      do 1 iz=1,nz
 1      sum=sum+z(iz)*conjg(w(iz))
      dotp=sum/nz
      end

c                                                                               
















