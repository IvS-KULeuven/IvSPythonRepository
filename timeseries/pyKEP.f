c     ---------------------------------------------------------------
c     GLS version 2.0.10
c     Calculates the generalised Lomb-Scargle periodogram
c                                  (Zechmeister & Kuerster 2009)
c     This program searches for periods with help of
c     the best sine fit and the best Keplerian fit.
c     by Mathias Zechmeister 08.03.2009  mathiaszechmeister@googlemail.com
c     2009/05/01: improved e-T0 grid
c     ---------------------------------------------------------------

      SUBROUTINE SineFit(A,B,off,pow,powLS)
c     fit sine wave y=A*cosx+B*sinx+off
c     A,B,off - fit parameter
      DOUBLE PRECISION pow,powLS,A,B,off, CC,SS,CS,C,S,YC,YS,Y,YY,D
      DOUBLE PRECISION v(10000),z(10000),ww(10000),cosx,sinx
      COMMON /Messwerte/ v,z,ww,YY,Y,N
      CC =0.
      CS =0.
      C  =0.
      S  =0.
      YC =0.
      YS =0.
      do 20 i=1,N                        ! calculate the sums
         cosx = cos(v(i))
         sinx = sin(v(i))
         CC = CC + ww(i) * cosx**2
         CS = CS + ww(i) * cosx*sinx
         C  = C  + ww(i) * cosx
         S  = S  + ww(i) * sinx
         YC = YC + z(i) * cosx
         YS = YS + z(i) * sinx
 20   continue
      SS = 1. - CC
      D  = CC*SS - CS*CS
      powLS=(SS*YC**2/D+CC*YS**2/D-2*CS*YC*YS/D)/ YY         ! Lomb-Scargle power
      CC = CC - C*C
      SS = SS - S*S
      CS = CS - C*S
      D  = CC*SS - CS*CS
      A = (YC*SS-YS*CS) / D
      B = (YS*CC-YC*CS) / D
      off = Y-A*C-B*S
c      pow= (A*YC+B*YS)/ YY
      pow=(SS*YC**2/D+CC*YS**2/D-2*CS*YC*YS/D)/ YY           ! GLS power
      END

      SUBROUTINE Spectral_Window(v,N,WS)
      INTEGER i,N
      DOUBLE PRECISION WS,WC,v(10000)
      WC=0.
      WS=0.
      do 45 i=1,N
         WC=WC + cos(v(i))
         WS=WS + sin(v(i))
 45   continue
      WS = (WC**2+WS**2) / N**2
      return
      END

c                   c! calculation of the true anomaly
      subroutine WaAn(M,e,zwischen)
      DOUBLE PRECISION M,Mn,e,F,Fn,zwischen
      INTEGER i,itermax
      PARAMETER(itermax=8)
      Fn = M+e*sin(M)+e**2/2*sin(2*M)      ! Initial value
c     iterative solving of the transcendent Kepler's equation
      do 50 i=1,itermax
         F  = Fn
         Mn = F-e*sin(F)
         Fn = F+(M-Mn)/(1-e*cos(F))
         if (abs((Fn-F)/F) .lt. 0.00001) goto 51
 50   continue
 51   zwischen = 2.0*ATAN(sqrt((1.+e)/(1.-e))*tan(Fn/2))
      return
      END

      SUBROUTINE Phases(JD,twopiFreq,N,phase)
      DOUBLE PRECISION JD(*),twopiFreq,phase(*)
      INTEGER i,N
      do 60 i=1,N
         phase(i) = JD(i) * twopiFreq
 60   continue
      END


c     ------------------------------------------------------------------
      subroutine kepler(JD_,RV_,RVerr_,N_,Fbeg,Fend,step,wexp,emin,emax,
     $estep,x0min,x0max,f1,s1,p1,l1,s2,maxstep,k2)
      implicit none
      INTEGER i,j,N,wexp,N_, index_bn
      INTEGER maxstep
c      ! wexp = weightening exponent (default Chi2)
c     !!! weightening of errors: wexp=0 (variance), wexp=2 (Chi2) !!!
      LOGICAL exists
      CHARACTER TV*40, fittype*8
      DOUBLE PRECISION twopi,G,AU,Msun,Mjup, fa
      PARAMETER(twopi=2.*3.141592654,G=6.6726E-11,AU=1.496E11)
      PARAMETER(Msun=1.989E30, Mjup = Msun/1047.39)
      INTEGER, PARAMETER :: Nmax=10000  ! max_. data points
      DOUBLE PRECISION JD_(N_),RV_(N_),RVerr_(N_)
      DOUBLE PRECISION f1(maxstep),s1(maxstep),p1(maxstep),l1(maxstep)
      DOUBLE PRECISION s2(maxstep)
      DOUBLE PRECISION k2(6)
      DOUBLE PRECISION JD(Nmax),RV(Nmax),RVerr(Nmax),phase(Nmax)
      DOUBLE PRECISION t(Nmax),JDmin,JDmax
      DOUBLE PRECISION PKep,Freq,Fbeg,Fend,step
      DOUBLE PRECISION ww(Nmax),v(Nmax),wy(Nmax),RVmean,m33,YY
      DOUBLE PRECISION pow,powLS,powSin,powKep,SW, z,wA
      DOUBLE PRECISION A,B,C,CBest, Amp,ph, K,RV0,x0,e,w,xx0,ee
      DOUBLE PRECISION A1sinK,A1sinS,AsinK,AsinS,powKe,zwischen
      DOUBLE PRECISION mass,mass2S,mass2K, dummy,FAP,M,prob
      DOUBLE PRECISION emin,emax,estep                 ! default_ values
      DOUBLE PRECISION x0min,x0max,x0step                 ! default_ values
      COMMON /Messwerte/ v,wy,ww,YY,RVmean,N
c      DATA RVmean,m33,YY, mass2S,mass2K, powSin,powKep /7* 0./
c     set all defaults to zero
      do 99 i=1,Nmax
         ww(i) = 0
         v(i) = 0
c         wy = 0
 99   continue
      YY = 0
      RVmean = 0
      N = 0
      m33 = 0
      mass2S = 0
      mass2K = 0
      powSin = 0
      powKep = 0



      fa=(86400./twopi/G/Msun)**(1./3.)           ! fa= 4.69679026E-06
c     read parameters
      x0min =x0min *twopi/360.                  ! convert to rad
      x0max =x0max *twopi/360.

      if (wexp.eq.2) fittype="chi2"
      if (wexp.eq.0) fittype="variance"
      do 70 i=1,N_
         JD(i) = JD_(i)
         RV(i) = RV_(i)
         RVerr(i)=RVerr_(i)
         JDmin =MIN(JDmin,JD(i))
         JDmax =MAX(JDmax,JD(i))
         N =i
         ww(i)=(1./Rverr(i))**wexp
         m33  = m33 + ww(i)
 70   continue
      do 71 i=1,N
         t(i) = JD(i)-JDmin
         ww(i) = ww(i) / m33               ! normalize weights
         RVmean = RVmean + RV(i)*ww(i)     ! weighted mean
 71   continue
      do 72 i=1,N
         wy(i)=RV(i)-RVmean                ! subtract mean
         YY   =YY +wy(i)**2*ww(i)          ! sum for chi2 above the mean
         wy(i)=wy(i)*ww(i)
 72   continue

c     fit Keplerian orbits RV=A*cosv+B*sinv+C
      j=1
      do 90 Freq=Fbeg,Fend, step
         call Phases(t,twopi*Freq,N,phase)
         powKe=0.
         do 10 ee=emin,emax,estep
            x0step=twopi/int(twopi*ee/estep+1)                       ! polar grid e-x0
            do 20 xx0=x0min,x0max,x0step
c              with known x0 and e parameter A, B und C can be minimized,
c              because the true anomaly is now calcuable:
               do 30 i=1,N
                  call WaAn(phase(i)-xx0,ee,zwischen)
                  v(i)= zwischen
 30            continue
               call SineFit(A,B,C,pow,powLS)
               powKe=MAX (pow, powKe)
               if (pow .GT. powKep) then
                  powKep=pow
                  x0  =xx0
                  e   =ee
                  w   =MOD(-DATAN2(B,A)+twopi, twopi)
                  K   =-B/sin(w)
                  RV0 =C-A*e
                  PKep=1./Freq
               endif
 20         continue
 10      continue
         f1(j) = Freq
         s2(j) = powKe
         j = j+1
 90   continue
      k2(1) = 1./PKep
      k2(2) = x0
      k2(3) = e
      k2(4) = w
      k2(5) = K
      k2(6) = RV0

      END
      program dummy
      end
