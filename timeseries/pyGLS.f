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
      DOUBLE PRECISION v(500000),z(500000),ww(500000),cosx,sinx
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
      DOUBLE PRECISION WS,WC,v(500000)
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
      subroutine gls(JD_,RV_,RVerr_,N_,Fbeg,Fend,step,wexp,
     $f1,s1,p1,l1,maxstep)
      implicit none
      INTEGER i,j,Nmax,N,wexp,N_, index
      INTEGER maxstep
c      ! wexp = weightening exponent (default Chi2)
c     !!! weightening of errors: wexp=0 (variance), wexp=2 (Chi2) !!!
      LOGICAL exists
      CHARACTER TV*40, fittype*8/"unknown"/
      DOUBLE PRECISION twopi
      PARAMETER(twopi=2.*3.141592654)
      PARAMETER(Nmax=500000)  ! max. data points
      DOUBLE PRECISION JD_(N_),RV_(N_),RVerr_(N_)
      DOUBLE PRECISION f1(maxstep),s1(maxstep),p1(maxstep),l1(maxstep)
      DOUBLE PRECISION JD(Nmax),RV(Nmax),RVerr(Nmax),phase(Nmax)
      DOUBLE PRECISION t(Nmax),JDmin/1.E38/,JDmax/-1.E38/
      DOUBLE PRECISION PSin,Freq,Fbeg,Fend,step
      DOUBLE PRECISION ww(Nmax),v(Nmax),wy(Nmax),RVmean,m33,YY
      DOUBLE PRECISION pow,powLS,powSin,SW
      DOUBLE PRECISION A,B,C,CBest, Amp,ph
      DOUBLE PRECISION dummy,FAP,M,prob
      COMMON /Messwerte/ v,wy,ww,YY,RVmean,N
c      DATA RVmean,m33,YY,powSin /4* 0./
      
c     SET default values to some stuff
      do 90 i=1,Nmax
         v(i) = 0.
         wy(i) = 0.
         ww(i) = 0.
 90   continue
c      write(*,*) 'before',RVMean
      RVmean = 0
      m33 = 0
      YY = 0
      powSin = 0
c      write(*,*) 'after',RVmean
c     read parameters

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

c     Fit a sine function
      i=1
      do 80 Freq=Fbeg,Fend, step
         call Phases(t,twopi*Freq,N,v)
         call SineFit(A,B,C,pow, powLS)
         call Spectral_Window(v,N,SW)
         f1(i) = Freq
         s1(i) = pow
         p1(i) = SW
         l1(i) = powLS
         i = i+1
         if (pow .GT. powSin) then
            powSin=pow
            ph   =MOD( DATAN2(A,B)+twopi, twopi)
            Amp  =A/sin(ph)
            CBest=C
            PSin =1./Freq
         endif
 80   continue

      prob=(1.-powSin)**((N-3)/2)              ! spectral false alarm probability
      M   = (JDmax-JDmin) * abs(Fend-Fbeg)     ! Number of independent frequencies
      FAP = M * prob
      if( FAP .gt. 0.01 ) FAP = 1. - (1.-prob)**M

      END