C FILE: PYDFT.F
      SUBROUTINE FT(XX,TSAM,NN,WZ,NFREQ,SI,LFREQ,T0,MM,DF,FTRX,FTIX,O,W)
	  
	  DIMENSION XX(NN), TSAM(NN)
	  integer NN
	  real*8 WZ
	  integer NFREQ
      real*8 si
      integer LFREQ
      real*8 T0,DF
      integer MM
	  DIMENSION FTRX(MM), FTIX(MM), O(MM),W(MM)
	  REAL*8 WTAN, SSIN
	  COMPLEX WORK

Cf2py intent(in) XX
Cf2py intent(in) TSAM
Cf2py intent(in) NN
Cf2py intent(in) WZ
Cf2py intent(in) NFREQ
Cf2py intent(in) SI
Cf2py intent(in) LFREQ
Cf2py intent(in) T0
Cf2py intent(in) MM
Cf2py intent(in) DF
Cf2py intent(in) FTRX
Cf2py intent(in) FTIX
Cf2py intent(in) O
Cf2py intent(in) W
Cf2py intent(out) FTRX
Cf2py intent(out) FTIX
Cf2py intent(out) O
Cf2py intent(out) W

	  TOL1 = 1.0 E -04
	  TOL2 = 1.0 E -08
	  WUSE = WZ
	  FNN = FLOAT( NN )
	  CONST1 = 1.0 / SQRT(2.0)
	  CONST2 = SI * CONST1
	  SUMT = 0.0
	  SUMX = 0.0
	  DO 100 I=1,NN
	     SUMT = SUMT + TSAM( I )
		 SUMX = SUMX + XX( I )
100   CONTINUE
	  ISTOP = NFREQ
	  
	  
	  
	  
	  TAU0 = SUMT / FNN         ! LIMIT
	  CSUM = FNN
	  SSUM = 0.0
	  FTRX(1) = SUMX / SQRT( FNN )
	  FTIX(1) = 0.0
	  WDEL = DF
	  WRUN = WUSE
      O(1) = 0.
      W(1) = FTRX(1)
	  II = 2
	  
	  
150   CONTINUE
      CSUM = 0.0
	  SSUM = 0.0
	  SUMTC = 0.0
	  SUMTS = 0.0
	  
	  DO 190 I = 1,NN
	     TTT = TSAM( I )
		 ARG1 = 2.0 * WRUN * TTT
		 ARG = FOLD( ARG1 )
		 TCOS = COS ( ARG )
		 TSIN = SIN ( ARG )
		 CSUM = CSUM + TCOS
		 SSUM = SSUM + TSIN
190   CONTINUE
      
	  WATAN = ATAN2( SSUM, CSUM )
	  IF(ABS(SSUM).GT.TOL1 .OR. ABS(CSUM).GT.TOL1) GOTO 200
	  WATAN = ATAN2( -SUMTC, SUMTS)
200   CONTINUE
      
	  WTAU = 0.5 * WTAN
	  WTNEW = WTAU
	  SUMR = 0.0
	  SUMI = 0.0
	  SCOS2 = 0.0
	  SSIN2 = 0.0
	  CROSS = 0.0
	  
	  DO 440 I = 1,NN
	  
	     TIM = TSAM(I)
		 ARG1 = WRUN * TIM - WTNEW
		 ARG = FOLD( ARG1 )
		 TCOS = COS(ARG)
		 TSIN = SIN(ARG)
		 
		 CROSS = CROSS + TIM * TCOS * TSIN
		 SCOS2 = SCOS2 + TCOS * TCOS
		 SSIN2 = SSIN2 + TSIN * TSIN
		 
		 XD = XX(I)
		 SUMR = SUMR + XD * TCOS
		 SUMI = SUMI + XD * TSIN
		 
440      CONTINUE
      
	  FTRD = CONST1 * SUMR / SQRT(SCOS2)
	  IF ( SSIN .LE. TOL1 ) GOTO 450
	  FTID = CONST2 * SUMI / SQRT(SSIN2)
	  GOTO 460

450   CONTINUE
      FTID = CONST2 * SUMX / SQRT( FNN )
	  IF ( ABS(CROSS) .GT. TOL2 ) FTID = 0.0

460   CONTINUE
      
	  PHASE1 = WTNEW - WRUN * T0
	  PHASE = FOLD( PHASE1 )
	  WORK = CMPLX( FTRD, FTID)*CEXP( CMPLX(0.0, PHASE ) )
	  FTRX(II) = REAL( WORK )
	  FTIX(II) = AIMAG( WORK )
	  O(II) = WRUN
	  W(II) = WORK
	  II = II + 1
	  WRUN = WRUN + WDEL
	  IF( II .LE. ISTOP )GOTO 150
	  
	  
	  
	  
	  IF( 2 * NFREQ .GT. LFREQ ) GOTO 999
	  I1 = NFREQ + 1
	  
	  DO 320 I= I1,LFREQ
	     FTRX(I) = 0.0
		 FTIX(I) = 0.0
320      CONTINUE
      
	  NSTOP = LFREQ / 2
	  DO 340 I=2, NSTOP
	     IPUT = LFREQ -I + 2
		 FTRX(IPUT) =  FTRX(I)
		 FTIX(IPUT) = -FTIX(I)
340      CONTINUE
465   CONTINUE
      
	  RETURN
999   CONTINUE
      END
	  
	  
	  FUNCTION FOLD( ARG )
      
	  PI = 3.1415926535898
	  ARGMAX = 8000.0 * PI
	  FOLD = ARG
10    CONTINUE
      IF( FOLD .LE. ARGMAX) GOTO 20
	     FOLD = FOLD - ARGMAX
		 GOTO 10
20    CONTINUE
      IF( FOLD .GT. -ARGMAX) GOTO 30
	     FOLD = FOLD + ARGMAX
		 GOTO 20
30    CONTINUE
      RETURN
	  END