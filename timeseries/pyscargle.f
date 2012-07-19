C FILE: PYSCARGLE.F
      SUBROUTINE SCAR2 (N,X,T,F0,NF,DF,F1,S1,SS,SC,SS2,SC2)

C     Computation of Scargles periodogram without explicit tau
C     calculation, with iteration (Method Cuypers)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(N),T(N)
      DIMENSION F1(NF),S1(NF)
      DIMENSION SS(NF),SC(NF),SS2(NF),SC2(NF)
      DATA TWOPI,DTWO,DNUL/6.28318530717959D0,2.0D0,0.0D0/
Cf2py intent(in) N
Cf2py intent(in) X
Cf2py intent(in) T
Cf2py intent(in) F0
Cf2py intent(in) NF
Cf2py intent(in) DF
Cf2py intent(in) F1
Cf2py intent(in) S1
Cf2py intent(in) SS
Cf2py intent(in) SC
Cf2py intent(in) SS2
Cf2py intent(in) SC2
Cf2py intent(out) S1
Cf2py intent(out) F1

      F = F0
      TPF = TWOPI * F
      TDF = TWOPI * DF
      TN = DFLOAT(N)
      TNSQ=TN*TN

      DO 20 K=1,NF
         SS(K)  = DNUL
         SC(K)  = DNUL
         SS2(K) = DNUL
         SC2(K) = DNUL
20       CONTINUE

      DO 40 I=1,N

         A = T(I)
         AF0 = DMOD(A*TPF,TWOPI)
         S0 = DSIN(AF0)
         C0 = DCOS(AF0)
         S20 = DTWO * S0 * C0
         C20 = C0 * C0 - S0 * S0

         ADF = DMOD(A*TDF,TWOPI)
         SDF = DSIN(ADF)
         CDF = DCOS(ADF)
         S2DF = DTWO * SDF * CDF
         C2DF = CDF * CDF - SDF * SDF
         XI=X(I)
         C0X = C0 * XI
         S0X = S0 * XI
         DO 30 K=1,NF
            SS(K) = SS(K) + S0X
            SC(K) = SC(K) + C0X
            CTX=C0X
            C0X = CTX * CDF - S0X * SDF
            S0X = S0X * CDF + CTX * SDF
            SS2(K) = SS2(K) + S20
            SC2(K) = SC2(K) + C20
            C2T = C20
            C20 = C2T * C2DF - S20 * S2DF
            S20 = S20 * C2DF + C2T * S2DF
30          CONTINUE

40       CONTINUE 

      DO 50 K=1,NF
         SSK  =  SS(K)
         SS2K = SS2(K)
         SCK  =  SC(K)
         SC2K = SC2(K)
         
         F1(K)=F
         S1(K)=(SCK*SCK*(TN-SC2K)+SSK*SSK*(TN+SC2K)-DTWO*SSK*SCK*SS2K)/
     &      (TNSQ-SC2K*SC2K-SS2K*SS2K)
         F=F+DF
                  
50    CONTINUE
 
      RETURN
      END

      SUBROUTINE JUSTEL(N,X,T,F0,NF,DF,NB,NC,XVAR,XX,F1,S1)

C     BEREKENING VAN HET PERIODOGRAM VOLGENS JURKEVICH
C                                         OF STELLINGWERF
 
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),T(N),F1(NF),S1(NF)
      DIMENSION NBIN(1000),XB(1000),NH(1000),XH(1000)
Cf2py intent(in) N
Cf2py intent(in) X
Cf2py intent(in) T
Cf2py intent(in) F0
Cf2py intent(in) NF
Cf2py intent(in) DF
Cf2py intent(in) NB
Cf2py intent(in) NC
Cf2py intent(in) XVAR
Cf2py intent(in) XX
Cf2py intent(in) F1
Cf2py intent(in) S1
Cf2py intent(out) F1
Cf2py intent(out) S1
      F=F0
       
C-------------------------------------------------------------------
C     NC = 1 : DE METHODE VAN JURKEVICH
C-------------------------------------------------------------------
 
      IF (NC.LE.1) THEN
 
         DF1=DFLOAT(N-NB)
 
         DO 16 K=1,NF
            DO 10 J=1,NB
               NBIN(J)=0
               XB(J)=0.0
   10       CONTINUE
 
C     BEREKENING VAN DE BININDEX VOOR ELKE WAARNEMING
 
            DO 12 I=1,N
 
C      DE PHASE PH VOOR EEN TIJDSTIP T(I)
 
               PH=DMOD(T(I)*F,1D0)
 
C      HET BINNUMMER VOOR DE IDE WAARNEMING
 
               JB=INT(PH*NB)+1
 
C      HET AANTAL WAARNEMINGEN IN DE JBDE BIN
 
               NBIN(JB)=NBIN(JB)+1
 
C      DE SOM VAN DE WAARNEMINGEN IN DE JBDE BIN
 
               XB(JB)=XB(JB)+X(I)
 
   12       CONTINUE
 
            DFRE=DF1
            VM=XX
 
            DO 14 J=1,NB
               NBJ=NBIN(J)
 
C     DE SOM VAN DE VARIANTIES BINNEN DE NIET-LEDIGE BINS
 
               IF (NBJ.NE.0) THEN
                  VM=VM-XB(J)*XB(J)/NBJ
               ELSE
 
C     HET AANTAL WAARNEMINGEN MIN HET AANTAL NIET-LEDIGE BINS
 
                 DFRE=DFRE+1.0
               END IF
 
   14       CONTINUE
 
C     XVAR: DE SCHATTING VAN DE VARIANTIE, BEREKEND IN HET HOOFDPROGR.
 
C     OPVULLEN VAN DE MATRIX MET FREQUENTIE EN TESTSTATISTIEK
 
            F1(K) = F
            S1(K) = VM/DFRE/XVAR
            F=F+DF
             
   16    CONTINUE
 
      ELSE
 
C---------------------------------------------------------------------
C     NC>1 : DE METHODE VAN STELLINGWERF
C---------------------------------------------------------------------
 
C     OPMERKING : EEN BIN-COVER IS HET FASE-INTERVAL/(NB * NC)
C                 EEN BIN IS HET FASE-INTERVAL/NB
 
         NBC=NB*NC
         NC1=NC-1
         NBC1=NBC-NC1
         NBC2=NBC1+1
         DFC=DFLOAT(N*NC-NBC)
         XXNC=XX*NC

         F=F0
          
         DO 40 K=1,NF
            DO 20 J=1,NBC
               NBIN(J)=0
               XB(J)=0.0
   20       CONTINUE
 
            DO 22 I=1,N
 
C     DE FASE PH VOOR EEN TIJDSTIP T(I)
 
               PH=DMOD(T(I)*F,1D0)
 
C     HET BINNUMMER VAN DE IDE WAARNEMING
 
               JB=INT(PH*NBC)+1
 
C     HET AANTAL WAARNEMINGEN IN DE JBDE BIN-COVER
 
               NBIN(JB)=NBIN(JB)+1
 
C     DE SOM VAN HET AANTAL WAARNEMINGEN IN DE JBDE BIN-COVER
 
               XB(JB)=XB(JB)+X(I)
 
   22       CONTINUE
 
C     XH : MATRIX MET ALS ELEMENTEN DE SOM VAN DE WAARNEMINGEN IN
C          EENZELFDE BIN-COVER
C     NH : MATRIX MET ALS ELEMENTEN HET AANTAL WAARNEMINGEN IN EEN-
C          ZELFDE BIN-COVER
 
            DO 24 L=1,NC1
               XH(L)=XB(L)
               NH(L)=NBIN(L)
   24       CONTINUE
 
C    BEREKENING VOOR ELK VAN DE NC COVERS : HET AANTAL WAARNEMINGEN
C    EN DE SOM VAN DE WAARNEMINGEN PER BIN
 
            DO 28 J=1,NBC1
               J1=J+1
               JNC1=J+NC1
               DO 26 L=J1,JNC1
                  XB(J)=XB(J)+XB(L)
                  NBIN(J)=NBIN(J)+NBIN(L)
   26          CONTINUE
   28       CONTINUE
 
C      ALS ER AAN HET EINDE VAN HET FASE-INTERVAL ONVOLDOENDE BIN-
C      COVERS ZIJN VOOR EEN VOLLEDIGE BIN MOETEN ER VOORAAN NOG NC - J
C      BINCOVERS BIJGENOMEN WORDEN
 
            DO 34 J=NBC2,NBC
               J1=J+1
               DO 30 L=J1,NBC
                  XB(J)=XB(J)+XB(L)
                  NBIN(J)=NBIN(J)+NBIN(L)
   30          CONTINUE
               J2=J-NBC1
               DO 32 L=1,J2
                  XB(J)=XB(J)+XH(L)
                  NBIN(J)=NBIN(J)+NH(L)
   32          CONTINUE
   34       CONTINUE
 
            DFRE=DFC
            VM=XXNC
 
            DO 36 J=1,NBC
               NBJ=NBIN(J)
 
C     DE SOM VAN DE VARIANTIES IN DE NIET-LEDIGE BINS
 
               IF (NBJ.NE.0) THEN
                  VM=VM-XB(J)*XB(J)/NBJ
               ELSE
                  DFRE=DFRE+1.0
               END IF
 
   36       CONTINUE
 
C     OPVULLEN VAN DE MATRIX MET FREQUENTIE F EN TESTSTATISTIEK S
 
            F1(K) = F
            S1(K) = VM/DFRE/XVAR
            F=F+DF
   40    CONTINUE
 
      END IF
 
      RETURN
      END

      SUBROUTINE JUSTEL2(N,X,T,F0,NF,DF,NB,NC,XVAR,XX,DD,F1,S1)

C     BEREKENING VAN HET PERIODOGRAM VOLGENS JURKEVICH
C                                         OF STELLINGWERF
 
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),T(N),F1(NF),S1(NF)
      DIMENSION NBIN(1000),XB(1000),NH(1000),XH(1000)
Cf2py intent(in) N
Cf2py intent(in) X
Cf2py intent(in) T
Cf2py intent(in) F0
Cf2py intent(in) NF
Cf2py intent(in) DF
Cf2py intent(in) NB
Cf2py intent(in) NC
Cf2py intent(in) XVAR
Cf2py intent(in) XX
Cf2py intent(in) DD
Cf2py intent(in) F1
Cf2py intent(in) S1
Cf2py intent(out) F1
Cf2py intent(out) S1
      F=F0
       
C-------------------------------------------------------------------
C     NC = 1 : DE METHODE VAN JURKEVICH
C-------------------------------------------------------------------
 
      IF (NC.LE.1) THEN
 
         DF1=DFLOAT(N-NB)
 
         DO 16 K=1,NF
            DO 10 J=1,NB
               NBIN(J)=0
               XB(J)=0.0
   10       CONTINUE
 
C     BEREKENING VAN DE BININDEX VOOR ELKE WAARNEMING
 
            DO 12 I=1,N
 
C      DE PHASE PH VOOR EEN TIJDSTIP T(I)
 
               PH=DMOD(T(I)*F+DD/2D0*T(I)*T(I),1D0)
 
C      HET BINNUMMER VOOR DE IDE WAARNEMING
 
               JB=INT(PH*NB)+1
 
C      HET AANTAL WAARNEMINGEN IN DE JBDE BIN
 
               NBIN(JB)=NBIN(JB)+1
 
C      DE SOM VAN DE WAARNEMINGEN IN DE JBDE BIN
 
               XB(JB)=XB(JB)+X(I)
 
   12       CONTINUE
 
            DFRE=DF1
            VM=XX
 
            DO 14 J=1,NB
               NBJ=NBIN(J)
 
C     DE SOM VAN DE VARIANTIES BINNEN DE NIET-LEDIGE BINS
 
               IF (NBJ.NE.0) THEN
                  VM=VM-XB(J)*XB(J)/NBJ
               ELSE
 
C     HET AANTAL WAARNEMINGEN MIN HET AANTAL NIET-LEDIGE BINS
 
                 DFRE=DFRE+1.0
               END IF
 
   14       CONTINUE
 
C     XVAR: DE SCHATTING VAN DE VARIANTIE, BEREKEND IN HET HOOFDPROGR.
 
C     OPVULLEN VAN DE MATRIX MET FREQUENTIE EN TESTSTATISTIEK
 
            F1(K) = F
            S1(K) = VM/DFRE/XVAR
            F=F+DF
             
   16    CONTINUE
 
      ELSE
 
C---------------------------------------------------------------------
C     NC>1 : DE METHODE VAN STELLINGWERF
C---------------------------------------------------------------------
 
C     OPMERKING : EEN BIN-COVER IS HET FASE-INTERVAL/(NB * NC)
C                 EEN BIN IS HET FASE-INTERVAL/NB
 
         NBC=NB*NC
         NC1=NC-1
         NBC1=NBC-NC1
         NBC2=NBC1+1
         DFC=DFLOAT(N*NC-NBC)
         XXNC=XX*NC

         F=F0
          
         DO 40 K=1,NF
            DO 20 J=1,NBC
               NBIN(J)=0
               XB(J)=0.0
   20       CONTINUE
 
            DO 22 I=1,N
 
C     DE FASE PH VOOR EEN TIJDSTIP T(I)
 
               PH=DMOD(T(I)*F+DD/2D0*T(I)*T(I),1D0)
 
C     HET BINNUMMER VAN DE IDE WAARNEMING
 
               JB=INT(PH*NBC)+1
 
C     HET AANTAL WAARNEMINGEN IN DE JBDE BIN-COVER
 
               NBIN(JB)=NBIN(JB)+1
 
C     DE SOM VAN HET AANTAL WAARNEMINGEN IN DE JBDE BIN-COVER
 
               XB(JB)=XB(JB)+X(I)
 
   22       CONTINUE
 
C     XH : MATRIX MET ALS ELEMENTEN DE SOM VAN DE WAARNEMINGEN IN
C          EENZELFDE BIN-COVER
C     NH : MATRIX MET ALS ELEMENTEN HET AANTAL WAARNEMINGEN IN EEN-
C          ZELFDE BIN-COVER
 
            DO 24 L=1,NC1
               XH(L)=XB(L)
               NH(L)=NBIN(L)
   24       CONTINUE
 
C    BEREKENING VOOR ELK VAN DE NC COVERS : HET AANTAL WAARNEMINGEN
C    EN DE SOM VAN DE WAARNEMINGEN PER BIN
 
            DO 28 J=1,NBC1
               J1=J+1
               JNC1=J+NC1
               DO 26 L=J1,JNC1
                  XB(J)=XB(J)+XB(L)
                  NBIN(J)=NBIN(J)+NBIN(L)
   26          CONTINUE
   28       CONTINUE
 
C      ALS ER AAN HET EINDE VAN HET FASE-INTERVAL ONVOLDOENDE BIN-
C      COVERS ZIJN VOOR EEN VOLLEDIGE BIN MOETEN ER VOORAAN NOG NC - J
C      BINCOVERS BIJGENOMEN WORDEN
 
            DO 34 J=NBC2,NBC
               J1=J+1
               DO 30 L=J1,NBC
                  XB(J)=XB(J)+XB(L)
                  NBIN(J)=NBIN(J)+NBIN(L)
   30          CONTINUE
               J2=J-NBC1
               DO 32 L=1,J2
                  XB(J)=XB(J)+XH(L)
                  NBIN(J)=NBIN(J)+NH(L)
   32          CONTINUE
   34       CONTINUE
 
            DFRE=DFC
            VM=XXNC
 
            DO 36 J=1,NBC
               NBJ=NBIN(J)
 
C     DE SOM VAN DE VARIANTIES IN DE NIET-LEDIGE BINS
 
               IF (NBJ.NE.0) THEN
                  VM=VM-XB(J)*XB(J)/NBJ
               ELSE
                  DFRE=DFRE+1.0
               END IF
 
   36       CONTINUE
 
C     OPVULLEN VAN DE MATRIX MET FREQUENTIE F EN TESTSTATISTIEK S
 
            F1(K) = F
            S1(K) = VM/DFRE/XVAR
            F=F+DF
   40    CONTINUE
 
      END IF
 
      RETURN
      END


      SUBROUTINE SCAR3 (N,X,T,F0,NF,DF,F1,S1,SS,SC,SS2,SC2,W)

C     Computation of Scargles periodogram without explicit tau
C     calculation, with iteration (Method Cuypers)
C     Weighted version!
C     SC2 = R2, SS2 = 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(N),T(N)
      DIMENSION F1(NF),S1(NF)
      DIMENSION SS(NF),SC(NF),SS2(NF),SC2(NF), W(N)
      DATA TWOPI,DTWO,DNUL/6.28318530717959D0,2.0D0,0.0D0/
Cf2py intent(in) N
Cf2py intent(in) X
Cf2py intent(in) T
Cf2py intent(in) F0
Cf2py intent(in) NF
Cf2py intent(in) DF
Cf2py intent(in) F1
Cf2py intent(in) S1
Cf2py intent(in) SS
Cf2py intent(in) SC
Cf2py intent(in) SS2
Cf2py intent(in) SC2
Cf2py intent(in) W
Cf2py intent(out) S1
Cf2py intent(out) F1

      F = F0
      TPF = TWOPI * F
      TDF = TWOPI * DF
      TN = DFLOAT(N)
      TNSQ=TN*TN

      DO 70 K=1,NF
         SS(K)  = DNUL
         SC(K)  = DNUL
         SS2(K) = DNUL
         SC2(K) = DNUL
70       CONTINUE

      DO 72 I=1,N

         A = T(I)
         AF0 = DMOD(A*TPF,TWOPI)
         S0 = DSIN(AF0)
         C0 = DCOS(AF0)
         S20 = DTWO * S0 * C0
         C20 = C0 * C0 - S0 * S0

         ADF = DMOD(A*TDF,TWOPI)
         SDF = DSIN(ADF)
         CDF = DCOS(ADF)
         S2DF = DTWO * SDF * CDF
         C2DF = CDF * CDF - SDF * SDF
         XI=X(I)
         WI=W(I)
         C0X = C0 * XI
         S0X = S0 * XI
         DO 72 K=1,NF
            SS(K) = SS(K) +  WI * S0X
            SC(K) = SC(K) +  WI * C0X
            CTX=C0X
            C0X = CTX * CDF - S0X * SDF
            S0X = S0X * CDF + CTX * SDF
            SS2(K) = SS2(K) +  WI * S20
            SC2(K) = SC2(K) +  WI * C20
            C2T = C20
            C20 = C2T * C2DF - S20 * S2DF
            S20 = S20 * C2DF + C2T * S2DF
73          CONTINUE

72       CONTINUE 

      DO 74 K=1,NF
         SSK  =  SS(K)
         SS2K = SS2(K)
         SCK  =  SC(K)
         SC2K = SC2(K)
         
         F1(K)=F
         S1(K)=(SCK*SCK*(TN-SC2K)+SSK*SSK*(TN+SC2K)-DTWO*SSK*SCK*SS2K)/
     &      (TNSQ-SC2K*SC2K-SS2K*SS2K)
         F=F+DF
                  
74    CONTINUE
 
      RETURN
      END
      program dummy
      end
