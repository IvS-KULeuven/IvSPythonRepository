      subroutine deeming1(T,X,F0,DF,N,NF,FR,S1)
Cf2py intent(in) T
Cf2py intent(in) X
Cf2py intent(in) F0
Cf2py intent(in) DF
Cf2py intent(in) N
Cf2py intent(in) NF
Cf2py intent(out) FR
Cf2py intent(out) S1
      implicit real*8 (a-h,o-z)
      dimension T(N),X(N)
      dimension FR(NF),S1(NF)
      DOUBLE PRECISION TPF,F0,DF
      DOUBLE PRECISION COX,SOX,TF
      DOUBLE PRECISION S0,C0
      DATA TWOPI/6.28318530717959D0/
      
      TPF = F0*TWOPI
      do 11 K=1,NF
      C0X = 0.
      S0X = 0.
      do 12 I=1,N
      TF = T(I)*TPF	 
      S0 = SIN(TF)
      C0 = COS(TF)
      C0X = C0X + C0 * X(I)
      S0X = S0X + S0 * X(I)
12    end do
       
      S1(K) = C0X*C0X + S0X*S0X
      FR(K) = TPF/TWOPI
      TPF=TPF+DF*TWOPI
11    end do
      RETURN
      END


      subroutine deeming2(T,X,F0,DF,N,NF,FR,S1)
Cf2py intent(in) T
Cf2py intent(in) X
Cf2py intent(in) F0
Cf2py intent(in) DF
Cf2py intent(in) N
Cf2py intent(in) NF
Cf2py intent(out) FR
Cf2py intent(out) S1
      implicit real*8 (a-h,o-z)
      dimension T(N),X(N)
      dimension FR(NF),S1(NF)
      dimension SS(NF),SC(NF)
      DOUBLE PRECISION TPF,F0,DF,TDF
      DOUBLE PRECISION COX,SOX,TF
      DOUBLE PRECISION S0,C0
      DOUBLE PRECISION A,AF0,ADF
      DOUBLE PRECISION SDF,CDF,XI,CTX
      DOUBLE PRECISION SSK,SCK
      DATA TWOPI/6.28318530717959D0/

      TPF = F0*TWOPI
      TDF = DF*TWOPI
      
      do 13 K=1,NF
         SS(K)  = 0.
         SC(K)  = 0.
13    end do

      do 14 I=1,N
         A = T(I)
         
         AF0 = A*TPF
         
         S0 = SIN(AF0)
         C0 = COS(AF0)

         ADF = A*TDF
         
         SDF = SIN(ADF)
         CDF = COS(ADF)

         XI=X(I)
         C0X = C0 * XI
         S0X = S0 * XI
         
         do 15, K=1,NF
            SS(K) = SS(K) + S0X
            SC(K) = SC(K) + C0X
            CTX=C0X
            C0X = CTX * CDF - S0X * SDF
            S0X = S0X * CDF + CTX * SDF
15       end do
14    end do
      TPF = F0
      do K=1,NF
         SSK  =  SS(K)
         SCK  =  SC(K)
         
         S1(K)=SCK*SCK+SSK*SSK
         FR(K)=TPF
         TPF = TPF + DF
      end do
      RETURN
      END
