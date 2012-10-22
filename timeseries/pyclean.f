      SUBROUTINE main_clean(t,bin,nt,hifreq,nf,gain,niter,
     &  nbins,startfreq,endfreq,f,wpow,wpha)
Cf2py intent(in) t
Cf2py intent(in) bin
Cf2py intent(in) nt
Cf2py intent(in) hifreq
Cf2py intent(in) nf
Cf2py intent(in) gain
Cf2py intent(in) niter
Cf2py intent(in) nbins
Cf2py intent(in) startfreq
Cf2py intent(in) endfreq
Cf2py intent(out) f
Cf2py intent(out) wpow
Cf2py intent(out) wpha


c t: times of observation
c bin: observed signal
c nt: number of observations
c hifreq: highest frequency
c nf: number of frequencies
c f: frequency domain for output final spectrum
c wpow: power spectrum for output final spectrum
c wpha: phase diagram, I assume
c niter: the number of CLEAN iterations
c gain: gain value

      implicit real*8 (a-h,o-z)
      dimension t(nt), bin(nt), f(nf), wpow(nf), wpha(nf)
      dimension startfreq(*), endfreq(*)

      dFREQ = hifreq/nf
c INITIALIZE CLEAN
 0010 call fourier(t,bin,nt,hifreq,nf,f,wpow,wpha,jmax,gain,niter,
     &     0,startfreq,endfreq,1)
C COMPUTE CLEAN
      call fourier(t,bin,nt,hifreq,nf,f,wpow,wpha,jmax,gain,niter,
     &     nbins,startfreq,endfreq,3)

      return
      end


      subroutine fourier(timesin,datain,ndata,hifreq,nfreq,freqout,
     &     power,phase,jmax,gain,niter,ncleanbins,startfreq,endfreq,
     &     mode)
C  Produces the dirty spectrum and the spectral window for a time
C  series contained in  timesin  and  datain.
c  The time average and the data mean are removed.

c Array POWER returns depending on MODE:
c Mode=0 :  power and phase diagram of window function
c Mode=1 :  power and phase diagram of window function and CLEAN init
c Mode=2 :  power and phase diagram of dirty spectrum
c Mode=3 :  power and phase diagram of CLEANed spectrum

      implicit real*8 (a-h,o-z)
      dimension timesin(ndata), datain(ndata)
      dimension freqout(nfreq), power(nfreq), phase(nfreq)
      dimension startfreq(*), endfreq(*)

c maxfiles: maxnr datapoints, maxfreqs: maxnr frequency bins
c maxfbins: maxnr clean bins
      parameter (maxfiles=1000000)
      parameter (maxfreqs=1000000)
      parameter (maxfbins=10)
      dimension  times(maxfiles), data(maxfiles)
      dimension  freq(0:2*maxfreqs), ones(maxfiles)
      dimension  isf(maxfbins), ief(maxfbins)
      complex*16 d(0:maxfreqs), w(0:2*maxfreqs)
      complex*16 dfour

      save freq, w, times, tmean


      dfreq= hifreq/nfreq         ! frequency element width
      twopi = 6.28318530717958
      pi=twopi/2.0d0

c Find the mean time and subtract it. Fill the ONES array
      if (mode.le.1) then
         tmean = 0.0d0
         do i=1,ndata
            tmean = tmean + timesin(i)
         enddo
         tmean = tmean / ndata
         do i=1,ndata
            times(i) = timesin(i) - tmean
         enddo
         do i=1,ndata
            ones(i)=1.0d0
         enddo
      endif

c Find the data mean and subtract it
      dmean = 0.0d0
      do i=1,ndata
         dmean = dmean + datain(i)
      enddo
      dmean = dmean / ndata
      do i=1,ndata
         data(i) = datain(i) - dmean
      enddo

c Calculate the Fourier transform
      do i= 0,nfreq
         if (mode.le.1) then
c            if (mod(i,200).eq.0) write(*,'(i5,$)') i
            freq(i)= dfreq*i
            w(i)= dfour(ndata,times,ones,freq(i),i,maxfreqs)
         else
            d(i)= dfour(ndata,times,data,freq(i),i,maxfreqs)
         endif
      enddo

c Complete the spectral window
      if (mode.eq.1) then
         do i= nfreq+1,2*nfreq
c            if (mod(i,200).eq.0) write(*,'(i5,$)') i
            freq(i)= dfreq*i
            w(i)= dfour(ndata,times,ones,freq(i),i,maxfreqs)
         enddo
         call clean(niter,gain,w,d,nfreq,ncleanbins,isf,ief,0)
      else
         call findrange(nfreq,freq,ncleanbins,startfreq,endfreq,isf,ief)
         call clean(niter,gain,w,d,nfreq,ncleanbins,isf,ief,1)
      endif

c Return power or phase spectrum
      powmax=0.0
      do i=1,nfreq
         if (mode.le.1) then
c            a=real(w(i))
c            b=imag(w(i))
            a=dble(w(i))
            b=dimag(w(i))
         elseif (mode.ge.2) then
            a=dble(d(i))
            b=dimag(d(i))
         endif

         if (a.eq.0.0) then
            phi=0.0
         else
            phi=atan(b/a)
         endif
         if ((a.lt.0) .and. (b.gt.0)) phi=phi+pi
         if ((a.lt.0) .and. (b.lt.0)) phi=phi-pi
         if (phi.lt.0) phi=phi+twopi

         power(i)=a**2 + b**2   ! spectral window
         phase(i)=phi           ! phase diagram 0 - 2pi

c NEW!! ******
c Convert phases such that they are the same as for
c  fitting sinusoids, with phase relative to t=0.
c        phase(i)=dmod((phi + pi/2 - twopi*tmean*freq(i)),twopi)

         if (mode.le.1) then
            freqout(i)=freq(i)
         elseif (mode.ge.2) then
            if (power(i).gt.powmax) then
               powmax=power(i)
               jmax=i
            endif
         endif
      enddo
      return
      end


      subroutine findrange(nf,freq,ncleanbins,startfreq,endfreq,isf,ief)
      implicit real*8 (a-h,o-z)
      dimension freq(0:nf), startfreq(*), endfreq(*), isf(*), ief(*)

      ndata = 0
      if (ncleanbins.eq.0) return

      do j=1,ncleanbins
         isf(j)=0
         ief(j)=0
      enddo
      
      do i=1,nf
       do j=1,ncleanbins
        if ((freq(i).ge.startfreq(j)).and.(freq(i).le.endfreq(j)))then
           ndata = ndata + 1
           if (i.eq.1) then
              isf(j)=1
           elseif (freq(i-1).lt.startfreq(j)) then
              isf(j)=i
           endif

           if (i.eq.nf) then
              ief(j)=nf
           elseif (freq(i+1).gt.endfreq(j)) then
              ief(j)=i
           endif
        endif
       enddo
      enddo

      if (ndata.eq.0) stop 'No points in Clean range ...'
      do j=1,ncleanbins
         if ((isf(j).eq.0) .and .(ief(j).eq.0))
     &        stop 'No points in one or more Clean range(s) ...'
      enddo

      end




      subroutine clean(ncl,gain,w,s,nfreq,ncleanbins,isf,ief,mode)
C------------------------------------------------------------------------------
C  Deconvolves the spectral window in WFILE from the dirty spectrum in DFILE,
C  by using an iterative, one-dimensional, Hoegbom CLEAN algorithm.
C  The resulting clean spectrum is in SFILE, the residual and CLEAN component
C  spectra are found in RFILE & CFILE, respectively.  The input parameters
C  include the filenames, number of CLEANs (iterations), and the gain per iter.
C     To produce the clean spectrum, a Gaussian beam is used which has been
C  fitted to the HWHM of the primary peak in the spectral window; the phases
C  are determined by the offset of the mean time TMEAN from the "origin" TZERO.
C    Since all spectra from real data are Hermitian, we define only the non-
C  negative frequencies.  The negative frequency elements are recovered by
C  the use of the function CVAL, which returns the complex conjugate for
C  negative frequencies, and zero for frequencies outside the defined range.
C------------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c     maxfreqs maximum # ft bins
      parameter (maxfreqs=1000000)
      complex*16 w(0:2*maxfreqs), r(0:maxfreqs)
      complex*16 b(0:maxfreqs), c(0:maxfreqs), s(0:maxfreqs), cc, alpha1
      dimension isf(*), ief(*)
      integer fillb

      save b, mb

      ms=nfreq
      mw=2*ms

C  fit a restoring beam B(0:MB) to the spectral window peak
      if (mode.eq.0) then
         dPHASE= 0.0d0
         HWIDTH= HWHM(MW,W)                  ! HWHM of W (index units)
         MB= FILLB(maxfreqs,B,HWIDTH,dPHASE) ! fill the restoring beam
         return
      endif

      do i=0,maxfreqs
         c(i) = (0.0d0, 0.0d0)
         r(i) = s(i)
      enddo

C  CLEAN the residual spectrum, storing the components in C(0:MS)
c     WRITE(*,'(A,$)') 'CLEANing ... '
      DO ICL=1,NCL
c        L= MAXLOC(MS,R)           ! element location of max. peak in R
         L= johnMAXLOC(MS,R,ncleanbins,isf,ief)
c! estimate the component to remove
         CC= GAIN*ALPHA1(L,MS,R,W) 
c! subtract the component from R
         CALL SUBCMP(MS,R,W,L,CC)
c! store the component in C
         C(L)= C(L) + CC
      ENDDO

C  Generate the clean spectrum S(0:MS)
c ! convolve C with B to form S
      CALL CONVOLV(MS,S,MS,C,MB,B)
c ! then add R to S
      CALL ADD(MS,S,MS,S,MS,R)    

c     WRITE(*,'(A,$)') 'done !'
      END


      INTEGER FUNCTION johnMAXLOC(M,ARRAY,ncleanbins,isf,ief)
C-----------------------------------------------------------------------------
C  returns the array index of the ARRAY(0:M) element with max. CDABS value.
C-----------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      COMPLEX*16 ARRAY(0:M)
      dimension isf(*), ief(*)

C  find the max. location
      if (ncleanbins.eq.0) then
         johnMAXLOC=MAXLOC2(m,array)
         RETURN
      endif

      LMAX= 0
      AMAX= CDABS(ARRAY(0))
      do j=1,ncleanbins
         DO I=isf(j),ief(j)
            ARRAYI= CDABS(ARRAY(I))
            if (ARRAYI.GT.AMAX) THEN
               AMAX= ARRAYI
               LMAX= I
            ENDIF
         ENDDO
      enddo
      johnMAXLOC= LMAX
      RETURN
      END


      SUBROUTINE CONVOLV(M,A,M1,A1,M2,A2)
C-----------------------------------------------------------------------------
C  Convolves complex arrays A1(0:M1) and A2(0:M2) to form A(0:M).
C-----------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      COMPLEX*16 A(0:M), A1(0:M1), A2(0:M2), CVAL
 
C  Convolve A1 with A2 to form A
      DO J=0,M
c! reset A(J)
         A(J)= (0.0d0,0.0d0)          
         DO I=-M2,M2
            A(J)= A(J)+CVAL(A2,I,M2)*CVAL(A1,J-I,M1)
         ENDDO
      ENDDO
      RETURN
      END


      SUBROUTINE ADD(M,A,M1,A1,M2,A2)
C-----------------------------------------------------------------------------
C  Adds complex arrays A1(0:M1) and A2(0:M2) to form A(0:M).
C-----------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      COMPLEX*16 A(0:M), A1(0:M1), A2(0:M2), CVAL
 
C  Add A1 to A2, forming A
      DO J=0,M
         A(J)= CVAL(A1,J,M1) + CVAL(A2,J,M2)
      ENDDO
      RETURN
      END


      REAL*8 FUNCTION HWHM(M,ARRAY)
C-----------------------------------------------------------------------------
C  Finds the half-width,half-maximum of the CDABS of ARRAY(0:M)
C  This is done by linear interpolation between successive elements of ARRAY
C  Returns the HWHM in the units of the array elements.
C-----------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      COMPLEX*16 ARRAY(0:M)
      REAL*8 HALFMX, HALFWD, LAST, CURRENT, XLAST, SLOPE

C  Initial conditions for search
      HALFMX= CDABS(ARRAY(0))/2.0d0 ! half-maximum (assume max at I=0)
      HALFWD= 0.0d0                 ! half-width not yet found
 
C  loop through until less than HALFMX
      DO I=1,M
c ! current array value
         CURRENT= CDABS(ARRAY(I))            
         if (CURRENT.LT.HALFMX) THEN
c last array value // last x coordinate // inverse slope between them
            LAST= CDABS(ARRAY(I-1))           
            XLAST= dFLOAT(I-1)                
            SLOPE= 1.0d0/(CURRENT-LAST)       
c interpolate to halfmx // pop out loop
            HALFWD= XLAST+SLOPE*(HALFMX-LAST)
            GOTO 10                          
         ENDIF
      ENDDO

C  return with the result
 0010 if (HALFWD.LE.0.0d0) WRITE(*,*) '*NOT FIND HALF-WIDTH,HALF-MAX*'
      HWHM= HALFWD
      RETURN
      END


      integer FUNCTION FILLB(MB,B,HWIDTH,PINCR)
C-----------------------------------------------------------------------------
C  Fills the restoring beam B(0:MB) with a gaussian of half-width HWIDTH
C  (in Index units).  The complex phase is zero for I=0 and increases
C  by PINCR for each element.  Returns the max. filled element.
C-----------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      COMPLEX*16 B(0:MB)
 
C  Set up some initial conditions
      DO I=0,MB
         B(I)= (0.0d0,0.0d0)
      ENDDO

C  if Gaussian, fit a Gaussian beam
C    ---Calculate SIGMA and the normalization constant
c sigma in element units // normalizn constant // maximum filled element
      SIGMA= HWIDTH/SQRT(2.0d0*LOG(2.0d0))
      CONST= 1.0d0/(2.0d0*SIGMA*SIGMA)    
      MFILL= INT(5*SIGMA)+1               
      if (MFILL.GT.MB) MFILL=MB
C    ---Fill B with the gaussian
      DO I=0,MFILL
         X= dFLOAT(I)
         GAUSS= EXP(-CONST*X*X)
         B(I)= dCMPLX(GAUSS,0.0d0)
      ENDDO
C    ---Include the phase information
      DO I=1,MFILL
         PHASEI= dFLOAT(I)*PINCR
         B(I)= B(I)*dCMPLX(COS(PHASEI),SIN(PHASEI))
      ENDDO

C  return the maximum filled element
      FILLB= MFILL
      RETURN
      END

      complex*16 function dfour(n,times,data,freq,ifreq,np)
C----------------------------------------------------------------------------
C  Returns the Fourier transform of the time series specified by TIMES(1:N)
C  and DATA(1:N), evaluated at the frequency FREQ.
C  The FT is normalized to have the data mean at DC
C----------------------------------------------------------------------------

      implicit real*8 (a-h,o-z) 
      dimension data(n), times(n)
      complex*16 ft

      twopi = 6.28318530717958

C  Evaluate FT at FREQ...
      FT= (0.d0, 0.d0)
      DO K=1,N
         PHASE= -TWOPI*FREQ*TIMES(K)
         COSPHASE = COS(PHASE)
         SINPHASE = SIN(PHASE)
         FT= FT + DATA(K)*DCMPLX(COSPHASE,SINPHASE)
      ENDDO
      
C  return with FT properly normalized
      DFOUR= FT/N
      RETURN
      END

      COMPLEX*16 FUNCTION ALPHA1(L,M,SPEC,WIND)
C-----------------------------------------------------------------------------
C  Returns an estimate for the component A, responsible for SPEC(L) through
C  the relation:
C                  SPEC(L) = A + (A*)WIND(2L)
C  SPEC is defined (0:M), and WIND is defined (0:2M)
C-----------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      COMPLEX*16 SPEC(0:M), WIND(0:2*M), WIN2L, CVAL

c allowed error in wnorm
      ERROR=0.0001d0  
 
C  Find the (L,-L) components which produce SPEC(L) through WIND
c                               ! (L,-L) interference
      WIN2L= WIND(2*L) 
      WNORM= 1.0d0-CDABS(WIN2L)**2
c ! avoid singularities
      if (WNORM.LT.ERROR) ALPHA1= 0.5d0*SPEC(L) 
      if (WNORM.GE.ERROR) ALPHA1= (SPEC(L)-WIN2L*CVAL(SPEC,-L,M))/WNORM
 
C  return with the estimate of A
      RETURN
      END


      SUBROUTINE SUBCMP(M,SPEC,WIND,L,CPOS)
C-----------------------------------------------------------------------------
C     Subtracts the complex component CPOS at array location L, and
C  its complex conjugate, CNEG, at -L from the spectrum SPEC(0:M).
C  The spectral window WIND(0:2M) is matched to the component and
C  subtracted from the entire spectral array.
C-----------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      COMPLEX*16 SPEC(0:M), WIND(0:2*M), CPOS, CNEG, CVAL
 
C  specify the -L component
      CNEG= CONJG(CPOS)
 
C  Remove effects of both +L and -L components !     (from L)      (from -L)
      DO I=0,M
         SPEC(I)= SPEC(I)-CPOS*CVAL(WIND,I-L,2*M)-CNEG*WIND(I+L)
      ENDDO  

C  return to the program
      RETURN
      END


      COMPLEX*16 FUNCTION CVAL(ARRAY,I,MARR)
C----------------------------------------------------------------------------
C  Returns the "value" of the complex ARRAY at the "location" I
C  If I > 0 then it is equivalent to ARRAY(I), but if I < 0,
C  it returns the complex conjugate of ARRAY(-I).  ARRAY must be
C  indexed from 0 to at least ABS(I) and COMPLEX.
C  if ABS(I) > MARR (max elt in array), CVAL=0.
C----------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      COMPLEX*16 ARRAY(0:MARR)

C  set value to conjugate of absolute location if necessary
      LOC= ABS(I)
      if (LOC.GT.MARR)THEN
         CVAL= (0.0d0,0.0d0)
         RETURN
      ENDIF
      if (I.GE.0) CVAL= ARRAY(LOC)
      if (I.LT.0) CVAL= CONJG(ARRAY(LOC))
      RETURN
      END

      INTEGER FUNCTION MAXLOC2(M,ARRAY)
C-----------------------------------------------------------------------------
C  returns the array index of the ARRAY(0:M) element with max. CDABS value.
C-----------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      COMPLEX*16 ARRAY(0:M)
 
C  find the max. location
      LMAX= 0
      AMAX= CDABS(ARRAY(0))
      DO I=1,M
         ARRAYI= CDABS(ARRAY(I))
         if (ARRAYI.GT.AMAX) THEN
            AMAX= ARRAYI
            LMAX= I
         ENDIF
      ENDDO
      MAXLOC= LMAX
      RETURN
      END
      program dummy
      end
