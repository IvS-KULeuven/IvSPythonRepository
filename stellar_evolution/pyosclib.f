c   osclib.f   2006-09-01

c   This file is part of the source code of OSC version 37. It consists of the
c   subroutines of OSC that you can link to the standard user interface (source
c   code in osc.f) or to your own Fortran code.

c   The gravitational constant is defined in the procedures readmod and osc06a
c   and the solar mass is defined in osc06a with the best values known by 2006:
c   G=6.6742d-8 and Msol=1.9884d33.

c===============================================================================

      blockdata oscdata
      implicit none
      include 'osc.i'
      data vers/37/
      data lmod,lres,losc,lwriteosc,lscreen/5*.false./
      data ll/0/,boundary/1/
      end

c===============================================================================

      subroutine version(v)
Cf2py intent(in) v
Cf2py intent(out) v
      implicit none
      include 'osc.i'
      integer v
  901 format('OSC   version ',i3)
      v=vers
      if(lscreen)write(6,901)vers
      return
      end

c===============================================================================

      subroutine setscreen(i)
      implicit none
      integer i
      include 'osc.i'
      lscreen=i.ne.0
      return
      end

c===============================================================================

c   setboundary
c   
c   This procedure sets the type of surface boundary condition to be applied to
c   delta P
c   --- in
c   b = 1 for our usual boundary condition derived in the case of a model with
c         zero surface pressure (the default value)
c       2 to impose delta P = 0 at the last point of the model
c   --- out
c   rc = 0 if the boundary condition has been successfully changed
c        1 if an illegal choice has been made ; in this case, the previously
c          chosen boundary condition is retained

      subroutine setboundary(b,rc)
      implicit none
      include 'osc.i'
      integer b,rc
  901 format('*** error: illegal value ',i3,
     .' for boundary condition')
      if(b.eq.1.or.b.eq.2)then
        boundary=b
        rc=0
      else
        rc=1
        if(lscreen)write(6,901)b
      endif
      return
      end

c===============================================================================

c     readmod
c     OSC06
c     Lecture du mod�e,
c
      subroutine readmod(modname,ir,dstc,dsts,st,out_r,out_m)
Cf2py intent(in) modname
Cf2py intent(out) ir
Cf2py intent(out) dstc
Cf2py intent(out) dsts
Cf2py intent(out) st
Cf2py intent(out) out_r
Cf2py intent(out) out_m
c     subroutine osc06(modname,ir)
c     from input file: descript,imax,radius,amass
c     we need: imaxmod,am,ray,gmsr3,tdyn,rhom,c3,c4
c     x,q,st,lmod,lres,ires,nres
      implicit double precision(a-h,o-z)
      real out_r,out_m
      parameter(G=6.6742d-8,damsol=1.9884d33,pi=3.14159265358979d0)
      include 'osc.i'
      dimension dstc(18),dsts(24),st(6,np0max)
      character modname*(*)
  904 format('*** error : imax =',i6,' too large')
  905 format('*** error opening ',a)
  906 format('*** error reading ',a)
      lmod=.false.
      lres=.false.
      open(11,file=modname,form='unformatted',status='old',err=510)
      read(11,err=500)descript
      read(11)imax,radius,amass
      write(6,*)'descript'
      write(6,*)descript
      write(6,*)'imax,radius,amass'
      write(6,*)imax,radius,amass
      out_r = radius
      out_m = amass
      if(imax.gt.np0max)then
        if(lscreen)write(6,904)imax
        ir=1
        return
      endif
      imaxmod=imax
      do 220 i=1,imax
  220 read(11,err=500)(st(k,i),k=1,5)
      close(11)
      am=amass
      ray=radius
      if(radius.eq.1d0.and.amass.eq.1d0)then
        gmsr3=1d0
      else
        gmsr3=G*amass/radius**3
      endif
      tdyn=1d0/sqrt(gmsr3)
      rhom=3d0*amass/4/pi/radius**3
c.... Si le point en surface a une pression nulle, on le supprime
      if(st(3,imax).eq.0d0)imax=imax-1
      c3=radius/G/amass
      c4=4d0*pi*radius**3/amass
      do i=1,imax
        x=st(1,i)/radius
        q=st(2,i)/amass
        st(6,i)=0d0
        st(3,i)=c3*st(3,i)/st(4,i)
        st(4,i)=c4*st(4,i)
        st(1,i)=x
        if(x.eq.0d0)then
          st(2,i)=st(4,i)/3d0
        else
          st(2,i)=q/x**3
        endif
      enddo
      call osc06a(imax,st,ir,G,damsol)
      if(ir.ne.0)return
      lmod=.true.
      lres=.true.
      ires=0
      nres=0
c     call osc26(ires,nres,ir)
      call mesh(ires,nres,ir)
      return
  500 ir=2
      close(11)
      l=len_trim(modname)
      if(lscreen)write(6,906)modname(1:l)
      return
  510 ir=3
      l=len_trim(modname)
      if(lscreen)write(6,905)modname(1:l)
      end


c
c     PYTHON READMOD
c
      subroutine readmodpy(st,ir,imax_,radius_,amass_,modname)
c     subroutine osc06(modname,ir)
c     from input file: descript,imax,radius,amass
c     we need: imaxmod,am,ray,gmsr3,tdyn,rhom,c3,c4
c     x,q,st,lmod,lres,ires,nres
      implicit double precision(a-h,o-z)
      parameter(G=6.6742d-8,damsol=1.9884d33,pi=3.14159265358979d0)
      include 'osc.i'
      dimension dstc(18),dsts(24),st(6,np0max)
      character modname*(*)
  904 format('*** error : imax =',i6,' too large')
  905 format('*** error opening ',a)
  906 format('*** error reading ',a)
      lmod=.false.
      lres=.false.
      imax = imax_
      radius = radius_
      amass = amass_
      write(6,*)'reading python input1',imax,radius,amass
      if(imax.gt.np0max)then
        if(lscreen)write(6,904)imax
        ir=1
        return
      endif
      imaxmod=imax
      am=amass
      ray=radius
      if(radius.eq.1d0.and.amass.eq.1d0)then
        gmsr3=1d0
      else
        gmsr3=G*amass/radius**3
      endif
      write(6,*)'reading python input1',imaxmod,ray,am
      tdyn=1d0/sqrt(gmsr3)
      rhom=3d0*amass/4/pi/radius**3
c.... Si le point en surface a une pression nulle, on le supprime
      if(st(3,imax).eq.0d0)imax=imax-1
      c3=radius/G/amass
      c4=4d0*pi*radius**3/amass
c      write(6,*)'reading python input1'
      call osc06a(imax,st,ir,G,damsol)
      if(ir.ne.0)return
      lmod=.true.
      lres=.true.
      ires=0
      nres=0
c      call osc26(ires,nres,ir)
      call mesh(ires,nres,ir)
c      write(6,*)'reading python input1'
      return
  500 ir=2
      l=len_trim(modname)
      if(lscreen)write(6,906)modname(1:l)
      return
  510 ir=3
      l=len_trim(modname)
      if(lscreen)write(6,905)modname(1:l)
      end


c
c     PYTHON READMOD WITH FREE GRAV
c
      subroutine readmodpyg(st,ir,imax_,radius_,amass_,modname,G,damsol)
c     subroutine osc06(modname,ir)
c     from input file: descript,imax,radius,amass
c     we need: imaxmod,am,ray,gmsr3,tdyn,rhom,c3,c4
c     x,q,st,lmod,lres,ires,nres
      implicit double precision(a-h,o-z)
      parameter(pi=3.14159265358979d0)
      include 'osc.i'
      dimension dstc(18),dsts(24),st(6,np0max)
      character modname*(*)
  904 format('*** error : imax =',i6,' too large')
  905 format('*** error opening ',a)
  906 format('*** error reading ',a)
      lmod=.false.
      lres=.false.
      imax = imax_
      radius = radius_
      amass = amass_
c      write(6,*)'reading python input1',imax,radius,amass
      if(imax.gt.np0max)then
        if(lscreen)write(6,904)imax
        ir=1
        return
      endif
      imaxmod=imax
      am=amass
      ray=radius
      if(radius.eq.1d0.and.amass.eq.1d0)then
        gmsr3=1d0
      else
        gmsr3=G*amass/radius**3
      endif
c      write(6,*)'reading python input1',imaxmod,ray,am
      tdyn=1d0/sqrt(gmsr3)
      rhom=3d0*amass/4/pi/radius**3
c.... Si le point en surface a une pression nulle, on le supprime
      if(st(3,imax).eq.0d0)imax=imax-1
      c3=radius/G/amass
      c4=4d0*pi*radius**3/amass
c      write(6,*)'reading python input1'
      call osc06a(imax,st,ir,G,damsol)
      if(ir.ne.0)return
      lmod=.true.
      lres=.true.
      ires=0
      nres=0
c      call osc26(ires,nres,ir)
      call mesh(ires,nres,ir)
c      write(6,*)'reading python input1'
      return
  500 ir=2
      l=len_trim(modname)
      if(lscreen)write(6,906)modname(1:l)
      return
  510 ir=3
      l=len_trim(modname)
      if(lscreen)write(6,905)modname(1:l)
      end

      subroutine getmodel(st)
Cf2py intent(in) st
Cf2py intent(out) st
c     Return model currently stored
      implicit double precision(a-h,o-z)
      include 'osc.i'     
      dimension st(6,np0max)
      j = 0
      k = 0
      do k=1,6
        do j=1,np0max
          st(k,j) = star(k,j)
        enddo
      enddo
      return
      end

c
c     OSC06A
c
      subroutine osc06a(imax,st,ir,g,damsol)
      implicit double precision(a-h,o-z)
      include 'osc.i'
      dimension st(6,*)
      parameter(d1max=5d-3,d4max=1d-1)
      dimension astar(6,np0max),dstar(6,np0max)
      integer istar(np0max)
  904 format(' *** error : x(j) < x(j-1) for j = ',i4)
  906 format(' *** error : too many points')
  907 format(' *** error : null pressure at the surface')
      if(st(3,imax).eq.0.)then
      if(lscreen)write(6,907)
      ir=4
      return
      endif
c     reset astar first!!! <--- modification by Pieter Degroote
      do 998 j=1,np0max
      do 999 l=1,6
      astar(l,j) = 0d0
      dstar(l,j) = 0d0
      star(l,j) = 0d0
  999 continue
  998 continue
 
      

c   imax = nombre de points du mod�e
c   st(1)=x, st(2)=q/x**3, st(3)=RP/GMrho, st(4)=4piR3rho/M, st(5)=Gamma1
c   Dans la section qui suit, on �imine les points de multiplicit�sup�ieure �c   2 et le mod�e est recopi�dans astar. Les variables 4 et 5 sont stock�s
c   sous leurs formes logarithmiques. Le nombre de points du mod�e est alors
c   kmax
      k=0
      do 140 j=1,imax
      if(j.eq.1)goto 110
      if(st(1,j).lt.st(1,j-1))goto 810
      if(st(1,j).ne.st(1,j-1))goto 120
  110 if(j.eq.imax)goto 140
      if(st(1,j).eq.st(1,j+1))goto 140
  120 k=k+1
      do 130 l=1,5
      if(l.ne.3.and.l.ne.4)astar(l,k)=st(l,j)
  130 if(l.eq.3.or.l.eq.4)astar(l,k)=log(st(l,j))
      if(k.eq.1)astar(2,1)=st(4,1)/3.d0
      astar(6,k)=st(2,j)/st(3,j)/st(5,j)
  140 continue
      kmax=k

c   On note dans le vecteur istar le statut des points du mod�e. En d�inissant
c   un zone comme l'ensemble des points entre deux points doubles (le premier
c   point et le dernier point du mod�e sont �alement des limites de zone), on
c   attribue un indice -1 �la limite inf�ieure d'une zone, +1 �sa limite
c   sup�ieure et 0 �un point interne.
      do 230 k=1,kmax
      if(k.eq.1)goto 210
      if(astar(1,k).eq.astar(1,k-1))goto 210
      if(k.eq.kmax)goto 220
      if(astar(1,k).eq.astar(1,k+1))goto 220
      istar(k)=0
      goto 230
  210 istar(k)=-1
      goto 230
  220 istar(k)=1
  230 continue

c   On calcule les d�iv�s de astar par rapport �x. Elles sont stock�s dans
c   dstar.

      do k=1,kmax
        dstar(1,k)=1d0
        if(astar(1,k).eq.0d0)then
          do l=2,5
            dstar(l,k)=0d0
          enddo
        else
          k1=k-1
          if(istar(k).eq.-1)k1=k+2
          k2=k+1
          if(istar(k).eq.1)k2=k-2
          dx1=astar(1,k)-astar(1,k1)
          dx2=astar(1,k2)-astar(1,k)
          a1=dx2/(dx1+dx2)
          a2=1d0-a1
          do l=2,5
            dy1=(astar(l,k)-astar(l,k1))/dx1
            dy2=(astar(l,k2)-astar(l,k))/dx2
            dstar(l,k)=a1*dy1+a2*dy2
          enddo
        endif
      enddo

c   On calcule �pr�ent astar(6) = RA/x (A=dlnrho/dr-dlnP/dr/Gamma1) et sa
c   d�iv� dstar(6)

      do k=1,kmax
        if(astar(1,k).eq.0d0)cycle
        astar(6,k)=astar(6,k)+dstar(4,k)/astar(1,k)
      enddo
      if(astar(1,1).eq.0d0)then
        astar(6,1)=(astar(1,3)**2*astar(6,2)-astar(1,2)**2*astar(6,3))
     .  /(astar(1,3)**2-astar(1,2)**2)
      endif

      do k=1,kmax
        if(astar(1,k).eq.0d0)then
          dstar(6,k)=0d0
        else
          k1=k-1
          if(istar(k).eq.-1)k1=k+2
          k2=k+1
          if(istar(k).eq.1)k2=k-2
          dx1=astar(1,k)-astar(1,k1)
          dx2=astar(1,k2)-astar(1,k)
          a1=dx2/(dx1+dx2)
          a2=1d0-a1
          dstar(6,k)=a1*dy1+a2*dy2
        endif
      enddo

c   On ajoute des points l�o c'est vraiment n�essaire et on revient aux
c   variables sans logarithmes. Le mod�e est alors d�rit par la variable star0
c   et il poss�e jmax points.
c   star0(1)=x, star0(2)=q/x**3, star0(3)=RP/GMrho, star0(4)=4piR3rho/M,
c   star0(5)=Gamma1 et star0(6)=RA/x
      i=0
      do 540 k=1,kmax
      if(k.eq.1)goto 520
      k1=k-1
      if(astar(1,k).eq.astar(1,k1))goto 520
      ni1=dabs(astar(1,k)-astar(1,k1))/d1max+1
      ni4=dabs(astar(4,k)-astar(4,k1))/d4max+1
      ni=max0(ni1,ni4)
      if(ni.lt.2)goto 520
      x=astar(1,k)
      x1=astar(1,k1)
      dx=x-x1
      do 510 m=2,ni
      if(i.ge.np0max)goto 830
      i=i+1
      xi=dfloat(m-1)/dfloat(ni)
      xx=x1+xi*dx
      star0(1,i)=xx
      star0(7,i)=1.d0
      phi0=(1.d0-xi)**2*(1.d0+2.d0*xi)
      psi0=(3.d0-2.d0*xi)*xi**2
      phi1=xi*(1.d0-xi)**2*dx
      psi1=-xi**2*(1.d0-xi)*dx
      dphi0=-6.d0*xi*(1.d0-xi)/dx
      dpsi0=6.d0*xi*(1.d0-xi)/dx
      dphi1=(1.d0-xi)*(1.d0-3.d0*xi)
      dpsi1=xi*(3.d0*xi-2.d0)
      do 510 l=2,6
      yy=phi0*astar(l,k1)+psi0*astar(l,k)
     *+phi1*dstar(l,k1)+psi1*dstar(l,k)
      dyy=dphi0*astar(l,k1)+dpsi0*astar(l,k)
     *+dphi1*dstar(l,k1)+dpsi1*dstar(l,k)
      if(l.eq.3.or.l.eq.4)then
      star0(l,i)=dexp(yy)
      star0(l+6,i)=dyy*star0(l,i)
      else
      star0(l,i)=yy
      star0(l+6,i)=dyy
      endif
  510 continue
  520 if(i.ge.np0max)goto 830
      i=i+1
      do 530 l=1,6
      if(l.eq.3.or.l.eq.4)then
      star0(l,i)=dexp(astar(l,k))
      star0(l+6,i)=dstar(l,k)*star0(l,i)
      else
      star0(l,i)=astar(l,k)
      star0(l+6,i)=dstar(l,k)
      endif
  530 continue
  540 continue
      jmax=i

c   Le vecteur st1 donne les coefficients du d�eloppement des variables du
c   mod�e en s�ie de puissances de x au premier point et jusqu'�l'ordre 2. Si
c   le premier point n'est pas au centre on n'a pas besoin des coefficients
c   d'ordre 2.
      do 610 l=1,12
  610 st1(l)=star0(l,1)
      st1(13)=0.d0
      if(star0(1,1).eq.0)then
        x2=star0(1,2)**2
        do l=14,18
          m=l-12
          st1(l)=(star0(m,2)-star0(m,1))/x2
        enddo
      else
        do l=14,18
          st1(l)=0d0
        enddo
      endif

c   Le vecteur st1 donnait les coefficients du d�eloppement des variables du
c   mod�e en s�ie de puissances de x au dernier point et jusqu'�l'ordre 3. La
c   version actuelle n'utilise plus les coefficients des ordres 2 et 3.
      do 630 l=1,12
      st2(l)=star0(l,jmax)
  630 st2(l+12)=0.d0

      call asymp
      ir=0
      return

  810 if(lscreen)write(6,904)j
      ir=4
      return
  830 if(lscreen)write(6,906)
      ir=4
      return
      end
c
c     ASYMP
c     Calcule des grandeurs utiles dans l'�aluation asymptotique
c     des fr�uences.
c
      subroutine asymp
      implicit double precision(a-h,o-z)
      include 'osc.i'
      logical lp,lm,lp1,lm1
      asp=0.d0
      asg1=0.d0
      asg2=0.d0
      do 150 j=1,jmax
      x=star0(1,j)
      p=1.d0/dsqrt(star0(3,j)*star0(5,j))
      lp=.false.
      lm=.false.
      if(star0(6,j).lt.0.d0)lp=.true.
      if(star0(6,j).gt.0.d0)lm=.true.
      g=dsqrt(star0(2,j)*dabs(star0(6,j)))
      if(j.eq.1)goto 110
      if(x1.eq.x)goto 110
      h=x-x1
      h2=h/2.d0
      asp=asp+(p+p1)*h2
      if(lp.and.lp1)asg1=asg1+(g+g1)*h2
      if(lp.and..not.lp1)asg1=asg1+2.d0*g**3*h/3.d0/(g**2+g1**2)
      if(.not.lp.and.lp1)asg1=asg1+2.d0*g1**3*h/3.d0/(g**2+g1**2)
      if(lm.and.lm1)asg2=asg2+(g+g1)*h2
      if(lm.and..not.lm1)asg2=asg2+2.d0*g**3*h/3.d0/(g**2+g1**2)
      if(.not.lm.and.lm1)asg2=asg2+2.d0*g1**3*h/3.d0/(g**2+g1**2)
  110 star0(13,j)=asp
      star0(14,j)=asg1
      star0(15,j)=asg2
      x1=x
      p1=p
      g1=g
      lp1=lp
      lm1=lm
  150 continue
      return
      end

c===============================================================================

c   oscfile
c   OSC07
c   Opening and closing OSC files

c     subroutine osc07(oscname,ir)
      subroutine oscfile(oscname,ir)
      implicit none
      include 'osc.i'
      integer ir,ios
      character oscname*(*)
      if(lwriteosc)then
        close(21)
        lwriteosc=.false.
      endif
      if(oscname.eq.' ')return
      open(21,file=oscname,iostat=ios,form='unformatted')
      if(ios.eq.0)then
        ir=0
        lwriteosc=.true.
      else
        ir=1
        lwriteosc=.false.
      endif
      return
      end

c===============================================================================

c     OSC08
c     Ecriture d'un mode

      subroutine osc08
      implicit double precision(a-h,o-z)
      include 'osc.i'
      if(.not.losc)return
      write(21)descript(1)
      write(21)imaxmod,ray,am
      write(21)np,nres,ires,ll,nz,mode,om,beta,sigma,freq,periode,
     .ev,xm,del,modeLee,parity,boundary
      if(ll.le.0)then
        kmax=2
      else
        kmax=4
      endif
      do 110 j=1,np
  110 write(21)(star(k,j),k=1,5),(sol(k,j),k=1,kmax)
      return
      end


      subroutine osc08py(star_,sol_)
Cf2py intent(out) sol_
Cf2py intent(out) star_
      implicit double precision(a-h,o-z)
      real*8 sol_(4,12001)
      real*8 star_(6,12001)
      include 'osc.i'
      if(ll.le.0)then
        jmax=2
      else
        jmax=4
      endif
      do j=1,jmax
      do k=1,12001
      sol_(j,k) = sol(j,k)
      enddo
      enddo
      do j=1,6
      do k=1,12001
      star_(j,k) = star(j,k)
      enddo
      enddo
      return
      end
      

c===============================================================================
c
c     OSC11
c     Triangularisation d'une matrice bloc-diagonale
c
      subroutine osc11(n0,n,ip,a,io,d,id)
      implicit double precision(a-h,o-z)
      dimension a(1),io(1)
      nn=n+n
      nn1=nn+1
      ip1=ip+1
      jmax1=n
      jmax2=nn
      kmax=n+ip
      ii=n+1
      mm=1
      d=1.d0
      id=0
      do 113 i0=1,n0
      if(i0.ne.n0)goto 101
      kmax=n
      jmax2=n
  101 jmax=jmax1
      inc1=1
      inc2=nn
      jm=jmax
      km=kmax
      do 112 i=1,n
      if(i.ne.ip1)goto 102
      ii=ii-n
      jmax=jmax2
      inc1=nn
      inc2=1
      jm=kmax
      km=jmax
102   i1=i+1
      io1=0
      if(i1.gt.jm)goto 105
      r=dabs(a(ii))
      jj=ii+inc1
      do 104 j=i1,jm
      s=dabs(a(jj))
      if(r.ge.s)goto 104
      io1=j-i
      r=s
  104 jj=jj+inc1
  105 io(mm)=io1
      if(io1.eq.0)goto 108
      d=-d
      jj=ii+io1*inc1
      kk=ii
      ll=jj
      do 107 k=i,km
      if(k.ne.ip1.or.i.gt.ip)goto 106
      kk=kk-n
      ll=ll-n
  106 r=a(kk)
      a(kk)=a(ll)
      a(ll)=r
      kk=kk+inc2
  107 ll=ll+inc2
  108 r=a(ii)
      d=d*r
      if(d.eq.0.d0)goto 120
      if(dabs(d).ge.1.d1)then
  121 id=id+1
      d=d/1.d1
      if(dabs(d).ge.1.d1)goto 121
      goto 120
      endif
  122 if(dabs(d).lt.1.d0)then
      id=id-1
      d=d*1.d1
      if(dabs(d).lt.1.d0)goto 122
      endif
  120 if(r.eq.0.d0)goto 111
      if(i1.gt.kmax)goto 111
      jj=ii+nn
      do 110 k=i1,kmax
      if(k.eq.ip1)jj=jj-n
      s=a(jj)/r
      a(jj)=s
      if(i1.gt.jmax)goto 110
      kk=ii+1
      ll=jj+1
      do 109 j=i1,jmax
      a(ll)=a(ll)-a(kk)*s
      kk=kk+1
  109 ll=ll+1
  110 jj=jj+nn
  111 ii=ii+nn1
  112 mm=mm+1
      if(ip.eq.n)ii=ii-n
113   continue
      return
      end
c
c     OSC12
c     Solution de Ax = b (apr� osc11)
c
      subroutine osc12(n0,n,ip,a,io,b)
      implicit double precision(a-h,o-z)
      dimension a(1),io(1),b(1)
      n01=n0+1
      n1=n+1
      nn=n+n
      nn1=nn+1
      ip1=ip+1
      kmax=n+ip
      mm=1
      ii=n1
      do 205 i0=1,n0
      if(i0.eq.n0)kmax=n
      do 204 i=1,n
      if(i.eq.ip1)ii=ii-n
      i1=i+1
      if(i.le.ip)goto 201
      io1=io(mm)
      if(io1.eq.0)goto 201
      jj=mm+io1
      r=b(mm)
      b(mm)=b(jj)
      b(jj)=r
201   if(i1.gt.kmax)goto 203
      s=b(mm)
      kk=ii+nn
      ll=mm+1
      do 202 k=i1,kmax
      if(k.eq.ip1)kk=kk-n
      b(ll)=b(ll)-a(kk)*s
      kk=kk+nn
202   ll=ll+1
203   ii=ii+nn1
204   mm=mm+1
      if(ip.eq.n)ii=ii-n
205   continue
      jmax1=n
      jmax2=n
      mm=n0*n
      ii=mm*nn-n
      do 211 i0v=1,n0
      i0=n01-i0v
      if(i0v.eq.2)jmax2=nn
      jmax=jmax2
      do 209 iv=1,n
      i=n1-iv
      if(i.ne.ip)goto 210
      jmax=jmax1
      ii=ii+n
210   i1=i+1
      s=b(mm)
      if(i1.gt.jmax)goto 207
      jj=ii+1
      ll=mm+1
      do 206 j=i1,jmax
      s=s-a(jj)*b(ll)
      jj=jj+1
206   ll=ll+1
207   b(mm)=s/a(ii)
      if(i.gt.ip)goto 208
      io1=io(mm)
      if(io1.eq.0)goto 208
      jj=mm+io1
      r=b(mm)
      b(mm)=b(jj)
      b(jj)=r
208   ii=ii-nn1
209   mm=mm-1
      if(ip.eq.0)ii=ii+n
211   continue
      return
      end
c
c     OSC13
c     Calcul de Ax
c
      subroutine osc13(n0,n,ip,a,b)
      implicit double precision(a-h,o-z)
      dimension a(1),b(1),bb(8)
      nn=2*n
      iq=n-ip
      ii=1
      l=1
      mm=1
      do 210 m=1,nn
  210 bb(m)=0.d0
      do 240 i0=1,n0+1
      if(i0.le.n0)then
      do 220 m=1,n
      bb(m)=bb(m+n)
      bb(m+n)=b(mm)
  220 mm=mm+1
      else
      do 225 m=1,n
      bb(m)=bb(m+n)
  225 bb(m+n)=0.d0
      endif
      if(i0.eq.1)imax=ip
      if(i0.gt.1.and.i0.le.n0)imax=n
      if(i0.gt.n0)imax=iq
      do 240 i=1,imax
      s=0.d0
      do 230 k=1,nn
      s=s+a(ii)*bb(k)
  230 ii=ii+1
      b(l)=s
  240 l=l+1
      return
      end
c
c     OSC14
c     Triangularisation d'une matrice
c
      subroutine osc14(n,a,aa,ioc,d)
      implicit double precision(a-h,o-z)
      dimension a(1),ioc(1),aa(1)
      n2=n**2
      do 111 i=1,n2
  111 aa(i)=a(i)
      d=1.d0
      np1=n+1
      ii=1
      do 106 i=1,n
      l=i
      if(i.eq.n)goto 103
      ip1=i+1
      r=dabs(aa(ii))
      jj=ii+n
      do 101 j=ip1,n
      s=dabs(aa(jj))
      if(r.ge.s)goto 101
      r=s
      l=j
      ll=jj
  101 jj=jj+n
      if(i.eq.l)goto 103
      d=-d
      jj=ii
      do 102 j=i,n
      r=aa(jj)
      aa(jj)=aa(ll)
      aa(ll)=r
      jj=jj+1
  102 ll=ll+1
  103 r=aa(ii)
      ioc(i)=l
      d=d*r
      if(d.eq.0.d0)return
      if(i.eq.n)goto 106
      jj=ii+n
      do 105 j=ip1,n
      s=aa(jj)/r
      aa(jj)=s
      kk=ii+1
      ll=jj+1
      do 104 k=ip1,n
      aa(ll)=aa(ll)-s*aa(kk)
      kk=kk+1
  104 ll=ll+1
  105 jj=jj+n
  106 ii=ii+np1
      return
      end
c
c     OSC15
c     Solution de Ax = b (apr� osc14)
c
      subroutine osc15(n,a,ioc,b,x)
      implicit double precision(a-h,o-z)
      dimension a(1),ioc(1),b(1),x(1)
      do 101 i=1,n
  101 x(i)=b(i)
      np1=n+1
      ii=1
      do 108 i=1,n
      ip1=i+1
      r=x(i)/a(ii)
      x(i)=r
      if(i.eq.n)goto 108
      jj=ii+1
      do 107 j=ip1,n
      x(j)=x(j)-a(jj)*r
  107 jj=jj+1
      ii=ii+np1
  108 continue
      ii=ii-np1
      do 110 l=2,n
      i=np1-l
      ip1=i+1
      r=x(i)
      jj=ii+n
      do 109 j=ip1,n
      r=r-a(jj)*x(j)
  109 jj=jj+n
      x(i)=r
      k=ioc(i)
      if(k.eq.i)goto 110
      x(i)=x(k)
      x(k)=r
  110 ii=ii-np1
      return
      end
c
c     OSC16
c     Multiplication de deux matrices
c
      subroutine osc16(n,a,b,c)
      implicit double precision(a-h,o-z)
      dimension a(n,n),b(n,n),c(n,n)
      do 102 il=1,n
      do 102 ic=1,n
      s=0.d0
      do 101 ik=1,n
  101 s=s+a(il,ik)*b(ik,ic)
  102 c(il,ic)=s
      return
      end
c
c     OSC17
c     Calcul de A**(-1)*B (apr� osc14)
c
      subroutine osc17(n,a,ioc,b,c)
      implicit double precision(a-h,o-z)
      dimension a(n,n),ioc(n),b(n,n),c(n,n)
      do 110 i=1,n
  110 call osc15(n,a,ioc,b(1,i),c(1,i))
      return
      end
c
c     OSC19
c     Interpolation
c
      subroutine osc19(x1,y1,dy1,x2,y2,dy2,x,y,dy)
      implicit real*8(a-h,o-z)
      phi0(u)=(2.d0*u-3.d0)*u**2+1.d0
      phi1(u)=((u-2.d0)*u+1.d0)*u
      psi0(u)=phi0(1.d0-u)
      psi1(u)=-phi1(1.d0-u)
      dphi0(u)=-6.d0*u*(1.d0-u)
      dphi1(u)=(3.d0*u-4.d0)*u+1.d0
      dpsi0(u)=-dphi0(1.d0-u)
      dpsi1(u)=dphi1(1.d0-u)
      h=x2-x1
      xi=(x-x1)/h
      y=y1*phi0(xi)+y2*psi0(xi)+h*(dy1*phi1(xi)+dy2*psi1(xi))
      dy=(y1*dphi0(xi)+y2*dpsi0(xi))/h+dy1*dphi1(xi)+dy2*dpsi1(xi)
      return
      end

c     setdegree
c     OSC24

c     subroutine osc24(l)
      subroutine setdegree(l)
Cf2py intent(in) l
      implicit real*8(a-h,o-z)
      include 'osc.i'
      ll=l
      return
      end

c     computemode
c     OSC25
c     Calcul d'un mode dont on conna� omega approch�
c     subroutine osc25(im,om0,ir)
      subroutine computemode(im,om0,ir)
Cf2py intent(in) im
Cf2py intent(in) om0
Cf2py intent(in) ir
      implicit real*8(a-h,o-z)
      include 'osc.i'
  906 format('*** error: undefined model')
  910 format('*** error: undefined mesh')
      losc=.false.
      if(.not.lmod)then
        ir=3
        if(lscreen)write(6,906)
        return
      endif
      if(.not.lres)then
        ir=4
        if(lscreen)write(6,910)
        return
      endif
      if(ll.eq.0)then
        call osc35(im,om0,ir)
      else
        call osc45(im,om0,ir)
      endif
      if(ir.ne.0)return
      losc=.true.
      if(lwriteosc)call osc08
      return
      end

c     mesh
c     OSC26
c     Calcul d'un r�eau et du mod�e sur ce r�eau
c     ires=0 : m�e r�eau que le mod�e
c          1 : r�eau adapt�au calcul de modes p
c          2 : r�eau adapt�au calcul de modes g
c          3 : r�eau adapt�au calcul des deux types de modes

      subroutine mesh(ires1,nres1,ir)
Cf2py intent(in) ires1
Cf2py intent(in) nres1
Cf2py intent(in) ir
c     subroutine osc26(ires1,nres1,ir)
      implicit double precision(a-h,o-z)
      include 'osc.i'
  901 format('   *** error : too many points')
  902 format('   ',i5,' points')
      ires=ires1
      nres=nres1
      if(ires.eq.0.or.nres.le.0)then
      ires=0
      np=jmax
      do 210 j=1,jmax
      do 210 l=1,12
  210 star(l,j)=star0(l,j)
      if(lscreen)write(6,902)np
      lres=.true.
      ir=0
      return
      endif
      if(nres.gt.npmax)then
      if(lscreen)write(6,901)
      lres=.false.
      ir=1
      return
      endif
      p0=star0(13,jmax)/nres
      g0=(star0(14,jmax)+star0(15,jmax))/nres
      k=0
      do 150 j=1,jmax
      if(j.eq.1)goto 130
      j1=j-1
      n1=(star0(13,j)-star0(13,j1))/p0+1
      n2=(star0(14,j)+star0(15,j)-star0(14,j1)-star0(15,j1))/g0+1
      if(ires.eq.1)n=n1
      if(ires.eq.2)n=n2
      if(ires.eq.3)n=max0(n1,n2)
      if(n.eq.1)goto 130
      xx1=star0(1,j1)
      xx=star0(1,j)
      do 120 l=1,n-1
      if(k.ge.npmax)then
      if(lscreen)write(6,901)
      lres=.false.
      ir=1
      return
      endif
      k=k+1
      x=((n-l)*xx1+l*xx)/n
      star(1,k)=x
      star(7,k)=1.d0
      do 110 m=2,6
      mm=m+6
  110 call osc19(xx1,star0(m,j1),star0(mm,j1),xx,star0(m,j),star0(mm,j),
     *x,star(m,k),star(mm,k))
  120 continue
  130 if(k.ge.npmax)then
      if(lscreen)write(6,901)
      ir=1
      lres=.false.
      return
      endif
      k=k+1
      do 140 m=1,12
  140 star(m,k)=star0(m,j)
  150 continue
      np=k
      if(lscreen)write(6,902)np
      ir=0
      lres=.true.
      return
      end

c     scanomega
c     OSC27
      subroutine scanomega(om1,om2,n,ind,m,omarr,k,ir)
Cf2py intent(in) om1
Cf2py intent(in) om2
Cf2py intent(in) n
Cf2py intent(in) ind
Cf2py intent(in) m
Cf2py intent(in,out) omarr
Cf2py intent(in,out) k
Cf2py intent(in) ir
      implicit real*8(a-h,o-z)
      include 'osc.i'
      real*8 omarr(*)
  906 format('*** error: undefined model')
  910 format('*** error: undefined mesh')
      if(.not.lmod)then
        ir=3
        if(lscreen)write(6,906)
        k=0
        return
      endif
      if(.not.lres)then
        ir=4
        if(lscreen)write(6,910)
        k=0
        return
      endif
      if(ll.eq.0)then
        call osc37(om1,om2,n,ind,ir,m,omarr,k)
      else
        call osc47(om1,om2,n,ind,ir,m,omarr,k)
      endif
      return
      end
c
c     OSC31
c     El�ents de matrice en un point ordinaire, l=0
c
      subroutine osc31(j,star,parm,aa,bb,daa,dbb,ind)
      implicit double precision(a-h,o-z)
      dimension star(12,1),a0(2,2),a1(2,2),da0(2,2),
     *da1(2,2),s(12),aa(2,2),bb(2,2),daa(2,2),dbb(2,2),t(2,2),
     *u(2,2)
      do 110 k=1,12
  110 s(k)=star(k,j)
      a0(1,1)=-3.d0/s(1)
      a0(1,2)=-1.d0/s(1)/s(5)
      a0(2,1)=(4.d0*s(2)+parm)*s(1)/s(3)
      a0(2,2)=s(1)*s(2)/s(3)
      do 112 k=1,2
      do 112 l=1,2
  112 aa(k,l)=a0(k,l)
      if(ind.eq.1)return
      a1(1,1)=3.d0/s(1)**2
      a1(1,2)=(1.d0/s(1)+s(11)/s(5))/s(1)/s(5)
      a1(2,1)=s(1)/s(3)*((4.d0*s(2)+parm)*(1.d0/s(1)
     *-s(9)/s(3))+4.d0*s(8))
      a1(2,2)=s(1)*s(2)/s(3)*(1.d0/s(1)+s(8)/s(2)-s(9)/s(3))
      call osc16(2,a0,a0,t)
      do 115 k=1,2
      do 115 l=1,2
  115 bb(k,l)=(t(k,l)+a1(k,l))/2.d0
      if(ind.eq.2)return
      do 120 k=1,2
      do 120 l=1,2
      da0(k,l)=0.d0
  120 da1(k,l)=0.d0
      da0(2,1)=s(1)/s(3)
      da1(2,1)=s(1)/s(3)*(1.d0/s(1)-s(9)/s(3))
      do 130 k=1,2
      do 130 l=1,2
  130 daa(k,l)=da0(k,l)
      call osc16(2,a0,da0,t)
      call osc16(2,da0,a0,u)
      do 150 k=1,2
      do 150 l=1,2
  150 dbb(k,l)=(t(k,l)+u(k,l)+da1(k,l))/2.d0
      return
      end
c
c     OSC32
c     El�ents de matrice au centre, l=0
c
      subroutine osc32(s,parm,cl,bb,dcl,dbb,ind)
      implicit double precision(a-h,o-z)
      dimension s(18),a0(2,2),a2(2,2),da0(2,2),da2(2,2),cl(2),
     *bb(2,2),dbb(2,2),t(2,2),ioc(2),u(2,2),dcl(2)
      a0(1,1)=-3.d0
      a0(1,2)=-1.d0/s(5)
      a0(2,1)=0.d0
      a0(2,2)=0.d0
      a2(1,1)=0.d0
      a2(1,2)=s(17)/s(5)**2
      a2(2,1)=(4.d0*s(2)+parm)/s(3)
      a2(2,2)=s(2)/s(3)
      do 105 l=1,2
      cl(l)=a0(1,l)
      do 105 k=1,2
      u(k,l)=-a0(k,l)
  105 if(k.eq.l)u(k,l)=u(k,l)+2.d0
      call osc14(2,u,t,ioc,det)
      call osc17(2,t,ioc,a2,bb)
      if(ind.eq.2)return
      do 110 k=1,2
      do 110 l=1,2
      da0(k,l)=0.d0
  110 da2(k,l)=0.d0
      da2(2,1)=1.d0/s(3)
      do 120 k=1,2
  120 dcl(k)=da0(1,k)
      call osc16(2,da0,bb,u)
      do 170 k=1,2
      do 170 l=1,2
  170 u(k,l)=u(k,l)+da2(k,l)
      call osc17(2,t,ioc,u,dbb)
      return
      end
c
c     OSC33
c     El�ents de matrice en surface, l=0
c
      subroutine osc33(s,parm,cl,aa,bb,dcl,daa,dbb,ind)
      implicit double precision(a-h,o-z)
      include 'osc.i'
      dimension s(24),cl(2),aa(2,2),bb(2,2),daa(2,2),dbb(2,2),dcl(2)
  901 format(' error: the boundary condition parameter has an',
     .' illegal value')
      if(boundary.eq.1)then
        call osc331(s,parm,cl,aa,bb,dcl,daa,dbb,ind)
      elseif(boundary.eq.2)then
        call osc332(s,parm,cl,aa,bb,dcl,daa,dbb,ind)
      else
        write(6,901)
        call exit(1)
      endif
      return
      end
c
c     OSC331
c     El�ents de matrice en surface, l=0
c     La condition limite sur delta P est celle qu'il convient d'utiliser pour
c     les mod�es de temp�ature nulle
c
      subroutine osc331(s,parm,cl,aa,bb,dcl,daa,dbb,ind)
      implicit double precision(a-h,o-z)
      dimension s(24),a0(2,2),a1(2,2),a2(2,2),da0(2,2),da1(2,2),
     *da2(2,2),cl(2),aa(2,2),bb(2,2),daa(2,2),dbb(2,2),t(2,2),
     *u(2,2),ioc(2),v(2,2),w(2,2),dcl(2)
      if(s(3).eq.0.d0)then
      a0(1,1)=0.d0
      a0(1,2)=0.d0
      a0(2,1)=(4.d0*s(2)+parm)*s(1)/s(9)
      a0(2,2)=s(1)*s(2)/s(9)
      a1(1,1)=-3.d0/s(1)
      a1(1,2)=-1.d0/s(1)/s(5)
      a1(2,1)=s(1)/s(9)*((4.d0*s(2)+parm)*(1.d0/s(1)-s(15)/s(9))
     *+4.d0*s(8))
      a1(2,2)=s(1)*s(2)/s(9)*(1.d0/s(1)+s(8)/s(2)-s(15)/s(9))
      do 103 l=1,2
      cl(l)=a0(2,l)
      do 103 k=1,2
      u(k,l)=-a0(k,l)
  103 if(k.eq.l)u(k,l)=u(k,l)+1.d0
      call osc14(2,u,t,ioc,det)
      call osc17(2,t,ioc,a1,aa)
      if(ind.eq.1)return
      a2(1,1)=3.d0/s(1)**2
      a2(1,2)=(1.d0/s(1)+s(11)/s(5))/s(1)/s(5)
      a2(2,1)=s(1)/s(9)*((4.d0*s(2)+parm)*((s(15)/s(9))**2
     *-s(21)/s(9)-s(15)/s(1)/s(9))+4.d0*s(8)*(1.d0/s(1)-s(15)/s(9))
     *+4.d0*s(14))
      a2(2,2)=s(1)*s(2)/s(9)*(s(14)/s(2)+(s(15)/s(9))**2-s(21)/s(9)
     *+s(8)/s(1)/s(2)-s(15)/s(1)/s(9)-s(8)*s(15)/s(2)/s(9))
      do 106 l=1,2
      do 106 k=1,2
      u(k,l)=-a0(k,l)
  106 if(k.eq.l)u(k,l)=u(k,l)+2.d0
      call osc14(2,u,t,ioc,det)
      call osc16(2,a1,aa,u)
      do 180 k=1,2
      do 180 l=1,2
  180 u(k,l)=u(k,l)+a2(k,l)
      call osc17(2,t,ioc,u,bb)
      if(ind.eq.2)return
      do 110 k=1,2
      do 110 l=1,2
      da0(k,l)=0.d0
      da1(k,l)=0.d0
  110 da2(k,l)=0.d0
      da0(2,1)=s(1)/s(9)
      da1(2,1)=s(1)/s(9)*(1.d0/s(1)-s(15)/s(9))
      da2(2,1)=s(1)/s(9)*((s(15)/s(9))**2-s(21)/s(9)-s(15)/s(1)/s(9))
      do 120 k=1,2
  120 dcl(k)=da0(2,k)
      call osc16(2,da0,aa,u)
      do 150 k=1,2
      do 150 l=1,2
  150 u(k,l)=u(k,l)+da1(k,l)
      call osc17(2,t,ioc,u,daa)
      call osc16(2,da1,aa,u)
      call osc16(2,a1,daa,v)
      call osc16(2,da0,bb,w)
      do 190 k=1,2
      do 190 l=1,2
  190 u(k,l)=da2(k,l)+u(k,l)+v(k,l)+w(k,l)
      call osc17(2,t,ioc,u,dbb)
      return
      else
      call osc31(1,s,parm,aa,bb,daa,dbb,ind)
      cl(1)=4.d0+parm/s(2)
      cl(2)=1.d0
      if(ind.eq.2)return
      dcl(1)=1.d0/s(2)
      dcl(2)=0.d0
      return
      endif
      end
c
c     OSC332
c     El�ents de matrice en surface, l=0
c     Condition limite delta P = 0
c
      subroutine osc332(s,parm,cl,aa,bb,dcl,daa,dbb,ind)
      implicit double precision(a-h,o-z)
      dimension s(24),a0(2,2),a1(2,2),a2(2,2),da0(2,2),da1(2,2),
     *da2(2,2),cl(2),aa(2,2),bb(2,2),daa(2,2),dbb(2,2),t(2,2),
     *u(2,2),ioc(2),v(2,2),w(2,2),dcl(2)
      if(s(3).eq.0.d0)then
      a0(1,1)=0.d0
      a0(1,2)=0.d0
      a0(2,1)=(4.d0*s(2)+parm)*s(1)/s(9)
      a0(2,2)=s(1)*s(2)/s(9)
      a1(1,1)=-3.d0/s(1)
      a1(1,2)=-1.d0/s(1)/s(5)
      a1(2,1)=s(1)/s(9)*((4.d0*s(2)+parm)*(1.d0/s(1)-s(15)/s(9))
     *+4.d0*s(8))
      a1(2,2)=s(1)*s(2)/s(9)*(1.d0/s(1)+s(8)/s(2)-s(15)/s(9))
      do 103 l=1,2
      cl(l)=a0(2,l)
      do 103 k=1,2
      u(k,l)=-a0(k,l)
  103 if(k.eq.l)u(k,l)=u(k,l)+1.d0
      call osc14(2,u,t,ioc,det)
      call osc17(2,t,ioc,a1,aa)
      if(ind.eq.1)return
      a2(1,1)=3.d0/s(1)**2
      a2(1,2)=(1.d0/s(1)+s(11)/s(5))/s(1)/s(5)
      a2(2,1)=s(1)/s(9)*((4.d0*s(2)+parm)*((s(15)/s(9))**2
     *-s(21)/s(9)-s(15)/s(1)/s(9))+4.d0*s(8)*(1.d0/s(1)-s(15)/s(9))
     *+4.d0*s(14))
      a2(2,2)=s(1)*s(2)/s(9)*(s(14)/s(2)+(s(15)/s(9))**2-s(21)/s(9)
     *+s(8)/s(1)/s(2)-s(15)/s(1)/s(9)-s(8)*s(15)/s(2)/s(9))
      do 106 l=1,2
      do 106 k=1,2
      u(k,l)=-a0(k,l)
  106 if(k.eq.l)u(k,l)=u(k,l)+2.d0
      call osc14(2,u,t,ioc,det)
      call osc16(2,a1,aa,u)
      do 180 k=1,2
      do 180 l=1,2
  180 u(k,l)=u(k,l)+a2(k,l)
      call osc17(2,t,ioc,u,bb)
      if(ind.eq.2)return
      do 110 k=1,2
      do 110 l=1,2
      da0(k,l)=0.d0
      da1(k,l)=0.d0
  110 da2(k,l)=0.d0
      da0(2,1)=s(1)/s(9)
      da1(2,1)=s(1)/s(9)*(1.d0/s(1)-s(15)/s(9))
      da2(2,1)=s(1)/s(9)*((s(15)/s(9))**2-s(21)/s(9)-s(15)/s(1)/s(9))
      do 120 k=1,2
  120 dcl(k)=da0(2,k)
      call osc16(2,da0,aa,u)
      do 150 k=1,2
      do 150 l=1,2
  150 u(k,l)=u(k,l)+da1(k,l)
      call osc17(2,t,ioc,u,daa)
      call osc16(2,da1,aa,u)
      call osc16(2,a1,daa,v)
      call osc16(2,da0,bb,w)
      do 190 k=1,2
      do 190 l=1,2
  190 u(k,l)=da2(k,l)+u(k,l)+v(k,l)+w(k,l)
      call osc17(2,t,ioc,u,dbb)
      return
      else
      call osc31(1,s,parm,aa,bb,daa,dbb,ind)
      cl(1)=0d0
      cl(2)=1.d0
      if(ind.eq.2)return
      dcl(1)=0d0
      dcl(2)=0.d0
      return
      endif
      end
c
c     OSC34
c     Matrice du probl�e radial
c
      subroutine osc34(np,star,st1,st2,parm,aco,daco,ind)
      implicit double precision(a-h,o-z)
      dimension star(12,np),st1(18),st2(24),aco(1),daco(1),
     *cl(2),aa(2,2),bb(2,2),dcl(2),daa(2,2),dbb(2,2),c(2,2),
     *d(2,2),dc(2,2),dd(2,2)
      if(ind.eq.2)goto 100
      kk=0
      do 260 j=1,np
      if(j.eq.1)then
      if(star(1,j).eq.0d0)then
      call osc32(st1,parm,cl,bb,dcl,dbb,2)
      do 210 k=1,2
      do 210 l=1,2
  210 aa(k,l)=0.d0
      else
      call osc31(j,star,parm,aa,bb,daa,dbb,2)
      cl(1)=1d0
      cl(2)=0d0
      endif
      do 220 k=1,2
      kk=kk+1
  220 aco(kk)=0.d0
      do 230 k=1,2
      kk=kk+1
  230 aco(kk)=cl(k)
      goto 250
      endif
      if(j.eq.np)then
      call osc33(st2,parm,cl,aa,bb,dcl,daa,dbb,2)
      goto 240
      endif
      call osc31(j,star,parm,aa,bb,daa,dbb,2)
  240 h=star(1,j)-star(1,j-1)
      h2=h/2.d0
      h6=h**2/6.d0
      do 270 k=1,2
      do 270 l=1,2
      d(k,l)=-h2*aa(k,l)+h6*bb(k,l)
  270 if(k.eq.l)d(k,l)=d(k,l)+1.d0
      do 280 k=1,2
      do 275 l=1,2
      kk=kk+1
  275 aco(kk)=-c(k,l)
      do 280 l=1,2
      kk=kk+1
  280 aco(kk)=d(k,l)
      if(j.eq.np)then
      do 290 k=1,2
      kk=kk+1
  290 aco(kk)=cl(k)
      do 295 k=1,2
      kk=kk+1
  295 aco(kk)=0.d0
      goto 260
      endif
  250 h=star(1,j+1)-star(1,j)
      h2=h/2.d0
      h6=h**2/6.d0
      do 255 k=1,2
      do 255 l=1,2
      c(k,l)=h2*aa(k,l)+h6*bb(k,l)
  255 if(k.eq.l)c(k,l)=c(k,l)+1.d0
  260 continue
      return
  100 kk=0
      do 160 j=1,np
      if(j.eq.1)then
      if(star(1,j).eq.0d0)then
      call osc32(st1,parm,cl,bb,dcl,dbb,3)
      do 103 k=1,2
      do 103 l=1,2
      aa(k,l)=0.d0
  103 daa(k,l)=0.d0
      else
      call osc31(j,star,parm,aa,bb,daa,dbb,3)
      cl(1)=1d0
      cl(2)=0d0
      dcl(1)=0d0
      dcl(2)=0d0
      endif
      do 101 k=1,2
      kk=kk+1
      aco(kk)=0.d0
  101 daco(kk)=0.d0
      do 102 k=1,2
      kk=kk+1
      aco(kk)=cl(k)
  102 daco(kk)=dcl(k)
      goto 150
      endif
      if(j.eq.np)then
      call osc33(st2,parm,cl,aa,bb,dcl,daa,dbb,3)
      goto 140
      endif
      call osc31(j,star,parm,aa,bb,daa,dbb,3)
  140 h=star(1,j)-star(1,j-1)
      h2=h/2.d0
      h6=h**2/6.d0
      do 142 k=1,2
      do 141 l=1,2
      d(k,l)=-h2*aa(k,l)+h6*bb(k,l)
  141 dd(k,l)=-h2*daa(k,l)+h6*dbb(k,l)
  142 d(k,k)=d(k,k)+1.d0
      do 145 k=1,2
      do 143 l=1,2
      kk=kk+1
      aco(kk)=-c(k,l)
  143 daco(kk)=-dc(k,l)
      do 144 l=1,2
      kk=kk+1
      aco(kk)=d(k,l)
  144 daco(kk)=dd(k,l)
  145 continue
      if(j.eq.np)then
      do 146 k=1,2
      kk=kk+1
      aco(kk)=cl(k)
  146 daco(kk)=dcl(k)
      do 147 k=1,2
      kk=kk+1
      aco(kk)=0.d0
  147 daco(kk)=0.d0
      goto 160
      endif
  150 h=star(1,j+1)-star(1,j)
      h2=h/2.d0
      h6=h**2/6.d0
      do 152 k=1,2
      do 151 l=1,2
      c(k,l)=h2*aa(k,l)+h6*bb(k,l)
  151 dc(k,l)=h2*daa(k,l)+h6*dbb(k,l)
  152 c(k,k)=c(k,k)+1.d0
  160 continue
      return
      end
c
c     OSC35
c     Calcul d'un mode radial
c
      subroutine osc35(im,om0,ir)
      implicit double precision(a-h,o-z)
      include 'osc.i'
      common io(2*npmax),x(2*npmax),x1(2*npmax),aco(8*npmax),
     $daco(8*npmax)
      parameter(deuxpi=6.28318530717959d0)
      character c
      data e1,e2,e3/1.d-8,1.d-2,1.d-4/
  901 format(3x,'it= 0',3x,'omega=',1pd12.5)
  905 format(3x,a1,1x,'it=',i2,3x,'omega=',1pd12.5,3x,'corr=',d12.5)
  907 format('   l=',i3,'   omega=',1pd15.8,'   sigma=',d15.8/
     *'   nz=',i4,'   frequency=',d15.8,'   period=',
     *d15.8/'   rotational splitting beta=',d15.8,'   parity=',i1/
     *'   ev=',d15.8,'   xm=',d15.8,'   delta=',d15.8/
     *'   mode=',i5,',',i5)
  911 format('   *** error : too many iterations')
      if(lscreen)write(6,901)om0
      call flush(6)
      om=om0
      om2=om**2
      if(om.lt.0.d0)om2=-om2
      iss=0
      np2=2*np
      js=np2-1
      do 1305 j=1,np2
      x(j)=1d0
 1305 x1(j)=x(j)
 1310 omega2=om2
      call osc34(np,star,st1,st2,omega2,aco,daco,2)
      call osc11(np,2,1,aco,io,det,idet)
      c='*'
      is=0
 1320 is=is+1
      iss=iss+1
      if(iss.gt.20)goto 1370
      call osc13(np,2,1,daco,x)
      call osc12(np,2,1,aco,io,x)
      s1=0.d0
      s2=0.d0
      xx=x(np)
      do 1322 j=1,np2
      x(j)=x(j)/xx
      s1=s1+x1(j)*x(j)
      s2=s2+x(j)**2
 1322 x1(j)=x(j)
      dom2=-s1/s2/xx
      om2a=om2
      om2=omega2+dom2
      oma=om
      om=dsqrt(dabs(om2))
      if(om2.lt.0.d0)om=-om
      corr=om-oma
      if(lscreen)write(6,905)c,iss,om,corr
      call flush(6)
      c=' '
c     if(iss.ge.2.and.dabs((om2-om2a)/omega2).le.e1.and.
c    *dabs(dom2/om2).le.e3)goto 1350
      if(iss.ge.2.and.dabs((om2-om2a)/omega2).le.e1)goto 1350
      if(im.eq.1.and.is.ge.2.and.dabs((om2-om2a)/dom2).le.e2)
     *goto 1310
      goto 1320
 1350 continue
      omeg2=om**2
      sigma=om/tdyn
      freq=sigma/deuxpi
      periode=1d0/freq
      call osc36(np,star,st2,om,x,nz,mode,xm,del,ev,beta,sol)
      modeLee=mode
      call osc38
      if(lscreen)write(6,907)0,om,sigma,nz,freq,periode,beta,parity,
     .ev,xm,del,mode,modeLee
      call flush(6)
      ir=0
      return
 1370 if(lscreen)write(6,911)
      ir=2
      return
      end
c
c     OSC36
c     Analyse d'un mode radial
c
      subroutine osc36(np,star,st2,omega,x,nz,mode,xm,del,ev,beta,sol)
      implicit double precision(a-h,o-z)
      dimension star(12,1),sol(9,1),st2(24),cl(2),aa(2,2),
     *bb(2,2),dcl(2),daa(2,2),dbb(2,2),x(2,1)
      ev=1.d0
      eh=0.d0
      beta=0d0
      nz=0
      a2=0.d0
      sx=0.d0
      sx2=0.d0
      omega2=omega**2
      if(omega.lt.0.d0)omega2=-omega2
      do 120 j=1,np
      sol(1,j)=x(1,j)
      sol(2,j)=x(2,j)
      sol(3,j)=0.d0
      sol(4,j)=0.d0
      sol(7,j)=0.d0
      sol(8,j)=0.d0
      sol(9,j)=0.d0
c     if(j.eq.1)then
      if(star(1,j).eq.0d0)then
      sol(5,j)=0.d0
      sol(6,j)=0.d0
      else
      if(j.eq.np)then
      call osc33(st2,omega2,cl,aa,bb,dcl,daa,dbb,1)
      else
      call osc31(j,star,omega2,aa,bb,daa,dbb,1)
      endif
      do 105 k=1,2
      s=0.d0
      do 103 m=1,2
  103 s=s+aa(k,m)*sol(m,j)
  105 sol(k+4,j)=s
      endif
      xx=star(1,j)
      xx0=xx**5
      xx1=xx0*xx
      xx2=xx1*xx
      y=sol(1,j)
      f=star(4,j)*y**2
      if(j.eq.1)goto 110
      if(y*y1.lt.0.d0.or.y.eq.0.d0)then
      nz=nz+1
      endif
      ff=f+f1
      a2=a2+ff*(xx0-xx01)
      sx=sx+ff*(xx1-xx11)
      sx2=sx2+ff*(xx2-xx21)
  110 xx01=xx0
      xx11=xx1
      xx21=xx2
      y1=y
  120 f1=f
      mode=nz+1
      a2=a2/1.d1
      sx=sx/a2/12.d0
      sx2=sx2/a2/14.d0
      xm=sx
      del=2.d0*dsqrt(sx2-sx**2)
      a=dsqrt(a2)
      if(sol(1,1).lt.0.d0)a=-a
      do 130 j=1,np
      do 125 l=1,2
  125 sol(l,j)=sol(l,j)/a
      do 130 l=5,6
  130 sol(l,j)=sol(l,j)/a
      return
      end
c
c     OSC37
c     Balayage en fr�uence, l=0
c     On stocke dans le vecteur omarr les m premi�es valeurs propres approch�s
c     obtenues au cours du balayage. Au retour m = le nombre de valeurs stock�s
c     dans le vecteur omarr.
c
      subroutine osc37(om1,om2,n,ind,ir,m,omarr,k)
      implicit double precision(a-h,o-z)
      include 'osc.i'
      real*8 omarr(*)
      common aco(8*npmax),daco(8*npmax),io(2*npmax)
  901 format(' error : om1 and om2 have opposite signs')
  902 format(' l=  0   omega=',1pd12.5,3x,'det=',0pf8.5,' (',i3,')')
  903 format(26x,'interpolated value : omega=',1pd12.5)
      if(om1*om2.le.0.d0)then
      ir=2
      if(lscreen)write(6,901)
      return
      endif
      k=0
      do 200 i=1,n+1
      goto(110,120),ind
  110 om=((n+1-i)*om1+(i-1)*om2)/n
      goto 130
  120 om=n/((n+1-i)/om1+(i-1)/om2)
  130 omega2=om**2
      if(om.lt.0.d0)omega2=-omega2
      call osc34(np,star,st1,st2,omega2,aco,daco,1)
      call osc11(np,2,1,aco,io,det,idet)
      if(i.gt.1)then
      sdet=dsign(1.d0,det)
      sdeto=dsign(1.d0,deto)
      if(sdet*sdeto.le.0d0)then
      det1=det
      deto1=deto
      if(idet.gt.ideto)deto1=deto*1.d1**(ideto-idet)
      if(idet.lt.ideto)det1=det*1.d1**(idet-ideto)
      omi=(det1*omo-deto1*om)/(det1-deto1)
      if(lscreen)write(6,903)omi
      call flush(6)
      if(k.lt.m)then
        k=k+1
        omarr(k)=omi
      endif
      endif
      endif
      if(lscreen)write(6,902)om,det,idet
      call flush(6)
      omo=om
      deto=det
      ideto=idet
  200 continue
      ir=0
      return
      end

c===============================================================================

c   OSC38
c   Sign of the derivative of the determinant with respect to omega2

      subroutine osc38
      implicit none
      include 'osc.i'
      integer io,idet1,idet2
      real*8 omega21,omega22,aco,daco,det1,det2
      common io(2*npmax),aco(8*npmax),daco(8*npmax)
      omega21=omeg2*(1d0-1d-6)
      call osc34(np,star,st1,st2,omega21,aco,daco,1)
      call osc11(np,2,1,aco,io,det1,idet1)
      omega22=omeg2*(1d0+1d-6)
      call osc34(np,star,st1,st2,omega22,aco,daco,1)
      call osc11(np,2,1,aco,io,det2,idet2)
      if(det1-det2.lt.0d0)then
        parity=1
      else
        parity=0
      endif
      return
      end

c===============================================================================
c
c     OSC41
c     El�ents de matrice en un point ordinaire, l>0
c
      subroutine osc41(j,star,ll,om2,aa,bb,daa,dbb,ind)
      implicit double precision(a-h,o-z)
      dimension star(12,1),a0(4,4),a1(4,4),da0(4,4),da1(4,4),
     *s(12),aa(4,4),bb(4,4),daa(4,4),dbb(4,4),t(4,4),u(4,4),v(4,4)
      do 110 k=1,12
  110 s(k)=star(k,j)
      al=ll
      al1=al+1.d0
      alo=al*al1/om2
      dalo=-alo/om2
      a0(1,1)=(alo*s(2)-al1)/s(1)
      a0(1,2)=alo*s(3)/s(1)-s(1)/s(5)
      a0(1,3)=alo/s(1)
      a0(1,4)=0.d0
      a0(2,1)=(om2+(4.d0-alo*s(2))*s(2))/s(1)/s(3)
      a0(2,2)=s(1)*s(2)/s(3)-(al+alo*s(2))/s(1)
      a0(2,3)=-alo*s(2)/s(1)/s(3)
      a0(2,4)=-1.d0/s(1)/s(3)
      a0(3,1)=-s(4)/s(1)
      a0(3,2)=0.d0
      a0(3,3)=-al/s(1)
      a0(3,4)=1.d0/s(1)
      a0(4,1)=alo*s(2)*s(4)/s(1)
      a0(4,2)=alo*s(3)*s(4)/s(1)
      a0(4,3)=(al*al1+alo*s(4))/s(1)
      a0(4,4)=-al1/s(1)
      do 115 k=1,4
      do 115 l=1,4
  115 aa(k,l)=a0(k,l)
      if(ind.eq.1)return
      a1(1,1)=(alo*s(8)-(alo*s(2)-al1)/s(1))/s(1)
      a1(1,2)=alo*s(3)/s(1)*(s(9)/s(3)-1.d0/s(1))-s(1)/s(5)
     **(1.d0/s(1)-s(11)/s(5))
      a1(1,3)=-alo/s(1)**2
      a1(1,4)=0.d0
      a1(2,1)=((4.d0-2.d0*alo*s(2))*s(8)-(om2+4.d0*s(2)
     *-alo*s(2)**2)*(1.d0/s(1)+s(9)/s(3)))/s(1)/s(3)
      a1(2,2)=s(1)*s(2)/s(3)*(1.d0/s(1)+s(8)/s(2)-s(9)/s(3))
     *+((al+alo*s(2))/s(1)-alo*s(8))/s(1)
      a1(2,3)=-alo*s(2)/s(1)/s(3)*(s(8)/s(2)-1.d0/s(1)-s(9)/s(3))
      a1(2,4)=1.d0/s(1)/s(3)*(1.d0/s(1)+s(9)/s(3))
      a1(3,1)=-s(4)/s(1)*(s(10)/s(4)-1.d0/s(1))
      a1(3,2)=0.d0
      a1(3,3)=al/s(1)**2
      a1(3,4)=-1.d0/s(1)**2
      a1(4,1)=alo*s(2)*s(4)/s(1)*(s(8)/s(2)+s(10)/s(4)-1.d0/s(1))
      a1(4,2)=alo*s(3)*s(4)/s(1)*(s(9)/s(3)+s(10)/s(4)-1.d0/s(1))
      a1(4,3)=(alo*s(10)-(al*al1+alo*s(4))/s(1))/s(1)
      a1(4,4)=al1/s(1)**2
      call osc16(4,a0,a0,t)
      do 120 k=1,4
      do 120 l=1,4
  120 bb(k,l)=(a1(k,l)+t(k,l))/2.d0
      if(ind.eq.2)return
      da0(1,1)=s(2)*dalo/s(1)
      da0(1,2)=dalo*s(3)/s(1)
      da0(1,3)=dalo/s(1)
      da0(1,4)=0.d0
      da0(2,1)=(1.d0-dalo*s(2)**2)/s(1)/s(3)
      da0(2,2)=-dalo*s(2)/s(1)
      da0(2,3)=-dalo*s(2)/s(1)/s(3)
      da0(2,4)=0.d0
      da0(3,1)=0.d0
      da0(3,2)=0.d0
      da0(3,3)=0.d0
      da0(3,4)=0.d0
      da0(4,1)=dalo*s(2)*s(4)/s(1)
      da0(4,2)=dalo*s(3)*s(4)/s(1)
      da0(4,3)=dalo*s(4)/s(1)
      da0(4,4)=0.d0
      da1(1,1)=dalo*(s(8)-s(2)/s(1))/s(1)
      da1(1,2)=dalo*s(3)/s(1)*(s(9)/s(3)-1.d0/s(1))
      da1(1,3)=-dalo/s(1)**2
      da1(1,4)=0.d0
      da1(2,1)=((dalo*s(2)**2-1.d0)*(1.d0/s(1)+s(9)/s(3))-2.d0
     **dalo*s(2)*s(8))/s(1)/s(3)
      da1(2,2)=dalo*(s(2)/s(1)-s(8))/s(1)
      da1(2,3)=-dalo*s(2)/s(1)/s(3)*(s(8)/s(2)-1.d0/s(1)-s(9)/s(3))
      da1(2,4)=0.d0
      da1(3,1)=0.d0
      da1(3,2)=0.d0
      da1(3,3)=0.d0
      da1(3,4)=0.d0
      da1(4,1)=dalo*s(2)*s(4)/s(1)*(s(8)/s(2)+s(10)/s(4)-1.d0/s(1))
      da1(4,2)=dalo*s(3)*s(4)/s(1)*(s(9)/s(3)+s(10)/s(4)-1.d0/s(1))
      da1(4,3)=dalo*(s(10)-s(4)/s(1))/s(1)
      da1(4,4)=0.d0
      call osc16(4,a0,da0,u)
      call osc16(4,da0,a0,v)
      do 130 k=1,4
      do 130 l=1,4
      daa(k,l)=da0(k,l)
  130 dbb(k,l)=(da1(k,l)+u(k,l)+v(k,l))/2.d0
      return
      end
c
c     OSC42
c     El�ents de matrice au centre, l>0
c
      subroutine osc42(s,ll,om2,cl,bb,dcl,dbb,ind)
      implicit double precision(a-h,o-z)
      dimension s(18),a0(4,4),a2(4,4),da0(4,4),da2(4,4),ioc(4),
     *bb(4,4),dbb(4,4),t(4,4),u(4,4),cl(2,4),dcl(2,4)
      al=ll
      al1=al+1.d0
      alo=al*al1/om2
      dalo=-alo/om2
      a0(1,1)=-al1+alo*s(2)
      a0(1,2)=alo*s(3)
      a0(1,3)=alo
      a0(1,4)=0.d0
      a0(2,1)=(om2+(4.d0-alo*s(2))*s(2))/s(3)
      a0(2,2)=-(al+alo*s(2))
      a0(2,3)=-alo*s(2)/s(3)
      a0(2,4)=-1.d0/s(3)
      a0(3,1)=-s(4)
      a0(3,2)=0.d0
      a0(3,3)=-al
      a0(3,4)=1.d0
      a0(4,1)=alo*s(2)*s(4)
      a0(4,2)=alo*s(3)*s(4)
      a0(4,3)=alo*s(4)+al*al1
      a0(4,4)=-al1
      a2(1,1)=alo*s(14)
      a2(1,2)=alo*s(15)-1.d0/s(5)
      a2(1,3)=0.d0
      a2(1,4)=0.d0
      a2(2,1)=((4.d0-2.d0*alo*s(2))*s(14)-(om2+
     *(4.d0-alo*s(2))*s(2))*s(15)/s(3))/s(3)
      a2(2,2)=-alo*s(14)+s(2)/s(3)
      a2(2,3)=-alo*s(2)/s(3)*(s(14)/s(2)-s(15)/s(3))
      a2(2,4)=-s(15)/s(3)**2
      a2(3,1)=-s(16)
      a2(3,2)=0.d0
      a2(3,3)=0.d0
      a2(3,4)=0.d0
      a2(4,1)=alo*(s(14)*s(4)+s(2)*s(16))
      a2(4,2)=alo*(s(15)*s(4)+s(3)*s(16))
      a2(4,3)=alo*s(16)
      a2(4,4)=0.d0
      do 105 l=1,4
      cl(1,l)=a0(1,l)
      cl(2,l)=a0(3,l)
      do 105 k=1,4
      u(k,l)=-a0(k,l)
  105 if(k.eq.l)u(k,l)=u(k,l)+2.d0
      call osc14(4,u,t,ioc,det)
      call osc17(4,t,ioc,a2,bb)
      if(ind.eq.2)return
      da0(1,1)=dalo*s(2)
      da0(1,2)=dalo*s(3)
      da0(1,3)=dalo
      da0(1,4)=0.d0
      da0(2,1)=(1.d0-dalo*s(2)**2)/s(3)
      da0(2,2)=-dalo*s(2)
      da0(2,3)=-dalo*s(2)/s(3)
      da0(2,4)=0.d0
      da0(3,1)=0.d0
      da0(3,2)=0.d0
      da0(3,3)=0.d0
      da0(3,4)=0.d0
      da0(4,1)=dalo*s(2)*s(4)
      da0(4,2)=dalo*s(3)*s(4)
      da0(4,3)=dalo*s(4)
      da0(4,4)=0.d0
      da2(1,1)=dalo*s(14)
      da2(1,2)=dalo*s(15)
      da2(1,3)=0.d0
      da2(1,4)=0.d0
      da2(2,1)=-(2.d0*dalo*s(2)*s(14)+(1.d0-dalo*s(2)**2)
     **s(15)/s(3))/s(3)
      da2(2,2)=-dalo*s(14)
      da2(2,3)=-dalo*s(2)/s(3)*(s(14)/s(2)-s(15)/s(3))
      da2(2,4)=0.d0
      da2(3,1)=0.d0
      da2(3,2)=0.d0
      da2(3,3)=0.d0
      da2(3,4)=0.d0
      da2(4,1)=dalo*(s(14)*s(4)+s(2)*s(16))
      da2(4,2)=dalo*(s(15)*s(4)+s(3)*s(16))
      da2(4,3)=dalo*s(16)
      da2(4,4)=0.d0
      call osc16(4,da0,bb,u)
      do 130 k=1,4
      do 130 l=1,4
  130 u(k,l)=u(k,l)+da2(k,l)
      call osc17(4,t,ioc,u,dbb)
      do 140 l=1,4
      dcl(1,l)=da0(1,l)
  140 dcl(2,l)=da0(3,l)
      return
      end
c
c     OSC43
c     El�ents de matrice en surface
c
      subroutine osc43(s,l,om2,cl,aa,bb,dcl,daa,dbb,ind)
      implicit double precision(a-h,o-z)
      include 'osc.i'
      dimension s(24),cl(2,4),aa(4,4),bb(4,4),dcl(2,4),daa(4,4),dbb(4,4)
  901 format(' error: the boundary condition parameter has an',
     .' illegal value')
      if(boundary.eq.1)then
        call osc431(s,l,om2,cl,aa,bb,dcl,daa,dbb,ind)
      elseif(boundary.eq.2)then
        call osc432(s,l,om2,cl,aa,bb,dcl,daa,dbb,ind)
      else
        write(6,901)
        call exit(1)
      endif
      return
      end
c
c     OSC431
c     El�ents de matrice en surface
c     La condition limite sur delta P est celle qu'il convient d'utiliser pour
c     les mod�es de temp�ature nulle
c
      subroutine osc431(s,ll,om2,cl,aa,bb,dcl,daa,dbb,ind)
      implicit double precision(a-h,o-z)
      dimension s(24),a0(4,4),a1(4,4),a2(4,4),da0(4,4),da1(4,4),
     *da2(4,4),cl(2,4),aa(4,4),bb(4,4),dcl(2,4),daa(4,4),dbb(4,4),
     *ioc(4),t(4,4),u(4,4),v(4,4),w(4,4)
      al=ll
      al1=al+1.d0
      alo=al*al1/om2
      dalo=-alo/om2
      if(s(3).eq.0.d0)then
      do 110 k=1,4
      do 110 l=1,4
  110 a0(k,l)=0.d0
      a0(2,1)=(om2+(4.d0-alo*s(2))*s(2))/s(1)/s(9)
      a0(2,2)=s(1)*s(2)/s(9)
      a0(2,3)=-alo*s(2)/s(1)/s(9)
      a0(2,4)=-1.d0/s(1)/s(9)
      a1(1,1)=(alo*s(2)-al1)/s(1)
      a1(1,2)=-s(1)/s(5)
      a1(1,3)=alo/s(1)
      a1(1,4)=0.d0
      a1(2,1)=(2.d0*s(8)*(2.d0-alo*s(2))-(om2+(4.d0-alo*s(2))*s(2))
     **(1.d0/s(1)+s(15)/s(9)))/s(1)/s(9)
      a1(2,2)=(s(2)+s(1)*s(8)-s(1)*s(2)*s(15)/s(9))/s(9)
     *-(al+alo*s(2))/s(1)
      a1(2,3)=-alo*(s(8)-s(2)/s(1)-s(2)*s(15)/s(9))/s(1)/s(9)
      a1(2,4)=(1.d0/s(1)+s(15)/s(9))/s(1)/s(9)
      a1(3,1)=-s(4)/s(1)
      a1(3,2)=0.d0
      a1(3,3)=-al/s(1)
      a1(3,4)=1.d0/s(1)
      a1(4,1)=alo*s(2)*s(4)/s(1)
      a1(4,2)=0.d0
      a1(4,3)=(al*al1+alo*s(4))/s(1)
      a1(4,4)=-al1/s(1)
      do 113 l=1,4
      cl(1,l)=a0(2,l)
      cl(2,l)=0.d0
      do 113 k=1,4
      u(k,l)=-a0(k,l)
  113 if(k.eq.l)u(k,l)=u(k,l)+1.d0
      cl(2,3)=al1
      cl(2,4)=1.d0
      call osc14(4,u,t,ioc,d)
      call osc17(4,t,ioc,a1,aa)
      if(ind.eq.1)return
      a2(1,1)=(alo*s(8)-(alo*s(2)-al1)/s(1))/s(1)
      a2(1,2)=(s(1)*s(11)/s(5)-1.d0)/s(5)+alo*s(9)/s(1)
      a2(1,3)=-alo/s(1)**2
      a2(1,4)=0.d0
      a2(2,1)=((om2+(4.d0-alo*s(2))*s(2))*((s(15)/s(9))**2-s(21)/s(9)
     *+(s(15)/s(9)+1.d0/s(1))/s(1))-2.d0*s(8)*(2.d0-alo*s(2))*
     *(1.d0/s(1)+s(15)/s(9))+4.d0*s(14)-alo*(2.d0*s(2)*s(14)
     *+s(8)**2))/s(1)/s(9)
      a2(2,2)=(s(1)*s(14)+s(1)*s(2)*((s(15)/s(9))**2-s(21)/s(9))
     *-s(1)*s(8)*s(15)/s(9)-s(2)*s(15)/s(9)+s(8))/s(9)
     *-(alo*s(8)-(al+alo*s(2))/s(1))/s(1)
      a2(2,3)=-alo/s(1)/s(9)*(s(14)+s(2)/s(1)**2+s(2)*((s(15)/s(9))**2
     *-s(21)/s(9))+s(2)*s(15)/s(1)/s(9)-s(8)*s(15)/s(9)-s(8)/s(1))
      a2(2,4)=-(1.d0/s(1)**2+(s(15)/s(9))**2-s(21)/s(9)+s(15)/s(1)
     */s(9))/s(1)/s(9)
      a2(3,1)=(s(4)/s(1)-s(10))/s(1)
      a2(3,2)=0.d0
      a2(3,3)=al/s(1)**2
      a2(3,4)=-1.d0/s(1)**2
      a2(4,1)=alo/s(1)*(s(4)*s(8)+s(2)*s(10)-s(2)*s(4)/s(1))
      a2(4,2)=alo*s(4)*s(9)/s(1)
      a2(4,3)=(alo*s(10)-(al*al1+alo*s(4))/s(1))/s(1)
      a2(4,4)=al1/s(1)**2
      do 116 l=1,4
      do 116 k=1,4
      u(k,l)=-a0(k,l)
  116 if(k.eq.l)u(k,l)=u(k,l)+2.d0
      call osc14(4,u,t,ioc,d)
      call osc16(4,a1,aa,u)
      do 210 k=1,4
      do 210 l=1,4
  210 u(k,l)=u(k,l)+a2(k,l)
      call osc17(4,t,ioc,u,bb)
      if(ind.eq.2)return
      do 120 k=1,4
      do 120 l=1,4
  120 da0(k,l)=0.d0
      da0(2,1)=(1.d0-dalo*s(2)**2)/s(1)/s(9)
      da0(2,3)=-dalo*s(2)/s(1)/s(9)
      da1(1,1)=dalo*s(2)/s(1)
      da1(1,2)=0.d0
      da1(1,3)=dalo/s(1)
      da1(1,4)=0.d0
      da1(2,1)=(-2.d0*s(2)*s(8)*dalo-(1.d0-dalo*s(2)**2)*
     *(1.d0/s(1)+s(15)/s(9)))/s(1)/s(9)
      da1(2,2)=-dalo*s(2)/s(1)
      da1(2,3)=-dalo*(s(8)-s(2)/s(1)-s(2)*s(15)/s(9))/s(1)/s(9)
      da1(2,4)=0.d0
      do 130 l=1,4
  130 da1(3,l)=0.d0
      da1(4,1)=dalo*s(2)*s(4)/s(1)
      da1(4,2)=0.d0
      da1(4,3)=dalo*s(4)/s(1)
      da1(4,4)=0.d0
      da2(1,1)=dalo*(s(8)-s(2)/s(1))/s(1)
      da2(1,2)=dalo*s(9)/s(1)
      da2(1,3)=-dalo/s(1)**2
      da2(1,4)=0.d0
      da2(2,1)=((1.d0-dalo*s(2)**2)*((s(15)/s(9))**2-s(21)/s(9)
     *+(s(15)/s(9)+1.d0/s(1))/s(1))+2.d0*s(2)*s(8)*dalo*(1.d0/s(1)
     *+s(15)/s(9))-dalo*(2.d0*s(2)*s(14)+s(8)**2))/s(1)/s(9)
      da2(2,2)=dalo*(s(2)/s(1)-s(8))/s(1)
      da2(2,3)=-dalo/s(1)/s(9)*(s(14)+s(2)/s(1)**2+s(2)*((s(15)
     */s(9))**2-s(21)/s(9))+s(2)*s(15)/s(1)/s(9)-s(8)*s(15)/s(9)
     *-s(8)/s(1))
      da2(2,4)=0.d0
      do 140 l=1,4
  140 da2(3,l)=0.d0
      da2(4,1)=dalo/s(1)*(s(4)*s(8)+s(2)*s(10)-s(2)*s(4)/s(1))
      da2(4,2)=dalo*s(4)*s(9)/s(1)
      da2(4,3)=dalo*(s(10)-s(4)/s(1))/s(1)
      da2(4,4)=0.d0
      do 150 l=1,4
      dcl(1,l)=da0(2,l)
  150 dcl(2,l)=0.d0
      call osc16(4,da0,aa,u)
      do 180 k=1,4
      do 180 l=1,4
  180 u(k,l)=u(k,l)+da1(k,l)
      call osc17(4,t,ioc,u,daa)
      call osc16(4,da1,aa,u)
      call osc16(4,a1,daa,v)
      call osc16(4,da0,bb,w)
      do 220 k=1,4
      do 220 l=1,4
  220 u(k,l)=da2(k,l)+u(k,l)+v(k,l)+w(k,l)
      call osc17(4,t,ioc,u,dbb)
      return
      else
      call osc41(1,s,ll,om2,aa,bb,daa,dbb,ind)
      cl(1,1)=(om2+4.d0*s(2)-alo*s(2)**2)/s(1)
      cl(1,2)=s(1)*s(2)
      cl(1,3)=-alo*s(2)/s(1)
      cl(1,4)=-1.d0/s(1)
      cl(2,1)=0.d0
      cl(2,2)=0.d0
      cl(2,3)=al1
      cl(2,4)=1.d0
      if(ind.eq.2)return
      dcl(1,1)=(1.d0-dalo*s(2)**2)/s(1)
      dcl(1,2)=0.d0
      dcl(1,3)=-dalo*s(2)/s(1)
      dcl(1,4)=0.d0
      do 310 k=1,4
  310 dcl(2,k)=0.d0
      return
      endif
      end
c
c     OSC432
c     El�ents de matrice en surface
c     Condition limite delta P = 0
c
      subroutine osc432(s,ll,om2,cl,aa,bb,dcl,daa,dbb,ind)
      implicit double precision(a-h,o-z)
      dimension s(24),a0(4,4),a1(4,4),a2(4,4),da0(4,4),da1(4,4),
     *da2(4,4),cl(2,4),aa(4,4),bb(4,4),dcl(2,4),daa(4,4),dbb(4,4),
     *ioc(4),t(4,4),u(4,4),v(4,4),w(4,4)
      al=ll
      al1=al+1.d0
      alo=al*al1/om2
      dalo=-alo/om2
      if(s(3).eq.0.d0)then
      do 110 k=1,4
      do 110 l=1,4
  110 a0(k,l)=0.d0
      a0(2,1)=(om2+(4.d0-alo*s(2))*s(2))/s(1)/s(9)
      a0(2,2)=s(1)*s(2)/s(9)
      a0(2,3)=-alo*s(2)/s(1)/s(9)
      a0(2,4)=-1.d0/s(1)/s(9)
      a1(1,1)=(alo*s(2)-al1)/s(1)
      a1(1,2)=-s(1)/s(5)
      a1(1,3)=alo/s(1)
      a1(1,4)=0.d0
      a1(2,1)=(2.d0*s(8)*(2.d0-alo*s(2))-(om2+(4.d0-alo*s(2))*s(2))
     **(1.d0/s(1)+s(15)/s(9)))/s(1)/s(9)
      a1(2,2)=(s(2)+s(1)*s(8)-s(1)*s(2)*s(15)/s(9))/s(9)
     *-(al+alo*s(2))/s(1)
      a1(2,3)=-alo*(s(8)-s(2)/s(1)-s(2)*s(15)/s(9))/s(1)/s(9)
      a1(2,4)=(1.d0/s(1)+s(15)/s(9))/s(1)/s(9)
      a1(3,1)=-s(4)/s(1)
      a1(3,2)=0.d0
      a1(3,3)=-al/s(1)
      a1(3,4)=1.d0/s(1)
      a1(4,1)=alo*s(2)*s(4)/s(1)
      a1(4,2)=0.d0
      a1(4,3)=(al*al1+alo*s(4))/s(1)
      a1(4,4)=-al1/s(1)
      do 113 l=1,4
      cl(1,l)=a0(2,l)
      cl(2,l)=0.d0
      do 113 k=1,4
      u(k,l)=-a0(k,l)
  113 if(k.eq.l)u(k,l)=u(k,l)+1.d0
      cl(2,3)=al1
      cl(2,4)=1.d0
      call osc14(4,u,t,ioc,d)
      call osc17(4,t,ioc,a1,aa)
      if(ind.eq.1)return
      a2(1,1)=(alo*s(8)-(alo*s(2)-al1)/s(1))/s(1)
      a2(1,2)=(s(1)*s(11)/s(5)-1.d0)/s(5)+alo*s(9)/s(1)
      a2(1,3)=-alo/s(1)**2
      a2(1,4)=0.d0
      a2(2,1)=((om2+(4.d0-alo*s(2))*s(2))*((s(15)/s(9))**2-s(21)/s(9)
     *+(s(15)/s(9)+1.d0/s(1))/s(1))-2.d0*s(8)*(2.d0-alo*s(2))*
     *(1.d0/s(1)+s(15)/s(9))+4.d0*s(14)-alo*(2.d0*s(2)*s(14)
     *+s(8)**2))/s(1)/s(9)
      a2(2,2)=(s(1)*s(14)+s(1)*s(2)*((s(15)/s(9))**2-s(21)/s(9))
     *-s(1)*s(8)*s(15)/s(9)-s(2)*s(15)/s(9)+s(8))/s(9)
     *-(alo*s(8)-(al+alo*s(2))/s(1))/s(1)
      a2(2,3)=-alo/s(1)/s(9)*(s(14)+s(2)/s(1)**2+s(2)*((s(15)/s(9))**2
     *-s(21)/s(9))+s(2)*s(15)/s(1)/s(9)-s(8)*s(15)/s(9)-s(8)/s(1))
      a2(2,4)=-(1.d0/s(1)**2+(s(15)/s(9))**2-s(21)/s(9)+s(15)/s(1)
     */s(9))/s(1)/s(9)
      a2(3,1)=(s(4)/s(1)-s(10))/s(1)
      a2(3,2)=0.d0
      a2(3,3)=al/s(1)**2
      a2(3,4)=-1.d0/s(1)**2
      a2(4,1)=alo/s(1)*(s(4)*s(8)+s(2)*s(10)-s(2)*s(4)/s(1))
      a2(4,2)=alo*s(4)*s(9)/s(1)
      a2(4,3)=(alo*s(10)-(al*al1+alo*s(4))/s(1))/s(1)
      a2(4,4)=al1/s(1)**2
      do 116 l=1,4
      do 116 k=1,4
      u(k,l)=-a0(k,l)
  116 if(k.eq.l)u(k,l)=u(k,l)+2.d0
      call osc14(4,u,t,ioc,d)
      call osc16(4,a1,aa,u)
      do 210 k=1,4
      do 210 l=1,4
  210 u(k,l)=u(k,l)+a2(k,l)
      call osc17(4,t,ioc,u,bb)
      if(ind.eq.2)return
      do 120 k=1,4
      do 120 l=1,4
  120 da0(k,l)=0.d0
      da0(2,1)=(1.d0-dalo*s(2)**2)/s(1)/s(9)
      da0(2,3)=-dalo*s(2)/s(1)/s(9)
      da1(1,1)=dalo*s(2)/s(1)
      da1(1,2)=0.d0
      da1(1,3)=dalo/s(1)
      da1(1,4)=0.d0
      da1(2,1)=(-2.d0*s(2)*s(8)*dalo-(1.d0-dalo*s(2)**2)*
     *(1.d0/s(1)+s(15)/s(9)))/s(1)/s(9)
      da1(2,2)=-dalo*s(2)/s(1)
      da1(2,3)=-dalo*(s(8)-s(2)/s(1)-s(2)*s(15)/s(9))/s(1)/s(9)
      da1(2,4)=0.d0
      do 130 l=1,4
  130 da1(3,l)=0.d0
      da1(4,1)=dalo*s(2)*s(4)/s(1)
      da1(4,2)=0.d0
      da1(4,3)=dalo*s(4)/s(1)
      da1(4,4)=0.d0
      da2(1,1)=dalo*(s(8)-s(2)/s(1))/s(1)
      da2(1,2)=dalo*s(9)/s(1)
      da2(1,3)=-dalo/s(1)**2
      da2(1,4)=0.d0
      da2(2,1)=((1.d0-dalo*s(2)**2)*((s(15)/s(9))**2-s(21)/s(9)
     *+(s(15)/s(9)+1.d0/s(1))/s(1))+2.d0*s(2)*s(8)*dalo*(1.d0/s(1)
     *+s(15)/s(9))-dalo*(2.d0*s(2)*s(14)+s(8)**2))/s(1)/s(9)
      da2(2,2)=dalo*(s(2)/s(1)-s(8))/s(1)
      da2(2,3)=-dalo/s(1)/s(9)*(s(14)+s(2)/s(1)**2+s(2)*((s(15)
     */s(9))**2-s(21)/s(9))+s(2)*s(15)/s(1)/s(9)-s(8)*s(15)/s(9)
     *-s(8)/s(1))
      da2(2,4)=0.d0
      do 140 l=1,4
  140 da2(3,l)=0.d0
      da2(4,1)=dalo/s(1)*(s(4)*s(8)+s(2)*s(10)-s(2)*s(4)/s(1))
      da2(4,2)=dalo*s(4)*s(9)/s(1)
      da2(4,3)=dalo*(s(10)-s(4)/s(1))/s(1)
      da2(4,4)=0.d0
      do 150 l=1,4
      dcl(1,l)=da0(2,l)
  150 dcl(2,l)=0.d0
      call osc16(4,da0,aa,u)
      do 180 k=1,4
      do 180 l=1,4
  180 u(k,l)=u(k,l)+da1(k,l)
      call osc17(4,t,ioc,u,daa)
      call osc16(4,da1,aa,u)
      call osc16(4,a1,daa,v)
      call osc16(4,da0,bb,w)
      do 220 k=1,4
      do 220 l=1,4
  220 u(k,l)=da2(k,l)+u(k,l)+v(k,l)+w(k,l)
      call osc17(4,t,ioc,u,dbb)
      return
      else
      call osc41(1,s,ll,om2,aa,bb,daa,dbb,ind)
      cl(1,1)=0d0
      cl(1,2)=1d0
      cl(1,3)=0d0
      cl(1,4)=0d0
      cl(2,1)=0.d0
      cl(2,2)=0.d0
      cl(2,3)=al1
      cl(2,4)=1.d0
      if(ind.eq.2)return
      dcl(1,1)=0d0
      dcl(1,2)=0.d0
      dcl(1,3)=0d0
      dcl(1,4)=0.d0
      do 310 k=1,4
  310 dcl(2,k)=0.d0
      return
      endif
      end
c
c     OSC44
c     Matrice du probl�e non radial
c
      subroutine osc44(np,star,st1,st2,ll,om2,aco,daco,ind)
      implicit double precision(a-h,o-z)
      dimension star(12,np),st1(18),st2(24),aco(1),daco(1),
     *cl(2,4),aa(4,4),bb(4,4),dcl(2,4),daa(4,4),dbb(4,4),c(4,4),
     *d(4,4),dc(4,4),dd(4,4)
      if(ind.eq.2)goto 100
      kk=0
      do 260 j=1,np
      if(j.eq.1)then
      if(star(1,j).eq.0d0)then
      call osc42(st1,ll,om2,cl,bb,dcl,dbb,2)
      do 210 k=1,4
      do 210 l=1,4
  210 aa(k,l)=0.d0
      else
      call osc41(j,star,ll,om2,aa,bb,daa,dbb,2)
      do l=1,2
      do k=1,4
      cl(l,k)=0d0
      enddo
      enddo
      cl(1,1)=1d0
      cl(2,3)=-ll
      cl(2,4)=1d0
      endif
      do 225 k=1,2
      do 220 l=1,4
      kk=kk+1
  220 aco(kk)=0.d0
      do 225 l=1,4
      kk=kk+1
  225 aco(kk)=cl(k,l)
      goto 250
      endif
      if(j.eq.np)then
      call osc43(st2,ll,om2,cl,aa,bb,dcl,daa,dbb,2)
      goto 240
      endif
      call osc41(j,star,ll,om2,aa,bb,daa,dbb,2)
  240 h=star(1,j)-star(1,j-1)
      h2=h/2.d0
      h6=h**2/6.d0
      do 270 k=1,4
      do 270 l=1,4
      d(k,l)=-h2*aa(k,l)+h6*bb(k,l)
  270 if(k.eq.l)d(k,l)=d(k,l)+1.d0
      do 280 k=1,4
      do 275 l=1,4
      kk=kk+1
  275 aco(kk)=-c(k,l)
      do 280 l=1,4
      kk=kk+1
  280 aco(kk)=d(k,l)
      if(j.eq.np)then
      do 290 k=1,2
      do 285 l=1,4
      kk=kk+1
  285 aco(kk)=cl(k,l)
      do 290 l=1,4
      kk=kk+1
  290 aco(kk)=0.d0
      goto 260
      endif
  250 h=star(1,j+1)-star(1,j)
      h2=h/2.d0
      h6=h**2/6.d0
      do 255 k=1,4
      do 255 l=1,4
      c(k,l)=h2*aa(k,l)+h6*bb(k,l)
  255 if(k.eq.l)c(k,l)=c(k,l)+1.d0
  260 continue
      return
  100 kk=0
      do 160 j=1,np
      if(j.eq.1)then
      if(star(1,j).eq.0d0)then
      call osc42(st1,ll,om2,cl,bb,dcl,dbb,3)
      do 103 k=1,4
      do 103 l=1,4
      aa(k,l)=0.d0
  103 daa(k,l)=0.d0
      else
      call osc41(j,star,ll,om2,aa,bb,daa,dbb,3)
      do l=1,2
      do k=1,4
      cl(l,k)=0d0
      dcl(l,k)=0d0
      enddo
      enddo
      cl(1,1)=1d0
      cl(2,3)=-ll
      cl(2,4)=1d0
      endif
      do 102 k=1,2
      do 101 l=1,4
      kk=kk+1
      aco(kk)=0.d0
  101 daco(kk)=0.d0
      do 102 l=1,4
      kk=kk+1
      aco(kk)=cl(k,l)
  102 daco(kk)=dcl(k,l)
      goto 150
      endif
      if(j.eq.np)then
      call osc43(st2,ll,om2,cl,aa,bb,dcl,daa,dbb,3)
      goto 140
      endif
      call osc41(j,star,ll,om2,aa,bb,daa,dbb,3)
  140 h=star(1,j)-star(1,j-1)
      h2=h/2.d0
      h6=h**2/6.d0
      do 142 k=1,4
      do 141 l=1,4
      d(k,l)=-h2*aa(k,l)+h6*bb(k,l)
  141 dd(k,l)=-h2*daa(k,l)+h6*dbb(k,l)
  142 d(k,k)=d(k,k)+1.d0
      do 145 k=1,4
      do 143 l=1,4
      kk=kk+1
      aco(kk)=-c(k,l)
  143 daco(kk)=-dc(k,l)
      do 144 l=1,4
      kk=kk+1
      aco(kk)=d(k,l)
  144 daco(kk)=dd(k,l)
  145 continue
      if(j.eq.np)then
      do 147 k=1,2
      do 146 l=1,4
      kk=kk+1
      aco(kk)=cl(k,l)
  146 daco(kk)=dcl(k,l)
      do 147 l=1,4
      kk=kk+1
      aco(kk)=0.d0
  147 daco(kk)=0.d0
      goto 160
      endif
  150 h=star(1,j+1)-star(1,j)
      h2=h/2.d0
      h6=h**2/6.d0
      do 152 k=1,4
      do 151 l=1,4
      c(k,l)=h2*aa(k,l)+h6*bb(k,l)
  151 dc(k,l)=h2*daa(k,l)+h6*dbb(k,l)
  152 c(k,k)=c(k,k)+1.d0
  160 continue
      return
      end
c
c     OSC45
c     Calcul d'un mode non radial
c
      subroutine osc45(im,om0,ir)
      implicit double precision(a-h,o-z)
      include 'osc.i'
      common io(4*npmax),x(4*npmax),x1(4*npmax),aco(32*npmax),
     $daco(32*npmax)
      parameter(deuxpi=6.28318530717959d0)
      character c
      data e1,e2,e3/1.d-8,1.d-2,1.d-4/
  901 format(3x,'it= 0',3x,'omega=',1pd12.5)
  905 format(3x,a1,1x,'it=',i2,3x,'omega=',1pd12.5,3x,'corr=',d12.5)
  907 format('   l=',i3,'   omega=',1pd15.8,'   sigma=',d15.8/
     *'   nz=',i4,'   frequency=',d15.8,'   period=',
     *d15.8/'   rotational splitting beta=',d15.8,'   parity=',i1/
     *'   ev=',d15.8,'   xm=',d15.8,'   delta=',d15.8/
     *'   mode=',i5,',',i5)
  911 format('   *** error : too many iterations')
      if(lscreen)write(6,901)om0
      call flush(6)
      om=om0
      om2=om**2
      if(om.lt.0.d0)om2=-om2
      iss=0
      np4=4*np
      js=np4-3
      do 1305 j=1,np4
      x(j)=1d0
 1305 x1(j)=x(j)
 1310 omega2=om2
      call osc44(np,star,st1,st2,ll,omega2,aco,daco,2)
      call osc11(np,4,2,aco,io,det,idet)
      c='*'
      is=0
 1320 is=is+1
      iss=iss+1
      if(iss.gt.20)goto 1370
      call osc13(np,4,2,daco,x)
      call osc12(np,4,2,aco,io,x)
      s1=0.d0
      s2=0.d0
      xx=x(np)
      do 1322 j=1,np4
      x(j)=x(j)/xx
      s1=s1+x1(j)*x(j)
      s2=s2+x(j)**2
 1322 x1(j)=x(j)
      dom2=-s1/s2/xx
      om2a=om2
      om2=omega2+dom2
      oma=om
      om=dsqrt(dabs(om2))
      if(om2.lt.0.d0)om=-om
      corr=om-oma
      if(lscreen)write(6,905)c,iss,om,corr
      call flush(6)
      c=' '
c     if(iss.ge.2.and.dabs((om2-om2a)/omega2).le.e1.and.
c    *dabs(dom2/om2).le.e3)goto 1350
      if(iss.ge.2.and.dabs((om2-om2a)/omega2).le.e1)goto 1350
      if(im.eq.1.and.is.ge.2.and.dabs((om2-om2a)/dom2).le.e2)
     *goto 1310
      goto 1320
 1350 continue
      omeg2=om**2
      sigma=om/tdyn
      freq=sigma/deuxpi
      periode=1d0/freq
      call osc46(np,star,st2,ll,om,x,nz,mode,xm,del,ev,beta,sol)
      call osc48
      if(lscreen)write(6,907)ll,om,sigma,nz,freq,periode,beta,parity,
     .ev,xm,del,mode,modeLee
      call flush(6)
      ir=0
      return
 1370 if(lscreen)write(6,911)
      ir=2
      return
      end
c
c     OSC46
c     Analyse d'un mode non radial
c
      subroutine osc46(np,star,st2,l,omega,x,nz,mode,xm,del,ev,beta,sol)
      implicit double precision(a-h,o-z)
      dimension star(12,1),sol(9,1),st2(24),cl(2,4),aa(4,4),
     *bb(4,4),dcl(2,4),daa(4,4),dbb(4,4),x(4,1)
      nz=0
      mode=0
      erot=0d0
      ev=0.d0
      eh=0.d0
      sx=0.d0
      sx2=0.d0
      omega2=omega**2
      if(omega.lt.0.d0)omega2=-omega2
      a2l1=2*l+1
      c=dsqrt(dfloat(l*(l+1)))
      do 120 j=1,np
      do 210 k=1,4
  210 sol(k,j)=x(k,j)
c     if(j.eq.1)then
      if(star(1,j).eq.0d0)then
      do 102 m=5,8
  102 sol(m,j)=0.d0
      else
      if(j.eq.np)then
      call osc43(st2,l,omega2,cl,aa,bb,dcl,daa,dbb,1)
      else
      call osc41(j,star,l,omega2,aa,bb,daa,dbb,1)
      endif
      do 105 k=1,4
      s=0.d0
      do 103 m=1,4
  103 s=s+aa(k,m)*sol(m,j)
  105 sol(k+4,j)=s
      endif
      xx=star(1,j)
      xx0=xx**a2l1
      xx1=xx0*xx
      xx2=xx1*xx
      r=sol(1,j)
      w=(star(2,j)*sol(1,j)+star(3,j)*sol(2,j)+sol(3,j))/omega2
      t=c*w
      frot=star(4,j)*(w**2+2d0*r*w)
      fv=star(4,j)*r**2
      fh=star(4,j)*t**2
      f=fv+fh
      sol(9,j)=f-frot
      if(j.eq.1)goto 110
      if(r*r1.lt.0.d0.or.r.eq.0.d0)then
      nz=nz+1
      mode=mode-1
      if((r-r1)*t.lt.0.d0)mode=mode+2
      endif
      ffrot=frot1+frot
      ffv=fv+fv1
      ffh=fh+fh1
      ff=ffv+ffh
      erot=erot+ffrot*(xx0-xx01)
      ev=ev+ffv*(xx0-xx01)
      eh=eh+ffh*(xx0-xx01)
      sx=sx+ff*(xx1-xx11)
      sx2=sx2+ff*(xx2-xx21)
  110 r1=r
      w1=w
      t1=t
      xx01=xx0
      xx11=xx1
      xx21=xx2
      frot1=frot
      fv1=fv
  120 fh1=fh
      if(omega.lt.0.d0)mode=mode-1
      a2=ev+eh
      erot=erot/a2
      ev=ev/a2
      eh=eh/a2
      a2=a2/2.d0/a2l1
      sx=sx/2.d0/(2*l+2)/a2
      sx2=sx2/2.d0/(2*l+3)/a2
      xm=sx
      del=2.d0*dsqrt(sx2-sx**2)
      beta=1d0-erot
      a=dsqrt(a2)
      if(sol(1,1).lt.0.d0)a=-a
      do j=1,np
        do k=1,8
          sol(k,j)=sol(k,j)/a
        enddo
        sol(9,j)=sol(9,j)/a2
      enddo
      return
      end
c
c     OSC47
c     Balayage en fr�uence, l>0
c     On stocke dans le vecteur omarr les m premi�es valeurs propres approch�s
c     obtenues au cours du balayage. Au retour m = le nombre de valeurs stock�s
c     dans le vecteur omarr.
c
      subroutine osc47(om1,om2,n,ind,ir,m,omarr,k)
      implicit double precision(a-h,o-z)
      include 'osc.i'
      real*8 omarr(*)
      common aco(32*npmax),daco(32*npmax),io(4*npmax)
  901 format(' error : om1 and om2 have opposite signs')
  902 format(' l=',i3,3x,'omega=',1pd12.5,3x,'det=',0pf8.5,' (',i3,')')
  903 format(26x,'interpolated value : omega=',1pd12.5)
      if(om1*om2.le.0.d0)then
      if(lscreen)write(6,901)
      ir=2
      return
      endif
      k=0
      do 200 i=1,n+1
      goto(110,120),ind
  110 om=((n+1-i)*om1+(i-1)*om2)/n
      goto 130
  120 om=n/((n+1-i)/om1+(i-1)/om2)
  130 omega2=om**2
      if(om.lt.0.d0)omega2=-omega2
      call osc44(np,star,st1,st2,ll,omega2,aco,daco,1)
      call osc11(np,4,2,aco,io,det,idet)
      if(i.gt.1)then
      sdet=dsign(1.d0,det)
      sdeto=dsign(1.d0,deto)
      if(sdet*sdeto.le.0.d0)then
      det1=det
      deto1=deto
      if(idet.gt.ideto)deto1=deto*1.d1**(ideto-idet)
      if(idet.lt.ideto)det1=det*1.d1**(idet-ideto)
      omi=(det1*omo-deto1*om)/(det1-deto1)
      if(lscreen)write(6,903)omi
      call flush(6)
      if(k.lt.m)then
        k=k+1
        omarr(k)=omi
      endif
      endif
      endif
      if(lscreen)write(6,902)ll,om,det,idet
      call flush(6)
      omo=om
      deto=det
      ideto=idet
  200 continue
      ir=0
      return
      end

c===============================================================================

c   OSC48
c   Sign of the derivative of the determinant with respect to omega2 gives us
c   parity. We compute also the mode number according to Lee prescription.

      subroutine osc48
      implicit none
      include 'osc.i'
      integer io,idet1,idet2,j
      real*8 omega21,omega22,aco,daco,det1,det2,a,b,c,d,aa,bb,cc,dd
      common aco(32*npmax),daco(32*npmax),io(4*npmax)
      omega21=omeg2*(1d0-1d-6)
      call osc44(np,star,st1,st2,ll,omega21,aco,daco,1)
      call osc11(np,4,2,aco,io,det1,idet1)
      omega22=omeg2*(1d0+1d-6)
      call osc44(np,star,st1,st2,ll,omega22,aco,daco,1)
      call osc11(np,4,2,aco,io,det2,idet2)
      if(det1-det2.lt.0d0)then
        parity=1
      else
        parity=0
      endif
      a=star(2,1)*sol(1,1)+sol(3,1)
      b=star(2,1)*sol(1,1)+star(3,1)*sol(2,1)+sol(3,1)
      modeLee=0
      do j=2,np
        aa=a
        bb=b
        a=star(2,j)*sol(1,j)+sol(3,j)
        b=star(2,j)*sol(1,j)+star(3,j)*sol(2,j)+sol(3,j)
        if(a*aa.lt.0d0.or.a.eq.0d0)then
          modeLee=modeLee-1
          if((a-aa)*(b+bb).lt.0d0)modeLee=modeLee+2
        endif
      enddo
      if(ll.eq.1.and.modeLee.ge.0)modeLee=modeLee+1
      return
      end

c===============================================================================

c   modinfo fournit des informations sur le mod�e lu

      subroutine modinfo(mass,radius,taudyn,npts,rc)
Cf2py intent(out) mass
Cf2py intent(out) radius
Cf2py intent(out) taudyn
Cf2py intent(out) npts
Cf2py intent(out) rc
      implicit none
      include 'osc.i'
      integer npts,rc
      real*8 mass,radius,taudyn
      if(lmod.and.lres)then
        mass=am
        radius=ray
        taudyn=tdyn
        npts=np
        rc=0
      else
        mass=0d0
        radius=0d0
        taudyn=0d0
        npts=0
        rc=1
      endif
      return
      end

c===============================================================================

c   oscinfo fournit des informations sur l'oscillation calcul�

      subroutine oscinfo(l,nz_,mode_,modeLee_,parity_,omega,sigma_,
     .beta_,ev_,xm_,delta_,boundary_,rc)
Cf2py intent(out) l
Cf2py intent(out) nz_
Cf2py intent(out) mode_
Cf2py intent(out) modeLee_
Cf2py intent(out) parity_
Cf2py intent(in) omega
Cf2py intent(out) sigma_
Cf2py intent(out) beta_
Cf2py intent(out) ev_
Cf2py intent(out) xm_
Cf2py intent(out) delta_
Cf2py intent(out) boundary_
Cf2py intent(in) rc
      implicit none
      include 'osc.i'
      integer l,nz_,mode_,modeLee_,parity_,rc,boundary_
      real*8 omega,sigma_,beta_,ev_,xm_,delta_
      if(losc)then
        l=ll
        nz_=nz
        mode_=mode
        modeLee_=modeLee
        parity_=parity
        omega=om
        sigma_=sigma
        beta_=beta
        ev_=ev
        xm_=xm
        delta_=del
        boundary_=boundary
        rc=0
      else
        l=0
        nz_=0
        mode_=0
        modeLee_=0
        parity_=0
        omega=0d0
        sigma_=0d0
        beta_=0d0
        ev_=0d0
        xm_=0d0
        delta_=0d0
        boundary_=0
        rc=1
      endif
      return
      end

c===============================================================================
