c $Header: /cvs/ppgen/src/dnlcc7LAPACK.f,v 1.2 2006/08/30 21:22:56 jelen Exp $
cmf 21-06-03 explicit initialization of a0 etc. to zero 
cmf          modified iteration procedure of a0, see a0inc
c**********************************************************************
c
c density for nonlinear pseudocore correction
c
c input
c mmax......working space dimension
c rnlc......matching radius for pseudocore
c r().......radial mesh
c dps().....pseudo valence density 
c dc()......full core density
c dcp().....	""		1st derivative
c dcpp()....	""		2nd	""	
c io .......unit for monitoring output
c
c output
c dc() = modified for r < rnlc
c dc() = full core density r > rnlc
c etc. for derivatives
c
c density is integrated over solid angle
c
c currently: x^6 polynomial 
c zero 1st and 2nd derivative of density at origin
c continuity up to 3rd derivative enforced
c monotonously decaying 
c
c Martin Fuchs, FHI der MPG, Berlin, 12-1995
c**********************************************************************
      subroutine dnlcc7(mmax,rnlc,r,dps,dc,dcp,dcpp,iu)
c
      implicit none
      include  'parameter.h'
      integer  iu,mmax,i,inlc,it,itmx,mstep
      integer  ipvt(4)
      real*8   rnlc,a0,a1,a2,a3,a4,a5,a6,x,fmom,f0,f1,f2,sc,scp,scpp
      real*8   al,a0inc,t0,dummy,dummy1
      real*8   r(mmax),dps(mmax),dc(mmax),dcp(mmax),dcpp(mmax),dd(mx)
      real*8   b(4),a(4,4)
      parameter (itmx=50,a0inc=1.25d0)
      external fmom

      f0(x)=a0+a1*x+a2*x*x+a3*x*x*x+a4*x**4+a5*x**5+a6*x**6
      f1(x)=a1+2.d0*a2*x+3.d0*a3*x*x+4.d0*a4*x**3
     1     +5.d0*a5*x**4+6.d0*a6*x**5
      f2(x)=2.d0*a2+6.d0*a3*x+12.d0*a4*x**2+20.d0*a5*x**3
     1     +30.d0*a6*x**4
c
c initialize
      a0=0.d0
      a1=0.d0
      a2=0.d0
      a3=0.d0
      a4=0.d0
      a5=0.d0
      a6=0.d0
c
c pseudocore cutoff radius
      do i=mmax,1,-1
        if(r(i) .le. rnlc)  goto 30
      enddo
      write(ie,*) '& dnlcc7 --- stop: pseudocore radius not found'
      write(ie,*) '& dnlcc7 --- rnlc',rnlc
      write(ie,*) '& dnlcc7 --- r   ',r(i),i
      stop 
  30  inlc=i
      rnlc=r(inlc)

c pseudocore function coefficients
c loop: want (i)  monotonously decaying pseudocore for gga
c            (ii) small core charge
c start dc(r=0)=a0 [violation of (i)] & gradually increase until (i) & (ii) are met
      a0=dc(inlc)
      do it=1,itmx
!       write(ie,*) '& dnlcc7 --- a0',a0/dc(inlc),f0(r(i))/dc(inlc),i

cmf replace
c       a0=1.8*f0(r(i))
        a0=a0*a0inc

        a1=0.d0
        a2=0.d0
        a3=(dcpp(inlc+1)-dcpp(inlc-1))/(r(inlc+1)-r(inlc-1))/rnlc
        a(1,1)=1.d0
        a(1,2)=4.d0*rnlc
        a(1,3)=10.d0*rnlc**2
        a(1,4)=20.d0*rnlc**3
        a(2,1)=3.d0
        a(2,2)=6.d0*rnlc
        a(2,3)=10.d0*rnlc**2
        a(2,4)=15.d0*rnlc**3
        a(3,1)=1.d0
        a(3,2)=rnlc
        a(3,3)=rnlc**2
        a(3,4)=rnlc**3
        a(4,1)=3.d0
        a(4,2)=4*rnlc
        a(4,3)=5*rnlc**2
        a(4,4)=6*rnlc**3

        b(1)=a3/6.d0
        b(2)=dcpp(inlc)/(2*rnlc)-a2/rnlc
        b(3)=(dc(inlc)-a0-a2*rnlc**2)/rnlc**3
        b(4)=dcp(inlc)/rnlc**2-2*a2/rnlc
 
c numerical recipies
c       call ludcmp(a,4,4,ipvt,dummy)
c       call lubksb(a,4,4,ipvt,b)

c essl
c       call dgef(a,4,4,ipvt)
c       call dges(a,4,4,ipvt,b,0)

c lapack
        call dgesv(4,1,a,4,ipvt,b,4,i)
        if(i .ne. 0) then
          if(i .lt. 0) write(ie,*) '& dnlcc7 - stop: singular matrix'      
          if(i .gt. 0) write(ie,*) '& dnlcc7 - stop: bad input'      
          stop
        endif

        a3=b(1)
        a4=b(2)
        a5=b(3)
        a6=b(4)

c debug
!       write(ie,*) '& dnlcc7 --- rnlc = ',r(inlc)
!       write(ie,*) '& dnlcc7 --- a0   = ',a0
!       write(ie,*) '& dnlcc7 --- a1   = ',a1
!       write(ie,*) '& dnlcc7 --- a2   = ',a2
!       write(ie,*) '& dnlcc7 --- a3   = ',a3
!       write(ie,*) '& dnlcc7 --- a4   = ',a4
!       write(ie,*) '& dnlcc7 --- a5   = ',a5
!       write(ie,*) '& dnlcc7 --- a6   = ',a6

        dummy=f0(r(inlc))
        do i=inlc,2,-1
cmf replace
c         if(f0(r(i)) .gt. a0) goto 45
          dummy1=f0(r(i-1))
          if(dummy1 .lt. dummy) goto 45
          dummy=dummy1
        enddo
  45    continue
      if(i .eq. 1) goto 50
      enddo
      write(ie,*) '& dnlcc7 - stop: pseudocore convergence error',it 
      stop
  50  continue

      sc=dc(inlc)
      scp=dcp(inlc)
      scpp=dcpp(inlc)

c save core density beyond matching radius
      do i=1,inlc
        dd(i)=dc(i)
      enddo
c
      do i=1,inlc
        dd(i)=dc(i)
        dc(i)=f0(r(i))
        dcp(i)=f1(r(i))
        dcpp(i)=f2(r(i))
      enddo

c debug
!     open(iu)
!     write(iu,*) 'Improved Polynomial Pseudocore -- dnlcc7'
!     write(iu,*) 'Pseudocore matching at i',inlc,' r ',r(inlc)
!     write(iu,*) 'Model parameters:'
!     write(iu,*) 'a0 ',a0
!     write(iu,*) 'a1 ',a1
!     write(iu,*) 'a2 ',a2
!     write(iu,*) 'a3 ',a3
!     write(iu,*) 'a4 ',a4
!     write(iu,*) 'a5 ',a5
!     write(iu,*) 'a6 ',a6
c
!c check matching of density and derivatives
!     write(iu,*) 'Checking values at pseudocore radius:'
!     write(iu,*) sc,dc(inlc),'density'
!     write(iu,*) scp,dcp(inlc),'1st derivative'
!     write(iu,*) scpp,dcpp(inlc),'2nd derivative'
!
!c sum rule checks
!     al=0.1d0*log(r(11)/r(1))
!     t0=fmom(0,mmax,al,1.d0,r,dc)
!     write(iu,*) t0
!    &,		 'total pseudo charge'
!     write(iu,*) -fmom(1,mmax,al,1.d0,r,dcp)/(3*t0)
!    &,		 '1st derivative test, should be = 1'
!c
!     write(iu,*) fmom(2,mmax,al,1.d0,r,dcpp)/(12*t0)
!    &,	 '2nd derivative test, should be = 1 (test may be inaccurate)'
!c
!     do i=1,inlc+10
!       write(iu,'(a1,4(1x,e12.6))') 'p',r(i),dps(i),dc(i),dd(i)
!     enddo
!     close(iu)
c
      return
      end
