c qlqgrid.f. This program really just interpolates the
c files. We do a nasty trick here though. We change separable
c pseudopotential to U*R instead of U*u(r).
c See error correction Oct. 7, 1998 below.
C
C THIS PROGRAM SETS UP THE ONE-DIMENSIONAL GRID QLG.DAT
c
c The program computes:
c QL(g) = integral(0,infinity) r**2 dr jL(gr) QL(r)
c
c where QL(r) is the Kleinman-bylander seperable potetnial for angular
c momentum L.
c
c This program does one QL at a time.
c
c The results QL(g) are written to a file containing nq points.
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        REAL*8 JL

        logical exst

	parameter (nrmax=10001, nqmax=10001)
	CHARACTER*30 filein
	dimension r(nrmax),qlr(nrmax),rLOG(nrmax),qlrLOG(nrmax),
     1	simpson(nrmax),q(nqmax),qlq(nqmax)
	dimension cl(0:5)
c ===============================================================
C constants
        PI=4.D0*DATAN(1.D0)
	EQ2=14.39975D0
	abohr=0.529177D0
	hartree=EQ2/ABOHR
c ===============================================================
	write(*,*)'  '
        WRITE(*,*)'                      WELCOME'
        WRITE(*,*)' PROGRAM QLG grid: Bessel transform of NON-LOCAL', 
     1  ' separable QL potential'
        WRITE(*,*)'  '
        WRITE(*,*)' '
C ===========================================================
c
c First we read in qlclATOMIC.dat and change it to qlcl (eV) units.
	open(unit=30,file='qlclATOMIC.dat',status='old')
	do L=0,5
	read(30,*,end=314)cl(L)
	write(*,*)' cl(L)=',cl(L),' 1/Hartree'
	LOK=L
	end do
314	continue
	write(*,*)'  LOK=',LOK
	close(unit=30)
	write(*,*)' We read in to L=',LOK
	write(*,*)' Now convert to 1/eV and write to qlcl.dat'
	open(unit=30,file='qlcl.dat',status='unknown')
	do L=0,LOK
	cl(L)=cl(L)/hartree
	write(*,*)' cl(L)=',cl(L)
	write(30,818)cl(L),L
818	format(f14.7,'                         L=',i3)
	end do
	close(unit=30)

	
C READ IN QL PSEUDOPOTENTIAL. Insert file ppkb0.dat, or ppkb1.dat , or ...
	write(*,*)' '
c	write(*,*)' insert the angular momentum L (0,1,2,..)'
c	read(*,*)L
	do 999 L=0,LOK
	write(*,*)'  '
	write(*,*)' Begin Next L value.......................'
	write(*,*)' L=',L
	write(*,*)'  '
	write(*,*)' Use one of the QL data files output from gncpp'
	write(*,*)' the files are ppkb0.dat or ppkb1.dat or ....'
c
        write(*,*)' UPDATE: This program was updated June 25, 1992'
        write(*,*)' We now read in eV-A units of everything from'
        write(*,*)' gncpp.f --- No conversions need to be made here!'
c
	filein='xxxxxx'
	if(L.eq.0)filein='vkb0.dat'
	if(L.eq.1)filein='vkb1.dat'
	if(L.eq.2)filein='vkb2.dat'
	if(L.eq.3)filein='vkb3.dat'
	write(*,89)filein
89	format('  ',a30)
	write(*,*)'  '
	write(*,*)' we read in ql(r) in atomic units.'
	write(*,*)' We now convert to eV,A.'
        inquire(file=filein,exist=exst)
         if(exst.eqv..true.)then
	open(unit=51,file=filein,status='old')
	else
	write(*,889)filein
889	format(' File not existing:',a30)
	write(*,*)' Skip to next L'
	go to 999
	end if
	
	
	
c Read in data and determine the number of points.
c NOTE: The non-local potential is in Hartree and Abohr units.
c
c We start off with the first point r=0
	rlog(1)=0.
	qlrLOG(1)=0.0
	ic=1

	do 10 i=2,nrmax
	read(51,*,end=11)rLOG(i),qlrLOG(i)
	if(i.eq.2.or.i.eq.3)write(*,*)' rLOG(i),qlrLOG=',rLOG(i),qlrLOG(i)
c ******* error correction Oct. 9, 1998.
c JJDong Oct 09, 1998
c	rLOG(i)=rLOG(i)*abohr
c	qlrLOG(i)=qlrLOG(i)*hartree/sqrt(abohr)
c JJDong, I think the original is V_kb*U, so I changed it to V_kb*U/r,
c which is V_kb*R(r)
 	if (rLOG(i).le.1.E-6) then
	 qlrLOG(i)=0.0
	else
	 qlrLOG(i)=qlrLOG(i)/rLOG(i)
	endif
c ******** end error correction.
	rLOG(i)=rLOG(i)*abohr
	qlrLOG(i)=qlrLOG(i)*hartree/sqrt(abohr)/abohr
	ic=ic+1
10	continue
11	continue
	close(unit=51)

	nr=ic
	write(*,*)'  '
	write(*,*)' read in nr points=',nr

c Now we interpolate to put on a linear grid
	rmax=rLOG(ic)
	write(*,*)' rmax=',rmax,' in Ang units.'

	if(L.eq.0)open(unit=30,file='vkblinear0.dat',status='unknown')
	if(L.eq.1)open(unit=30,file='vkblinear1.dat',status='unknown')
	if(L.eq.2)open(unit=30,file='vkblinear2.dat',status='unknown')
	if(L.eq.3)open(unit=30,file='vkblinear3.dat',status='unknown')
	if(L.eq.4)open(unit=30,file='vkblinear4.dat',status='unknown')
	do i=1,nrmax-1
	x=0.9999999d0*rmax*(float(i-1))/float(nrmax-1-1)
	norder=5
        call polyint(rLOG,qlrLOG,nrmax,nr,norder,x,y,dy)
	r(i)=x
	qlr(i)=y
	write(30,*)x,y
	end do
	nr=nrmax-1
	close(unit=30)
c
c JJDong
c setup here not to do FF transform
	iff=0
	if (iff.eq.0) then
	print *, 'L =', L,  '  JJ: NO FF Transform now'
	goto 999
	endif
c JJDong
c
        write(*,*)' caution ****** **** caution'
        write(*,*)' we read in ql(r)*r not ql(r)'
        do 1010 i=1,nr
        if(r(i).gt.0.0001)qlr(i)=qlr(i)/r(i)
1010      continue


	write(*,*)'  '
c we use nrmax-1 since for simpsons rule we need an odd number of
c points and sometimes we pad the data with an extra points of zero.
	if(nr.gt.nrmax-1)stop' too many points!'
	rmin=r(1)
	rmax=r(nr)
	dr=r(2)-r(1)
c the points are assumed to be equally spaced. lets check!
	if(dabs(  ((rmax-rmin)/dfloat(nr-1)) -dr ).gt.0.000001)
     1	stop' points are not equally space!'
c ===============================================================
c we now determine how many q points and how large is qmax.
	write(*,*)'  '
	write(*,*)' we now determine the q-grid. We need to decide'
	write(*,*)' how many q points and how large qmax is.'
	write(*,*)'  '
	write(*,*)'  insert # of q points. 1001 is a good choice'
	nq=3001
c	read(*,*)nq
	write(*,*)' nq=',nq
        WRITE(*,*)'  '
	write(*,*)' Determine qmax. use hbar**2 qmax**2 /2m =ecut'
	write(*,*)' to determine ecut.'
        WRITE(*,*)' INSERT ENERGY CUTOFF IN EV.'
c	READ(*,*)ECUT
	ecut=2000.00
        WRITE(*,*)' ECUT=',ECUT
C FACTR= HBAR**2/2M=7.62/2. EV*ANG**2
        FACTR=7.62D0/2.
        qmax=DSQRT( ECUT/FACTR)
        WRITE(*,*)'  '
        WRITE(*,*)' qmax (Angstrom)=',qmax
        WRITE(*,*)'  '
        WRITE(*,*)' Please be patient. This will take 2 minutes. '
	qmin=0.0
c ===============================================================
c set up factors for simpson's rule.
	do 22 ir=1,nr
	simpson(ir)=dr*4.d0/3.d0
	if(mod(ir,2).eq.1)simpson(ir)=dr*2.d0/3.d0
22	continue
	simpson(1)=dr*1.d0/3.d0
c nr must be odd for simpsons rule to work. If it is not odd, then we
c add one more point withthe function being zero!
	if(mod(nr,2).eq.0)then
	nr=nr+1
	r(nr)=rmax+dr
	qlr(nr)=0.0d0
	end if
	simpson(nr)=1.d0/3.d0
c ===============================================================
c loop over all q values:
	do 100 iq=1,nq
	q(iq)=qmin+dfloat(iq-1)*(qmax-qmin)/dfloat(nq-1)
	qlq(iq)=0.0d0
c ===============================================================
c no the integral over all space r**2 dr jl(gr) QL(r)
	do 200 ir=1,nr
	x=q(iq)*r(ir)
	arg=(r(ir)**2)*jl(L,x)*qlr(ir)
	qlq(iq)=qlq(iq)+simpson(ir)*arg
200	continue
c ===============================================================
c write out every hundreth one to the screen.
	if(mod(iq,500).eq.1)write(*,81)iq,q(iq),qlq(iq)
81	format('  iq=',i5,'   q=',f9.4,' qlq(iq)=',e16.6)
100	continue
c end of loop over all q values.
c ===============================================================
	write(*,*)'  '
	write(*,*)'  write out results to qlq.dat '
	write(*,*)' '
	if(L.eq.0)write (*,*)' Actually qlq_0.dat'
	if(L.eq.1)write (*,*)' Actually qlq_1.dat'
	if(L.eq.2)write (*,*)' Actually qlq_2.dat'
	if(L.eq.3)write(*,*)' Actually qlq_3.dat'
	write(*,*)'  '
        write(*,*)' remember: q in 1/anstrom, QL(q) in eV*sqrt(A**3)'

	open(unit=30,file='qlq_0.dat',status='unknown')
	close(unit=30)
	open(unit=30,file='qlq_1.dat',status='unknown')
	close(unit=30)
	open(unit=30,file='qlq_2.dat',status='unknown')
	close(unit=30)
	if(L.eq.0)open(unit=30,file='qlq_0.dat',status='unknown')
	if(L.eq.1)open(unit=30,file='qlq_1.dat',status='unknown')
	if(L.eq.2)open(unit=30,file='qlq_2.dat',status='unknown')
	if(L.eq.3)open(unit=30,file='qlq_3.dat',status='unknown')
	do 40 iq=1,nq
	write(30,41)q(iq),qlq(iq)
41	format(2e18.9)
40	continue
	close(unit=30)
c ===============================================================
999	continue
	stop
	end
c ===============================================================
c
c
C SPHERICAL BESSEL FUNCTION
C
        FUNCTION JL(L,X)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        REAL*8 JL
        IF(L.EQ.0)GO TO 100
        IF(L.EQ.1)GO TO 101
        IF(L.EQ.2)GO TO 102
        IF(L.EQ.3)GO TO 103
        IF(L.EQ.4)GO TO 104
        JL=0.0D0
        RETURN
C
100     CONTINUE
C
C  L=0  J0=SINX/X
C
	XP=0.001
        IF(X.LT.XP)GO TO 1000
        JL=DSIN(X)/X
        RETURN
1000    X2=X*X
	X4=X2*X2
	JL=(1.D0-X2/6.D0+X4/120.D0)
        RETURN
101     CONTINUE
C
C L=1  J1 = SINX/X**2-COSX/X
C
	XP=0.1D0
        IF(X.LT.XP)GO TO 1001
        JL=(DSIN(X)/X-COS(X))/X
        RETURN
1001	X2=X*X
	X4=X2*X2
	X6=X2*X4
	X8=X2*X6
        JL=1.0D0-X2/10.0D0+X4/280.0D0-X6/15120.0D0
     1 +X8/1330560.0D0
	JL=(X/3.0D0)*JL
        RETURN
102     CONTINUE
C
C L=2  J2=(3/X**3-1/X)*SINX-3*COSX/X*X
C
	XP=0.5D0
        IF(X.LT.XP)GO TO 1002
        JL=( (3.D0/(X*X)-1.D0)*DSIN(X)-3.D0*DCOS(X)/X)/X
        RETURN
1002    X2=X*X
	X4=X2*X2
	X6=X2*X4
	X8=X2*X6
	X10=X2*X8
	JL=1.0D0-X2/14.0D0+X4/504.0D0-X6/33264.0D0
     1 +X8/3459456.0D0-X10/518918400.0D0
	JL=(X2/15.0D0)*JL
        RETURN
103	CONTINUE
C
C L=3  J3=(15/X**4-6/X**2)*SINX-(15/X**3-1/X)*COSX
C
	XP=0.5D0
        IF(X.LT.XP)GO TO 1003
	X2=X*X
	X3=X2*X
	X4=X2*X2
	JL=(15.0D0/X4-6.0D0/X2)*DSIN(X)-(15.0D0/X3-
     1 1.0D0/X)*DCOS(X)
	RETURN
1003	X2=X*X
	X3=X2*X
	X4=X2*X2
	X6=X2*X4
	X8=X2*X6
	JL=1.0D0-X2/18.0D0+X4/792.0D0-X6/61776.0D0
     1 +X8/7413120.0D0
	JL=(X3/105.0D0)*JL
	RETURN
104	CONTINUE
C
C L=4  J4=(105/X**5-45/X**3+1/X)*SINX-(105/X**4-10/X**2)*COSX
C
	XP=0.5D0
        IF(X.LT.XP)GO TO 1004
	X2=X*X
	X3=X2*X
	X4=X2*X2
	X5=X3*X2
	JL=(105.0D0/X5-45.0D0/X3+1.0D0/X)*DSIN(X)-
     1  (105.0D0/X4-10.0D0/X2)*DCOS(X)
	RETURN
1004	X2=X*X
	X4=X2*X2
	X6=X2*X4
	X8=X2*X6
	JL=1.0D0-X2/22.0D0+X4/1144.0D0-X6/102960.0D0
     1 +X8/14002560.0D0
	JL=(X4/945.0D0)*JL
	RETURN
        END
