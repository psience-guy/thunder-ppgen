C getlinear.f
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

	integer    nrmax, nr, np, nlmax,lmax
	parameter (nrmax=10001,nlmax=4)
	integer    i,j,k
	CHARACTER*30 foutl(0:nlmax),fpsil(0:nlmax),fvppl(0:nlmax) 
	dimension rLOG(nrmax),qlrLOG(nrmax),psiLOG(nrmax)
	dimension vl(nrmax)
	character*72	junk
	character*6	fouthd
	character*9	fvpphd,fpsihd
	character*20    flocal,psilocal
	parameter (fouthd='pseudo',fvpphd='vpplinear')
	parameter (fpsihd='psilinear')
	character*2	c
c ===============================================================
C constants
        PI=4.D0*DATAN(1.D0)
	EQ2=14.39975D0
	abohr=0.529177D0
	hartree=EQ2/ABOHR
c ===============================================================
c jjdong
	print *, ' This code is to transform the output of psgen into'
	print *, ' a linear grid format.'
	write(*,*)' Ready to open ncpp.ini'
	open (9, file='ncpp.ini', status='old')
!        read (9,910) a1,i1,i2,i3,a2
        read (9,*) a1,i1,i2,i3,a2
        write (*,910) a1,i1,i2,i3,a2
910	 format(f6.2,3(1x,i3),1x,f6.2)
	 do i=1, i1+i2
	  read (9,920) i3,i4,a2
	  write (*,920) i3,i4,a2
	 enddo
920	format(2(1x,i3),1x,f6.2)	
	read (9,930) lmax, c
930	format(i1,a2)
	print *, 'Lmax =', lmax
	print *, 'PP =', c
	close(9)
c Get lmax and l_loc from lmaxNlloc.dat
        write(*,*)' get lmax and L_loc from lmaxNlloc.dat'
	lmaxold=lmax
        open(unit=77,file='lmaxNlloc.dat',status='old')
        read(77,*)lmax
        read(77,*)l_loc
        close(unit=77)
	write(*,*)' lmax=',lmax,' L_loc=',L_loc
	if(lmax.ne.lmaxold)stop' lmax does not agree.'

	do i=0, lmax
C	 write(junk, '(a)') 'junk'
	 write(fvppl(i),'(a,i1.1,''.dat'')') fvpphd,i
	 write(fpsil(i),'(a,i1.1,''.dat'')') fpsihd,i
	 Write(foutl(i),'(a,i1.1,''.dat'')') fouthd,i
	 print *, i, ' ', foutl(i), fvppl(i), fpsil(i)
c	 pause
	enddo
cofs	 fvppl(lmax)='vpplinear_loc.dat'
	flocal='vpplinear_loc.dat'
	psilocal='psilinear_loc.dat'

c	print *, '	First: local part ', fvppl(lmax)
	print *, '	First: local part ', flocal
	print *,'  read from foutl(l_loc)=',foutl(l_loc)
        write(*,*)' Ready to open 21'
	open (21, file=foutl(l_loc), status='old')
	 read (21,'(a)') junk
c We start off with the first point r=0
	 rlog(1)=0.0d0
	 qlrLOG(1)=0.0d0
	 vl(1)=0.0d0
	 nr=1
	 do i=2,nrmax 
	  read (21,*, end=800) j,rlog(i),psilog(i),qlrlog(i)
	  nr=nr+1
	  vl(i)=qlrlog(i)
	  rlog(i)=rlog(i)*abohr
	  qlrlog(i)=qlrlog(i)*hartree
	  psilog(i)=psilog(i)/sqrt(abohr)
	 enddo
800	continue
	close(21, status='keep')
	  qlrlog(1)=2.0d0*qlrlog(2)-qlrlog(3)

	np=nr-1	
	write(*,*)'  '
	write(*,*)' read in nr points=',nr

c Now we interpolate to put on a linear grid
	rmax=rLOG(nr)
	write(*,*)' rmax=',rmax,' in Ang units.'
        write(*,*)' Ready to open 22'
c	open (22, file=fvppl(lmax), status='unknown')
	open (22, file=flocal, status='unknown')
	 do i=1,nrmax-1
	  x=0.9999999d0*rmax*(float(i-1))/float(nrmax-1-1)
	  norder=5
          call polyint(rLOG,qlrLOG,nrmax,nr,norder,x,y,dy)
	  write(22,*)x,y
	 enddo
	close(22, status='keep')
        write(*,*)' Ready to open 23'
c	open (23, file=fpsil(lmax), status='unknown')
	open (23, file=psilocal, status='unknown')

	 do i=1,nrmax-1
	  x=0.9999999d0*rmax*(float(i-1))/float(nrmax-1-1)
	  norder=5
          call polyint(rLOG,psiLOG,nrmax,nr,norder,x,y,dy)
	  write(23,*)x,y
	 enddo
	nr=nrmax-1
	close(23, status='keep')

	write(*,*)' '
	print *, '	Second: non-local part '
	write(*,*)'  '
!
!
!
cOFS	do k=0, lmax-1
c OFS Note that we include the L value for the
c local potential also.
c The result for that L value will be ZERO. Thus we pad the data.
	do k=0, lmax
!
!
	write(*,*)' '
	write(*,*)' L value=',k
	if(k.eq.L_loc)write(*,*)' This is local. We set vnl=0'
	print *, '		',fvppl(k)
        write(*,*)' Ready to open New 21'

	open (21, file=foutl(k), status='old')

	 read (21,'(a)') junk
c We start off with the first point r=0
	 rlog(1)=0.0d0
	 qlrLOG(1)=0.0d0
	 nr=1
	 do i=2, nrmax 
	  read (21,*, end=810) j,rlog(i),psilog(i),qlrlog(i)
	  nr=nr+1
!
!         The non-local potential now will be
!         v(nonlocal) = v(individual) - v(local)
!
!         Later: v(individual) = v(nonlocal) + v(local)
!
!         Since v(individual) and v(local) should have
!         the same long range behavior, v(nonlocal)
!         should have short range character.
!
!
	  qlrlog(i)=qlrlog(i)-vl(i)
!
!
	  rlog(i)=rlog(i)*abohr
	  qlrlog(i)=qlrlog(i)*hartree
	  psilog(i)=psilog(i)/sqrt(abohr)
	 enddo
810	continue
	close(21, status='keep')
	qlrlog(1)=2.0d0*qlrlog(2)-qlrlog(3)

	np=nr-1
	write(*,*)'  '
	write(*,*)' read in nr points=',nr

c Now we interpolate to put on a linear grid
	rmax=rLOG(nr)
	write(*,*)' rmax=',rmax,' in Ang units.'
        write(*,*)' Ready to open 12'

	open (12, file=fvppl(k), status='unknown')

	 do i=1,nrmax-1
	  x=0.9999999d0*rmax*(float(i-1))/float(nrmax-1-1)
	  norder=5
          call polyint(rLOG,qlrLOG,nrmax,nr,norder,x,y,dy)
	  write(12,*)x,y
	 enddo
	close(12, status='keep')

	open (13, file=fpsil(k), status='unknown')

	 do i=1,nrmax-1
	  x=0.9999999d0*rmax*(float(i-1))/float(nrmax-1-1)
	  norder=5
          call polyint(rLOG,psiLOG,nrmax,nr,norder,x,y,dy)
	  write(13,*)x,y
	 enddo
	nr=nrmax-1
	close(13, status='keep')
	enddo
	write(*,*)'  '
	write(*,*)' Done with linear. Bye now.'
	write(*,*)'  '

	end
