c polyint.f
c
* Official Cactus-98 Package Feb. 26, 1998
* Subsection: begin.
c
c                        Professor Otto F. Sankey
c                        Department of Physics and Astronomy
c                        Box 871504
c                        Arizona State University
c                        Tempe, AZ 85287-1504
c                        (602) 965-4334
c                        (602) 965-7954 FAX
c                        
c otto.sankey@asu.edu
c 
c ==================================================================
c polyint.f  A subroutine for interoplating (and extrapolating) a set
c data using a polynominal. From Numerical recipes, p. 82. Note:
c In Num. Rec. they call this program polint.
c
c In this subroutine we put in the entire array of 
c xarray(npoints),yarray(npoints). 
c We also put in a value of x and the output
c is the interpolated (or extrapolated) value of y and an error estimate dy.
c The polynomial interpolation is to a n-order polynomial.
c In the original polint.f (of numerical recipes), you did not put in the
c entire array, but only the array near the point x. This is a pain in the
c rear. We only do this thing once are twice, so who cares about an extra milisecond
c of computation time
c
c The program first finds the value of k, where k is defined by the fact that
c x lies between xarray(k) and xarray(k+1). If x<xarray(1) then we use the
c the first npoints of the array for our polynomial. If x>xarray(npoints), then
c we use the last n points of the array. Also if k is near the ends of the
c array, then in order to have n points of information, we use the first (or
c last) n points of xarray.
c ========================================================================
c test.f Remove the c's to test the program
c	implicit double precision (a-h,o-z)
c	dimension xxa(100),yya(100),xa(5),ya(5)
c We create an array of x and y. 
c We have two examples. A linear x, and a nonlinear x.
c In both cases the y function is a gaussian.
c	do i=1,100
c
c
c The linear x array.
c	xxa(i)=i
c end of linear x array.
c The non-linear x array.
c 	ff=exp( float(i)/100.)-exp(1./100.)
c 	ff=ff/(exp(1.)-1.)
c 	xxa(i)=1.+(100.-1.)*ff
c End of non-linear x array
c
c
c 	yya(i)=dexp( -(xxa(i)-50.0d0)**2/(25.0d0**2) )
c 	end do
c We will use a 5'th order polynomial.
c 	n=5
c 	ndim=100
c 	npoints=100
c 	open(unit=30,file='junk1.dat',status='unknown')
c 	open(unit=31,file='junk2.dat',status='unknown')
c Now lets check the interpolator (extrapolator) by computing from
c -50 to 150 not at integers values of x but at an integer + 1/sqrt(2).
c
c 	do i=-50,155
c 	x=dfloat(i)+1.d0/dsqrt(2.0d0)
c
c 	call polyint(xxa,yya,ndim,npoints,n,x,y,dy)
c
c 	write(*,*)'  '
c 	if(i.ge.1.and.i.le.100)write(*,1)i,xxa(i),yya(i)
c 	write(*,2)i,x,y,dy
c 1	format(' i=',i3,' true x,y=',2e16.6)
c 2	format(' i=',i3,' fit  x,y=',2e16.6,' err=',e16.6)
c 	if(i.ge.1.and.i.le.100)write(30,*)xxa(i),yya(i)
c 	write(31,*)x,y
c 	end do
c 	close(unit=30)
c 	close(unit=31)
c 	stop
c 	end
c ========================================================================
c polyint.f  A subroutine for interoplating (and extrapolating) a set
c data using a polynominal. From Numerical recipes, p. 82.
c Otto F. Sankey, Dept. of Physics, Araizona State University, tempe AZ 85287
c 602 965-4334, otto.sankey@asu.edu
c Dec. 10, 1991.
c
	subroutine polyint(xarray,yarray,ndim,npoints,n,x,y,dy)
c Input:
c xarray(npoints),yarray(npoints) array of x,y values. 
c         Note: The xarray(i) do not have to be on a regular grid, 
c         but xarray(i+1) must always be greater than xarray(i). 
c ndim -- dimension of xarray,yarray exactly as defined in calling program.
c npoints -- number of points in xarray,yarray. Of course npoints.le.ndim.
c n -- the order of the polynomial we use to interpolate y at the value of x.
c x -- The input value of x in which we wish the interpolated value of y.
c
c Output:
c y -- The interpolated value of y at the input point x.
c dy -- An error estimate for y.
c
	parameter (nmax=10)
	implicit double precision (a-h,o-z)
	dimension xarray(ndim),yarray(ndim)
	dimension xa(nmax),ya(nmax),c(nmax),d(nmax)
c ========================================================================
	if(npoints.gt.ndim)stop' too many points in polyint!'
	if(n.gt.nmax)stop' nmax not big enough in polyint! sorry!'
c ===================================================================
c Determine the value of k. 
c The input point x lies between xarray(k) and xarray(k+1).
	k=0
	if(x.lt.xarray(1))then
	k=1
	write(*,*)' caution k=1 in polyint'
	end if
	if(x.gt.xarray(npoints))then
	k=npoints
	write(*,*)' caution, k=npoints in polyint'
	end if
	if(k.ne.0)go to 900
c We find the point k as x is between xarray(k) and xarray(k+1).
	do 800 i=1,npoints
	if(xarray(i).gt.x)go to 810
	k=i
800	continue
810	continue
c
900	continue
c ===================================================================
c Now fill up the xa and ya arrays which go up to n.
c
	nover2=n/2
	if(mod(n,2).eq.0)nover2=nover2-1
c Suppose x deep in the array.
	if(k-nover2.lt.1)go to 100
	if(k-nover2+n-1.gt.npoints)go to 200
c For x somewhere deep in xarray, 
c Then for n=odd,
c xa(1)=xarray(k-n/2), xa(2)=xarray(k+1-n/2), ... xa(n)=xarray(k+n-n/2)
c Then for n=even,
c xa(1)=xarray(k+1-n/2), xa(2)=xarray(k+2-n/2), ... xa(n)=xarray(k+1+n-n/2)
	do 23 i=1,n
	xa(i)=xarray(k-nover2+(i-1))
	ya(i)=yarray(k-nover2+(i-1))
23	continue
	go to 27
c k-nover2.lt.1
100	continue
	ishift=1-(k-nover2)
	do 123 i=1,n
	xa(i)=xarray(k-nover2+(i-1)+ishift)
	ya(i)=yarray(k-nover2+(i-1)+ishift)
123	continue
	go to 27

c k-nover2.lt.1
200	continue
	ishift=-((k-nover2+n-1)-npoints)
	do 223 i=1,n
	xa(i)=xarray(k-nover2+(i-1)+ishift)
	ya(i)=yarray(k-nover2+(i-1)+ishift)
223	continue
c
27	continue
c ===================================================================
c The rest of the program is right out of numerical recipes, p. 82. For comments
c on the program on this part, see numerical recipes.
	ns=1
	dif=dabs(x-xa(1))
	do 11 i=1,n
	dift=dabs(x-xa(i))
	if(dift.lt.dif)then
		ns=i
		dif=dift
	end if
	c(i)=ya(i)
	d(i)=ya(i)
11	continue
	y=ya(ns)
	ns=ns-1
	do 13 m=1,n-1
	do 12 i=1,n-m
	ho=xa(i)-x
	hp=xa(i+m)-x
	w=c(i+1)-d(i)
	den=ho-hp
	if(den.eq.0.d0)stop' two xa"s are equal in polyint!!'
	den=w/den
	d(i)=hp*den
	c(i)=ho*den
12	continue
	if(2*ns.lt.n-m)then
	dy=c(ns+1)
	else
	dy=d(ns)
	ns=ns-1
	end if
	y=y+dy
13	continue
	return
	end
c ==================== Last line of subs.f ===========================
