c mk_pseudoFile.f. Make the final pseudopotential file.
c Version June 5, 2000.
c =====================================================
c OFS June 5, 2000. A major change was made in that I
c altered the program so that it wrote out the shells
c of the pseudopotential skipping l_loc. Here l_loc is
c the local potential. Before, this program assumed that l_loc
c was the highest Lvalue. This causes grief for transition
c metals.
c =====================================================
c
C This code is to save Pseuodo Potential data into one date file.
C Jianjun Dong, 06-26-98, 07-14-98
C Some part of code is adapted from Otto Sankey's vlocal2erf.f
C Before running this code, you must have following input data files:
C 1) vpplinear_loc.dat
C 2) vpplinear0.dat, vpplinear1.dat, ... etc.
C 3) qlcl.dat (containing cl(l) )
C 4) vkblinear0.dat, vknlinear1.dat, ... etc.
C and you have to check the cutoff radius of Vpp and Vkb
C To compile: f77 -o wtvpp.x wtvpp.f

	integer		nz, iheader, lmax, nlmax, nrmax, i,j,k
	parameter      (nlmax=4, iheader=14, nrmax=10999)
	dimension	lsshPP(nlmax)
	character*2	cotto
	character*70	cheader(iheader)
	character*30	fout, fvpp, fvkb
c	real		rc_vpp(0:nlmax), rc_vkb(0:nlmax),cl(0:nlmax)
	real		cl(0:nlmax)
	real		r(nrmax), v(nrmax)
	real		vshort(nrmax), verf(nrmax)
	real		ave, nearzero
	parameter      (nearzero=1.0e-6)
	integer		ncsh,nvsh,ntypeexc
        real            exmix,fnumv,fnumc
	real		znuc,rc
        real            rc_max
	integer Lshell(nlmax)
c ============================================
	write(*,*)'  '
	write(*,*)' Welcome to mk_pseudoFile.x'
	write(*,*)' This writes out the final XXX.pp file.'
	write(*,*)' E.g. 014.pp etc.'
	write(*,*)'  '
c First read lmax AND l_loc AND the Lvalues of the NL potential
c from lmaxNlloc.dat
        write(*,*)' Read lmax,L_loc, and L-values from lmaxNlloc.dat'
        open(unit=77,file='lmaxNlloc.dat',status='old')
        read(77,*) lmax
        read(77,*) L_loc
c We now read the Lvalues that we actually have a non-local pseudo
c potential for.
c We assume that the L values go in order, S,P,D,..
c We then figure out the allowed l values by eliminating l_loc.
c For example, if lmax=2 (SPD) and l_loc=2, then L=0,1.
c For example, if lmax=2 (SPD) and l_loc=1, then L=0,2.
c For example, if lmax=2 (SPD) and l_loc=0, then L=1,2.
	nsshPP=(lmax+1)-1
	if(lmax.gt.nlmax)stop' lmax too big.. Redimension nlmax'
        do 1492  issh=1,nsshPP
        read(77,*)lsshPP(issh)
1492	continue
        close(unit=77)
c The number of pseudopotential shells is nsshPP, and
c their L values are lsshPP(1:nsshPP)


C Reading
	open (30, file='ncpp.ini', status='old')
c
c JUERGEN begin
c
c       read (30,'(f6.2,3(1x,i3),1x,f6.2,f6.2)')znuc,ncsh,nvsh,ntypeexc,rc
        read (30,*) znuc,ncsh,nvsh,ntypeexc,rc
        exmix = 0
c
c JUERGEN end
c
c JUERGEN begin
c
c       charges fnumc, fnumv have to be determined
c
        fnumc = 0.0
        nz=int(znuc+0.5)
        fnumc = 0.0
        do i=1,ncsh
!          read (30,'(2(1x,i3),1x,f6.2)')j,k,rc
          read (30,*)j,k,rc
          fnumc = fnumc + rc 
        enddo
        fnumv = 0.0
        do i=ncsh+1,ncsh+nvsh
!          read (30,'(2(1x,i3),1x,f6.2)')j,k,rc
          read (30,*) j,k,rc
          write (*,*) j,k,rc
          fnumv = fnumv + rc 
        enddo
        read (30,'(i1,a2)') lmax,cotto
	close(30, status='keep')
	if (lmax.le.0.or.lmax.gt.nlmax) then
	 print *, ' and should be within 0,', nlmax
	 stop 'error'
	endif
c
c JUERGEN end
c
C header: information about this PseudoPotential'
	 open (30, file='ncpp.dat', status='old')
	 do i=1, 70
	  cheader(1)(i:i)='!'
	 enddo
c
c
         do i=2, iheader-1
	  read (30,'(a)') cheader(i)
	 enddo
	 do i=1, 70
           cheader(iheader)(i:i)='!'
	 enddo
	 close(30, status='keep')

C Writing
	write(*,*)' The OUTPUT file:'
	write(*,*)' It will have a name like 014.pp (Si), 003.pp (Li)'
	write (fout, '(i3.3,''.pp'')') nz
	print *, ' The output file is:', fout
c
	open (40, file=fout, status='unknown') 
        do i=1, iheader
          write(40,'(a)') cheader(i)
	enddo
c
c
c New fixup. Here is what used to be there.
cofs	write(40, *) lmax
c ***********ERROR This was fixed June 5, 2000.
COFSc I want to put the "shell" structure in. lmax means L=0,1,... lmax-1.
COFSc So lmax=2 means L=0,1.
COFS	nshells=Lmax
COFS	do i=1,nshells
COFSc So lmax=2 means L=0,1.
COFS	Lshell(i)=i-1
COFS	end do
COFSc Now write out nshells
COFSc               L(1),L(2),...,L(nshells)
	nshells=nsshPP
	do i=1,nshells
	Lshell(i)=lsshPP(i)
	end do

! JPL/RBE 1999 - Write out iexc which is the exchange-correlation approximation
! that was used to generate this pseudopotential data.
c
c JUERGEN begin
c
        IF (ntypeexc .NE. 12) THEN
          write (40, 906) ntypeexc,0.0
        ELSE
          write (40,906) ntypeexc,exmix
        END IF
c
905       format(i5,'                                   ! iexc')
906     format(i5,1x,f6.2,25x,'   ! iexc , exmix')
c
c
c JURGEN end
c
        
	write(40,900)nshells
        write (*,*) 'nshells = ', nshells
	if(nshells.eq.0)then
	write(*,*)' There is no pseudopotential. nshells=',nshells
	stop' There is no pseudopotential.'
	end if
	if(nshells.eq.1)write(40,901)(Lshell(i),i=1,nshells)
	if(nshells.eq.2)write(40,902)(Lshell(i),i=1,nshells)
	if(nshells.eq.3)write(40,903)(Lshell(i),i=1,nshells)
	if(nshells.eq.4)write(40,904)(Lshell(i),i=1,nshells)
900	format(i5,'                                   ! nshells')
901	format(i5,'                                   ! L values')
902	format(2i5,'                              ! L values')
903	format(3i5,'                         ! L values')
904	format(4i5,'                    ! L values')

200	print *, ' reading V_local (r):'
c
c ==========================================================
c June 5, 2000 OFS This program has a funny problem.
c For some reason it asks about nz=1 or 2. What the heck for.
c	if (nz.eq.1.or.nz.eq.2)then
c	write(*,*)'  June 5, 2000 OFS'
c	write(*,*)' This program has a funny problem.'
c	write(*,*)' For some reason it asks about nz=1 '
c	write(*,*)' or 2. What the heck for.'
c	write(*,*)' We need to fix that.'
c	stop' June 5, 2000 OFS Fix mk_pseudoFile.f '
c	end if
c =========================================================
        if (nz .ne. 1 .and. nz .ne. 2) then
	 open(unit=30,file='vpplinear_loc.dat',status='old')
        end if
	ic=0

c e2=e^2 eV Angstrom
	e2=14.39975

	do i=1,nrmax
         if (nz .ne. 1 .and. nz .ne. 2) then
	  read(30,*,end=210)r(i),v(i)
         else
          r(i) = 0.005d0*float(i - 1)
          if (r(i) .le. 1.0d-4) then
           v(i) = - float(nz)*e2/1.0d-04
          else 
           v(i) = - float(nz)*e2/r(i)
          end if
         end if

c        rv=-r*v/e2. For large r this should be z.
c        rv=-r(i)*v(i)/e2
c
	 ic=ic+1
	 if(ic.gt.nrmax)stop' ic.gt.nrmax redimenison'
c
	end do
c
210	continue
	close(unit=30)
c ic is how many points we read in.
	write(*,*)' ic=',ic
c ====================================
c We now fit V_local to -Z * e^2 * erf( alpha * r)/r.
c For large r, this becomes -Z e^2/r.
c For small r, this become -Z *e^2 2/sqrt(pi) *alpha
c We choose alpha so that when V*r/e^2 is safely constant,
c then that means erf is safely 1. For erf to be safely 1,
c means x in erf(x)=about=4. 
c So let r_safe be the value where v(r_safe) is very very close
c to -Ze^2/r.
c So we have alpha*r_safe=4. Then alpha=4/r_safe.
c 
c We now look throught the data and find rsafe.
c
c       Juergen use a somewhat lower tolerance.
c       For exact exchange I observed energies
c       which were about 0.05 eV to small. 
c
c       tolerance=1.0e-05
c
        tolerance=1.0e-07
c
	npoint=ic-11
c
	isafe=-1
c
	do i=1,npoint

c         Add 10 points ahead of i. 
c         If the average of these 10 points ahead of i is very close to i,
c         then we have found i_safe corresponding to r_safe.
c
          avg=0.
c
c JUERGEN begin 
c 
c         We do not get reasonable values for z when dealing
c         with exact exchange potentials. We therefore
c         determine z from the data of ncpp.ini
c
          do iadd=1,10	
            avg=(-r(i+iadd)*v(i+iadd)/e2)+avg
          end do
          avg=avg/10.0
c
          point=-v(i)*r(i)/e2
          if(abs(point-avg).lt.tolerance)then
            isafe=i
            Z=point
            go to 220
          endif
        enddo
c
c
c JUERGEN end
c 
c
c
220	continue
c
	write(*,*)' isafe=',isafe
	write(*,*)' A very close approximation to Z is',Z
c
c JUERGEN begin
c
!       Z = znuc - fnumc
!       write(*,*) ' '
!       write(*,*) ' Z is reset to   ',Z
!       write(*,*) ' '
c
c JUERGEN end
c
	if(isafe.le.0)then
	write(*,*)' No isafe found *** error ****'
	end if
	if(isafe.le.0)stop' No isafe found ***errror****'
	rsafe=r(isafe)
	write(*,*)' rsafe=',rsafe
c Saftey is erf(4). 4=alpha*rsafe, or alpha=4./rsafe
	alpha=4./rsafe
	write(*,*)' alpha=',alpha
c Now calculate the shortrange potential.
c V_local = -Z e^2 * erf(alpha r)/r + V_local_shortrange
	pi=3.1415926

! JPL 1999 reevaluate v(1) from 2.0d0*v(2) - v(3)
        v(1) = 2.0d0*v(2) - v(3)
	Vshort(1)=v(1)+Z*e2*(2./sqrt(pi))*alpha
	verf(1)=-Z*e2*(2./sqrt(pi))*alpha
        if (nz .ne. 1 .and. nz .ne. 2) then
	 do i = 2, ic
	  arg=alpha*r(i)
	  Vshort(i)=v(i)-(-Z*e2*OFSerf(arg)/r(i))
	  Verf(i)=(-Z*e2*OFSerf(arg)/r(i))
   	 end do
        else
         v(2) = 2.0d0*v(3) - v(4)
         v(3) = 2.0d0*v(4) - v(5)
         v(4) = 2.0d0*v(5) - v(6)
         v(5) = 2.0d0*v(6) - v(7)
         Vshort(2)=v(2)+Z*e2*(2./sqrt(pi))*alpha
         Vshort(3)=v(3)+Z*e2*(2./sqrt(pi))*alpha
         Vshort(4)=v(4)+Z*e2*(2./sqrt(pi))*alpha
         Vshort(5)=v(5)+Z*e2*(2./sqrt(pi))*alpha
        end if

c Jel  
c find rc_max
        rc_max = 0.0d0
        do nssh = 1,nshells
         i = Lshell(nssh)
         write (fvpp,'(a,i1.1,''.dat'')') 'vpplinear',i
         ic=0
         if (nz .ne. 1 .and. nz .ne. 2) then
          open (30, file=fvpp, status='old')
          do j = 1, nrmax
           read (30, *, end=7777) r(j), v(j)
          enddo
         endif
7777     do j=1, nrmax-9, 10
          ave=0.0
          do k=0,9
           ave=ave+v(j+k)*v(j+k)
          enddo
          ave=sqrt(ave/10.0)
          if (ave.gt.nearzero) ic=j+9
         enddo
         if (ic .eq. 0) then
          ic = 5
          r(2) = r(1) + 0.005
          r(3) = r(2) + 0.005
          r(4) = r(3) + 0.005
          r(5) = r(4) + 0.005
         end if

         if ( r(ic) .gt. rc_max) rc_max = r(ic)
        enddo

        write (40,910) (znuc - fnumc)    
	write (40,911) alpha
        write (40,912) rc_max+0.01

910     format(f14.5,'                        ! Z val')
911     format(f14.8,'                        ! alpha')
912     format(f14.2,'                        ! Rcutoff')


        if (nz .eq. 1 .or. nz .eq. 2) isafe = isafe + 3
        write(40,*)isafe

c Vshort is the short-range potential. This we need.
c We need it for the neutral atom potential (We don't need it for the
c pseudopotential).
	do i=1,isafe
	write(40,*)r(i),Vshort(i)
	end do

300	print *, ' Reading and saving V_NL(l) (SR)'
c OFS June 5, 2000 Error ***fixup
cOFS	do i=0, lmax-1
	do 1110 nssh=1,nshells
	i=Lshell(nssh)
	 write (fvpp,'(a,i1.1,''.dat'')') 'vpplinear',i
	 do j=1, nrmax
	  v(j)=0.0
	 enddo
	 ic=0
         if (nz .ne. 1 .and. nz .ne. 2) then
	  open (30, file=fvpp, status='old')
	  do j=1, nrmax
	   read (30, *, end=301) r(j), v(j)
	  enddo
         end if
301	 do j=1, nrmax-9, 10
	  ave=0.0
	  do k=0,9
	   ave=ave+v(j+k)*v(j+k)
	  enddo
	   ave=sqrt(ave/10.0)
	  if (ave.gt.nearzero) ic=j+9   
	 enddo
         if (ic .eq. 0) then
          ic = 5
          r(2) = r(1) + 0.005
          r(3) = r(2) + 0.005
          r(4) = r(3) + 0.005
          r(5) = r(4) + 0.005
         end if
310	 write (40, 391) i, ic
	 write (*, 391) i, ic
391	 format(' L=',i1, ' Npoint=',i5)
	 close(30) 
	 do j=1, ic
	  write (40,*) r(j), v(j)
	 enddo
COFS	enddo
1110	continue

400	print *, ' Reading V^(KB) information'
        if (nz .ne. 1 .and. nz .ne. 2) then
	 open (30, file='qlcl.dat', status='old') 
c OFS June 5, 2000 Error ***fixup
COFS	 do i=0, lmax-1
c We read all(!!) cl(L) values, even the local 
c value (where cl=0). This is how qlqgrid.f writes them.
c Note that it is lmax, not lmax-1!! OFS June 5, 2000
	do i=0,lmax
	  read (30, 490) cl(i)
	enddo
490	format(f14.7)
	 close(30)
        end if


410	print *, ' Reading and saving V^KB(l) (SR)'
c OFS June 5, 2000 Error ***fixup
COFS	do i=0, lmax-1
        do 1120 nssh=1,nshells
        i=Lshell(nssh)
	 write (fvkb,'(a,i1.1,''.dat'')') 'vkblinear',i
	 do j=1, nrmax
	  v(j)=0.0
	 enddo
	 ic=0
         if (nz .ne. 1 .and. nz .ne. 2) then
	  open (30, file=fvkb, status='old')
	  do j=1, nrmax
	   read (30, *, end=411) r(j), v(j)
	  enddo
411	  do j=1, nrmax-9, 10
	   ave=0.0
	   do k=0,9
	    ave=ave+v(j+k)*v(j+k)
	   enddo
	    ave=sqrt(ave/10.)
	   if (ave.gt.nearzero) ic=j+9   
	  enddo
         else
          ic = 5
          r(2) = r(1) + 0.005
          r(3) = r(2) + 0.005
          r(4) = r(3) + 0.005
          r(5) = r(4) + 0.005
         end if 
420	 write (40, 491) i, ic, cl(i)
	 write (*, 491) i, ic, cl(i)
491	 format(' L=',i1, ' Npoint=',i5,' cl=', f14.7)
c Let me now check to make sure we have not gotten confused and written
c a cl=0.
	if(abs(cl(i)).lt.0.000001)then
	write(*,*)' cl=',cl(i),' i=',i
	write(*,*)' I"m afraid you wrote out the wrong one.'
	write(*,*)' Please check.'
	stop' please check cl in mk_pseudoFile.f'
	end if
	 close(30) 

! JPL 1999 Reevaluate v(1) from 2.0*v(2) - v(3)
         v(1) = 2.0d0*v(2) - v(3)
	 do j=1, ic
	  write (40,*) r(j), v(j)
	 enddo
COFS	enddo
1120	continue
	
	close(40)
	write(*,*)'  '
	write(*,*)' *******************************************'
	write(*,*)'  '
	write(*,*)' Done with mk_pseudoFile.x'
        write(*,*)' The OUTPUT file:'
        print *, ' The output file is:', fout
	write(*,*)'  '
	write(*,*)' *******************************************'

	end
c ofserf.f  O Sankey's erf for machines which don't
c have an erf.
c ===============================================
c
	FUNCTION OFSERF(X)
C SEE ABRAMOWITZ AND STEGUN.
c	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        Y=X
        IF(X.LT.0.)Y=-X
        IF(Y.LT.1.0E-7)GO TO 1
        P=0.3275911D0
c	P=0.3275911
        T=1.0/(1.0+P*Y)
        A1=0.254829592
        A2=-0.284496736
        A3=1.421413741
        A4=-1.453152027
        A5=1.061405429
        Z=1.0-(A1*T+A2*T*T+A3*(T**3)+A4*(T**4)+A5*(T**5))*EXP(-Y*Y)
        IF(X.LT.0.)Z=-Z
        GO TO 2
1       Z=0.0
2       CONTINUE
        OFSERF=Z
        RETURN
        END
c
	function OFSerfc(x)
	OFSERFC=1.0-OFSERF(X)
	RETURN
	END
