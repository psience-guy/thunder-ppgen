! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Motorola, Physical Sciences Research Labs - Alex Demkov
! University of Regensburg - Juergen Fritsch
! Universidad de Madrid - Jose Ortega

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola - Jun Wang
! Ohio State University - Dave Drabold

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in
! subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! initial.f90
! Program Description
! ===========================================================================
!       This program pepares the ncpp.ini file for running the pseudopotential 
! program. This automates the procedure considerably. Much of the necessary 
! data is in elements.input, which this program reads.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658 
! Provo, UT 841184602-4658
! FAX 801-422-2265
! Office telephone 801-422-7444
! ===========================================================================
! 
! Program Declaration
! ===========================================================================
        program pickanatom
        use precision
        implicit none

! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: max_z = 110
        integer, parameter :: nsh_max = 18

! Local Variable Declaration and Description
! ===========================================================================
        integer ielement
        integer ioption
        integer ipoint
        integer issh
        integer ncsh
        integer numlmax
        integer ntypepsp
        integer nvsh
        integer nznuc

        integer, dimension (max_z, nsh_max) :: lq
        integer, dimension (max_z) :: ncoresh
        integer, dimension (max_z, nsh_max) :: nq
        integer, dimension (max_z) :: nssh
        integer, dimension (max_z) :: ntexc
        integer, dimension (max_z) :: ntpsp
        integer, dimension (max_z) :: numL
        integer, dimension (max_z) :: nvalesh
        integer, dimension (max_z) :: nz

        real(kind=long) rnlc
        real(kind=long) bmix
        real(kind=long) exmix

        real(kind=long), dimension (max_z) :: rnlc_in
        real(kind=long), dimension (nsh_max) :: xocc
        real(kind=long), dimension (max_z, nsh_max) :: xocc_in

        character(len=1) answer
        character(len=50) bar
        character(len=1) copyright
        character(len=2) pptype

        character(len=2), dimension (max_z) :: element
        character(len=1), dimension (0:9) :: z

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! First we initialize everything.
        element = 'XX'

        do ielement = 1, max_z
         nz(ielement) = ielement
        end do 

! Initialize some symbols
        copyright = char(169)
        do ipoint = 0, 9
         z(ipoint) = char(48+ipoint)
        end do

! Read elements.input file.
! *****************************************************************************
        open (unit=50, file = '../ppprogs/element.input', status = 'old')
        write (*,*) ' Begin reading main data in element.input '
        write (*,*) '  '
        write (*,*) '  '
         
! 7 comment lines
        do ielement = 1, 7
         read(50,*)
        end do

! Now begin real(kind=long) data
        do ielement = 1, max_z
         read (50,*)
         read (50,*) nz(ielement), element(ielement)
         if (ielement .ne. nz(ielement)) stop ' element.input mixup.'
         read (50,*) ncoresh(ielement), nvalesh(ielement), numL(ielement)
         if (ncoresh(ielement) .lt. 0) then
          write (*,*) ' ielement, ncoresh = ', ielement, ncoresh(ielement)
          write (*,*) ' Number of core shell cannot be negative! '
          stop
         end if
         if (nvalesh(ielement) .le. 0) then
          write (*,*) ' ielement, nvalesh = ', ielement, nvalesh(ielement)
          write (*,*) ' Number of valence shell must be positive! '
          stop
         end if

         nssh(ielement) = ncoresh(ielement) + nvalesh(ielement)
         if (nssh(ielement) .gt. nsh_max) then
          write (*,*) ' ielement, nssh = ', ielement, nssh(ielement)
          write (*,*) ' nssh is too large!'
          stop
         end if

         if (numL(ielement) .le. 0 .and. element(ielement) .ne. 'XX') then
          write (*,*) ' ielement, numL = ', ielement, numL(ielement)
          write (*,*) ' Number of non-local orbitals must be positive! '
          stop
         end if

         do issh = 1, nssh(ielement)
          read (50,*) nq(ielement,issh), lq(ielement,issh),                  &
     &                xocc_in(ielement,issh)
         end do
         read (50,*) ntexc(ielement), ntpsp(ielement), rnlc_in(ielement)

        end do
        close (unit = 50)

! *****************************************************************************
! Set up the bar (horizontal line for optics of output file)
        do ipoint = 1, 60
         bar(ipoint:ipoint) = ' = '
        end do
        bar(25:36) = 'Fireball2000'

        write (*,100) bar
        write (*,*) '          *------------------------------------* '
        write (*,*) '          |     THIS CODE IS PROPRIETORY,      | '
        write (*,*) '          |      SEE COPYRIGHT INFORMATION!    | '
        write (*,*) '          *------------------------------------* '
        write (*,*) '  '
        write (*,*) '             Usable only with permission from '
        write (*,*) '           the Fireball2000 executive committee. '
        write (*,*) '        This program is NOT, under ANY circumstances,'
        write (*,*) '          to be transfered to an unauthorized user! '
        write (*,*) '  '
        write (*,100) bar

        write (*,*) '  '
        write (*,*) '    *-------------------------------------------------* '
        write (*,*) '    |                                                 | '
        write (*,*) '    |                    Welcome.                     | '
        write (*,*) '    |       The Fireball2000 team welcomes you.       | '
        write (*,*) '    |                                                 | '
        write (*,*) '    *-------------------------------------------------* '
        write (*,*) ' '
        write (*,*) ' This program helps you to initialize the '
        write (*,*) ' "pseudopotential" part of the Fireball2000',copyright
        write (*,*) ' package. '
        write (*,*) '  '
        write (*,*) ' You will be asked 4 questions. Hit return for defaults. '
        write (*,*) '  '
        write (*,*) ' Ready. '
        write (*,*) '           Set. '
        write (*,*) '                    Go. '
        write (*,*) '  '
        write (*,*) '  '

! Open output file.
        open (12, file = 'ncpp.ini', status = 'unknown')

! *****************************************************************************
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' =================== Question No. 1 =================== '
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' Which atom would you like to generate pseudopotential ?'
        write (*,*) ' Your choice:'
        write (*,*) '  '

! Write out the periodic table.
        call periodensystem

        write (*,*) ' ===> Input desired atomic number:'
        read (*,*) nznuc
        if (nznuc .lt. 0 .or. nznuc .gt. max_z) then
         write (*,*) ' This atom is not contained in the elements.input file. '	
         stop
        end if
        write (*,*) ' Your chosen element: ', element(nznuc)

        if (element(nznuc) .eq. 'XX' .or. element(nznuc) .eq.'xx') then
         write (*,*) ' Sorry, this element is not yet in element.input '
         write (*,*) ' Go add an entry to element.input.  It will take you '
         write (*,*) ' a few minutes.  You must put in some default values. '
         write (*,*) ' Please check other entries as examples. '
         write (*,*) '  '
         stop
        end if

! *****************************************************************************
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' =================== Question No. 2 =================== '
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' We suggest the following partition of core and valence '
        write (*,*) ' electron shells (This can of course change depending on '
        write (*,*) ' your choice): '
        write (*,*) '  '
        write (*,200) ncoresh(nznuc), nvalesh(nznuc)
        write (*,*) '  '

        write (*,*) '  '
        write (*,*) ' Use the default? Y/N '
        read (*,*) answer 

        if (answer .eq. 'Y' .or. answer .eq. 'y' .or. answer .eq. ' ') then
         write (*,*) '  '
         write (*,*) ' OK: We USE the default. '
         ncsh = ncoresh(nznuc)
         nvsh = nvalesh(nznuc)
        else 
         write (*,*) '  '
         write (*,*) ' OK: We change the defaults. '
         write (*,*) '  '
         write (*,*) ' Insert your own choices of ncsh and nvsh '
         read (*,*) ncsh, nvsh 
        end if

        if (ncsh .lt. 0) then
         write (*,*) ' ncsh = ', ncsh 
         write (*,*) ' Number of core shell cannot be negative!'
         stop
        end if

        if (nvsh .le. 0) then
         write (*,*) ' nvsh = ', nvsh 
         write (*,*) 'Number of valence shell must be positive!'
         stop
        end if
        if ((ncsh + nvsh) .gt. nsh_max) then
         write (*,*) ' nsh_max = ', nsh_max
         write (*,*) ' ncsh + nvsh = ', ncsh + nvsh
         write (*,*) ' nssh is too large!'
        end if

! *****************************************************************************
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' =================== Question No. 3 =================== '
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' Determine number of orbitals for non-local potential. '
        write (*,*) '  '

        write (*,*) ' Do you want to use (insert corresponding number): '
        write (*,*) ' [1] a s basis '
        write (*,*) ' [2] a sp3 basis '
        write (*,*) ' [3] a sp3d5 basis '
        write (*,*) ' [4] a sp3d5f14 basis? '
        write (*,*) '  '
        write (*,*) ' Use the default? Y/N '
        read (*,*) answer 

        if (answer .eq. 'Y' .or. answer .eq. 'y' .or. answer .eq. ' ') then
         write (*,*) '  '
         write (*,*) ' OK: We USE the default. '
         numlmax = numL(nznuc)
        else 
         write (*,*) '  '
         write (*,*) ' OK: We change the defaults. '
         write (*,*) '  '
         write (*,*) ' Insert 1, 2, 3, or 4 '
         read (*,*) ioption
         if (ioption .eq. 1) then 
          numlmax = 1 - 1
          write (*,*) ' You have chosen numL = 1 '
         else if (ioption .eq. 2) then
          numlmax = 2 - 1
          write (*,*) ' You have chosen numL = 2 '
         else if (ioption .eq. 3) then
          numlmax = 3 - 1
          write (*,*) ' You have chosen numL = 3 '
         else if (ioption .eq. 4) then
          numlmax = 4 - 1
          write (*,*) ' You have chosen numL = 4 '
         else  
          numlmax = numL(nznuc)
         end if
        end if 

! *****************************************************************************
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' =================== Question No. 4 =================== '
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' Determine occupation number of each shells. '
        write (*,*) '  '
        if (ncoresh(nznuc) .gt. 0) then
         write (*,*) ' For core shells, the default occupation is:'
         do issh = 1, ncoresh(nznuc)
          write (*,300) nq(nznuc,issh), lq(nznuc,issh), xocc_in(nznuc,issh)
         end do

         write (*,*) '  '
         write (*,*) ' Use the default? Y/N '
         read (*,*) answer 

         if (answer .eq. 'Y' .or. answer .eq. 'y' .or. answer .eq. ' ') then
          write (*,*) '  '
          write (*,*) ' OK: We USE the default. '
          write (*,*) '  '
	  do issh = 1, ncoresh(nznuc)
	   xocc(issh) = xocc_in(nznuc,issh)
	  end do
         else
          write (*,*) '  '
          write (*,*) ' OK: We change the defaults. '
          write (*,*) '  '
          write (*,*) ' Which exchange-correlation functional would you like? ' 
          do issh = 1, ncoresh(nznuc)
           write (*,301) nq(nznuc,issh), lq(nznuc,issh)
           read (*,*) xocc(issh)
          end do
         end if
         write (*,*) '  '
         do issh = 1, ncoresh(nznuc)
          write (*,300) nq(nznuc,issh), lq(nznuc,issh), xocc(issh)
         end do
        end if
 
        write (*,*) '  '
        write (*,*) ' For valence shells, the default occupation is:'	
        do issh = ncoresh(nznuc) + 1, nssh(nznuc)
         write (*,300) nq(nznuc,issh), lq(nznuc,issh), xocc_in(nznuc,issh)
        end do

        write (*,*) '  '
        write (*,*) ' Use the default? Y/N '
        read (*,*) answer 

        if (answer .eq. 'Y' .or. answer .eq. 'y' .or. answer .eq. ' ') then
         write (*,*) '  '
         write (*,*) ' OK: We USE the default. '
         write (*,*) '  '
         do issh = ncoresh(nznuc) + 1, nssh(nznuc)
          xocc(issh) = xocc_in(nznuc,issh)
         end do
        else
         write (*,*) '  '
         write (*,*) ' OK: We change the defaults. '
         write (*,*) '  '
         do issh = ncoresh(nznuc) + 1, nssh(nznuc)
          write (*,301) nq(nznuc,issh), lq(nznuc,issh)
          read (*,*) xocc(issh)
         end do
        end if
        write(*,*)'  '
        do issh = ncoresh(nznuc) + 1, nssh(nznuc)
         write (*,300) nq(nznuc,issh), lq(nznuc,issh), xocc(issh)
        end do

! *****************************************************************************
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' =================== Question No. 5 =================== '
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' Which exchange-correlation functional do you want to '
        write (*,*) ' use? The standard one we have used is the Ceperely- '
        write (*,*) ' Alder form as parameterized by Perdew-Zunger '
        write (*,*) ' (ioption = 3). Here are the different options: '
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' 1   LDA    Wigner '
        write (*,*) ' 2   LDA    Hedin/Lundqvist '
        write (*,*) ' 3   LDA    Ceperley/Alder Perdew/Zunger (1980) '
        write (*,*) ' 4   GGA    Perdew/Wang (1991) '
        write (*,*) ' 5   GGA    Becke (1988) X, Perdew (1986) '
        write (*,*) ' 6   GGA    Perdew/Burke/Ernzerhof (1996) '
        write (*,*) ' 7   LDA    Zhao/Parr '
        write (*,*) ' 8   LDA    Ceperley/Alder Perdew/Wang (1991) '
        write (*,*) ' 9   GGA    Becke (1988) X, Lee/Yang/Parr (1988) '
        write (*,*) ' 10  GGA    Perdew/Wang (1991) X, Lee/Yang/Parr (1988) '
        write (*,*) ' 11  LSDA   Volko/Wisk/Nusair (1980) '
        write (*,*) ' 12  B3LYP  Mix exact exchange and BLYP '
        write (*,*) '  '
        write (*,*) ' Use the default? Y/N '
        read (*,*) answer 

        if (answer .eq. 'Y' .or. answer .eq. 'y' .or. answer .eq. ' ') then
         write (*,*) '  '
         write (*,*) ' OK: We USE the default. '
         write (*,*) '  '
         ioption = 3 
         write (*,*) ' You have chosen option = ', ioption
        else
         write (*,*) '  '
         write (*,*) ' OK: We change the defaults. '
         write (*,*) '  '
         write (*,*) ' Which exchange-correlation functional would you like? ' 
         read (*,*) ioption
         write (*,*) ' You have chosen option = ', ioption
        end if

        bmix = 0.5            
        if (ioption .eq. 12) then
         write (*,*) '  '
         write (*,*) ' For the B3LYP case, you have the choice for the degree '
         write (*,*) ' degree of the mixing. The default is 0.3 of exact '
         write (*,*) ' exchange and 0.7 of gradient exchange. '   
         write (*,*) ' Use the default? Y/N'
         read (*,*) answer 
         if (answer .eq. 'Y' .or. answer .eq. 'y' .or. answer .eq.' ') then
          write (*,*) '  '
          write (*,*) ' OK: We USE the defaults.'
          write (*,*) '  '
          exmix = 0.3            
         else
          write (*,*) '  '
          write (*,*) ' OK: We change the defaults. '
          write (*,*) '  '
          write (*,*) ' What fraction of exact exchange for the mixing would '
          write (*,*) ' you like? '
          read (*,*) exmix
         end if
 
         write (*,*) '  '
         write (*,*) ' The convergence of the all electron calculation can be '
         write (*,*) ' smaller.'
         write (*,*) ' Default for the potential mixing is 0.5 '  
         write (*,*) ' Use the default? Y/N'
         read (*,*) answer
         if (answer .eq. 'Y' .or. answer .eq. 'y' .or. answer .eq.' ') then
          write (*,*) '  '
          write (*,*) ' OK: We USE the defaults.'
          write (*,*) '  '
         else
          write (*,*) '  '
          write (*,*) ' OK: We change the defaults. '
          write (*,*) '  '
          write (*,*) ' What mixing fraction would you prefer? ' 
          read (*,*) bmix
         end if
        end if 

! *****************************************************************************
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' =================== Question No. 6 =================== '
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' Determine the type of pseudopotential function.  Please '
        write (*,*) ' check manual for avaiable choices.  If you are not sure '
        write (*,*) ' or do not care, hit return to use default value. '
        write (*,*) ' Here are the choices of pseudopotential scheme: '
        write (*,*) ' (1) h  Hamann type; (2) t  Troullier-Martins type. '
        write (*,*) '  '
        write (*,*) ' Use the default? Y/N '
        read (*,*) answer 

        if (answer .eq. 'Y' .or. answer .eq. 'y' .or. answer .eq. ' ') then
         write (*,*) '  '
         write (*,*) ' OK: We USE the default. '
         write (*,*) '  '
         ntypepsp = ntpsp(nznuc)
        else
         write (*,*) '  '
         write (*,*) ' OK: We change the defaults.'
         write (*,*) '  '
         write (*,*) ' Input your choice of pseudopotential: '
         write (*,*) ' Do this only when you know what you are doing!'
         read (*,*) ntypepsp
        end if

! *****************************************************************************
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' =================== Question No. 7 =================== '
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' Determine the Rnlc.  Please check manual for the meaning.'
        write (*,*) ' If you are not sure or do not care, hit return to use '
        write (*,*) ' default value.'
        write (*,*) '  '
        write (*,*) ' Use the default? Y/N '
        read (*,*) answer 

        if (answer .eq. 'Y' .or. answer .eq. 'y' .or. answer .eq. ' ') then
         write (*,*) '  '
         write (*,*) ' OK: We USE the default. '
         rnlc = rnlc_in(nznuc)
        else
         write (*,*) '  '
         write (*,*) ' OK: We change the defaults.'
         write (*,*) '  '
         write (*,*) ' Input your choice of Rnlc: '
         write (*,*) ' Do this only if you know what you are doing!'
         read (*,*)  rnlc
        end if

! *****************************************************************************
! Write to output file - ncpp.ini
        write (12,400) 1.0d0*nznuc, ncsh, nvsh, ioption, rnlc, exmix, bmix
        do issh = 1, ncsh + nvsh
         write (12,401) nq(nznuc,issh), lq(nznuc,issh), xocc(issh)
        end do
        if (ntypepsp .eq. 2) then
         pptype = ' t'
        else
         pptype = ' h'
        end if
        write (12,402) numlmax, pptype

        write (*,*) '  '
        write (*,*) ' Thank you for your input.'
        write (*,*) '  '
        write (*,*) ' You have created: ncpp.ini '
        write (*,*) ' This file will be used to generate the pseudopotential. '
        write (*,*) ' pseudopotential.'
        write (*,*) '  '

        close (unit = 12)

! Format Statements
! ===========================================================================
100     format (a50)
200	format (' Number of core and valence electron shells: ', i3, 2x, i3)
300	format (2x, 'n = ', i3, 2x, ' l = ', i3, 2x, ' xocc = ', f6.2)
301	format (2x, 'n = ', i3, ' and l = ', i3,                             &
     &          ' Now input occupation number!')
400     format (f6.2, 3(1x,i3), 1x, 3(f6.2))
401     format (2(1x,i3), 1x, f6.2)
402     format (i1, a2)

        end 
