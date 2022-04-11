! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
! Universidad de Madrid - Jose Ortega

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang 
! Ohio University - Dave Drabold

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in
! subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! periodensystem.f90
! Program Description
! ===========================================================================
!       Prints out a periodic table
! 
! Program Declaration
! ===========================================================================
        subroutine periodensystem
      	implicit none

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
	write (*,*) '  '
	write (*,*) ' ------------------------------------------------------- '
	write (*,*) ' |1 |                                               |2 | '
	write (*,*) ' |H |                                               |He| '
	write (*,*) ' ------------------------------------------------------- '
	write (*,*) ' |3 |4 |                             |5 |6 |7 |8 |9 |10| '
	write (*,*) ' |Li|Be|                             |B |C |N |O |F |Ne| '
	write (*,*) ' ------------------------------------------------------- '
	write (*,*) ' |11|12|                             |13|14|15|16|17|18| '
	write (*,*) ' |Na|Mg|                             |Al|Si|P |S |Cl|Ar| '
	write (*,*) ' ------------------------------------------------------- '
	write (*,*) ' |19|20|21|22|23|24|25|26|27|28|29|30|31|32|33|34|35|36| '
	write (*,*) ' |K |Ca|Sc|Ti|V |Cr|Mn|Fe|Co|Ni|Cu|Zn|Ga|Ge|As|Se|Br|Kr| '
	write (*,*) ' ------------------------------------------------------- '
	write (*,*) ' |37|38|39|40|41|42|43|44|45|46|47|48|49|50|51|52|53|54| '
	write (*,*) ' |Rb|Sr|Y |Zr|Nb|Mo|Tc|Ru|Rh|Pd|Ag|Cd|In|Sn|Sb|Te|I |Xe| '
	write (*,*) ' ------------------------------------------------------- '
	write (*,*) ' |55|56|57|72|73|74|75|76|77|78|79|80|81|82|83|84|85|86| '
	write (*,*) ' |Cs|Ba|La|Hf|Ta|W |Re|Os|Ir|Pt|Au|Hg|Tl|Pb|Bi|Po|At|Rn| '
	write (*,*) ' ------------------------------------------------------- '
	write (*,*) ' |87|88|89| '
	write (*,*) ' |Fr|Ra|Ac|  ------------------------------------------- '
	write (*,*) ' ----------  |58|59|60|61|62|63|64|65|66|67|68|69|70|71| '
	write (*,*) '             |Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu| '
	write (*,*) '             ------------------------------------------- '
	write (*,*) '             |90|91|92|83|94| '
	write (*,*) '             |Th|Pa|U |Np|Pu| '
	write (*,*) '             ---------------- '
	write (*,*) '  '

! Format Statements
! ===========================================================================

        return 
        end 
