# This is the compiler options for LINUX only
F77 = ifort  -no-ipo
F90 = ifort  -no-ipo
FFLAGS = -g 

MODULES = precision.o constants.o
#XC = get_uxc.o ggacrad.o ggaxrad.o lsdavwn.o


OBJECTS1 = qlqgrid.o polyint.o

OBJECTS2 = localqgrid.o polyint.o

OBJECTS3 = linear.o polyint.o

OBJECTS4 = mk_pseudoFile.o 

all: qlqgrid.x localqgrid.x linear.x pickanatom.x mk_pseudoFile.x


qlqgrid.x : ${OBJECTS1} 
	${F77} ${OBJECTS1} -o qlqgrid.x

localqgrid.x : ${OBJECTS2} 
	${F77} ${OBJECTS2} -o localqgrid.x

linear.x : ${OBJECTS3} 
	${F77} ${OBJECTS3} -o linear.x

IMODULES = precision.o
INITIAL = pickanatom.o periodensystem.o
pickanatom.x: $(IMODULES) $(INITIAL)
	$(F90) $(IMODULES) $(INITIAL) -o pickanatom.x

mk_pseudoFile.x : ${OBJECTS4} 
	${F77} ${OBJECTS4} -o mk_pseudoFile.x

%.o :: %.f
	${F77} -c $<

.PHONY: clean veryclean

clean:
	rm -f -r core *.o .nfs* rii_files *.x.ip*  .#* ldtmp* *.mod *.vo

veryclean: clean
	rm -f *.x *.mod

# *****************************************************************************
# junk objects - just to get things to work
# *****************************************************************************
qlqgrid.o : qlqgrid.f
	$(F77) $(FFLAGS) -c qlqgrid.f
polyint.o : polyint.f
	$(F77) $(FFLAGS) -c polyint.f
localqgrid.o : localqgrid.f
	$(F77) $(FFLAGS) -c localqgrid.f
linear.o : linear.f
	$(F77) $(FFLAGS) -c linear.f
mk_pseudoFile.o : mk_pseudoFile.f
	$(F77) $(FFLAGS) -c mk_pseudoFile.f
   
 
# *****************************************************************************
# Modules
# *****************************************************************************
precision.o : MODULES/precision.f90
	$(F90) $(FFLAGS) -c MODULES/precision.f90
constants.o : MODULES/constants.f90
	$(F90) $(FFLAGS) -c MODULES/constants.f90
 
 
# *****************************************************************************
# pickanatom.x objects
# *****************************************************************************
pickanatom.o : pickanatom.f90
	$(F90) $(FFLAGS) -c pickanatom.f90
periodensystem.o : periodensystem.f90
	$(F90) $(FFLAGS) -c periodensystem.f90
