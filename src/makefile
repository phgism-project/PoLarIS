PHG_HOME  = ../../
PHG_LIB   =${PHG_HOME}/src/libphg

FILES = Makefile functions.h poisson.c simplest.c .cvsignore *-plot.sh jobs.sh \
	maxwell*.c heat.c non-smooth.c eigen.c elastic.c navier-stokes.[ch] \
	navier-stokes-unsteady.c

#default: test-grid
default: ins-flow
#default: test-solverT

all: lib poisson simplest maxwell maxwell-complex heat non-smooth eigen \
	elastic maxwell-eigen navier-stokes

examples.zip: ${FILES}
	@zip -9 -u -y $@ $^

clean:
	-/bin/rm -f *.o

distclean: clean
	-/bin/rm -f *.o ins-flow

matclean:
	-/bin/rm -f *.m *.m.dat

lib:
	@(cd ../src; $(MAKE))


include ${PHG_HOME}/Makefile.inc


EXTRAFLAGS=  -DDUMP_MAT_VEC=0


###########
#   INS   #
###########


ins-flow: ${PHG_LIB}${LIB_SUFFIX} ins-main.o ins-solver.o ins-utils.o save-load.o ins-system.o temp-solver.o ice-grid.o slip-bdry.o  moving-mesh.o  layers.o fv_vert.o core-shallow.o core-FO.o  netcdf-utils.o ins-bc.o 
	@echo "Linking ins-flow"
	${LINKER} ${LDFLAGS} -o $@ $^ ${LIBS}
#ins-mg.o multi-grid.o multi-grid-smoother.o mg-line-smoother.o




fv_vert.o: fv_vert.c
	@echo "compiling fv_vert.o ..."
	${CC} ${CFLAGS} ${CPPFLAGS} -O0 -c fv_vert.c 
