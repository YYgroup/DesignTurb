# Makefile for compiling with personal F90 and 
# To compile and link write 'make turb' or simply 'make'
#
#
#
SHELL = /bin/bash
FFLAG = -O4 -w 

IDIR  = -I/WORK2/pku_yyg/shenwy/opt/fftw-3.3.8/include
LDIR  = -L/WORK2/pku_yyg/shenwy/opt/fftw-3.3.8/lib
FCOMP = mpif90 -c ${FFLAG} 

LINK  = mpif90
LIBS  = -lfftw3 -lm
OBJ   = tube.o subroutine.o

.SUFFIXES: .o .f90

.f90.o:
	${FCOMP} $*.f90

comu:  ${OBJ}
	${LINK} -o tube ${OBJ} ${LDIR} ${LIBS} 

clean:
	rm -f *.o *~ tube *.out out

