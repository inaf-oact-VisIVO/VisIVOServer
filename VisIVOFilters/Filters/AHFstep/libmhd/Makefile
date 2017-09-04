#**
# Makefile for libmhd.a, a library for all things related to (magneto-)hydrodynamics
#
# version 0.0 (09.09.2008), Alexander Knebe (aknebe@aip.de)
#**

# If ${OPTIMIZE}, ${CCFLAGS} and/or ${DEFINEFLAGS} are set, they will be
# used for the compilation
CFLAGS+=${OPTIMIZE} ${CCFLAGS} ${DEFINEFLAGS}

# If a different MAKE is specified use that one, otherwise try your luck
# with just make; maybe it will work.
MAKE?=make

# The same for ar
AR?=ar

# PHONY targets
.PHONY: all clean

OBJS = magneto.o \
       hydro.o\
		 hydro_tests.o

# Catch empty make call
all:
	${MAKE} libmhd.a

# Cleaning up
clean:
	rm -f libmhd.a $(OBJS)

# The library itself
libmhd.a: $(OBJS)
	${AR} -r libmhd.a $(OBJS)

# The individual routines stuff
magneto.o: magneto.c
	$(CC) $(CFLAGS) -c -o magneto.o magneto.c

hydro.o: hydro.c
	$(CC) $(CFLAGS) -c -o hydro.o hydro.c

hydro_tests.o: hydro_tests.c
	$(CC) $(CFLAGS) -c -o hydro_tests.o hydro_tests.c


