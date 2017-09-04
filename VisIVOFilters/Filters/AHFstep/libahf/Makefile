#**
# Makefile for libahf.a, a library for all things AHF
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

OBJS = ahf_gridinfo.o \
       ahf_halos_sfc.o \
       ahf_halos.o

# Catch empty make call
all:
	${MAKE} libahf.a

# Cleaning up
clean:
	rm -f libahf.a $(OBJS)

# The library itself
libahf.a: $(OBJS)
	${AR} -r libahf.a $(OBJS)

# The individual routines
ahf_gridinfo.o: ahf_gridinfo.c
	$(CC) $(CFLAGS) -c -o ahf_gridinfo.o ahf_gridinfo.c

ahf_halos.o: ahf_halos.c
	$(CC) $(CFLAGS) -c -o ahf_halos.o ahf_halos.c

ahf_halos_sfc.o: ahf_halos_sfc.c ahf_halos_sfc.h
	$(CC) $(CFLAGS) -c -o ahf_halos_sfc.o ahf_halos_sfc.c
