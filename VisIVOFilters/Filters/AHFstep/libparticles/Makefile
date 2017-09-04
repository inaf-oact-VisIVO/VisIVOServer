#**
# Makefile for libparticles.a, a library for all things related to moving, assigning, etc. of particles
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

OBJS = 	density.o \
       	drift_pos.o\
	kick_mom.o\
	leavers.o\
	llmsort.o\
	lltools.o\
	relink.o \
	move_part.o

# Catch empty make call
all:
	${MAKE} libparticles.a

# Cleaning up
clean:
	rm -f libparticles.a $(OBJS)

# The library itself
libparticles.a: $(OBJS)
	${AR} -r libparticles.a $(OBJS)


# The individual routines stuff
density.o: density.c
	$(CC) $(CFLAGS) -c -o density.o density.c

drift_pos.o: drift_pos.c
	$(CC) $(CFLAGS) -c -o drift_pos.o drift_pos.c

kick_mom.o: kick_mom.c
	$(CC) $(CFLAGS) -c -o kick_mom.o kick_mom.c

leavers.o: leavers.c
	$(CC) $(CFLAGS) -c -o leavers.o leavers.c

llmsort.o: llmsort.c
	$(CC) $(CFLAGS) -c -o llmsort.o llmsort.c

lltools.o: lltools.c
	$(CC) $(CFLAGS) -c -o lltools.o lltools.c

relink.o: relink.c
	$(CC) $(CFLAGS) -c -o relink.o relink.c

move_part.o: move_part.c
	$(CC) $(CFLAGS) -c -o move_part.o move_part.c


