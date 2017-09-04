#**
# Makefile for libmhd.a, a library for all things related to (magneto-)hydrodynamics
#
# version 0.0 (04.09.2008), Alexander Knebe (aknebe@aip.de)
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

OBJS =	input.o \
	output.o\
	input_grid.o\
	output_grid.o\
	output_lc.o\
	io_logging.o\
	write_boundary.o \
       write_density.o \
       write_divB.o \
       write_forces.o \
       write_hydro.o \
       write_nodepart.o \
       write_nodes.o \
       write_positions.o \
       write_residual.o

# Catch empty make call
all:
	${MAKE} libio_serial.a

# Cleaning up
clean:
	rm -f libio_serial.a $(OBJS)

# The library itself
libio_serial.a: $(OBJS)
	${AR} -r libio_serial.a $(OBJS)

# The individual routines stuff
input.o: input.c
	$(CC) $(CFLAGS) -c -o input.o input.c

output.o: output.c
	$(CC) $(CFLAGS) -c -o output.o output.c

input_grid.o: input_grid.c
	$(CC) $(CFLAGS) -c -o input_grid.o input_grid.c

output_grid.o: output_grid.c
	$(CC) $(CFLAGS) -c -o output_grid.o output_grid.c

output_lc.o: output_lc.c
	$(CC) $(CFLAGS) -c -o output_lc.o output_lc.c

io_logging.o: io_logging.c
	$(CC) $(CFLAGS) -c -o io_logging.o io_logging.c

write_boundary.o:	write_boundary.c
	$(CC) $(CFLAGS) -c -o write_boundary.o write_boundary.c

write_density.o: write_density.c
	$(CC) $(CFLAGS) -c -o write_density.o write_density.c

write_divB.o: write_divB.c
	$(CC) $(CFLAGS) -c -o write_divB.o write_divB.c

write_forces.o: write_forces.c
	$(CC) $(CFLAGS) -c -o write_forces.o write_forces.c

write_hydro.o: write_hydro.c
	$(CC) $(CFLAGS) -c -o write_hydro.o write_hydro.c

write_nodepart.o: write_nodepart.c
	$(CC) $(CFLAGS) -c -o write_nodepart.o write_nodepart.c

write_nodes.o: write_nodes.c
	$(CC) $(CFLAGS) -c -o write_nodes.o write_nodes.c

write_positions.o: write_positions.c
	$(CC) $(CFLAGS) -c -o write_positions.o write_positions.c

write_residual.o: write_residual.c
	$(CC) $(CFLAGS) -c -o write_residual.o write_residual.c


