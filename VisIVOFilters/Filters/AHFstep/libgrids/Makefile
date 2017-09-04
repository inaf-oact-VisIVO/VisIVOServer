#**
# Makefile for libgrids.a, a library for all things related to the grid structures
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

OBJS = 	coarse_to_fine.o \
       	fine_to_coarse.o\
			generate_grids.o\
			get_nnodes.o\
			go_down.o\
			go_up.o\
			refine_grid.o\
			store_refs.o

# Catch empty make call
all:
	${MAKE} libgrids.a

# Cleaning up
clean:
	rm -f libgrids.a $(OBJS)

# The library itself
libgrids.a: $(OBJS)
	${AR} -r libgrids.a $(OBJS)


# The individual routines stuff
coarse_to_fine.o: coarse_to_fine.c
	$(CC) $(CFLAGS) -c -o coarse_to_fine.o coarse_to_fine.c

fine_to_coarse.o: fine_to_coarse.c
	$(CC) $(CFLAGS) -c -o fine_to_coarse.o fine_to_coarse.c

generate_grids.o: generate_grids.c
	$(CC) $(CFLAGS) -c -o generate_grids.o generate_grids.c

get_nnodes.o: get_nnodes.c
	$(CC) $(CFLAGS) -c -o get_nnodes.o get_nnodes.c

go_down.o: go_down.c
	$(CC) $(CFLAGS) -c -o go_down.o go_down.c

go_up.o: go_up.c
	$(CC) $(CFLAGS) -c -o go_up.o go_up.c

refine_grid.o: refine_grid.c
	$(CC) $(CFLAGS) -c -o refine_grid.o refine_grid.c

store_refs.o: store_refs.c
	$(CC) $(CFLAGS) -c -o store_refs.o store_refs.c



