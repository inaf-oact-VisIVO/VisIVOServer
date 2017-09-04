#******************************************************************#
#                                                                  #
# Makefile for              AMIGA                                  #
#                           =====                                  #
#                                                                  #
# version 1.0 (30/01/2005), Alexander Knebe                        #
#                                                                  #
# Lightcone output added by Enn Saar (02/10/2003)                  #
# AHF (AMIGA-Halo-Finder) added by Stuart Gill (30/01/2005)        #
#                                                                  #
#******************************************************************#

.PHONY: clean veryclean

#*******************************************************************#
AMIGA:			DEFINEFLAGS = ${MASTER_DEFINEFLAGS}
AMIGAHF:			DEFINEFLAGS = ${MASTER_DEFINEFLAGS} -DAHF
AHFstep:			DEFINEFLAGS = ${MASTER_DEFINEFLAGS} -DAHFstep -DAHF
CONVERT:			DEFINEFLAGS = ${MASTER_DEFINEFLAGS} -DCONVERT_TERM -DVERBOSE

export DEFINEFLAGS

# libraries in use
# ----------------
LIB_MATH 	=	m

# flags for compiling all individual routines
# -------------------------------------------
CFLAGS	=	${OPTIMIZE} ${CCFLAGS} ${DEFINEFLAGS}

# flags for linking all parts together
# ------------------------------------
LFLAGS	=	${LNFLAGS}

# all object files
# ----------------
OBJS	=	\
	comm.o \
	common.o\
	startrun.o\
	step.o

LIBS =	libahf/libahf.a \
			libgravity/libgravity.a \
			libgrids/libgrids.a \
			libio/libio.a \
			libio_serial/libio_serial.a \
			libmhd/libmhd.a \
			libparticles/libparticles.a \
			libsfc/libsfc.a \
			libutility/libutility.a

# what's to make ?
# ----------------
AMIGA: 		main.o ${OBJS} ${LIBS} 
	$(CC) ${PROFILE} -o AMIGA main.o ${OBJS} ${LIBS} $(LIBS) -l${LIB_MATH} ${LFLAGS} 

AHFstep:	main.o ${OBJS} ${LIBS}
	$(CC) ${PROFILE} -o AHFstep main.o  ${OBJS} ${LIBS} $(LIBS) -l${LIB_MATH} ${LFLAGS} 

AMIGAHF:	main.o ${OBJS} ${LIBS}
	$(CC) ${PROFILE} -o AMIGAHF main.o  ${OBJS} ${LIBS} $(LIBS) -l${LIB_MATH} ${LFLAGS} 

CONVERT:	main.o ${OBJS} ${LIBS}
	$(CC) ${PROFILE} -o CONVERT main.o  ${OBJS} ${LIBS} $(LIBS) -l${LIB_MATH} ${LFLAGS} 

clean:
	rm -f *.o AMIGA* AHF* CONVERT *.il
	cd libsfc ; ${MAKE} clean
	cd libio ; ${MAKE} clean
	cd libio_serial ; ${MAKE} clean
	cd libahf ; ${MAKE} clean
	cd libmhd ; ${MAKE} clean
	cd libgravity ; ${MAKE} clean
	cd libgrids ; ${MAKE} clean
	cd libparticles ; ${MAKE} clean
	cd libutility ; ${MAKE} clean

veryclean:
	rm -f *.o AMIGA* AHF* CONVERT *~ *~.* *.il
	cd libsfc ; ${MAKE} clean
	cd libio ; ${MAKE} clean
	cd libio_serial ; ${MAKE} clean
	cd libahf ; ${MAKE} clean
	cd libmhd ; ${MAKE} clean
	cd libgravity ; ${MAKE} clean
	cd libgrids ; ${MAKE} clean
	cd libparticles ; ${MAKE} clean
	cd libutility ; ${MAKE} clean


# dependencies for include files
# ------------------------------
common.h:	param.h\
				tdef.h
tdef.h: 		param.h



# how to make the libraries?
# ------------------------------
libsfc/libsfc.a:
	cd libsfc ; ${MAKE} libsfc.a

libio/libio.a:
	cd libio ; ${MAKE} libio.a

libio_serial/libio_serial.a:
	cd libio_serial ; ${MAKE} libio_serial.a

libahf/libahf.a:
	cd libahf ; ${MAKE} libahf.a

libmhd/libmhd.a:
	cd libmhd ; ${MAKE} libmhd.a

libgravity/libgravity.a:
	cd libgravity ; ${MAKE} libgravity.a

libgrids/libgrids.a:
	cd libgrids ; ${MAKE} libgrids.a

libparticles/libparticles.a:
	cd libparticles ; ${MAKE} libparticles.a

libutility/libutility.a:
	cd libutility ; ${MAKE} libutility.a



#dependencies for .o files
# ------------------------
comm.o:	comm.c
	${CC} ${CFLAGS} -c comm.c

main.o:			main.c
	${CC} ${CFLAGS} -c main.c
	
startrun.o:		startrun.c
	${CC} ${CFLAGS} -c startrun.c

step.o:			step.c
	${CC} ${CFLAGS} -c step.c
