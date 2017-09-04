
VERSION = 0.1

CC       = g++
cc	=gcc


#OPTIMIZE = -O3 -s -fomit-frame-pointer 
#CFLAGS   = $(DEFINES) $(OPTIMIZE)

LFLAGS   =   -L/opt/visivo/lib  -rdynamic -L/opt/vtk56/lib/vtk-5.6  -L/opt/cfitsio3_1_0/lib  -L/opt/hdf5-1.8.5-patch1/lib  -lVisIVOApiTest -lcfitsio -lcurl  -lssl -lcrypto -lidn -lrt -lhdf5 -lvtkRendering -lvtkHybrid -lvtkVolumeRendering -lvtkexoIIc -lvtkRendering -lvtkGraphics -lvtkverdict -lvtkImaging -lvtkftgl -lvtkfreetype -lXt -lSM -lICE -lX11 -lXext -lvtkIO -lvtkFiltering -lvtkCommon -lvtkDICOMParser -lvtkNetCDF -lvtkmetaio -lvtksqlite -lpthread -lvtkpng -lvtktiff -lvtkzlib -lvtkjpeg -lm -lvtkexpat -lvtksys -ldl /usr/local/lib64/libGL.so /usr/local/lib64/libOSMesa.so

LFLAGS   =   -L/opt/visivo/lib -lVisIVOApiTest -rdynamic -L/opt/vtk56/lib/vtk-5.6  -L/opt/cfitsio3_1_0/lib  -L/opt/hdf5-1.8.5-patch1/lib  -lcfitsio -lcurl  -lssl -lcrypto -lidn -lrt -lhdf5 -lvtkRendering -lvtkHybrid -lvtkVolumeRendering -lvtkexoIIc -lvtkRendering -lvtkGraphics -lvtkverdict -lvtkImaging -lvtkftgl -lvtkfreetype -lXt -lSM -lICE -lX11 -lXext -lvtkIO -lvtkFiltering -lvtkCommon -lvtkDICOMParser -lvtkNetCDF -lvtkmetaio -lvtksqlite -lpthread -lvtkpng -lvtktiff -lvtkzlib -lvtkjpeg -lm -lvtkexpat -lvtksys -ldl /usr/local/lib64/libGL.so /usr/local/lib64/libOSMesa.so

INCLUDE=-I/opt/visivo/include

PROGS    = myprogram.c

PROGS_O  =  myprogram.o
           
#$(PROGS_O):
#	$(cc) -c myprogram.c -o myprogram.o

all:
	$(cc) -c $(PROGS) $(INCLUDE) 
	$(CC) -o myprogram $(PROGS_O) $(CFLAGS) $(LFLAGS)  

clean           :  
	rm -f *.o 

