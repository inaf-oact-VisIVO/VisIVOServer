
VERSION = 0.1

CC      = g++
cc	= gcc


#OPTIMIZE = -O3 -s -fomit-frame-pointer 
#CFLAGS   = $(DEFINES) $(OPTIMIZE)


LFLAGS   =   -L/opt/lcg/lib64 -llcg_util -L/usr/local/VisIVOServerTools/lib -lVisIVOApiTest -rdynamic -L/opt/vtk-5.6.0/lib/vtk-5.6 -lcfitsio -L/opt/curl-7.19.4/lib -lcurl -lrt -lhdf5 -lvtkRendering -lvtkHybrid -lvtkVolumeRendering -lvtkexoIIc -lvtkRendering -lvtkGraphics -lvtkverdict -lvtkImaging -lvtkftgl -lvtkfreetype -lXt -lSM -lICE -lX11 -lXext -lvtkIO -lvtkFiltering -lvtkCommon -lvtkDICOMParser -lvtkNetCDF -lvtkmetaio -lvtksqlite -lpthread -lvtkpng -lvtktiff -lvtkzlib -lvtkjpeg -lm -lvtkexpat -lvtksys -ldl /opt/lib64/libGL.so /opt/lib64/libOSMesa.so

INCLUDE=-I/usr/local/VisIVOServerTools/include -I/opt/lcg/include

PROGS    = myprogram.c

PROGS_O  =  myprogram.o
           
#$(PROGS_O):
#	$(cc) -c myprogram.c -o myprogram.o

all:
	$(cc) -c $(PROGS) $(INCLUDE) 
	$(CC) -o myprogram $(PROGS_O) $(CFLAGS) $(LFLAGS)  

clean           :  
	rm -f *.o 

