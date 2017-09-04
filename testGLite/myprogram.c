#include "visivo.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define NB 16777
#define NVOL 262144

int main(int argc, char*argv[])
{
int errorCode;

char filename[256];

//*********************************
//*********************************
//********************************* VisIVOImporter
VisIVOImporter envVI1;

errorCode=VI_Init(&envVI1);
errorCode=VI_SetAtt(&envVI1,VI_SET_FFORMAT,"ascii");
errorCode=VI_SetAtt(&envVI1,VI_SET_FILEPATH,"lfn://grid/cometa/ube/mrvbt16.ascii");
//errorCode=VI_SetAtt(&envVI1,VI_SET_OUTFILEVBT,"mrvbt16.bin");
errorCode=VI_SetAtt(&envVI1, VI_SET_LFNOUT,"lfn://grid/cometa/ube/mrvbt16.bin");
errorCode=VI_SetAtt(&envVI1, VI_SET_VO,"cometa");


VI_Import(&envVI1);
//*********************
//*********************
//********************* VisIVOFilter
VisIVOFilter envVF1;

char operation[256];
strcpy(operation,"pointproperty");

errorCode=VF_Init(&envVF1);
errorCode=VF_SetAtt(&envVF1,VF_SET_OPERATION,operation);
errorCode=VF_SetAtt(&envVF1,VF_SET_FILEVBT,"lfn://grid/cometa/ube/mrvbt16.bin");
errorCode=VF_SetAtt(&envVF1,VF_SET_RESOLUTION,"32 32 32");
errorCode=VF_SetAtt(&envVF1,VF_SET_POINTCOLUMNS,"X Y Z");
errorCode=VF_SetAtt(&envVF1,VF_SET_APPEND,"");
errorCode=VF_SetAtt(&envVF1, VF_SET_OUTCOL,"density");
errorCode=VF_SetAtt(&envVF1, VF_SET_VO,"cometa");

VF_Filter(&envVF1);
//*********************************
//*********************************
//********************************* VisIVOViewer
VisIVOViewer envVV1;

errorCode=VV_Init(&envVV1);
errorCode=VV_SetAtt(&envVV1,VV_SET_FILEVBT,"lfn://grid/cometa/ube/mrvbt16.bin");
errorCode=VV_SetAtt(&envVV1,VV_SET_FIELD,"X Y Z");  // X Y Z in the VBT file
errorCode=VV_SetAtt(&envVV1,VV_SET_COLOR,"");
errorCode=VV_SetAtt(&envVV1,VV_SET_COLORSCALAR,"density");
//errorCode=VV_SetAtt(&envVV1,VV_SET_OUT,"VVmrvbt16_");
errorCode=VV_SetAtt(&envVV1,VV_SET_BOX,"");
errorCode=VV_SetAtt(&envVV1,VV_SET_AXES,"");
errorCode=VV_SetAtt(&envVV1,VV_SET_PALETTE,"");
//errorCode=VV_SetAtt(&envVV1,VV_SET_COLORRANGEFROM,"0");
//errorCode=VV_SetAtt(&envVV1,VV_SET_COLORRANGETO,"400");
errorCode=VV_SetAtt(&envVV1,VV_SET_OPACITY,"0.6");
errorCode=VV_SetAtt(&envVV1, VV_SET_VO,"cometa");
errorCode=VV_SetAtt(&envVV1, VV_SET_LFNOUT,"lfn://grid/cometa/ube/VVmrvbt16_");


//cycle File 

float camera[6];
camera[0]=0.0; //Az start
camera[1]=60.0; //Az End
camera[2]=0.0; //El start
camera[3]=10.0; //El End
camera[4]=1.0; //Zoom start
camera[5]=1.0; //Zoom end

int zoomend=1;
float zstepframe[1];
zstepframe[0]=0.1;
int framesec=5;
int length=2;

float campos[6];
float camfp[6];
float camroll[2];
camroll[0]=0;
camroll[1]=60;
int i;
for(i=0;i<6;i++) campos[i]=35;
for(i=0;i<6;i++) camfp[i]=NOSET_CAM;

errorCode==VV_SetCameraPath(&envVV1, 0, &camera[0], zoomend,zstepframe,framesec,length, campos, camfp, camroll, 0);


//*********************************
//*********************************
//*********************************

printf("\n Call VV_View \n");
VV_View(&envVV1);


printf("End\n");
return 0;
}
    
