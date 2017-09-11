/***************************************************************************
 *   Copyright (C) 2012 by Ugo Becciani   *
 *   ugo.becciani@oact.inaf.it   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#ifdef WIN32
    #include <time.h>
#endif
#include "vsmuportalop.h"
#include "vstable.h"

#include <vector>
#include <math.h>
#include "vsVector.h"
#include "vsLine.h"
#include "vsVoxel.h"


#ifndef VOLUME_PARAMETERS 
#define VOLUME_PARAMETERS //parametri del volume da visualizzare (spazio tra i piani interni del rivelatore) in cm
#define DIM_X dimX 
#define DIM_Y dimY
#define DIM_Z dimZ
#define DIM_VOXEL dimVox
#define DIM_VETTORE DIM_X*DIM_Y*DIM_Z/(DIM_VOXEL*DIM_VOXEL*DIM_VOXEL)
#endif

#define D_PIANI_FISICI 1 //distanza in cm tra i due piani fisici appartenenti allo stesso piano logico 
#define D_PIANI_MIN dMin//distanza in cm tra i piani 1 e 2 e tra i piani 3 e 4
#define D_PIANI_MAX dMax//distanza in cm tra i piani 2 e 3 in cm
#define D_PIANI_TOT D_PIANI_MAX+(D_PIANI_MIN*2) //distanza in cm tra il piano 1 e il piano 4



//prova
const unsigned int VSMuPortalOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSMuPortalOp::MIN_NUMBER_OF_ROW = 100;

//---------------------------------------------------------------------
VSMuPortalOp::VSMuPortalOp()
{
    m_fArray=NULL;
    m_outPointsfArray=NULL;
    m_outVolfArray=NULL;
    m_nOfRow=0;
    m_nOfCol=10;
    m_nOfColPoints=7;
    m_nOfColVol=7;
    m_volDim=60*30*30;
}
//---------------------------------------------------------------------


//---------------------------------------------------------------------
VSMuPortalOp::~VSMuPortalOp()
{	
    if(m_fArray!=NULL)
        for(unsigned int i=0;i<m_nOfCol;i++)
        {
            if(m_fArray[i] != NULL) delete [] m_fArray[i];
        }
	if(m_fArray!=NULL) delete [] m_fArray;
    
    if(m_outPointsfArray!=NULL)
        for(unsigned int i=0;i<m_nOfColPoints;i++)
        {
            if(m_outPointsfArray[i] != NULL) delete [] m_outPointsfArray[i];
        }
	if(m_outPointsfArray!=NULL) delete [] m_outPointsfArray;
    
    if(m_outVolfArray!=NULL)
        for(unsigned int i=0;i<m_nOfColVol;i++)
        {
            if(m_outVolfArray[i] != NULL) delete [] m_outVolfArray[i];
        }
	if(m_outVolfArray!=NULL) delete [] m_outVolfArray;
    
    
}
//---------------------------------------------------------------------

//---------------------------------------------------------------------
void VSMuPortalOp::printHelp()
//---------------------------------------------------------------------
{
/*
 USO: ./<nome eseguibile> <input data file name> <output scattering Vectors file name> <output voxels file name> <dim voxel>
    
    *****<input data file name>:
    
    file di dati simulati con Geant così generato:
    
    Struttura geometrica e materiali inseriti in GEANT:
    -8 Piani di rivelazione fisici (4 piani logici) con strip di scintillatore
    plastico (da x=-3 m a x=+3 m, da y=-1.5 m a y =+1.5 m) di spessore 1 cm
    -Supporto per ciascun piano di rivelazione logico, in Alluminio da 2mm di
    spessore e dimensioni 3 m x 6 m
    - Tetto e pavimento del container in ferro da 3 mm di spessore
    - Cubetto di Uranio da 10cmx10cmx10cm posto al centro (0,0,0) per il file
    simdata_run009.out
    - Struttura a "CT" per il file simdata_run011.out
    
    La struttura dei file dati comprende per ogni linea (cioè per ogni evento)
    10 variabili:
    N.evento
    X_A coordinata X piano fisico A
    Y_B coordinata Y piano fisico B
    X_C coordinata X piano fisico C
    Y_D coordinata Y piano fisico D
    X_E coordinata X piano fisico E
    Y_F coordinata Y piano fisico F
    X_G coordinata X piano fisico G
    Y_H coordinata Y piano fisico H
    P   modulo impulso muone iniziale
    
Note:
    a) Il N. evento non è progressivo, perché gli eventi presenti sul file
    sono gli eventi selezionati in base al criterio che ogni piano abbia 1 o 2
    hit (strip colpite). Gli eventi che in origine non avevano hit su qualcuno
    dei piani (o ne avevano troppi) sono stati scartati, per simulare quanto
    avverrà nella realtà.
    
    b) Le unità di misura e il formato delle variabili sono:
    N.evento:  intero
    Coordinate X,Y...  floating (cm)
    Impulso muone      floating (GeV/c)
    
    c) Le coordinate Z dei vari piani non sono state scritte sul file, in
    quanto costanti per tutti gli eventi. Esse sono
    Z per i piani A,B:    (+250 + 249)/2 = +249.5 cm
    Z per i piani C,D:    (+150 + 149)/2 = +149.5 cm
    Z per i piani E,F:    (-149 - 150)/2 = -149.5 cm
    Z per i piani G,H:    (-249 - 250)/2 = -249.5 cm
    
    d) Ognuno dei 4 punti nello spazio sarà dunque identificato dalle terne:
    (X_A,Y_B,Z_AB)
    (X_C,Y_D,Z_CD)
    (X_E,Y_F,Z_EF)
    (X_G,Y_H,Z_GH)
    
    *******<output scattering Vectors file name>
    file di output nel quale vengono salvati i
    punti di scattering con i rispettivi angoli di scattering
    
    *******<output voxels file name>
    file di output nel quale salvare la lista dei voxels che memorizza la sommatoria dei thetaQuadro
*/
    
    
    std::cout<<"Produce two tables. The first one contains scattering points and angles. The second is a volume: each mesh point contains the square sum of the scattering angle."<<std::endl;
    std::cout<<"The input file must be the output from the muportal importer with 10 column: Event number, X_A Y_B X_C Y_D X_E Y_F X_G Y_H (8 values coordinates in cm at the planes of the system), Energy pulse in GeV/C."<<std::endl<<std::endl;
    std::cout<<"Usage: VisIVOFilters --op poca [--resolution x_res y_res z_res] [--dimvox voxel_size] [--trackplanedist distance] [--innerdist distance]  [--outpoints points.bin] [--outvol vol.bin] [--history] [--historyfile filename.xml] [--help] [--file] inputFile.bin"<<std::endl<<std::endl;
    
    std::cout<<"Example: VisIVOFilters --op poca --resolution 600 300 300 --dimvox 10 --outpoints points.bin --outvol vol.bin --file inputFile.bin"<<std::endl;
    
    std::cout<<"Note:  "<<std::endl;
    std::cout<<"--resolution  3D mesh size in cm. Default value 600 300 300"<<std::endl;
    std::cout<<"--dimvox cubic voxel dimension in cm. Default value 10"<<std::endl;
    std::cout<<"--trackplanedist Distance in cm between planes 1 - 2  and planes  3 - 4. Default value 100"<<std::endl;
    std::cout<<"--innerdist. Distance in cm between planes 3 - 4. Default value 300 "<<std::endl;
    std::cout<<"--outpoints Name of the new table containing poca points."<<std::endl;
    std::cout<<"--outvol. Name of the new table containing the volume having the theta square value sum in each voxel."<<std::endl;
    std::cout<<"--history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
    std::cout<<"--historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;

    std::cout<<"--file Input table filename."<<std::endl;
    
    std::cout<<"--help produce this output."<<std::endl;
    
    return;
    
}

//---------------------------------------------------------------------
bool VSMuPortalOp::allocateArray()
//---------------------------------------------------------------------
{
	m_fArray=new  float*[m_nOfCol];
    m_outPointsfArray=new  float*[m_nOfColPoints];
    m_outVolfArray=new  float*[m_nOfColVol];
    
    for(unsigned int i=0;i<m_nOfCol;i++)
    {
            try
            {
                m_fArray[i]=new  float[m_tables[0]->getNumberOfRows()];
            }
            catch(std::bad_alloc &e)
            {
                m_fArray[i]=NULL;
            }
            
			if(m_fArray[i]==NULL)
			{
                    return false;
            }
	}
    for(unsigned int i=0;i<m_nOfColPoints;i++)
    {
        try
        {
            m_outPointsfArray[i]=new  float[m_tables[0]->getNumberOfRows()];
        }
        catch(std::bad_alloc &e)
        {
            m_outPointsfArray[i]=NULL;
        }
        
        if(m_outPointsfArray[i]==NULL)
        {
            return false;
        }
	}

    for(unsigned int i=0;i<m_nOfColVol;i++)
    {
        try
        {
            m_outVolfArray[i]=new  float[m_volDim];
        }
        catch(std::bad_alloc &e)
        {
            m_outVolfArray[i]=NULL;
        }
        
        if(m_outVolfArray[i]==NULL)
        {
            return false;
        }
	}

    
    return true;
}


//---------------------------------------------------------------------
bool VSMuPortalOp::execute()
//---------------------------------------------------------------------
{
    FILE *f, *fPoints, *fVoxel;
    
    int dimX=600, dimY=300, dimZ=300, dimVox=10, dMin=100, dMax=300;
    int i = 1;
    
    std::string fileNameOutputPoints;
    fileNameOutputPoints=getParameterAsString("outpoints");
    if(fileNameOutputPoints==""||fileNameOutputPoints=="unknown")
    {
        std::stringstream fileNameOutputSStream;
        std::string filenameInpuTable=m_tables[0]->getLocator();
        int len=filenameInpuTable.length();
        time_t rawtime;
        struct tm * timeinfo;
        char buffer [80];
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        strftime (buffer,80,"%Y%m%d%H%M",timeinfo);
        fileNameOutputSStream<<filenameInpuTable.substr(0, len-4)<<"_pocapointsop_"<<buffer<<".bin";
        fileNameOutputPoints=fileNameOutputSStream.str();
    }
    if(fileNameOutputPoints.find(".bin") == std::string::npos)
        fileNameOutputPoints.append(".bin");
    
    std::string fileNameOutputVol;
    fileNameOutputVol=getParameterAsString("outvol");
    if(fileNameOutputVol==""||fileNameOutputVol=="unknown")
    {
        std::stringstream fileNameOutputSStream;
        std::string filenameInpuTable=m_tables[0]->getLocator();
        int len=filenameInpuTable.length();
        time_t rawtime;
        struct tm * timeinfo;
        char buffer [80];
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        strftime (buffer,80,"%Y%m%d%H%M",timeinfo);
        fileNameOutputSStream<<filenameInpuTable.substr(0, len-4)<<"_pocavolop_"<<buffer<<".bin";
        fileNameOutputVol=fileNameOutputSStream.str();
    }
    if(fileNameOutputVol.find(".bin") == std::string::npos)
        fileNameOutputVol.append(".bin");

    if(fileNameOutputVol==fileNameOutputPoints)
    {
        std::cerr<<"outpoints and outvol are the same file. it is impossible to continue"<<std::endl;
        return 1;
        
    }

    std::stringstream temp1;
    int sampleDimensions[3];
    temp1.str(getParameterAsString("resolution"));
    int counterRes=0;
    
    while (!temp1.eof())
    {
        std::string paramField;
        temp1>>sampleDimensions[counterRes]; //set resolution
        if(sampleDimensions[counterRes]<0)
        {
            std::cerr<<"VSMuPortalOp: Invalid resolution is given"<<std::endl;
             return false;
        }
        counterRes++;
        if(counterRes==3)
            break;
    }
    if(counterRes>=1) dimX=sampleDimensions[0];
    if(counterRes>=2) dimY=sampleDimensions[1];
    if(counterRes==3) dimZ=sampleDimensions[2];
    if(isParameterPresent("dimvox")) dimVox=getParameterAsInt("dimvox");
    if(isParameterPresent("trackplanedist")) dMin=getParameterAsInt("trackplanedist");
    if(isParameterPresent("innerdist")) dMax=getParameterAsInt("innerdist");

    float c=0, energy;
	Vector p1, p2, p3, p4, scattPoint;
	double thetaDeg=0, thetaQuadro=0;
	int event=0, index=0, nEv=0, nPt=0, nEv_par=0;

    m_volDim=DIM_VETTORE;
    
    if(!allocateArray())
        return false;
    
    
// inizialize output points VBT  

    VSTable tablePoints;
    tablePoints.setLocator(fileNameOutputPoints);
    
#ifdef VSBIGENDIAN
	std::string endianism="big";
	
#else
	std::string endianism="little";
#endif
    
    tablePoints.setEndiannes(endianism);
    tablePoints.setType("float");
    tablePoints.addCol("id_Ev");
    tablePoints.addCol("X");
    tablePoints.addCol("Y");
    tablePoints.addCol("Z");
    tablePoints.addCol("theta");
    tablePoints.addCol("thetasqr");
    tablePoints.addCol("energy");

    // inizialize output Volume VBT  
    int nVoxelX=dimX/dimVox;
    int nVoxelY=dimY/dimVox;
    int nVoxelZ=dimZ/dimVox;
    VSTable tableVolume;
    tableVolume.setLocator(fileNameOutputVol);
    tableVolume.setIsVolume(true);
    tableVolume.setCellNumber(nVoxelX,nVoxelY,nVoxelZ);
    tableVolume.setCellSize(1,1,1); 
    tableVolume.setEndiannes(endianism);
    tableVolume.setType("float");
    tableVolume.setNumberOfRows(DIM_VETTORE);
    tableVolume.addCol("sumTheta");
    tableVolume.addCol("sumThetaQuadro");
    tableVolume.addCol("ThetaQuadroAvg");
    tableVolume.addCol("sigma");
    tableVolume.addCol("error");
    tableVolume.addCol("nPoints");
    tableVolume.writeHeader();  //overwrite if table exist!
    
    
    
	Voxel *voxelVector=new Voxel[DIM_VETTORE];
	for(int i=0; i<DIM_VETTORE; i++)
        voxelVector[i]=Voxel();
    /*	
     Z per i piani A,B:    (+250 + 249)/2 = +249.5 cm
     Z per i piani C,D:    (+150 + 149)/2 = +149.5 cm
     Z per i piani E,F:    (-149 - 150)/2 = -149.5 cm
     Z per i piani G,H:    (-249 - 250)/2 = -249.5 cm
     */	
    
    //#ifdef RIGGI
	//calcolo delle coordinate z dei piani che tiene conto dello spessore tra i piani fisici contigui (secondo prof. Riggi)
	p1.z=(double)(D_PIANI_TOT-D_PIANI_FISICI)/2;
	p2.z=(double)(D_PIANI_MAX-D_PIANI_FISICI)/2;	
	p3.z=-(double)(D_PIANI_MAX-D_PIANI_FISICI)/2;
	p4.z=-(double)(D_PIANI_TOT-D_PIANI_FISICI)/2;
    //#endif
	/* 
     //calcolo delle coordinate z dei piani che non tiene conto dello spessore tra i piani fisici contigui
     p1.z=(double)((D_PIANI_MAX/2)+D_PIANI_MIN);
     p2.z=(double)(D_PIANI_MAX/2);	
     p3.z=-(double)(D_PIANI_MAX/2);
     p4.z=-(double)((D_PIANI_MAX/2)+D_PIANI_MIN);
     */
    
    /*#ifdef SR
     //calcolo delle coordinate z dei piani che tiene conto delle simulazioni di SR
     p1.z=(double)10-150;
     p2.z=(double)110-150;
     p3.z=(double)430-150;	
     p4.z=(double)530-150;
     
     #endif*/
	/*char *res;
     char buf[200];
     
     
     while(1)
     {
     //printf("lettura\n");		
     res=fgets(buf, 200, f);	
     
     //if(res==NULL) break;
     printf("******fgets:%s ", buf);
     
     }*/
    m_nOfRow=m_tables[0]->getNumberOfRows();
    unsigned long long int fromRow=0;
    unsigned long long int toRow=m_nOfRow-1;
    unsigned int colList[10];
    for(int j=0;j<10;j++) colList[j]=j;
    m_tables[0]->getColumn(colList, m_nOfCol, fromRow, toRow, m_fArray);
    unsigned long long int numOfScattPoints=0;
    for(int k=0;k<m_nOfRow;k++)
    {

        event=m_fArray[0][k];
        p1.x=m_fArray[1][k];
        p1.y=m_fArray[2][k];
        p2.x=m_fArray[3][k];
        p2.y=m_fArray[4][k];
        p3.x=m_fArray[5][k];
        p3.y=m_fArray[6][k];
        p4.x=m_fArray[7][k];
        p4.y=m_fArray[8][k];
        energy=m_fArray[9][k];
        
        
        //TO DO: prima di invocare la Poca applicare una gaussiana ai dati letti (strip)
        
        /*double Line::scatteringAngle(Line *l)
         Vector Line::poca(Line *l)
         
         */
        
        
		//Poca Algorithm
		Line l1=Line(p1, p2);
		Line l2=Line(p3, p4);
		scattPoint=l1.poca(&l2);
		if(true) //se le rette non sono "quasi" parallele: (D<SMALL_NUM)----> esegue sempre!!!
		{
			
			thetaDeg=l1.scatteringAngle(&l2);			
			//thetaDeg=rad2deg(thetaRad);
			thetaQuadro=pow(thetaDeg, 2); //in deg
            
            
			index=scattPoint.Point2Voxel(DIM_X, DIM_Y, DIM_Z, DIM_VOXEL);
			if(index>=0 && index<DIM_VETTORE) //scrivo sul file le coordinate del pto di scattering solo se cade dentro il volume preso in considerazione e il thetaquadro
			{ 		
                m_outPointsfArray[0][numOfScattPoints]=event;
                m_outPointsfArray[1][numOfScattPoints]=scattPoint.x;
                m_outPointsfArray[2][numOfScattPoints]=scattPoint.y;
                m_outPointsfArray[3][numOfScattPoints]=scattPoint.z;
                m_outPointsfArray[4][numOfScattPoints]=thetaDeg;
                m_outPointsfArray[5][numOfScattPoints]=thetaQuadro;
                m_outPointsfArray[6][numOfScattPoints]=energy;

                numOfScattPoints++;
                
                
                /*
                 operazioni sul Vettore dei voxel
                 */
				voxelVector[index].nPts++;
				voxelVector[index].memTheta(thetaDeg);
				voxelVector[index].addTheta(thetaDeg);
				voxelVector[index].addThetaQuadro(thetaQuadro); 
            }
					
        }
        else //se gli eventi sono paralleli---> 
            nEv_par++;			

    }
    tablePoints.setNumberOfRows(numOfScattPoints);
    tablePoints.writeHeader();
    fromRow=0;
    toRow=numOfScattPoints-1;
    tablePoints.putColumn(colList,m_nOfColPoints,fromRow, toRow, m_outPointsfArray);
    for(int i=0; i<DIM_VETTORE; i++)
	{	
        
		voxelVector[i].setThetaM();
		voxelVector[i].setThetaQuadroM();
		voxelVector[i].setSigma();
		voxelVector[i].setError();
        m_outVolfArray[0][i]=voxelVector[i].getSumTheta();
        m_outVolfArray[1][i]=voxelVector[i].getThetaM();
        m_outVolfArray[2][i]=voxelVector[i].getSumThetaQuadro();
        m_outVolfArray[3][i]=voxelVector[i].getThetaQuadroM();
        m_outVolfArray[4][i]=voxelVector[i].getSigma();
        m_outVolfArray[5][i]=voxelVector[i].getError();
        m_outVolfArray[6][i]=voxelVector[i].nPts;
    
    }
    fromRow=0;
    toRow=DIM_VETTORE-1;
    tableVolume.putColumn(colList,m_nOfColVol,fromRow, toRow, m_outVolfArray);
	
	std::cout<<"\n\n***************\n";
	std::cout<<"Results for a disctance between planes (up-down): "<<D_PIANI_TOT<<" cm.\n";
	std::cout<< "Parallel tracks: "<<nEv_par<<".\n";
  	std::cout << "Run completed.\n";
    m_realOutFilename.push_back(fileNameOutputVol);
    m_realOutFilename.push_back(fileNameOutputPoints);
    
	return 0;
}

