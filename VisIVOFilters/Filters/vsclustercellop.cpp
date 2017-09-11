/***************************************************************************
 *   Copyright (C) 2008 by Ugo Becciani   *
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
#include "vsclustercellop.h"
#include "vspointdistributeop.h"
#include "vstable.h"
#include "VisIVOFiltersConfigure.h"

//prova
const unsigned int VSClusterCellOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSClusterCellOp::MIN_NUMBER_OF_ROW = 100;

//---------------------------------------------------------------------
VSClusterCellOp::VSClusterCellOp()
{
 m_fArray=NULL;
 m_grid=NULL;
 m_fresult=NULL;
 m_nOfRow=0;
 m_nOfCol=0;	

}
//---------------------------------------------------------------------


//---------------------------------------------------------------------
VSClusterCellOp::~VSClusterCellOp()
{	if(m_fArray!=NULL)
 	  for(unsigned int i=0;i<m_nOfCol;i++)
	  {
		if(m_fArray[i] != NULL) delete [] m_fArray[i];
	  }
	if(m_fArray!=NULL) delete [] m_fArray;
	if(m_fresult!=NULL) 
	{
		delete [] m_fresult[0];
		delete [] m_fresult;
	}
	if(m_grid!=NULL) 
	{
		delete [] m_grid[0];
		delete [] m_grid;
	}
}
//---------------------------------------------------------------------

//---------------------------------------------------------------------
void VSClusterCellOp::printHelp()
//---------------------------------------------------------------------
{
std::cout<<"Produce a new table or add a new field to the input table.The operation performs the following: 1) It creates a temporary  volume using a field distribution (CIC algorithm) on a regular mesh; 2) It computes, with the same CIC algorithm, the property for each data point, considering the cells where the point is spread on the volume; 3) It save the property in a new table or add the field to the original input table."<<std::endl<<std::endl;
std::cout<<"Usage: VisIVOFilters --op pointproperty  --resolution x_res y_res z_res --points x_col y_col z_col [--field column_name] [--constant value] [--append] [--out filename_out.bin] [--outcol col_name] [--periodic] [--help] [--file] inputFile.bin"<<std::endl<<std::endl;

std::cout<<"Example: VisIVOFilters --op pointproperty --resolution 16 16 16 --points X Y Z --field Mass  --append --outcol distribute --file inputFile.bin"<<std::endl;

std::cout<<"Note:  "<<std::endl;
std::cout<<"--resolution  3D mesh size."<<std::endl;
std::cout<<"--points Columns to be assumed for points coordinates."<<std::endl;
std::cout<<"--field Valid columns name list to be distributed in the grid."<<std::endl;
std::cout<<"--constant Assign a constant to all points to be distributed in the grid. Ignored if field option is given. Default value is a 1.0 for all points."<<std::endl;
std::cout<<"--append. No new table will be cretaed. The original table will have the new field. "<<std::endl;
std::cout<<"--out Name of the new table. Default name is given. Ignored if --append is specified."<<std::endl;
std::cout<<"--periodic. It specifies the box is periodic. Particles outside the box limits are considered inside on the other side."<<std::endl;
std::cout<<"--file Input table filename."<<std::endl;
std::cout<<"--outcol. Column name of the new field"<<std::endl;

std::cout<<"--help produce this output."<<std::endl;

return;

}

//---------------------------------------------------------------------
bool VSClusterCellOp::allocateArray()
//---------------------------------------------------------------------
{
	m_fresult=new float*[1];
	m_grid=new float*[1];
	m_fArray=new  float*[3];
	bool goodAllocation=false;
	while(!goodAllocation)
	{
		goodAllocation=true;
		for(unsigned int i=0;i<m_nOfCol;i++)
		{
try
{
			m_fArray[i]=new  float[m_nOfRow];
}
catch(std::bad_alloc &e)
{
	m_fArray[i]=NULL;
}

			if(m_fArray[i]==NULL) 
			{	
				goodAllocation=false;
				for(unsigned int j=0;j<i;j++) 
					delete [] m_fArray[j];
				if(m_nOfRow==MIN_NUMBER_OF_ROW)
				{ 
					delete [] m_fArray;
					m_fArray=NULL;
					return false;
				}
				m_nOfRow=m_nOfRow-MAX_NUMBER_TO_REDUCE_ROW;
				if(m_nOfRow<=MAX_NUMBER_TO_REDUCE_ROW) m_nOfRow=MIN_NUMBER_OF_ROW;
				break;
			}

		}
try
{
		m_fresult[0]=new float[m_nOfRow];
}	
catch(std::bad_alloc &e)
{
	m_fresult[0]=NULL;
}

		if(m_fresult==NULL) 
		{
			goodAllocation=false;
			for(unsigned int j=0;j<m_nOfCol;j++) 
				delete [] m_fArray[j];
			if(m_nOfRow==MIN_NUMBER_OF_ROW)
			{ 
				delete [] m_fArray;
				m_fArray=NULL;
				return false;
			}
			m_nOfRow=m_nOfRow-MAX_NUMBER_TO_REDUCE_ROW;
			if(m_nOfRow<=MAX_NUMBER_TO_REDUCE_ROW) 
				m_nOfRow=MIN_NUMBER_OF_ROW;
		}

		if(goodAllocation)
		{
try		
{
		  m_grid[0]=new  float[m_numNewPts];
}
catch(std::bad_alloc &e)
{
	m_grid=NULL;
}

		  if(m_grid[0]==NULL) 
		  {	
			goodAllocation=false;
			for(unsigned int j=0;j<m_nOfCol;j++) 
				delete [] m_fArray[j];
			delete [] m_fresult[0];

			if(m_numNewPts==MIN_NUMBER_OF_ROW)
			{ 
				delete [] m_fArray;
				m_fArray=NULL;
				delete [] m_grid;
				m_grid=NULL;
				delete [] m_fresult;
				m_fresult=NULL;
				return false;
			}
			m_nOfRow=m_nOfRow-MAX_NUMBER_TO_REDUCE_ROW;
			if(m_nOfRow<=MAX_NUMBER_TO_REDUCE_ROW) 
				m_nOfRow=MIN_NUMBER_OF_ROW;
			m_numNewPts=m_numNewPts-MAX_NUMBER_TO_REDUCE_ROW;
			if(m_numNewPts<=MAX_NUMBER_TO_REDUCE_ROW)
				m_numNewPts=MIN_NUMBER_OF_ROW;
		  }
		}
	}


return true;
}


//---------------------------------------------------------------------
bool VSClusterCellOp::execute()
//---------------------------------------------------------------------
{
//cannot be applied to volume tables
std::stringstream sstmp1;
std::string  randat;
time_t rawtime;
struct tm * timeinfo;
char buffer [80];
time ( &rawtime );
timeinfo = localtime ( &rawtime );
strftime (buffer,80,"%Y%m%d%H%M%S",timeinfo);
sstmp1<<"_"<<rand()<<buffer;  
randat=sstmp1.str();

bool periodic=false;
if(m_tables[0]->getIsVolume())
{
	std::cerr<<"VSClusterCellOp: cannot be applied to  volume."<<std::endl;
	return false;
}
std::string tempFilename;

/*** Execute the PointDistribute  OP **/
{
VSTable table(m_tables[0]->getLocator());

VSPointDistributeOp op;

std::string stmp;
stmp=getParameterAsString("resolution");
if(stmp==""||stmp=="unknown")
{
	std::cerr<<"VSClusterCellOp: Invalid resolution is given"<<std::endl;
	return false;
	
}
op.addParameter("resolution",stmp);
op.addParameter("cic",""); //Future: other options activable.
if(isParameterPresent("periodic"))
{
	periodic=true;
	op.addParameter("periodic",""); //Future: other options activable.
}
stmp.clear();
stmp=getParameterAsString("points");
if(stmp==""||stmp=="unknown")
{
	std::cerr<<"VSClusterCellOp: Invalid points is given"<<std::endl;
	return false;
	
}
op.addParameter("points",stmp);
stmp.clear();
stmp=getParameterAsString("field");
if(stmp!="" && stmp!="unknown")
	op.addParameter("field",stmp);
	
stmp.clear();
stmp=getParameterAsString("constant");
if(stmp!="" && stmp!="unknown")
	op.addParameter("constant",stmp);
	
stmp.clear();
stmp=getParameterAsString("out");
if(stmp=="" || stmp=="unknown")
{
		tempFilename="./_tempPDOp.bin"+randat;
}else
{
	tempFilename=stmp;
	tempFilename.append("_tempPDOp.bin"+randat);
}
op.addParameter("out",tempFilename);
op.addInput(&table);
op.execute();
float *ftmp;
ftmp=new float[3];
if(!op.getOrigin(ftmp))
{	
	std::cerr<<"VSClusterCellOp: Invalid Origin is given"<<std::endl;
	remove(tempFilename.c_str());
	tempFilename.append(".head");
	remove(tempFilename.c_str());

	return false;
}
for(int i=0;i<3;i++) m_origin[i]=ftmp[i];
if(!op.getSpacing(ftmp))
{	
	std::cerr<<"VSClusterCellOp: Invalid Origin is given"<<std::endl;
	remove(tempFilename.c_str());
	tempFilename.append(".head");
	remove(tempFilename.c_str());

	return false;
}
for(int i=0;i<3;i++) m_spacing[i]=ftmp[i];
} //block to destroy table and op
/*** END PointDistribute OP **/

bool append=true; 
if(!isParameterPresent("append")) append=false;
std::string fileNameOutput;
if(!append)
{
	fileNameOutput=getParameterAsString("out"); 
	if(fileNameOutput==""||fileNameOutput=="unknown")
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
  		fileNameOutputSStream<<filenameInpuTable.substr(0, len-4)<<"_pointpropertyop_"<<buffer<<".bin";
		fileNameOutput=fileNameOutputSStream.str();
	}
  	if(fileNameOutput.find(".bin") == std::string::npos)
	    		fileNameOutput.append(".bin");
} else
	fileNameOutput=m_tables[0]->getLocator();
	m_realOutFilename.push_back(fileNameOutput);

VSTable tableProperty;
tableProperty.setLocator(fileNameOutput);

#ifdef VSBIGENDIAN
	std::string endianism="big";
	
#else	
	std::string endianism="little";
#endif

tableProperty.setEndiannes(endianism);
tableProperty.setType("float");
tableProperty.setNumberOfRows(m_tables[0]->getNumberOfRows());
std::string colNameOutput;

if(getParameterAsString("outcol").empty() ||getParameterAsString("outcol")=="unknown")
	colNameOutput="__PointProperty_";
else
	colNameOutput=getParameterAsString("outcol");
int numOfCols=0;
if(append)
{
	for(unsigned int i=0;i<m_tables[0]->getNumberOfColumns();i++) 
		tableProperty.addCol(m_tables[0]->getColName(i));	
	numOfCols=m_tables[0]->getNumberOfColumns();
}
if(!tableProperty.addCol(colNameOutput))
{
	std::cerr<<"Error: Invalid or duplicate column name in existing table: "<<colNameOutput<<std::endl;
	remove(tempFilename.c_str());
	tempFilename.append(".head");
	remove(tempFilename.c_str());

	return false;
}	
numOfCols++;
tableProperty.writeHeader();  //overwrite if table exist!

unsigned long long int totEle;
unsigned long long int  fromRow,toRow,startCounter;
startCounter=0;
unsigned int colList[3],gridColList[1];
gridColList[0]=0;

int counterCols=0;
std::stringstream ssListparameters;
ssListparameters.str(getParameterAsString("points"));
while (!ssListparameters.eof())
{
	std::string paramField;
	ssListparameters>>paramField;
	if(m_tables[0] -> getColId(paramField)>=0)
	{
		colList[counterCols]=m_tables[0] -> getColId(paramField); 
		counterCols++;
		if(counterCols==3)
			break;
	}
}
if(counterCols !=3)
{
	std::cerr<<"VSPointDistributeOp: Invalid columns in --points argument is given"<<std::endl;
	remove(tempFilename.c_str());
	tempFilename.append(".head");
	remove(tempFilename.c_str());
	return false;
}

VSTable tableGrid(tempFilename);
m_numNewPts=tableGrid.getNumberOfRows();

int maxInt=getMaxNumberInt();
unsigned int maxEle;
unsigned long long int totRows=m_tables[0]->getNumberOfRows();
unsigned long long int nOfEle=totRows;

if(totRows>maxInt)
	maxEle=maxInt; 
else
	maxEle=totRows;

m_nOfRow=maxEle;
if(m_numNewPts>maxEle) m_numNewPts=maxEle;
m_nOfCol=3;
bool allocationArray=allocateArray();
maxEle=m_nOfRow;

if(m_fArray==NULL ||m_grid==NULL ||  !allocationArray )	
{
	std::cerr<<"Failed Array allocation. Mathematical Operation terminated"<<std::endl;
	remove(tempFilename.c_str());
	tempFilename.append(".head");
	remove(tempFilename.c_str());
	return false;
}
// load the grid
unsigned long long int gridIndex[2];  //limits of cashed grid
gridIndex[0]=0;
gridIndex[1]=m_numNewPts-1;
tableGrid.getColumn(gridColList,1,0,m_numNewPts-1,m_grid);

std::stringstream ssresolution;
int sampleDimensions[3];
ssresolution.str(getParameterAsString("resolution"));
int counterRes=0;

while (!ssresolution.eof())
{
	std::string paramField;
	ssresolution>>sampleDimensions[counterRes]; //set resolution
	if(sampleDimensions[counterRes]<0)
	{
		std::cerr<<"VSPointDistributeOp: Invalid resolution is given"<<std::endl;
		remove(tempFilename.c_str());
		tempFilename.append(".head");
		remove(tempFilename.c_str());
		return false;
	}
	counterRes++;
	if(counterRes==3)
		break;
}
int jkFactor = sampleDimensions[0]*sampleDimensions[1];
int jFactor = sampleDimensions[0];

totEle=totRows;
while(totEle!=0)
{
	fromRow=startCounter;
	toRow=fromRow+maxEle-1;
	if(toRow>totRows-1)toRow=totRows-1;
	m_tables[0]->getColumn(colList, m_nOfCol, fromRow, toRow, m_fArray);
 	int nCell=sampleDimensions[0]*sampleDimensions[1]*sampleDimensions[2];

	for (int ptId=0; ptId < toRow-fromRow+1; ptId++)
    	{
		m_fresult[0][ptId]=0.0000001; //avoid zeros for log scales 
      		float px[3];
		for (int j=0; j < 3; j++)
			px[j]=m_fArray[j][ptId];
 
    		float pos1 = (float) (px[0] - m_origin[0]) / m_spacing[0];
    		float pos2 = (float) (px[1] - m_origin[1]) / m_spacing[1];
    		float pos3 = (float) (px[2] - m_origin[2]) / m_spacing[2];
     		int i1 = floor(pos1);
/*       		if(fabs((pos1-i1))>0.5)
			i1++;*/
     		int i2 = floor(pos2);
/*       		if(fabs((pos2-i2))>0.5)
			i2++;*/
     		int i3 = floor(pos3);
/*       		if(fabs((pos3-i3))>0.5)
			i3++;*/
    		int i11=i1+1;
    		int i21=i2+1;
    		int i31=i3+1;
               	float dd1=pos1-(float)i1; 
               	float dd2=pos2-(float)i2; 
               	float dd3=pos3-(float)i3;
               	float de1=1.0-dd1;
               	float de2=1.0-dd2;
               	float de3=1.0-dd3;
		if(periodic)
		{
			while(i1<0)
  				i1=sampleDimensions[0]+i1;
			while(i2<0)
  				i2=sampleDimensions[1]+i2;
			while(i3<0)
  				i3=sampleDimensions[2]+i3;
			while(i1>sampleDimensions[0]-1)
  				i1=i1-sampleDimensions[0];
			while(i2>sampleDimensions[1]-1)
  				i2=i2-sampleDimensions[1];
			while(i3>sampleDimensions[2]-1)
  				i3=i3-sampleDimensions[2];
			while(i11<0)
  				i11=sampleDimensions[0]+i11;
			while(i21<0)
  				i21=sampleDimensions[1]+i21;
			while(i31<0)
  				i31=sampleDimensions[2]+i31;
			while(i11>sampleDimensions[0]-1)
  				i11=i11-sampleDimensions[0];
			while(i21>sampleDimensions[1]-1)
  				i21=i21-sampleDimensions[1];
			while(i31>sampleDimensions[2]-1)
  				i31=i31-sampleDimensions[2];
    		}

// calculate weights
			float d[8];
               		d[0]= de1*de2*de3;
               		d[1]= dd1*de2*de3;
               		d[2]= de1*dd2*de3;
               		d[3]= dd1*dd2*de3;
               		d[4]= de1*de2*dd3;
               		d[5]= dd1*de2*dd3;
               		d[6]= de1*dd2*dd3;
               		d[7]= dd1*dd2*dd3;

// linearize coordinates
			unsigned long long int ind[8];
               		ind[0] = i3*jkFactor + i2*jFactor + i1;
               		ind[1] = i3*jkFactor + i2*jFactor + i11;
               		ind[2] = i3*jkFactor + i21*jFactor + i1;
               		ind[3] = i3*jkFactor + i21*jFactor + i11;
               		ind[4] = i31*jkFactor + i2*jFactor + i1;
               		ind[5] = i31*jkFactor + i2*jFactor + i11;
               		ind[6] = i31*jkFactor + i21*jFactor + i1;
               		ind[7] = i31*jkFactor + i21*jFactor + i11;

		  for(int n=0;n<8;n++)
		  { 
		     if(ind[n]<0  || ind[n]>nCell-1)
				continue;
		   if(ind[n]<gridIndex[0] || ind[n]>gridIndex[1]) //ind1 NOT in cache!
		   {
			gridIndex[0]=ind[n];
			gridIndex[1]=gridIndex[0]+m_numNewPts-1;
			if(gridIndex[1]>tableGrid.getNumberOfRows())gridIndex[1]=tableGrid.getNumberOfRows()-1;
 			tableGrid.getColumn(gridColList,1,gridIndex[0],gridIndex[1],m_grid); 	
		   }
//		   int in1=ind[n]-gridIndex[0];
//		   float grv=m_grid[0][ind[n]-gridIndex[0]];
		   m_fresult[0][ptId]=m_fresult[0][ptId]+m_grid[0][ind[n]-gridIndex[0]]*d[n];
		  }
		
	}
	unsigned int resultColList[1];
	resultColList[0]=numOfCols-1;
	tableProperty.putColumn(resultColList, 1, fromRow, toRow, m_fresult);	
	startCounter=toRow+1;
	totEle=totEle-(toRow-fromRow+1);
	if(totEle<0) totEle=0;
}
remove(tempFilename.c_str());
tempFilename.append(".head");
remove(tempFilename.c_str());

return true;
}
