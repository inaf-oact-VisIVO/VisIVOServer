//
// C++ Implementation: vsgrid2pointdistr
//
// Description: 
//
//
// Author:  <>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <cmath>
#ifdef WIN32
	#include <time.h>
#endif
#include "vstable.h"
#include "vsstatisticop.h"


#include "vsgrid2pointdistr.h"

const unsigned int VSGrid2PointDistr::MIN_NUMBER_OF_ROW = 10000;
const unsigned int VSGrid2PointDistr::MAX_NUMBER_TO_REDUCE_ROW = 10000;

//---------------------------------------------------------------------
VSGrid2PointDistr::VSGrid2PointDistr()
{	
m_gridCellNumber=new unsigned int [3];
m_fArray=NULL;
m_grid=NULL;
m_nOfCol= 0;
m_nOfRow= 0;
m_nOfEle= 0;
  m_tsc=false;
  m_cic=true;
  m_ngp=false;
for(int i=0;i<3;i++)
{	
	m_spacing[i]=1.0;
	m_origin[i]=0.0;
//	m_gridCellNumber[i]=0;
}
return;
}
//---------------------------------------------------------------------


//---------------------------------------------------------------------
VSGrid2PointDistr::~VSGrid2PointDistr()
{	if(m_fArray!=NULL)
 	  for(unsigned int i=0;i<m_nOfCol;i++)
	  {
		if(m_fArray[i] != NULL) delete [] m_fArray[i];
	  }
	if(m_fArray!=NULL) delete [] m_fArray;
	if(m_grid!=NULL) 
	{
		delete [] m_grid[0];
		delete [] m_grid;
	}
//	delete [] m_gridCellNumber;
}
//---------------------------------------------------------------------
//---------------------------------------------------------------------
void VSGrid2PointDistr::printHelp()
{
std::cout<<"This filter produces a new table or add a new field to the input table.The operation performs the following: 1) It load a  volume (input volume data table) and a table with a  point distribution in the volume 2) It computes, using the CIC (default) or NGP or TSC algorithm, a value (density assumed) for each data point, considering the  cells value where the point is spread on the volume (cellValue*w_factor*cellVolume) 3) It save the property in a new table or add the field to the original input table."<<std::endl<<std::endl;
std::cout<<"Usage: VisIVOFilters --op grid2point   --points x_col y_col z_col [--field column_name] [--density] [--append] [--out filename_out.bin] [--outcol col_name] [--tsc] [--ngp] --volume inputVolmeData.bin [--gridOrigin xg0 xg1 xg2] [--gridSpacing sg0 sg1 sg2] [--box length] [--periodic] [--help] [--file] inputFile.bin"<<std::endl<<std::endl;

std::cout<<"Example: VisIVOFilters --op grid2point  --points X Y Z --field Mass  --append --outcol distribute --volume inputVolmeData.bin --file inputFile.bin"<<std::endl;

std::cout<<"Note:  "<<std::endl;
std::cout<<"--points Columns to be assumed for points coordinates."<<std::endl;
std::cout<<"--field Valid Volume Column Name. Default value is the first column name"<<std::endl;
std::cout<<"--density. CellVolume is not considered (cellVolume=1)"<<std::endl;
std::cout<<"--append. No new table will be created. The original table will have the new field. "<<std::endl;
std::cout<<"--out Name of the new table. Default name is given. Ignored if --append is specified."<<std::endl;
std::cout<<"--outcol. Column name of the new field"<<std::endl;
std::cout<<"--tsc. The TSC algorithm is adopted."<<std::endl;
std::cout<<"--ngp. The NGP algorithm is adopted."<<std::endl;
std::cout<<"--volume. Input data volume filename (a VisIVO Binary Table)."<<std::endl;
std::cout<<"--gridOrigin. It specifies the coordinate of the lower left corner of the grid. Default vaules are assumed from the box of inputFile.bin"<<std::endl;
std::cout<<"--gridSpacing. It specifies the length of each cell dimension in arbitray unit. This parameter is ignored if the box option is given. Default values are assumed from the box of inputFile.bin"<<std::endl;
std::cout<<"--box. It specifies the length of a box. Default value is assumed from the box of inputFile.bin if the gridSpacing option is not given"<<std::endl;
std::cout<<"--periodic. It specifies the grid is periodic. Particles outside the grid limits are considered inside on the other side."<<std::endl;

std::cout<<"--file Input table filename with point distribution."<<std::endl;

std::cout<<"--help produce this output."<<std::endl;
std::cout<<"WARNING: Cell geometry is considered only to compute the cellVolume value in this operation."<<std::endl;

return;


}
//---------------------------------------------------------------------
//---------------------------------------------------------------------
bool VSGrid2PointDistr::allocateArray()
{	
unsigned long long int tempLL=getMaxNumberInt();
if(((unsigned long long int)m_nOfEle*m_nOfCol)>tempLL) 
	m_nOfEle=(int)tempLL/m_nOfCol;

m_grid=new float*[1];
try {
	m_fArray= new float*[m_nOfCol];
}
catch(std::bad_alloc &e) {
	m_fArray=NULL;
}
if(m_fArray==NULL) 
	return false;

bool goodAllocation=false;
while(!goodAllocation){
	goodAllocation=true;

	for(unsigned int i=0;i<m_nOfCol;i++){
try{
		m_fArray[i]= new float[m_nOfEle];
}
catch(std::bad_alloc &e){
		m_fArray[i]= NULL;
}
		if(m_fArray[i]==NULL) {	
			goodAllocation=false;
			//DELETE MEMORY PREVIOUSLY ALLOCATED
			for(unsigned int j=0;j<i;j++) 
				delete [] m_fArray[j];
			if(m_nOfEle==MIN_NUMBER_OF_ROW)
			{ 
				delete [] m_fArray;
				m_fArray=NULL;
				return false;
			}//close if 
			m_nOfEle= m_nOfEle-MAX_NUMBER_TO_REDUCE_ROW;
			if(m_nOfEle<=MIN_NUMBER_OF_ROW)  m_nOfEle=MIN_NUMBER_OF_ROW;
			break;
		}//close if
	}//close for i
	if(goodAllocation)
	{
try		
{
	  	m_grid[0]=new  float[m_nOfGridEle];
}
catch(std::bad_alloc &e)
{
		m_grid[0]=NULL;
}
		  if(m_grid[0]==NULL) 
		  {	
			goodAllocation=false;
			for(unsigned int j=0;j<m_nOfCol;j++) 
				delete [] m_fArray[j];
			if(m_nOfGridEle==MIN_NUMBER_OF_ROW || m_nOfEle==MIN_NUMBER_OF_ROW)
			{ 
				delete [] m_fArray;
				m_fArray=NULL;
				delete [] m_grid;
				m_grid=NULL;
				return false;
			}
			m_nOfEle=m_nOfEle-MAX_NUMBER_TO_REDUCE_ROW;
			if(m_nOfEle<=MAX_NUMBER_TO_REDUCE_ROW)
				m_nOfEle=MIN_NUMBER_OF_ROW;
			m_nOfGridEle=m_nOfGridEle-MAX_NUMBER_TO_REDUCE_ROW;
			if(m_nOfGridEle<=MAX_NUMBER_TO_REDUCE_ROW) 
				m_nOfGridEle=MIN_NUMBER_OF_ROW;
		  }
	}//if good allocation
}//close while
	
for(int i=0;i<m_nOfEle;i++)
	m_fArray[m_nOfCol-1][i]=0.0;	
return true;
		
}
//---------------------------------------------------------------------
bool VSGrid2PointDistr::setOrigin(float x0, float y0,float z0)
{
	m_origin[0]=x0;
	m_origin[1]=y0;
	m_origin[2]=z0;
	return true;	
}
//---------------------------------------------------------------------
bool VSGrid2PointDistr::setSpacing(float xs, float ys,float zs)
{
	if(xs<=0. || ys <=0. || zs<=0.)
	{
		std::cerr<<"VSGrdid2Point. Invalid grid spacing."<<std::endl;
		return false;
	}
	m_spacing[0]=xs;
	m_spacing[1]=ys;
	m_spacing[2]=zs;
	return true;	
}
//---------------------------------------------------------------------
//---------------------------------------------------------------------
bool VSGrid2PointDistr::setOrigin()
{

if(getParameterAsString("gridOrigin").empty()||getParameterAsString("gridOrigin")=="unknown")
{	
	std::cerr<<"No valid input gridOrigin parameter."<<std::endl;
	return false;
}
std::stringstream ssOrigin; 
ssOrigin.str(getParameterAsString("gridOrigin"));
int count=0;
while (!ssOrigin.eof())
{
	ssOrigin>>m_origin[count];
	if(count==2)
		break;
	count++;
}

return true;	
}
//---------------------------------------------------------------------
bool VSGrid2PointDistr::setSpacing()
{
if(isParameterPresent("box"))
{
	float box=getParameterAsFloat("box");
	if(box<=0.0)	
	{	
		std::cerr<<"No valid input box parameter."<<std::endl;
		return false;
	}
	for(int i=0;i<3;i++)
	{
		if(m_gridCellNumber[i]<=0)
		{
			std::cerr<<"No valid grid cell number."<<std::endl;
			return false;
		}
		m_spacing[i]=box/m_gridCellNumber[i];
	}
	return true;
}
if(isParameterPresent("gridSpacing"))
{
	if(getParameterAsString("gridSpacing").empty()||getParameterAsString("gridSpacing")=="unknown")
	{	
		std::cerr<<"No valid input gridSpacing parameter. Operation aborted"<<std::endl;
		return false;
	}
	std::stringstream ssSpacing; 
	ssSpacing.str(getParameterAsString("gridSpacing"));
	int count=0;
	while (!ssSpacing.eof())
	{
		ssSpacing>>m_spacing[count];
		if(m_spacing[count] <=0.)
		{		
			std::cerr<<"No valid input gridSpacing parameter."<<std::endl;
			return false;
		}
		if(count==2)
			break;
		count++;
	}
}
return true;	
}
//---------------------------------------------------------------------

//---------------------------------------------------------------------
bool VSGrid2PointDistr::execute()
{
bool originGiven=false;
bool spacingGiven=false;
bool periodic=false;
if(isParameterPresent("periodic"))
	periodic=true;

if(!isParameterPresent("volume") || getParameterAsString("volume").empty() || getParameterAsString("volume")=="unknown" )
{
	std::cerr<<"VSGrid2Point: invalid  volume option."<<std::endl;
	return false;
}
std::string volumeFilename=getParameterAsString("volume");
if(volumeFilename.find(".bin") == std::string::npos) 
	volumeFilename.append(".bin");

VSTable VolumeTable(volumeFilename);

const float * gridCellSize;
//m_gridCellNumber=new unsigned int [3];
gridCellSize=new float [3];
m_gridCellNumber=VolumeTable.getCellNumber();
for(int i=0;i<3;i++)
	if(m_gridCellNumber[i]<=0)
	{
		std::cerr<<"No valid grid cell number."<<std::endl;
		return false;
	}

m_sampleDimensions[0]=m_gridCellNumber[0];
m_sampleDimensions[1]=m_gridCellNumber[1];
m_sampleDimensions[2]=m_gridCellNumber[2];
gridCellSize=VolumeTable.getCellSize();
if(!VolumeTable.tableExist()){
	std::cerr<<"Not existing Volume table is given"<<std::endl;
	return false;
}
if(!VolumeTable.getIsVolume()){
	std::cerr<<"No valid input Volume table is given"<<std::endl;
	return false;
}

if(isParameterPresent("gridOrigin"))
	originGiven=setOrigin();

if(isParameterPresent("box"))
	spacingGiven=setSpacing();

if(isParameterPresent("gridSpacing") && !spacingGiven)
{
	spacingGiven=setSpacing();
	if(!spacingGiven)
		return false;
}
   if(isParameterPresent("tsc"))
   {
	m_tsc=true;
	m_ngp=false;
	m_cic=false;
   }
   if(isParameterPresent("ngp"))
   {
	if(m_tsc)
		std::cerr<<"Ignored --tsc parameter"<<std::endl;
	m_tsc=false;
	m_cic=false;
	m_ngp=true;
   }


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
  		fileNameOutputSStream<<filenameInpuTable.substr(0, len-4)<<"_vsgrid2point_"<<buffer<<".bin";
		fileNameOutput=fileNameOutputSStream.str();
	}
  	if(fileNameOutput.find(".bin") == std::string::npos)
	    		fileNameOutput.append(".bin");
} else
	fileNameOutput=m_tables[0]->getLocator();

m_realOutFilename.push_back(fileNameOutput);


VSTable tableGrid2Point;
tableGrid2Point.setLocator(fileNameOutput);

#ifdef VSBIGENDIAN
	std::string endianism="big";
	
#else	
	std::string endianism="little";
#endif

tableGrid2Point.setEndiannes(endianism);
tableGrid2Point.setType("float");
tableGrid2Point.setNumberOfRows(m_tables[0]->getNumberOfRows());
std::string colNameOutput;

if(getParameterAsString("outcol").empty() ||getParameterAsString("outcol")=="unknown")
	colNameOutput="__Grid2Point_";
else
	colNameOutput=getParameterAsString("outcol");
int numOfCols=0;
if(append)
{
	for(unsigned int i=0;i<m_tables[0]->getNumberOfColumns();i++) 
		tableGrid2Point.addCol(m_tables[0]->getColName(i));	
	numOfCols=m_tables[0]->getNumberOfColumns();
}
if(!tableGrid2Point.addCol(colNameOutput))
{
	std::cerr<<"Error: Invalid or duplicate column name in existing table: "<<colNameOutput<<std::endl;
	return false;
}	
numOfCols++;
tableGrid2Point.writeHeader();  //overwrite if table exist!

unsigned long long int totEle;
unsigned long long int  fromRow,toRow,startCounter;
startCounter=0;
unsigned int colList[3];
int gridFieldList=0;
if(isParameterPresent("field"))
{
	if(getParameterAsString("field").empty() || getParameterAsString("field")=="unknown")
		gridFieldList=0;
	else
	{
		gridFieldList=VolumeTable.getColId( getParameterAsString("field"));
		if(gridFieldList==-1)
			gridFieldList=0;
	}
}

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
	std::cerr<<"VSGrid2Point: Invalid columns in --points argument is given"<<std::endl;
	return false;
}
m_nOfCol=counterCols+1;
m_nOfRow= m_tables[0]->getNumberOfRows();
m_nOfGridRow=VolumeTable.getNumberOfRows();
if(m_nOfRow>getMaxNumberInt()) 
	m_nOfEle= getMaxNumberInt();
else 
	m_nOfEle= m_nOfRow;
if(m_nOfGridRow>getMaxNumberInt()) 
	m_nOfGridEle= getMaxNumberInt();
else 
	m_nOfGridEle= m_nOfGridRow;
bool allocationArray=allocateArray();
if(!allocationArray)
{	
	std::cerr<<"Failed Array allocation. VSGrid2Point aborted"<<std::endl;
	return false;

}
/*** Execute the Statistic  OP **/
if(!spacingGiven || ! originGiven)
{
VSStatisticOp op1;
op1.addParameter("list",getParameterAsString("points"));
op1.addParameter("silent","");
op1.addInput(m_tables[0]);
op1.execute();
float max,min,avg;
for (unsigned int i=0;i<m_nOfCol-1;i++)
{
	float dummy;
	op1.getRange(i,max,min,avg,dummy);
	if(!originGiven) m_origin[i]=min;
	if(!spacingGiven) m_spacing[i]=(max-min)/m_gridCellNumber[i];
}	
}
/******************/
/***** Load Grid **/
unsigned long long int gridIndex[2];  //limits of cashed grid
gridIndex[0]=0;
gridIndex[1]=m_nOfGridEle-1;
unsigned int *gridColList;
gridColList=new unsigned int[1];
gridColList[0]=gridFieldList;
VolumeTable.getColumn(gridColList,1,0,m_nOfGridEle-1,m_grid);
/*******************/
unsigned long long int index=0;
float cellVolume=m_spacing[0]*m_spacing[1]*m_spacing[2];
if(isParameterPresent("density"))
	cellVolume=1.0;
while(index< m_nOfRow)
{	
	m_tables[0]->getColumn(colList,m_nOfCol-1,index,index+m_nOfEle-1,m_fArray);

        if(m_ngp)
	{
	  int jkFactor = m_gridCellNumber[0]*m_gridCellNumber[1];
          int jFactor = m_gridCellNumber[0];
	  int ind;
// NGP on points
  	 for (int ptId=0; ptId < m_nOfEle; ptId++)
    	 {
		m_fArray[m_nOfCol-1][ptId]=0.0000001;      		float px[3];
		px[0]=m_fArray[0][ptId];
		px[1]=m_fArray[1][ptId];
		px[2]=m_fArray[2][ptId];
 
    		float pos1 = (float) (px[0] - m_origin[0]) / m_spacing[0];
    		float pos2 = (float) (px[1] - m_origin[1]) / m_spacing[1];
    		float pos3 = (float) (px[2] - m_origin[2]) / m_spacing[2];
   
// find nearest grid points    
     		int i1 = floor(pos1);
       		if(fabs((pos1-i1))>0.5)
			i1++;
     		int i2 = floor(pos2);
       		if(fabs((pos2-i2))>0.5)
			i2++;
     		int i3 = floor(pos3);
       		if(fabs((pos3-i3))>0.5)
			i3++;
		if(periodic)
		{
			
			if(i1<0)
  				i1=m_sampleDimensions[0]+i1;
			if(i2<0)
  				i2=m_sampleDimensions[1]+i2;
			if(i3<0)
  				i3=m_sampleDimensions[2]+i3;
			if(i1>m_sampleDimensions[0]-1)
  				i1=i1-m_sampleDimensions[0];
			if(i2>m_sampleDimensions[1]-1)
  				i2=i2-m_sampleDimensions[1];
			if(i3>m_sampleDimensions[2]-1)
  				i3=i3-m_sampleDimensions[2];
		}
          	if (i1 <m_sampleDimensions[0] && i1 >= 0 && i2 < m_sampleDimensions[1] && i2 >= 0 && i3 < m_sampleDimensions[2] && i3 >= 0) 
                {
               		ind = i3*jkFactor + i2*jFactor + i1;


   
// calculate density

			 if(ind<gridIndex[0] || ind>gridIndex[1]) //ind1 NOT in cache!
			 {
				gridIndex[0]=ind;
				gridIndex[1]=gridIndex[0]+m_nOfGridEle-1;
				if(gridIndex[1]>VolumeTable.getNumberOfRows())gridIndex[1]=VolumeTable.getNumberOfRows()-1;
 				VolumeTable.getColumn(gridColList,1,gridIndex[0],gridIndex[1],m_grid);
			  }
		   	  m_fArray[m_nOfCol-1][ptId]+=m_grid[0][ind-gridIndex[0]]*cellVolume;
	    	} //if(i1..)
	  	m_fArray[0][ptId]=m_fArray[m_nOfCol-1][ptId]; //prepare array for write
         } //for(... ptId..)


	}
	if(m_cic)
	{
	  int jkFactor = m_gridCellNumber[0]*m_gridCellNumber[1];
          int jFactor = m_gridCellNumber[0];
 	   int nCell=m_sampleDimensions[0]*m_sampleDimensions[1]*m_sampleDimensions[2];
	   for (int ptId=0; ptId < m_nOfEle; ptId++)
    	   {
		m_fArray[m_nOfCol-1][ptId]=0.0000001; //avoid zeros for log scales 
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
			if(i1<0)
  				i1=m_sampleDimensions[0]+i1;
			if(i2<0)
  				i2=m_sampleDimensions[1]+i2;
			if(i3<0)
  				i3=m_sampleDimensions[2]+i3;
			if(i1>m_sampleDimensions[0]-1)
  				i1=i1-m_sampleDimensions[0];
			if(i2>m_sampleDimensions[1]-1)
  				i2=i2-m_sampleDimensions[1];
			if(i3>m_sampleDimensions[2]-1)
  				i3=i3-m_sampleDimensions[2];
			if(i11<0)
  				i11=m_sampleDimensions[0]+i11;
			if(i21<0)
  				i21=m_sampleDimensions[1]+i21;
			if(i31<0)
  				i31=m_sampleDimensions[2]+i31;
			if(i11>m_sampleDimensions[0]-1)
  				i11=i11-m_sampleDimensions[0];
			if(i21>m_sampleDimensions[1]-1)
  				i21=i21-m_sampleDimensions[1];
			if(i31>m_sampleDimensions[2]-1)
  				i31=i31-m_sampleDimensions[2];
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
		     if(ind[n]<0  || ind[n]>nCell)
				continue;
		     if(ind[n]<gridIndex[0] || ind[n]>gridIndex[1]) //ind1 NOT in cache!
		     {
			gridIndex[0]=ind[n];
			gridIndex[1]=gridIndex[0]+m_nOfGridEle-1;
			if(gridIndex[1]>VolumeTable.getNumberOfRows())gridIndex[1]=VolumeTable.getNumberOfRows()-1;
 			VolumeTable.getColumn(gridColList,1,gridIndex[0],gridIndex[1],m_grid); 	
		   	}
//		   int in1=ind[n]-gridIndex[0];
//		   float grv=m_grid[0][ind[n]-gridIndex[0]];

		   m_fArray[m_nOfCol-1][ptId]+=m_grid[0][ind[n]-gridIndex[0]]*d[n]*cellVolume;
		  }
	  	m_fArray[0][ptId]=m_fArray[m_nOfCol-1][ptId]; //prepare array for write
	  } // for

	}  //if(!tsc)
	if(m_tsc)
	{
  	   int nCell=m_sampleDimensions[0]*m_sampleDimensions[1]*m_sampleDimensions[2];
  	   for (int ptId=0; ptId < m_nOfEle; ptId++)
    	   {
		m_fArray[m_nOfCol-1][ptId]=0.0000001;      	        float px[3];
		int ind_x,ind_y,ind_z;       
		int xList[4],yList[4],zList[4];
		int nCell=m_gridCellNumber[0]*m_gridCellNumber[1]*m_gridCellNumber[2];
        	double node_x, node_y, node_z, dist_x, dist_y, dist_z, w_x, w_y, w_z, w;
		px[0]=m_fArray[0][ptId]-m_origin[0];
		if(periodic && px[0]<0.) px[0]+=m_sampleDimensions[0]*m_spacing[0];
		if(periodic && px[0]>=(m_sampleDimensions[0]*m_spacing[0])) px[0]-=m_sampleDimensions[0]*m_spacing[0];

		px[1]=m_fArray[1][ptId]-m_origin[1];
		if(periodic && px[1]<0.) px[1]+=m_sampleDimensions[1]*m_spacing[1];
		if(periodic && px[1]>=(m_sampleDimensions[1]*m_spacing[1])) px[1]-=m_sampleDimensions[1]*m_spacing[1];

		px[2]=m_fArray[2][ptId]-m_origin[2];
 		if(periodic && px[2]<0.) px[2]+=m_sampleDimensions[2]*m_spacing[2];
		if(periodic && px[2]>=(m_sampleDimensions[2]*m_spacing[2])) px[2]-=m_sampleDimensions[2]*m_spacing[2];
; 		
		ind_x = (int)floor(px[0]/m_spacing[0]);
		ind_y = (int)floor(px[1]/m_spacing[1]);
    		ind_z = (int)floor(px[2]/m_spacing[2]);
		for(int j=0;j<4;j++)
		{
			xList[j]=ind_x+j-1;    
			yList[j]=ind_y+j-1;    
			zList[j]=ind_z+j-1;
		}


		for(int ix = 0; ix <4;ix++)
		   for(int iy = 0; iy <4;iy++)
    			for(int iz = 0; iz <4;iz++)
			{
	        		dist_x = px[0]/m_spacing[0] - xList[ix];
                		dist_y = px[1]/m_spacing[1] - yList[iy];
                		dist_z = px[2]/m_spacing[2] - zList[iz];
                		if(dist_x<0) dist_x=-1*dist_x;
                		if(dist_y<0) dist_y=-1*dist_y;
                		if(dist_z<0) dist_z=-1*dist_z;

                		if(dist_x <= .5)   w_x = .75 - dist_x*dist_x;
                		else if (dist_x <= 1.5) w_x = .5*(1.5-dist_x)*(1.5-dist_x);
                		else w_x = 0;
                
                		if(dist_y <= .5)   w_y = .75 - dist_y*dist_y;
                		else if (dist_y <= 1.5) w_y = .5*(1.5-dist_y)*(1.5-dist_y);
                		else w_y = 0;
                
                		if(dist_z <= .5)   w_z = .75 - dist_z*dist_z;
                		else if (dist_z <= 1.5) w_z = .5*(1.5-dist_z)*(1.5-dist_z);
                		else w_z = 0;                
                
                		if( (w_x < 0. || w_y  < 0.) ||  w_z < 0.)
					std::cerr<<"WARNING: particle "<<ptId<<" xyzList "<<xList[ix]<<" "<< yList[iy]<<" "<<zList[iz]<<" AND w_x,y,z "<<w_x<<" "<<w_y<<" "<<w_z<<std::endl;
	
                		w = w_x * w_y * w_z;
                		if (w != 0)
                		{
					int facX= xList[ix];
					int facY= yList[iy];
					int facZ= zList[iz];
					if(w>1.01) 
					{
						std::cerr<<"Error 1 on tsc schema. operation Aborted."<<std::endl;
						return false;
					}
					if(periodic && facX<0)
						 facX+=m_sampleDimensions[0];
					if(periodic && facX>=m_sampleDimensions[0]) facX-=m_sampleDimensions[0];
					if(periodic && facY<0) facY+=m_sampleDimensions[1];
					if(periodic && facY>=m_sampleDimensions[1]) facY-=m_sampleDimensions[1];
					if(periodic && facZ<0) facZ+=m_sampleDimensions[2];
					if(periodic && facZ>=m_sampleDimensions[2]) facZ-=m_sampleDimensions[2];

 					if(facX<0  || facX>= m_sampleDimensions[0])
						continue;
					if(facY<0  || facY>= m_sampleDimensions[1])
						continue;
					if(facZ<0  || facZ>= m_sampleDimensions[2])
						continue;

                			int ind_test = facX + m_sampleDimensions[0]*facY + m_sampleDimensions[0]*m_sampleDimensions[1]*facZ;
					if(ind_test<0 || ind_test>(nCell-1))
					{
						std::cerr<<"Error on particle "<<ptId<<" xyzList "<<xList[ix]<<" "<< yList[iy]<<" "<<zList[iz]<<" AND w w_x,y,z "<<w<<" "<<w_x<<" "<<w_y<<" "<<w_z<<std::endl;
					}
					if(ind_test<gridIndex[0] || ind_test>gridIndex[1]) //ind1 NOT in cache!
		   			{
						gridIndex[0]=ind_test;
						gridIndex[1]=gridIndex[0]+m_nOfGridEle-1;
						if(gridIndex[1]>VolumeTable.getNumberOfRows())gridIndex[1]=VolumeTable.getNumberOfRows()-1;
 						VolumeTable.getColumn(gridColList,1,gridIndex[0],gridIndex[1],m_grid); 	
		   			}

					m_fArray[m_nOfCol-1][ptId]+=m_grid[0][ind_test-gridIndex[0]]*w*cellVolume; // TBV

                 		}
                 
           		}// closes for ix iy iz
	  m_fArray[0][ptId]=m_fArray[m_nOfCol-1][ptId]; //prepare array for write
	  }// for ptId
	} //if m_tsc
	unsigned int resultColList[1];
	resultColList[0]=tableGrid2Point.getNumberOfColumns()-1;
	tableGrid2Point.putColumn(resultColList,1,index,index+m_nOfEle-1, m_fArray);	

	index+= m_nOfEle;
	if(index+m_nOfEle>= m_nOfRow) 
		m_nOfEle= m_nOfRow-index;
}

return true; 
}
//---------------------------------------------------------------------


