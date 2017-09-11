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


#include <cstdlib>
#include <cstring>
#include <iostream>
/*#include <string>*/
#include <sstream>

#include <set>
#include <map>
#include <fstream>
#include <cmath>
#ifdef WIN32
	#include <time.h>
#endif
#include "vstable.h"
#include "vsmrcampos.h"
#include "VisIVOFiltersConfigure.h"
#include "visivoutils.h"
#include "vsstatisticop.h"



const unsigned int VSMRCamPos::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSMRCamPos::MIN_NUMBER_OF_ROW = 100;


//---------------------------------------------------------------------
VSMRCamPos::VSMRCamPos()
//---------------------------------------------------------------------
{
 m_fArray=NULL;
 m_fArrayWrite=NULL;
 m_nOfRow=0;
 m_nOfCol=0;	
}
//---------------------------------------------------------------------
VSMRCamPos::~VSMRCamPos()
//---------------------------------------------------------------------
{
	if(m_fArray!=NULL)
 		for(unsigned int i=0;i<m_nOfCol;i++)
		{
			if(m_fArray[i]!=NULL) delete [] m_fArray[i];
			if(m_fArrayWrite[i]!=NULL) delete [] m_fArrayWrite[i];
		}
	if(m_fArray!=NULL) delete [] m_fArray;
	if(m_fArrayWrite!=NULL) delete [] m_fArrayWrite;
}
//---------------------------------------------------------------------
bool VSMRCamPos::allocatefArray()
//---------------------------------------------------------------------
{
unsigned long long int tempLL=getMaxNumberInt();
if(((unsigned long long int)m_nOfRow*m_nOfCol)>tempLL) 
	m_nOfRow=(int)tempLL/m_nOfCol;

try
{
	m_fArray=new  float*[m_nOfCol];
}
catch(std::bad_alloc &e)
{
	m_fArray=NULL;
}

	if(m_fArray==NULL)
	{
		return false;
	}
try
{
	m_fArrayWrite=new  float*[m_nOfCol];
}
catch(std::bad_alloc &e)
{
	m_fArrayWrite=NULL;
}

	if(m_fArrayWrite==NULL)
	{
		return false;
	}

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
				{	
					delete [] m_fArray[j];
					delete [] m_fArrayWrite[j];
				}
				if(m_nOfRow==MIN_NUMBER_OF_ROW)
				{ 
					delete [] m_fArray;
					delete [] m_fArrayWrite;
					m_fArray=NULL;
					m_fArrayWrite=NULL;
					return false;
				}
				if(m_nOfRow<=MAX_NUMBER_TO_REDUCE_ROW)
					m_nOfRow=MIN_NUMBER_OF_ROW;
				else
					m_nOfRow=m_nOfRow-MAX_NUMBER_TO_REDUCE_ROW;
				break;
			}

try
{
			m_fArrayWrite[i]=new  float[m_nOfRow];
}
catch(std::bad_alloc &e)
{
	m_fArrayWrite[i]=NULL;
}

			if(m_fArrayWrite[i]==NULL) 
			{	
				goodAllocation=false;
				for(unsigned int j=0;j<i;j++)
				{ 	
					delete [] m_fArrayWrite[j];
					delete [] m_fArray[j];
				}
				delete [] m_fArray[i];
				if(m_nOfRow==MIN_NUMBER_OF_ROW)
				{ 
					delete [] m_fArray;
					delete [] m_fArrayWrite;
					m_fArray=NULL;
					m_fArrayWrite=NULL;
					return false;
				}
				if(m_nOfRow<=MAX_NUMBER_TO_REDUCE_ROW) 
					m_nOfRow=MIN_NUMBER_OF_ROW;
				else
					m_nOfRow=m_nOfRow-MAX_NUMBER_TO_REDUCE_ROW;
				break;
			}

		}
		
	}
	return true;
}
//---------------------------------------------------------------------
void VSMRCamPos::printHelp()
//---------------------------------------------------------------------
{
	std::cout<<"Create a  new table in Multi Resolution from an input table. Operation not yet allowed on volumes."<<std::endl<<std::endl;
	std::cout<<"Usage: VisIVOFilters --op mres --points x_col y_col z_col [--pos values]  [--geometry layer_file.txt] [--background vaule] [--out filename_out.bin] [--history] [--historyfile filename.xml] [--help] [--file] inputFile.bin"<<std::endl;

	std::cout<<"Example: VisIVOFilters --op mres --points X Y Z --pos 10.0 10.0 10.0  --geometry layer.txt --out pos_layer.bin --file pos.bin"<<std::endl;

	std::cout<<"Note: --help produce this output "<<std::endl;
	std::cout<<"--points columns to be assumed for points coordinates."<<std::endl;
	std::cout<<"--pos camera point coordinates. Default value is the center of the domain."<<std::endl;
	std::cout<<"--geometry file must contain a radius and a randomizator value: 1.0 all values included.";
	std::cout<<"0.1 means 1 per cent of value included in the layer. Three default layers are created"<<std::endl; 
	std::cout<<""<<std::endl; 
	std::cout<<"Examples:"<<std::endl; 
	std::cout<<"5.0	1.0"<<std::endl;
	std::cout<<"10.0 0.1"<<std::endl;
	std::cout<<"20.0 0.05"<<std::endl;
	std::cout<<"30.0 0.01"<<std::endl;
	std::cout<<"70.0 0.005"<<std::endl;
	std::cout<<"--background  a randomizator value for points outside the geometry.";
	std::cout<<" Default value is maximum 100000 values from input VBT"<<std::endl;
	std::cout<<"--out Name of the new table. Default name is given."<<std::endl;
    std::cout<<"--history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
    std::cout<<"--historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;
	std::cout<<"--help produce this output "<<std::endl;
	return;
}

//---------------------------------------------------------------------
bool VSMRCamPos::execute()
//---------------------------------------------------------------------
// totRows is the total rows of the out tables
// maxRows is the bigger number of rows among input tables
// 
// maxEle maximum number allowed to upload/download data from a generic table: It is the lower value between maxRows, maxInt and fArray allocation fArray[0][maxEle]
// 
// nOfEle total number of element to upload for the specific input table.
// 
// totEle is a counter: the number of elements that I have still to  upload for the specific table. Starts from nOfEle, its value decrese of maxEle each cycle
// 
{	
std::stringstream sstmp1;
std::string  randat;
time_t rawtime;
struct tm * timeinfo;
char buffer [80];
time ( &rawtime );
timeinfo = localtime ( &rawtime );
strftime (buffer,80,"%Y%m%d%H%M%S",timeinfo);
sstmp1<<"_"<<rand()<<buffer<<"_";  
randat=sstmp1.str();
bool defaultGeometry=false;
std::set<float>::iterator itset;

std::string filename;
unsigned long long int totRows,nOfCols=-1;
unsigned long long int maxRows=-1;
VSTable tableMR;

std::vector<int> colIdVector;
std::set<float> levelRow;

float pointsRange[3][2];
float interval[3];
float minInterval;
float camPos[3]; 
float background;
int eleLevNum[100];
for(int i=0;i<100;i++) eleLevNum[i]=0;

if(m_tables[0]->getIsVolume())
{
	 	std::cerr<<"vsmrcampos: operation not yet allowed on a volume table"<<std::endl; 
		return false;
}
background=100000.0/m_tables[0]->getNumberOfRows();
if(isParameterPresent("background"))
  background=getParameterAsFloat("background");
if(background>1.0) background=1.0;
  

if(!isParameterPresent("points"))
{
 	std::cerr<<"vsmrcampos: points option not given"<<std::endl; 
	return false;
  
} else{
  m_points=getParameterAsString("points");
  std::stringstream sstmp(m_points);
  for(int i=0;i<3;i++)
  {
    std::string col;
    sstmp>>col;
    if(m_tables[0]->getColId(col)==-1)
    {
 	std::cerr<<"vsmrcampos: invalid col name "<<col <<" in points option"<<std::endl; 
	return false;
    }
    colIdVector.push_back(m_tables[0]->getColId(col));
  }
}
{
  VSStatisticOp op;
  op.addParameter("list",m_points);
  op.addParameter("silent","");
//  std::stringstream tmp;
//  tmp<<numberOfBin;
//  op.addParameter("histogram",tmp.str());
  op.addInput(m_tables[0]);
  if(!op.execute())
  {
	std::cerr<<"Invalid statisticop. Operation aborted"<<std::endl;
	return false;
  }
  float dummy;
  for(int i=0;i<3;i++)
  {  
    op.getRange(i,pointsRange[i][0],pointsRange[i][1],dummy,dummy);
    interval[i]=pointsRange[i][0]-pointsRange[i][1];
  }
  minInterval=interval[0];
  for(int i=0;i<3;i++) if(minInterval<interval[i]) minInterval=interval[i];
  
//  histogramPointer=op.getHisto();
}
if(!isParameterPresent("pos") || getParameterAsString("pos").empty() || getParameterAsString("pos")=="unknown")
{
    camPos[0]=(pointsRange[0][0]-pointsRange[0][1])/2.0;
    camPos[1]=(pointsRange[1][0]-pointsRange[1][1])/2.0;
    camPos[2]=(pointsRange[2][0]-pointsRange[2][1])/2.0;
} else{

  std::stringstream sstmp(getParameterAsString("pos"));
  sstmp>>camPos[0];
  sstmp>>camPos[1];
  sstmp>>camPos[2];
}

if(getParameterAsString("geometry").empty() || getParameterAsString("geometry")=="unknown" )
{
// default geometry;
		float key, value;
		int fact=1;
		for(int i=0;i<3;i++)
		{
		  key=minInterval*0.01*fact;
		  value=1.0/fact;
		  fact=fact*10;
		  levelRow.insert(key);	
		  m_layers.insert(std::make_pair(key,value));
		}
} else{
	    filename=getParameterAsString("geometry");
	    std::ifstream fileInput(filename.c_str());
	    if(!fileInput)
	    {
		std::cerr<<"Cannot open geometry file"<<filename<<std::endl;
		return false;
	    }
	    while(!fileInput.eof()) 
	    {
	        float key;
		float value;
		fileInput >> key;
		fileInput >> value;
		levelRow.insert(key);		
		m_layers.insert(std::make_pair(key,value));
	    }

	    fileInput.close();	
}
if(levelRow.size()==0) return false;

//for (itset=levelRow.begin(); itset!=levelRow.end(); itset++)
//  std::clog<<*itset<<std::endl;

std::map<float,float>::iterator iter;
//for (iter =m_layers.begin(); iter!=m_layers.end(); iter++)
//  std::clog<<iter->first<<"  "<<iter->second<<std::endl;


int *levelNum,*levelNumSave;

levelNum=new int[levelRow.size()+1];//number of elements to discard for each level
levelNumSave=new int[levelRow.size()+1];//number of elements to discard for each level saved
int icount=0;
for (itset=levelRow.begin(); itset!=levelRow.end(); itset++)
{
//  std::clog<<*itset<<std::endl;
  iter=m_layers.find(*itset);
  levelNum[icount]=1.0/iter->second;
  levelNumSave[icount]=1.0/iter->second;
//  std::clog<<iter->second<<std::endl;
  icount++;
}
levelNum[levelRow.size()]=1.0/background;
levelNumSave[levelRow.size()]=1.0/background;

int maxInt=getMaxNumberInt();
unsigned int maxEle;
unsigned int numberOfTempTables=0;
unsigned int writeCounter=0;
totRows=m_tables[0]->getNumberOfRows();
unsigned long long int nOfEle=totRows;

if(totRows>maxInt)		
  maxEle=maxInt; 
else
  maxEle=totRows;
	
unsigned long long int fromRow, toRow, startCounter=0;

unsigned int *colList=NULL;
unsigned int nOfCol=m_tables[0]->getNumberOfColumns();
try
{
	colList=new unsigned int[nOfCol];
}
catch(std::bad_alloc &e)
{
	colList=NULL;
}

if(colList==NULL)
{
		std::cerr<<"Failed colList allocation. Select Field terminated"<<std::endl;
		return false;
}

m_nOfCol=(unsigned int) nOfCol;
m_nOfRow=maxEle;
bool allocationfArray=allocatefArray();
maxEle=m_nOfRow;
if(m_fArray==NULL ||m_fArrayWrite==NULL ||  !allocationfArray )	
{
		std::cerr<<"Failed fArray allocation. Select Field terminated"<<std::endl;
		delete [] colList;
		return false;
}

std::stringstream fileNameOutputSStream;
fileNameOutputSStream<<getParameterAsString("out");
std::string fileNameOutput;

if(fileNameOutputSStream.str()==""||fileNameOutputSStream.str()=="unknown")
{
		fileNameOutputSStream.str().erase(); //QUI verificare
  		std::string filenameInputTable=m_tables[0]->getLocator();
  		int len=filenameInputTable.length();
		time_t rawtime;
		struct tm * timeinfo;
		char buffer [80];
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		strftime (buffer,80,"%Y%m%d%H%M",timeinfo);
  		fileNameOutputSStream<<filenameInputTable.substr(0, len-4)<<"_mrcampos_"<<buffer<<".bin";  //QUI verificare
}

fileNameOutput=fileNameOutputSStream.str();
if(fileNameOutput.find(".bin") == std::string::npos)
	    fileNameOutput.append(".bin");
m_realOutFilename.push_back(fileNameOutput);
std::string tmpDir=getDir(fileNameOutput);
//Clean existing tab
remove(fileNameOutput.c_str());
tableMR.setLocator(fileNameOutput);
#ifdef VSBIGENDIAN
	std::string endianism="big";
	
#else	
	std::string endianism="little";
#endif

tableMR.setEndiannes(endianism);
tableMR.setType("float");
for(unsigned int k=0;k<m_tables[0]->getNumberOfColumns();k++)
		tableMR.addCol(m_tables[0]->getColName(k));	

for(int i=0;i<nOfCol;i++) colList[i]=i;

unsigned long long int totEle=nOfEle;
bool goodEle;

int *randomLevelNumber;
randomLevelNumber=new int[levelRow.size()+1]; //+1 for the background
for(int j=0;j<levelRow.size()+1;j++)
{
  int ck=rand();
  int num=ck%levelNumSave[j]-1;
  if(num<0) num=0;
//  std::clog<<num<<" " <<ck<<std::endl;
  randomLevelNumber[j]=num;
}

while(totEle!=0)
{
  fromRow=startCounter;
  toRow=fromRow+maxEle-1;
  if(toRow>totRows-1)toRow=totRows-1;
  m_tables[0]->getColumn(colList, nOfCol, fromRow, toRow, m_fArray);
  for(unsigned int i=0;i<toRow-fromRow+1;i++)
  {
	float dist=0.;
	float dist2=0.;
	int levelOfPoint=levelRow.size();  //this is an index for arraies
//	std::clog<<levelRow.size()<<std::endl;
	for(int j=0;j<3;j++)
	{
		int colId=colIdVector[j];
		dist2=dist2+(m_fArray[colId][i]-camPos[j])*(m_fArray[colId][i]-camPos[j]);
	}
	dist= sqrt(dist2); //QUI verifica
	std::set<float>::iterator iterSet;
	int j=0;
	for(iterSet=levelRow.begin(); iterSet!=levelRow.end();iterSet++)
	{  
	    
	    if(dist<=*iterSet) 
	    {  
	      levelOfPoint=j;
	      break;
	    }
	    j++;
	}
	if(writeCounter==maxEle ||writeCounter<0 )
	{
		std::stringstream convertStream;
		convertStream<<tmpDir<<"MRcampos_temp__"<<randat<<numberOfTempTables<<".bin";
		std::string filenameTemp;
		convertStream>>filenameTemp;
		numberOfTempTables++;
		VSTable tableTemp;
		tableTemp.setLocator(filenameTemp);
		tableTemp.setEndiannes(endianism);
 		tableTemp.setType("float");
		unsigned long long int wrNOfRow=writeCounter;
 		tableTemp.setNumberOfRows(wrNOfRow);

		for(unsigned int k=0;k<m_tables[0]->getNumberOfColumns();k++)
			tableTemp.addCol(m_tables[0]->getColName(k));

					tableTemp.writeTable(m_fArrayWrite);
					writeCounter=0;
	 }
//	 if(levelOfPoint<1) 
//	   std::clog<<levelOfPoint<<" "<<levelNum[levelOfPoint]<<" "<<randomLevelNumber[levelOfPoint]<<std::endl;
	 if((levelNum[levelOfPoint]-1)==randomLevelNumber[levelOfPoint])
	 {
	   for(unsigned int j=0;j<nOfCol;j++)
	        m_fArrayWrite[j][writeCounter]=m_fArray[j][i];
	    eleLevNum[levelOfPoint]++;
	    writeCounter++;
	 }
	 
	 levelNum[levelOfPoint]--;
	 if(levelNum[levelOfPoint]<=0)
	 {
	    int num=rand()%levelNumSave[levelOfPoint];
	    randomLevelNumber[levelOfPoint]=num;
	   levelNum[levelOfPoint]=levelNumSave[levelOfPoint];  
	 }
  }//for		
  startCounter=toRow+1;
  totEle=totEle-(toRow-fromRow+1);
}//while

if(numberOfTempTables==0) // No temp tables are written
{
    unsigned long long int wrNOfRow=writeCounter;
    tableMR.setNumberOfRows(wrNOfRow);
    tableMR.writeTable(m_fArrayWrite);
} else{
    unsigned long long int fromWRow, toWRow=-1;
    unsigned long long int wrNOfRow=(unsigned long long int) (maxEle)*numberOfTempTables+writeCounter;

    tableMR.setNumberOfRows(wrNOfRow);
    tableMR.writeTable();

    for(unsigned int i=0;i<numberOfTempTables;i++)
    {
	std::stringstream convertStream,convertStreamHead;
	convertStream<<tmpDir<<"MRcampos_temp__"<<randat<<i<<".bin";
	convertStreamHead<<tmpDir<<"MRcampos_temp__"<<randat<<i<<".bin.head";
	std::string filenameTemp;
	std::string filenameTempHead;
	convertStream>>filenameTemp;
	filenameTempHead=filenameTemp+".head";
	VSTable tableTemp(filenameTemp);
	fromWRow=toWRow+1;
	toWRow=fromWRow+maxEle-1;
	fromRow=0;
	toRow=maxEle-1;
	tableTemp.getColumn(colList,nOfCol, fromRow,  toRow, m_fArray);
	tableMR.putColumn(colList,nOfCol, fromWRow,  toWRow, m_fArray);
	remove(filenameTemp.c_str());
	remove(filenameTempHead.c_str());
    }
    if(writeCounter>0) //exist m_fArrayWrite partially filled
    {
	fromWRow=toWRow+1;
	toWRow=fromWRow+writeCounter-1;
	tableMR.putColumn(colList,nOfCol, fromWRow,  toWRow, m_fArrayWrite);
			
    }
}
for(int i=0;i<=levelRow.size();i++) std::cerr<<"Level "<<i<<" Elements: "<<eleLevNum[i]<<std::endl;
delete [] colList;		
return true;
}

