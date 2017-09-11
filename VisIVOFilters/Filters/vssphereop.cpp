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
/*#include <set>*/
#include <fstream>
#include <cmath>
#ifdef WIN32
	#include <time.h>
#endif
#include "vstable.h"
#include "vssphereop.h"
#include "VisIVOFiltersConfigure.h"
#include "visivoutils.h"

const unsigned int VSSphereOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSSphereOp::MIN_NUMBER_OF_ROW = 100;


//---------------------------------------------------------------------
VSSphereOp::VSSphereOp()
//---------------------------------------------------------------------
{
 m_fArray=NULL;
 m_fArrayWrite=NULL;
 m_nOfRow=0;
 m_nOfCol=0;	
}
//---------------------------------------------------------------------
VSSphereOp::~VSSphereOp()
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
bool VSSphereOp::allocatefArray()
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
void VSSphereOp::printHelp()
//---------------------------------------------------------------------
{
	std::cout<<"Create a  new table from an input table of a sub-box or of a sphere. Operation not allowed on volumes."<<std::endl<<std::endl;
	std::cout<<"Usage: VisIVOFilters --op extraction --geometry geometry_file [--out filename_out.bin] [--history] [--historyfile filename.xml] [--help] [--file] inputFile.bin"<<std::endl;

	std::cout<<"Example: VisIVOFilters --op extraction  --geometry geometry.txt --out pos_extracted.bin --file pos.bin"<<std::endl;

	std::cout<<"Note: --help produce this output "<<std::endl;
	std::cout<<"  	  --geometry file must contain three valid column names and a value for each column. The fourth field means the extraction mode: RADIUS: a sphere will be extracted, BOX and CORNER: a rectangular region will be extracted."<<std::endl; 
	std::cout<<""<<std::endl; 
	std::cout<<"A fourth filed named RADIUS  indicates the sphere radius and the columns values represent the sphere center"<<std::endl;
	std::cout<<"Examples:"<<std::endl; 
	std::cout<<"X	25.0"<<std::endl;
	std::cout<<"Y	25.0"<<std::endl;
	std::cout<<"Z	25.0"<<std::endl;
	std::cout<<"RADIUS	5.0"<<std::endl<<std::endl;
	std::cout<<"A fourth filed named BOX  indicates the half side measure and the columns values represent the rectangular center"<<std::endl;
	std::cout<<"Example:"<<std::endl; 
	std::cout<<"X	25.0"<<std::endl;
	std::cout<<"Y	25.0"<<std::endl;
	std::cout<<"Z	25.0"<<std::endl;
	std::cout<<"BOX	5.0"<<std::endl<<std::endl;
	std::cout<<"A fourth filed named CORNER  indicates the  side measure and the columns values represent the rectangular lower corner"<<std::endl;
	std::cout<<"Example:"<<std::endl; 
	std::cout<<"X	0.0"<<std::endl;
	std::cout<<"Y	0.0"<<std::endl;
	std::cout<<"Z	0.0"<<std::endl;
	std::cout<<"CORNER	10.0"<<std::endl;
	std::cout<<"--out Name of the new table. Default name is given."<<std::endl;
    std::cout<<"--history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
    std::cout<<"--historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;
	std::cout<<"--help produce this output "<<std::endl;
	return;
}

//---------------------------------------------------------------------
bool VSSphereOp::execute()
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

std::string filename;
	unsigned long long int totRows,nOfCols=-1;
	unsigned long long int maxRows=-1;
	VSTable tableSphere;
	float radius;
	std::vector<int> colIdVector;
	int extraction=-1;
	if(m_tables[0]->getIsVolume())
	{
	 	std::cerr<<"vsextractionop: operation not allowed: a volume is selected"<<std::endl; 
		return false;
	}
	if(getParameterAsString("geometry").empty() || getParameterAsString("geometry")=="unknown" )
	{
		std::cerr<<"vssextractionop: No file with sphere geometry is given"<<std::endl;
		return false;
	} else
	{
		filename=getParameterAsString("geometry");
	}

	std::ifstream fileInput(filename.c_str());
	if(!fileInput)
	{
		std::cerr<<"Cannot open geometry file"<<filename<<std::endl;
		return false;
	}
	while(!fileInput.eof()) 
	{
	        std::string colName;
		std::string sValue;
		float value;
		fileInput >> colName;
		fileInput >> sValue;
		if(colName=="RADIUS") 
		{	
			radius=atof(sValue.c_str());
			extraction=1;
		} 
		if(colName=="BOX") 
		{	
			radius=atof(sValue.c_str());
			extraction=2;
		} 
		if(colName=="CORNER") 
		{	
			radius=atof(sValue.c_str());
			extraction=3;
		} 
		if(m_tables[0]->getColId(colName)>=0)
		{
			value=atof(sValue.c_str());
			m_colVector.push_back(value);
			colIdVector.push_back(m_tables[0]->getColId(colName));
		}

		
	}
	if(extraction==-1)
	{		
		std::cerr<<"vssextractionop: Invalid Extraction mode. Only CORNER, BOX, RADIUS are allowed"<<std::endl;
		return false;
	}
	fileInput.close();	
	if(m_colVector.size()!=3) return false;
	if(extraction==-1) return false;
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
  		fileNameOutputSStream<<filenameInputTable.substr(0, len-4)<<"_extract_"<<buffer<<".bin";  //QUI verificare
	}

	fileNameOutput=fileNameOutputSStream.str();
  	if(fileNameOutput.find(".bin") == std::string::npos)
	    fileNameOutput.append(".bin");
	std::string tmpDir=getDir(fileNameOutput);
//Clean existing tab
	remove(fileNameOutput.c_str());
	tableSphere.setLocator(fileNameOutput);
#ifdef VSBIGENDIAN
	std::string endianism="big";
	
#else	
	std::string endianism="little";
#endif

	tableSphere.setEndiannes(endianism);
 	tableSphere.setType("float");
	for(unsigned int k=0;k<m_tables[0]->getNumberOfColumns();k++)
		tableSphere.addCol(m_tables[0]->getColName(k));	
// 

// 	tableSphere.


	for(int i=0;i<nOfCol;i++) colList[i]=i;

	unsigned long long int totEle=nOfEle;
	bool goodEle;


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

			switch(extraction)
			{
			  case 1:
			  {
			  /*** Sphere Extraction **/
			    goodEle=false;
			    for(int j=0;j<3;j++)
			    {
				int colId=colIdVector[j];
				dist2=dist2+(m_fArray[colId][i]-m_colVector[j])*(m_fArray[colId][i]-m_colVector[j]);
			    }
			    dist= sqrt(dist2); //QUI verifica
			    if(dist<=radius) 
				goodEle=true;
			    break;	
			  }
			

			  case 2:
			  {
			  /*** BOX Extraction **/
			    goodEle=true;
			    for(int j=0;j<3;j++)
			    {
				int colId=colIdVector[j];
				if(fabs(m_fArray[colId][i]-m_colVector[j]) > radius)
				{
					goodEle=false;
					break;
				}
			    }
			    break;	
			  }
			
			  case 3:
			  {
			  /*** CORNER Extraction **/
			    goodEle=true;
			    for(int j=0;j<3;j++)
			    {
				int colId=colIdVector[j];
				int distanza=m_fArray[colId][i]-m_colVector[j];
				if((m_fArray[colId][i]-m_colVector[j]) > radius || (m_fArray[colId][i]-m_colVector[j]) < 0.)
				{	
					goodEle=false;
					break;
				}
			    }
			    break;	
			  }
			}


			if(goodEle)
			{
				// write temporary table
				if(writeCounter==maxEle ||writeCounter<0 )
				{
					std::stringstream convertStream;
					convertStream<<tmpDir<<"Extract_temp__"<<randat<<numberOfTempTables<<".bin";
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
				for(unsigned int j=0;j<nOfCol;j++) m_fArrayWrite[j][writeCounter]=m_fArray[j][i];
					writeCounter++;
			}
		}		
		startCounter=toRow+1;
		totEle=totEle-(toRow-fromRow+1);
	}

	if(numberOfTempTables==0) // No temp tables are written
	{
		unsigned long long int wrNOfRow=writeCounter;
 		tableSphere.setNumberOfRows(wrNOfRow);
		tableSphere.writeTable(m_fArrayWrite);
	} else	
	{
		unsigned long long int fromWRow, toWRow=-1;
		unsigned long long int wrNOfRow=(unsigned long long int) (maxEle)*numberOfTempTables+writeCounter;

 		tableSphere.setNumberOfRows(wrNOfRow);
		tableSphere.writeTable();

		for(unsigned int i=0;i<numberOfTempTables;i++)
		{
			std::stringstream convertStream,convertStreamHead;
			convertStream<<tmpDir<<"Extract_temp__"<<randat<<i<<".bin";
			convertStreamHead<<tmpDir<<"Extract_temp__"<<randat<<i<<".bin.head";
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
			tableSphere.putColumn(colList,nOfCol, fromWRow,  toWRow, m_fArray);
			remove(filenameTemp.c_str());
			remove(filenameTempHead.c_str());
		}
		if(writeCounter>0) //exist m_fArrayWrite partially filled
		{
			fromWRow=toWRow+1;
			toWRow=fromWRow+writeCounter-1;
			tableSphere.putColumn(colList,nOfCol, fromWRow,  toWRow, m_fArrayWrite);
			
		}
	}
delete [] colList;		
return true;
}

