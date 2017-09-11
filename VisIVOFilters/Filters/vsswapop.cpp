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
#include <sstream>
#include <fstream>
#ifdef WIN32
	#include <time.h>
#endif
#include "vsswapop.h"
#include "vstable.h"
#include "visivoutils.h"

//prova
const unsigned int VSSwapOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSSwapOp::MIN_NUMBER_OF_ROW = 100;

//---------------------------------------------------------------------
VSSwapOp::VSSwapOp()
{
 m_fArray=NULL;
 m_colList=NULL;
 m_nOfRow=0;
 m_nOfCol=0;
 m_nOfEle=0;	

}
//---------------------------------------------------------------------


//---------------------------------------------------------------------
VSSwapOp::~VSSwapOp()
{
	if(m_fArray!=NULL)
 	  for(unsigned int i=0;i<m_nOfCol;i++)
	  {
		if(m_fArray[i] != NULL) delete [] m_fArray[i];
	  }
	if(m_fArray!=NULL) delete [] m_fArray;
	if(m_colList!=NULL) delete [] m_colList;
}
//---------------------------------------------------------------------

//---------------------------------------------------------------------
void VSSwapOp::printHelp()
//---------------------------------------------------------------------
{
std::cout<<"Endianism swap: little Endian to big Endian and viceversa of a VisIVO Binary Table."<<std::endl<<std::endl;
std::cout<<"Usage: VisIVOFilters --op swap   [--override] [--out filename_out.bin] [--history] [--historyfile filename.xml] [--help] [--file] inputFile.bin"<<std::endl<<std::endl;

std::cout<<"Example: VisIVOFilters --op swap  --override --file inputFile.bin"<<std::endl;

std::cout<<"Note:  "<<std::endl;
std::cout<<"--override The same input table is swapped."<<std::endl;
std::cout<<"--out Name of the new table. Default name is given. Ignored if --override is specified."<<std::endl;
std::cout<<"--history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
std::cout<<"--historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;

std::cout<<"--file Input table filename."<<std::endl;

std::cout<<"--help produce this output."<<std::endl;

return;

}

//---------------------------------------------------------------------
bool VSSwapOp::allocateArray()
//---------------------------------------------------------------------
{ 
unsigned long long int tempLL=getMaxNumberInt();
if(((unsigned long long int)m_nOfEle*m_nOfCol)>tempLL) 
	m_nOfEle=(int)tempLL/m_nOfCol;

	m_colList=new unsigned int[m_nOfCol];
try
{
	m_fArray=new  float*[m_nOfCol];
}
catch(std::bad_alloc &e)
{
	m_fArray=NULL;
}

	if(m_fArray==NULL)
		return false;

	bool goodAllocation=false;
	while(!goodAllocation)
	{
		goodAllocation=true;
		for(unsigned int i=0;i<m_nOfCol;i++)
		{
try
{
			m_fArray[i]=new  float[m_nOfEle];
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
				if(m_nOfEle==MIN_NUMBER_OF_ROW)
				{ 
					delete [] m_fArray;
					m_fArray=NULL;
					return false;
				}
				m_nOfEle=m_nOfEle-MAX_NUMBER_TO_REDUCE_ROW;
				if(m_nOfEle<=MAX_NUMBER_TO_REDUCE_ROW) m_nOfEle=MIN_NUMBER_OF_ROW;
				break;
			}

		}
		
	}
	return true;
}


//---------------------------------------------------------------------
bool VSSwapOp::execute()
//---------------------------------------------------------------------
{
bool override=false; 
if(isParameterPresent("override")) override=true;
std::string fileNameOutput;
if(!override)
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
  		fileNameOutputSStream<<filenameInpuTable.substr(0, len-4)<<"_swapEndianism_"<<buffer<<".bin";
		fileNameOutput=fileNameOutputSStream.str();
	}
  	if(fileNameOutput.find(".bin") == std::string::npos)
	    		fileNameOutput.append(".bin");
} else
	fileNameOutput=m_tables[0]->getLocator();

m_realOutFilename.push_back(fileNameOutput);



m_nOfRow=m_tables[0]->getNumberOfRows();
if(m_nOfRow>getMaxNumberInt()) 
	m_nOfEle= getMaxNumberInt();
else 
	m_nOfEle= m_nOfRow;

m_nOfCol=m_tables[0]->getNumberOfColumns();
bool allocationArray=allocateArray();
if(!allocationArray)
{	
	std::cerr<<"Failed Array allocation. Swap aborted"<<std::endl;
	return false;
}
for(int i=0;i<m_nOfCol;i++)
	m_colList[i]=i;
VSTable tableSwap;
tableSwap.setLocator(fileNameOutput);
std::string endianism=m_tables[0]->getEndiannes();
if(m_tables[0]->getIsVolume())
{
	tableSwap.setIsVolume(true);
	const unsigned int *temp;
	temp=new unsigned int[3];
	temp=m_tables[0]->getCellNumber();
	tableSwap.setCellNumber(temp[0],temp[1],temp[2]);
	const float *tempf;
	tempf=new float[3];
	tempf=m_tables[0]->getCellSize();
	tableSwap.setCellSize(tempf[0],tempf[1],tempf[2]);
}
std::string swapEndianism;
if(endianism=="big" || endianism=="b")
	swapEndianism="little";
if(endianism=="little" || endianism=="l")
	swapEndianism="big";
if(swapEndianism!="big" && swapEndianism!="little")
{
	std::cerr<<"Swap. Invalid input Endianism table. Operation aborted."<<std::endl;
	return false;
}
tableSwap.setEndiannes(swapEndianism);
std::string type=m_tables[0]->getType();
tableSwap.setType(type);
tableSwap.setNumberOfRows(m_tables[0]->getNumberOfRows());
std::string colNameOutput;
for(unsigned int i=0;i<m_tables[0]->getNumberOfColumns();i++) 
	tableSwap.addCol(m_tables[0]->getColName(i));	
tableSwap.writeHeader();  //overwrite if table exist!

unsigned long long int index=0;
while(index< m_nOfRow)
{	
	m_tables[0]->getColumn(m_colList,m_nOfCol,index,index+m_nOfEle-1,m_fArray);

#ifdef VSBIGENDIAN
  std::string endianismSS="big";
#else	
  std::string endianismSS="little";
#endif
	int comp=endianismSS.compare(swapEndianism);
	if(comp!=0)  // table is automatically swapped in the system native endianism from getColumn method
	{
	  if(tableSwap.getType()=="float" || tableSwap.getType()=="f")
	    for(int j=0;j<m_nOfCol;j++)
		for(int i=0;i<index+m_nOfEle;i++)
			m_fArray[j][i]=floatSwap((char *)(&m_fArray[j][i]));

	  if(tableSwap.getType()=="double" || tableSwap.getType()=="d")
	    for(int j=0;j<m_nOfCol;j++)
		for(int i=0;i<index+m_nOfEle;i++)
			m_fArray[j][i]=doubleSwap((char *)(&m_fArray[j][i]));
	}
	tableSwap.putColumn(m_colList,m_nOfCol,index,index+m_nOfEle-1,m_fArray);

	index+= m_nOfEle;
	if(index+m_nOfEle>= m_nOfRow) 
		m_nOfEle= m_nOfRow-index;
}

return true;
}