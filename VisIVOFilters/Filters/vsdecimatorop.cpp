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
//#include <cstdlib>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <set>
#include <math.h>
#include "vstable.h"
#include "vstableop.h"
#include "vsdecimatorop.h"
#include "visivoutils.h"

const unsigned int VSDecimatorTableOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSDecimatorTableOp::MIN_NUMBER_OF_ROW = 100;

//---------------------------------------------------------------------
VSDecimatorTableOp::VSDecimatorTableOp()
//---------------------------------------------------------------------
{
 m_fArray=NULL;
 m_f1Array=NULL;
 m_nOfRow=0;
 m_nOfCol=0;	
 m_numberOfExtractedValue=0;
 m_ranFromRow=0;
 m_ranToRow=0; 
 m_firstRanSet=true;	
 }
//---------------------------------------------------------------------
VSDecimatorTableOp::~VSDecimatorTableOp()
//---------------------------------------------------------------------
{
	if(m_fArray!=NULL)
 	     for(unsigned int i=0;i<m_nOfCol;i++)
			if(m_fArray[i]!=NULL) delete [] m_fArray[i];
	if(m_fArray!=NULL) delete [] m_fArray;
 	for(unsigned int i=0;i<m_nOfCol;i++)
			if(m_f1Array[i]!=NULL) delete [] m_f1Array[i];
	if(m_f1Array!=NULL) delete [] m_f1Array;
 
 }

//---------------------------------------------------------------------
bool VSDecimatorTableOp::allocatefArray()
//---------------------------------------------------------------------
{
if(m_nOfRow>getMaxNumberInt()) m_nOfRow=(unsigned long long int)getMaxNumberInt();
unsigned long long int tempLL=getMaxNumberInt();
if(((unsigned long long int)m_nOfRow*m_nOfCol)>tempLL) 
	m_nOfRow=(int)tempLL/m_nOfCol;


m_nOfDecimatedRow=m_nOfRow/(m_step+1) +1;
try
{
	m_fArray=new  float*[m_nOfCol];
	m_f1Array=new  float*[m_nOfCol];
}
catch(std::bad_alloc & e)
{
	m_fArray=NULL;
	m_f1Array=NULL;
}
	if(m_fArray==NULL || m_f1Array==NULL )
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
			m_fArray[i] = new float[m_nOfRow];
}
catch(std::bad_alloc & e)
{
	m_fArray[i]=NULL;
}
try
{
			m_f1Array[i] = new float[m_nOfDecimatedRow];
}
catch(std::bad_alloc & e)
{
	m_f1Array[i]=NULL;
}
 			if(m_fArray[i]==NULL ||  m_f1Array[i]==NULL ) 
 			{	
				goodAllocation=false;
				for(unsigned int j=0;j<i;j++) 
				{	
					delete [] m_fArray[j];
					delete [] m_f1Array[j];
				}
				if(m_nOfRow==MIN_NUMBER_OF_ROW)
				{ 
					delete [] m_fArray;
					delete [] m_f1Array;
					m_fArray=NULL;
					m_f1Array=NULL;
					return false;
				}
				if(m_nOfRow<=MAX_NUMBER_TO_REDUCE_ROW)
					m_nOfRow=MIN_NUMBER_OF_ROW;
				else
					m_nOfRow=m_nOfRow-MAX_NUMBER_TO_REDUCE_ROW;

				m_nOfDecimatedRow=m_nOfRow/(m_step+1) +1;
				break;
 			}
		}
	}

	computeDecimatedPoints();
	return true;
}
//---------------------------------------------------------------------
void VSDecimatorTableOp::printHelp()
//---------------------------------------------------------------------
{ 
std::cout<<"Create a decimator subset from the original data table."<<std::endl<<std::endl;
std::cout<<"Usage: VisIVOFilters --op decimator --skip step [--list parameters] [--iseed iseed] [--out filename_out.bin] [--help] [--file] inputFile.bin"<<std::endl;

std::cout<<"Example: VisIVOFilters --op decimator --skip 5 --iseed 1 --out filename_out.bin --file inputFile.bin"<<std::endl;

std::cout<<"Note: "<<std::endl;
std::cout<<"--skip An integer that represent the number of elements to skip"<<std::endl;
std::cout<<"--perc OBSOLETE: Percentage (from 0.0 to 50.0) of the input file obtained in the output file."<<std::endl;
std::cout<<"--list Valid columns names of the input table. Default: all columns are included"<<std::endl;
std::cout<<"--iseed Specify the seed for the random generation. Default value 0."<<std::endl;
std::cout<<"--out Output table filename. Default name is given."<<std::endl;
std::cout<<"--file  Input table filename."<<std::endl;
std::cout<<"--help produce this output "<<std::endl;

return;
 }
//---------------------------------------------------------------------
void VSDecimatorTableOp::computeDecimatedPoints()
//---------------------------------------------------------------------
// compute m_listElements
{
double param, fractpart, intpart;
unsigned long long int intpart1;
param = (double) m_tables[0]->getNumberOfRows();
param=param/m_nOfRow;
fractpart = modf (param , &intpart);
intpart1=(unsigned long long int) intpart;
m_listElements=ceil((double)m_nOfRow/(double)(m_step+1))*intpart1;
unsigned long long int residual=   m_tables[0]->getNumberOfRows()-intpart1*m_nOfRow;
m_listElements=m_listElements+ceil((double)residual/(double)(m_step+1));


}
//---------------------------------------------------------------------
bool VSDecimatorTableOp::subset(unsigned long long int fromRow,unsigned long long int toRow)
//---------------------------------------------------------------------
{

int index;
unsigned long long int numberOfElements=toRow-fromRow+1;

m_valueToExtract=0;
for(unsigned long long int i=0;i<numberOfElements;i=i+(m_step+1))
{
		for(int j=0;j< m_nOfCol;j++)
			m_f1Array[j][m_valueToExtract]=m_fArray[j][i];
		m_valueToExtract++;
}

m_numberOfExtractedValue=m_numberOfExtractedValue+m_valueToExtract;
if(m_numberOfExtractedValue>m_listElements)
	m_numberOfExtractedValue=m_listElements;
return true;
}
//---------------------------------------------------------------------
bool VSDecimatorTableOp::execute()
//---------------------------------------------------------------------
{
std::set<unsigned int> colNumberSet;
std::set<unsigned int>::iterator iter;
bool percUsed=false;

if(getParameterAsString("skip").empty() || getParameterAsString("skip")=="unknown")
{
	if(!(getParameterAsString("skip").empty() || getParameterAsString("skip")=="unknown"))
	{
		std::cerr<<"Obsolete option perc used"<<std::endl;
		percUsed=true;
	} else
	{
 		std::cerr<<"vsdecimartorableop: No skip value is given"<<std::endl;
 		return false;
	}
}


if(getParameterAsString("list").empty() || getParameterAsString("list")=="unknown" )
{
 for(unsigned int i = 0; i <m_tables[0]->getNumberOfColumns(); ++i)
  colNumberSet.insert(i);
} else
{
	std::stringstream ssListparameters;
	ssListparameters.str(getParameterAsString("list"));
	std::string paramField="";
	std::string paramFieldGlobal;
	while (!ssListparameters.eof())
	{
		ssListparameters>>paramField;
		if(m_tables[0] -> getColId(paramField)>=0)
			colNumberSet.insert(m_tables[0] -> getColId(paramField));

	}
}

m_nOfCol = colNumberSet.size();
/***** Allocate array  */
if(percUsed)
{
	float percentage=getParameterAsFloat("perc");

	if(percentage <0.0) percentage=0.0;
	if(percentage >95.0)
	{ 
  	std::cerr<<"Forced perc=95.0"<<std::endl;
  	percentage=95.0;
	}
	m_step=(int) (100.0/percentage) -1;
} else
	m_step=getParameterAsInt("skip");	

if(m_step<1) m_step=1;
m_nOfRow=m_tables[0]->getNumberOfRows();



unsigned int *colList=NULL;
try
{
colList = new unsigned int[m_nOfCol];
}
catch(std::bad_alloc & e)
{
	colList=NULL;
}

if(colList==NULL)
{
	std::cerr<<"Bad colList allocation. Decimator terminated"<<std::endl;
	return false;
}
int count=0;
for(iter = colNumberSet.begin(); iter != colNumberSet.end(); iter++) 	
{
	colList[count]=*iter;
	count++;
}
bool allocationfArray=allocatefArray();
unsigned int maxEle=(unsigned int) m_nOfRow;
if(m_fArray==NULL || !allocationfArray )	
{
	std::cerr<<"Failed fArray allocation. Select Field terminated"<<std::endl;
	delete [] colList;
	return false;
}

/*****END Allocate array  */


std::stringstream fileNameOutputSStream;
fileNameOutputSStream<<getParameterAsString("out");
std::string fileNameOutput;

if(fileNameOutputSStream.str()=="")
{
  std::string filenameInputTable=m_tables[0]->getLocator();
  int len=filenameInputTable.length();
  fileNameOutputSStream<<filenameInputTable.substr(0, len-4)<<"_samp_"<<getParameterAsString("perc")<<".bin";
}

fileNameOutput=fileNameOutputSStream.str();
if(fileNameOutput.find(".bin") == std::string::npos)
   fileNameOutput.append(".bin");

m_realOutFilename.push_back(fileNameOutput);

//Clean existing tab
remove(fileNameOutput.c_str());
/*** Write  new header SubTable ***/

VSTable percTable;
percTable.setLocator(fileNameOutput);
percTable.setType("float");
percTable.setNumberOfRows(m_listElements);
percTable.setIsVolume(m_tables[0]->getIsVolume());
if(m_tables[0]->getIsVolume())
{
   percTable.setCellNumber(m_tables[0]->getCellNumber()[0],
                           m_tables[0]->getCellNumber()[1],
                           m_tables[0]->getCellNumber()[2]);
   percTable.setCellSize(m_tables[0]->getCellSize()[0],
                           m_tables[0]->getCellSize()[1],
                           m_tables[0]->getCellSize()[2]);
}

percTable.setEndiannes(m_tables[0]->getEndiannes());

for(iter = colNumberSet.begin(); iter != colNumberSet.end(); iter++) 	
        percTable.addCol(m_tables[0]->getColName(*iter));


percTable.setNumberOfRows(m_listElements);
percTable.writeHeader();


/*** Start decimator and write binary file**/

unsigned long long int totRows, fromRow, toRow, startCounter=0;
totRows=m_tables[0]->getNumberOfRows();
bool firstin=true;
while(true)
{
	fromRow=startCounter;
	toRow=fromRow+maxEle-1;
	if(toRow>totRows-1)toRow=totRows-1;
	m_tables[0]->getColumn(colList, m_nOfCol, fromRow, toRow, m_fArray);
	bool goodRandomize=subset(fromRow,toRow);
	if(goodRandomize)
	{
		if(firstin)
		{
			m_ranFromRow=0;
			firstin=false;
		} else
			m_ranFromRow=m_ranToRow+1;
		m_ranToRow=m_numberOfExtractedValue-1;
		percTable.putColumn(colList, m_nOfCol, m_ranFromRow, m_ranToRow, m_f1Array);
	}
	else
	{
		std::cerr<<"Error: randomizer operation failed"<<std::endl;
		return false;
	}
	if(m_numberOfExtractedValue==m_listElements) break;
	startCounter=toRow+1;
	if(toRow==(totRows-1)) break;
}

delete [] colList;
return true;
}