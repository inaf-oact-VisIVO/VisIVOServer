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
#include <string>
#include <sstream>
#include <fstream>
#include <set>
#include "vstable.h"
#include "vstableop.h"
#include "vsrandomizertableop.h"
#include "visivoutils.h"

const unsigned int VSRandomizerTableOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSRandomizerTableOp::MIN_NUMBER_OF_ROW = 100;

//---------------------------------------------------------------------
VSRandomizerTableOp::VSRandomizerTableOp()
//---------------------------------------------------------------------
{
 m_fArray=NULL;
 m_tmpArray=NULL;
 m_nOfRow=0;
 m_nOfCol=0;	
 m_numberOfExtractedValue=0;
 m_ranFromRow=0;
 m_ranToRow=0; 
 m_firstRanSet=true;	
 }
//---------------------------------------------------------------------
VSRandomizerTableOp::~VSRandomizerTableOp()
//---------------------------------------------------------------------
{	if(m_fArray!=NULL)
 		for(unsigned int i=0;i<m_nOfCol;i++)
			if(m_fArray[i]!=NULL) delete [] m_fArray[i];
	if(m_fArray!=NULL) delete [] m_fArray;
	if(m_tmpArray!=NULL) delete [] m_tmpArray;
 
 }

//---------------------------------------------------------------------
bool VSRandomizerTableOp::allocatefArray()
//---------------------------------------------------------------------
{
unsigned long long int tempLL=getMaxNumberInt();
if(((unsigned long long int)m_nOfRow*m_nOfCol)>tempLL) 
	m_nOfRow=(int)tempLL/m_nOfCol;

try
{
	m_tmpArray=new  float[m_nOfCol];
}
catch(std::bad_alloc & e)
{
	m_tmpArray=NULL;
}
	if(m_tmpArray==NULL)
	{
		return false;
	}

try
{
	m_fArray=new  float*[m_nOfCol];
}
catch(std::bad_alloc & e)
{
	m_fArray=NULL;
}
	if(m_fArray==NULL)
	{
		return false;
	}
	if(m_nOfRow>getMaxNumberInt()) m_nOfRow=(unsigned long long int)getMaxNumberInt();

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
 			if(m_fArray[i]==NULL) 
 			{	
				goodAllocation=false;
				for(unsigned int j=0;j<i;j++) 
				{	
					delete [] m_fArray[j];
				}
				if(m_nOfRow==MIN_NUMBER_OF_ROW)
				{ 
					delete [] m_fArray;
					m_fArray=NULL;
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
void VSRandomizerTableOp::printHelp()
//---------------------------------------------------------------------
{ 
if(m_MpiRank==0)
{
std::cout<<"Create a random subset from the original data table."<<std::endl<<std::endl;
std::cout<<"Usage: VisIVOFilters --op randomizer --perc percentage [--field parameters] [--iseed iseed] [--out filename_out.bin] [--history] [--historyfile filename.xml] [--help] [--file] inputFile.bin"<<std::endl;

std::cout<<"Example: VisIVOFilters --op randomizer --perc 10.0 --iseed 1 --out filename_out.bin --file inputFile.bin"<<std::endl;

std::cout<<"Note: "<<std::endl;
std::cout<<"--perc Percentage (from 0.0 to 100.0) of the input file obtained in the output file."<<std::endl;
std::cout<<"--field Valid columns names of the input table. Default: all columns are included"<<std::endl;
std::cout<<"--iseed Specify the seed for the random generation. Default value 0."<<std::endl;
std::cout<<"--out Output table filename. Default name is given."<<std::endl;
std::cout<<"--file  Input table filename."<<std::endl;
std::cout<<"--history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
std::cout<<"--historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;
std::cout<<"--help produce this output "<<std::endl;

}
return;
}
//---------------------------------------------------------------------
bool VSRandomizerTableOp::subset(unsigned long long int fromRow,unsigned long long int toRow)
//---------------------------------------------------------------------
{
int iseed;
if(m_firstRanSet)
{	
	iseed=getParameterAsInt("iseed"); //QUI modifica
	m_firstRanSet=false;
}else
	iseed=m_lastRandomValue;

int iran,rran;
srand(iseed);
int index;
unsigned long long int numberOfElements=toRow-fromRow;
double tmp;
tmp=((double)(toRow-fromRow))/(double)(m_tables[0]->getNumberOfRows());
m_valueToExtract=tmp*m_listElement;
if(toRow==(m_tables[0]->getNumberOfRows()-1)) 
	m_valueToExtract=m_listElement-m_numberOfExtractedValue;
if((m_numberOfExtractedValue+m_valueToExtract)>m_listElement) 
	m_valueToExtract=m_listElement-m_numberOfExtractedValue;
for(unsigned long long int i=0;i<m_valueToExtract;i++)
{
	iran=rand();
        index= (int)((numberOfElements-i)*((float)(iran) /(float) RAND_MAX))+i;
	for(unsigned int j=0;j<m_nOfCol;j++)
		m_tmpArray[j]=m_fArray[j][i];
	for(unsigned int j=0;j<m_nOfCol;j++)
		m_fArray[j][i]=m_fArray[j][index];
	for(unsigned int j=0;j<m_nOfCol;j++)
		m_fArray[j][index]=m_tmpArray[j];

}
m_lastRandomValue=iran;
m_numberOfExtractedValue=m_numberOfExtractedValue+m_valueToExtract;
return true;
}
//---------------------------------------------------------------------
bool VSRandomizerTableOp::execute()
//---------------------------------------------------------------------
{
std::set<unsigned int> colNumberSet;
std::set<unsigned int>::iterator iter;

if(getParameterAsString("perc").empty() || getParameterAsString("perc")=="unknown")
{
std::cerr<<"vsrandomizertableop: No perc value is given"<<std::endl;
 return false;
}

if(isParameterPresent("field"))
{
  if(getParameterAsString("field").empty() || getParameterAsString("field")=="unknown" )
  {
    for(unsigned int i = 0; i <m_tables[0]->getNumberOfColumns(); ++i)
	colNumberSet.insert(i);
  } else
  {
	std::stringstream ssListparameters;
	ssListparameters.str(getParameterAsString("field"));
	std::string paramField="";
	std::string paramFieldGlobal;
	while (!ssListparameters.eof())
	{
		ssListparameters>>paramField;
		if(m_tables[0] -> getColId(paramField)>=0)
			colNumberSet.insert(m_tables[0] -> getColId(paramField));

	}
  }


} else  //obsolete --list
{

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

}

m_nOfCol = colNumberSet.size();
/***** Allocate array  */
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
	std::cerr<<"Bad colList allocation. Randomizer terminated"<<std::endl;
	return false;
}
int count=0;
for(iter = colNumberSet.begin(); iter != colNumberSet.end(); iter++) 	
{
	colList[count]=*iter;
	count++;
}
m_nOfRow=m_tables[0]->getNumberOfRows();;
bool allocationfArray=allocatefArray();
unsigned int maxEle=(unsigned int) m_nOfRow;
if(m_fArray==NULL || !allocationfArray )	
{
	std::cerr<<"Failed fArray allocation. Select Field terminated"<<std::endl;
	delete [] colList;
	return false;
}

/***** Allocate array  */
float percentage=getParameterAsFloat("perc");

if(percentage <0.0) percentage=0.0;
if(percentage >95.0)
{ 
  std::cerr<<"Forced perc=95.0"<<std::endl;
  percentage=95.0;
}


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
//Clean existing tab
remove(fileNameOutput.c_str());
/*** Write  new header SubTable ***/


m_realOutFilename.push_back(fileNameOutput);
m_listElement=(unsigned long long int) (m_tables[0]->getNumberOfRows()*percentage/100.);
VSTable percTable;
percTable.setLocator(fileNameOutput);
percTable.setType("float");
percTable.setNumberOfRows(m_listElement);
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
			firstin=false;
			m_ranFromRow=0;
		}
		else
			m_ranFromRow=m_ranToRow+1;
		m_ranToRow=m_ranFromRow+m_valueToExtract-1;
		percTable.putColumn(colList, m_nOfCol, m_ranFromRow, m_ranToRow, m_fArray);
	}
	else
	{
		std::cerr<<"Error: randomizer operation failed"<<std::endl;
		return false;
	}
	startCounter=toRow+1;
	if(toRow==(totRows-1)) break;
}

delete [] colList;
return true;
}