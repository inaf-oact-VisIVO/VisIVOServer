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
#include <math.h>
#ifdef WIN32
	#include <time.h>
#endif
#include "vstable.h"
#include "vssplittableop.h"
#include "vsstatisticop.h"
#include "vsselectfieldop.h"
#include "vsextractsubvolumeop.h"

const unsigned int VSSplitTableOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSSplitTableOp::MIN_NUMBER_OF_ROW = 100;


//---------------------------------------------------------------------
VSSplitTableOp::VSSplitTableOp()
//---------------------------------------------------------------------
{
 m_fArray=NULL;
 m_nOfRows=0;
 m_nOfCols=0;	
 m_numOfTables=0;
}
//---------------------------------------------------------------------
VSSplitTableOp::~VSSplitTableOp()
//---------------------------------------------------------------------
{
	if(m_fArray!=NULL)	
 	  for(unsigned int i=0;i<m_nOfCols;i++)
		{
			if(m_fArray[i]!=NULL) delete [] m_fArray[i];
		}
	if(m_fArray!=NULL) delete [] m_fArray;
}

//---------------------------------------------------------------------
bool VSSplitTableOp::allocatefArray()
//---------------------------------------------------------------------
{
unsigned long long int tempLL=getMaxNumberInt();
if(((unsigned long long int)m_nOfLocalRow*m_nOfCols)>tempLL) 
	m_nOfLocalRow=(int)tempLL/m_nOfCols;


try
{
	m_fArray=new  float*[m_nOfCols];
}
catch(std::bad_alloc & e)
{
	m_fArray=NULL;
}
	if(m_fArray==NULL)
	{
		return false;
	}

	bool goodAllocation=false;
	while(!goodAllocation)
	{
		goodAllocation=true;
		for(unsigned int i=0;i<m_nOfCols;i++)
		{
try
{
			m_fArray[i] = new float[m_nOfLocalRow];
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
				if(m_nOfLocalRow==MIN_NUMBER_OF_ROW)
				{ 
					delete [] m_fArray;
					m_fArray=NULL;
					return false;
				}
				if(m_nOfLocalRow<=MAX_NUMBER_TO_REDUCE_ROW)
					m_nOfLocalRow=MIN_NUMBER_OF_ROW;
				else
					m_nOfLocalRow=m_nOfLocalRow-MAX_NUMBER_TO_REDUCE_ROW;
				break;
 			}
		}
	}
	return true;
}


//---------------------------------------------------------------------
void VSSplitTableOp::printHelp()
//---------------------------------------------------------------------
{
	std::cout<<"Split an existing table into two or more tables, using a field that will be used to divide the table."<<std::endl<<std::endl;
	std::cout<<"Usage: VisIVOFilters --op splittable [--field column] [--volumesplit direction]  [--numoftables numberOfTable] [--maxsizetable MaxMbyteSize] [--hugesplit] [--out filename_out.bin] [--history] [--historyfile filename.xml] [--help] [--file] inputFile.bin"<<std::endl;

	std::cout<<"Example: VisIVOFilters --op splittable --field X --numoftables 10  --out filename_out --file inputFile.bin"<<std::endl<<std::endl;

	std::cout<<"This command will compute the X range (e.g range 0-100) and will produce 10 tables where the X range will be: 0-10, 10-20 ,20-30 and so on."<<std::endl<<std::endl;

	std::cout<<"Note:"<<std::endl;
	std::cout<<"--filed A valid column name along which the table will be split. Must be given to split a table."<<std::endl;
	std::cout<<"--volumesplit Direction (1,2 or3) along which the volume will be split. Must be given to split a volume"<<std::endl;
	std::cout<<"--numoftables gives the number of tables in which the input table will be split. It must be greater than 1."<<std::endl;
	std::cout<<"--maxsizetable indicates the maximum size of the split table. VisIVO Filter will compute how many tables will be created. This option is ignored if  --numoftable option is given. "<<std::endl;
	std::cout<<"--hugesplit. Must be given to force the generation of more than 100 tables from the input table, avoiding errors."<<std::endl;
	std::cout<<"--out Output prefix root table filename. A suffix _split_#.bin is added to each generated table,  to this prefix. Default name is given."<<std::endl;
    std::cout<<"--history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
    std::cout<<"--historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;
    std::cout<<"--file  Input table filename."<<std::endl;
	std::cout<<"--help produce this output "<<std::endl;


	return;
}
//---------------------------------------------------------------------
bool VSSplitTableOp::writeTable()
//---------------------------------------------------------------------
{
time_t rawtime;
struct tm * timeinfo;
char buffer [80];
time ( &rawtime );
timeinfo = localtime ( &rawtime );
strftime (buffer,80,"%Y%m%d%H%M",timeinfo);
std::stringstream  fileNameLimitsSStream;	
fileNameLimitsSStream<<"__VSsplitop_lim"<<buffer;  //QUI verificare
std::string fileLimit;
fileLimit=fileNameLimitsSStream.str();
std::ofstream limitFile(fileLimit.c_str());
limitFile<<m_field<<" "<<m_lowLimit<<" "<<m_upLimit<<std::endl;
limitFile.close();
VSSelectFieldOp op;
op.addInput(m_tables[0]);
std::stringstream locator;
locator<<m_fileNameOutput<<"_split_"<<m_tableNumber<<".bin";	
op.addParameter("out",locator.str());
op.addParameter("limits",fileLimit);
if(!op.execute())
{ 
	remove(fileLimit.c_str());
	return false;
}
m_realOutFilename.push_back(locator.str());
remove(fileLimit.c_str());
return true;
}
//---------------------------------------------------------------------
bool VSSplitTableOp::splitPoints()
//---------------------------------------------------------------------
{
int numberOfBin=10000;
unsigned long long int **histogramPointer;

float max,min;

VSStatisticOp op;
op.addParameter("list",m_field);
op.addParameter("silent","");
std::stringstream tmp;
tmp<<numberOfBin;
op.addParameter("histogram",tmp.str());
op.addInput(m_tables[0]);
if(!op.execute())
{
	std::cerr<<"Invalid statisticop. Operation aborted"<<std::endl;
	return false;
}
float dummy;
op.getRange(m_fieldId,max,min,dummy,dummy);
histogramPointer=op.getHisto();

unsigned long long int rowForTable=m_nOfRows/(unsigned long long int) m_numOfTables;
float deltaStep=(max-min)/numberOfBin;
m_lowLimit=min;
m_upLimit=min;
unsigned long long int tmp1=0;
unsigned long long int tmp2=0;
m_tableNumber=0;
for(int i=0;i<=numberOfBin;i++)
{
	if(m_tableNumber==m_numOfTables) break;
	tmp1+=histogramPointer[0][i];
	while(true)
	{
	 if(tmp1>=rowForTable)
	 {
		unsigned long long int diff=rowForTable-tmp2;
		float perc=((float)diff)/((float)histogramPointer[0][i]);
		m_upLimit+=deltaStep*perc;
		histogramPointer[0][i]-=diff;
		if(m_tableNumber==m_numOfTables-1) m_upLimit=max;
		if(!writeTable())
			return false;
		m_tableNumber++;
		if(m_tableNumber==m_numOfTables) break;
		m_lowLimit=m_upLimit+0.000000001;
	  	tmp1=histogramPointer[0][i];
		tmp2=0;
	 } else
	 {
		m_upLimit+=deltaStep;
		tmp2=tmp1;
		break; //next localHistogram value
	 }		 
	}
}

}
//---------------------------------------------------------------------
bool VSSplitTableOp::splitVolume()
//---------------------------------------------------------------------
{
const unsigned int *temp1;
unsigned int *temp;
const unsigned int *globalCells;
temp=new unsigned int[3];
temp1=m_tables[0]->getCellNumber();
for(int i=0;i<3;i++) temp[i]=temp1[i];
globalCells=m_tables[0]->getCellNumber();
if(temp[m_vDir-1]<m_numOfTables)
{
	std::cerr<<"Warning. Requested number of tables "<<m_numOfTables<<" is greather than the cell in the requested split direction "<<m_vDir<<" The number of tables will be "<<temp[m_vDir-1]<<std::endl; 
	m_numOfTables=temp[m_vDir-1];
}
int rowForTable=temp[m_vDir-1]/m_numOfTables;
temp[m_vDir-1]=rowForTable;
unsigned int countTotalRow=0;
for(unsigned int i=0;i<m_numOfTables;i++)
{
	std::stringstream globalCellsSS;
	globalCellsSS<<temp[0]<<" "<<temp[1]<<" "<<temp[2];
	std::stringstream startingCell;
	if(m_vDir==1) startingCell<<countTotalRow<<" 0 0";
	if(m_vDir==2) startingCell<<"0 "<<countTotalRow<<" 0";
	if(m_vDir==3) startingCell<<"0 0 "<<countTotalRow;
	VSExtractSubVolumeOp tempVolume;
	tempVolume.addParameter("startingcell",startingCell.str());
	tempVolume.addParameter("resolution",globalCellsSS.str());
	tempVolume.addParameter("out",globalCellsSS.str());
	std::stringstream locator;
	locator<<m_fileNameOutput<<"_split_"<<i<<".bin";	
	tempVolume.addParameter("out",locator.str());
	tempVolume.addInput(m_tables[0]);
	if(!tempVolume.execute())
	{
		std::cerr<<"Invalid execution subvolume. Operation Aborted"<<std::endl;
		delete [] temp;
		return false;
	}
	m_realOutFilename.push_back(locator.str());
	countTotalRow+=rowForTable;
	if(i==m_numOfTables-2)
		temp[m_vDir-1]=globalCells[m_vDir-1]-countTotalRow;
}
delete [] temp;
return true;
}
//---------------------------------------------------------------------
bool VSSplitTableOp::fastTableSplit()
//---------------------------------------------------------------------
{
//VSTable tabSplitted[m_numOfTables];
VSTable * tabSplitted = new VSTable[m_numOfTables];
	if(!allocatefArray())
{
  std::cerr<<"Not space available. Operation aborted"<<std::endl;
  return false;
}
unsigned int *colList=NULL;
try
{
colList = new unsigned int[m_nOfCols];
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
for(unsigned int i=0;i<m_nOfCols;i++) 
  colList[i]=i;
  unsigned long long int totRows, fromRow, toRow, startCounter=0;

for(unsigned int i=0;i<m_numOfTables;i++)
{
  unsigned long long int localRows, res;
  localRows=m_nOfRows/m_numOfTables;
  res=m_nOfRows-localRows*m_numOfTables;
  if(i<res) localRows++;
  
  std::stringstream locator;
  locator<<m_fileNameOutput<<"_split_"<<i<<".bin";	
  
  tabSplitted[i].setLocator(locator.str());
  tabSplitted[i].setType("float");
  tabSplitted[i].setNumberOfRows(localRows);
  tabSplitted[i].setEndiannes(m_tables[0]->getEndiannes());
  for(unsigned int j=0; j< m_nOfCols;j++)
      tabSplitted[i].addCol(m_tables[0]->getColName(j));
  tabSplitted[i].writeHeader();
  
  unsigned long long int localtoRows, localfromRow=0,writedRows=0, nOfLocalRow=m_nOfLocalRow;
  totRows=localRows;
  bool firstin=true;
  while(true)
  {
	fromRow=startCounter;
	toRow=fromRow+nOfLocalRow-1;
	if((toRow-fromRow)>totRows-1)toRow=totRows+fromRow-1;
	localtoRows=localfromRow+(toRow-fromRow);
	m_tables[0]->getColumn(colList, m_nOfCols, fromRow, toRow, m_fArray);
	tabSplitted[i].putColumn(colList, m_nOfCols, localfromRow, localtoRows, m_fArray);
	if(toRow==(m_nOfRows-1)) break;
	writedRows=writedRows+(localtoRows-localfromRow+1);
	startCounter=toRow+1;
	localfromRow=localtoRows+1;
	if(writedRows==localRows) 
	  break;
	if((localRows-writedRows)<m_nOfLocalRow) nOfLocalRow=localRows-writedRows;
  }
  if(toRow==(m_nOfRows-1)) break;


}


delete [] colList;
return true;

}
//---------------------------------------------------------------------
bool VSSplitTableOp::execute()
//---------------------------------------------------------------------
{
bool fastSplit=false;
if(!m_tables[0]->getIsVolume() && !isParameterPresent("field"))
{
//	std::cerr<<"Operation aborted. A valid field must be given."<<std::endl;
//	return false;
	fastSplit=true; //QUI
}
if(m_tables[0]->getIsVolume() && !isParameterPresent("volumesplit"))
{
	std::cerr<<"Operation aborted. A valid volumesplit must be given."<<std::endl;
	return false;
}
if(!m_tables[0]->getIsVolume() && isParameterPresent("field"))
{
	m_field=getParameterAsString("field");
	m_fieldId= m_tables[0]->getColId(m_field);;
	if(m_fieldId==-1)
	{
		std::cerr<<"Operation Aborted. Invalid field name "<<m_field<<std::endl;
		return false;
	} 
} else if(m_tables[0]->getIsVolume())
{
	m_vDir=getParameterAsInt("volumesplit");
	if(m_vDir<1 || m_vDir>3)
	{
		std::cerr<<"Operation Aborted. Invalid volumesplit option"<<m_vDir<<std::endl;
		return false;
	} 
}
float maxMbyteTable=0;	
if(isParameterPresent("numoftables"))
{
	m_numOfTables=getParameterAsInt("numoftables");
	if(m_numOfTables<=1)
	{
		std::cerr<<"An invalid number of tables is given: "<<m_numOfTables<<" .It must be given more than 1 tables."<<std::endl;
		return false;
	}
	if((m_numOfTables>100) && !isParameterPresent("hugesplit"))
	{
		std::cerr<<"An high  number of splitted tables is given: "<<m_numOfTables<<" .It must be forced with --hugesplit parameter."<<std::endl;
		return false;

	}
} else
{
	if(isParameterPresent("maxsizetable"))
	{
		maxMbyteTable=getParameterAsFloat("maxsizetable");
		if(maxMbyteTable<0)
		{
			std::cerr<<"An invalid maxsizetable parameter is given: "<<maxMbyteTable<<" .It must be given more than 0 MBytes."<<std::endl;
			return false;
			
		}
	}
}
if(maxMbyteTable>0)
{
	std::string filename=m_tables[0]->getLocator();
	std::ifstream fileInput(filename.c_str());
	fileInput.seekg (0, std::ios::end);
	unsigned long long int length=fileInput.tellg();
	length=length/(1024*1024); //MegaByte
	m_numOfTables=ceil((float)length/maxMbyteTable);
	fileInput.close();
	if((m_numOfTables>100) && !isParameterPresent("hugesplit"))
	{
		std::cerr<<"An high  number of splitted tables will occurr: "<<m_numOfTables<<" .It must be forced with --hugesplit parameter."<<std::endl;
		return false;
	}
}
if(m_numOfTables<=1)
{
	std::cerr<<"Invalid parameters. Operation aborted"<<std::endl;
	return false;
}
m_nOfRows=m_tables[0]->getNumberOfRows();
m_nOfCols=m_tables[0]->getNumberOfColumns();
m_nOfLocalRow=m_nOfRows/m_numOfTables + 1;
std::stringstream fileNameOutputSStream;
fileNameOutputSStream<<getParameterAsString("out");

if(fileNameOutputSStream.str()==""||fileNameOutputSStream.str()=="unknown")
{
	fileNameOutputSStream.str().erase(); //QUI verificare
	std::string filenameInputTable=m_tables[0]->getLocator();
 	int len=filenameInputTable.length();
	fileNameOutputSStream<<filenameInputTable.substr(0, len-4);  //QUI verificare
        m_fileNameOutput=fileNameOutputSStream.str();
} else
{
      std::string filenameOutTable=fileNameOutputSStream.str();
      if(filenameOutTable.find(".bin") != std::string::npos)
      {	
	int len=filenameOutTable.length();
	m_fileNameOutput=filenameOutTable.substr(0, len-4);
      } else
	m_fileNameOutput=filenameOutTable;
}

if(m_tables[0]->getIsVolume())
	splitVolume();
else if(fastSplit)
       fastTableSplit();
else
	splitPoints();

return true;

  
}

