//
// C++ Implementation: vsaddidop
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

#include "vsaddidop.h"
#include "vstable.h"
const unsigned int VSAddIdOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSAddIdOp::MIN_NUMBER_OF_ROW = 100;

//---------------------------------------------------------------------
VSAddIdOp::VSAddIdOp()
{
m_fArray=NULL;
m_nOfRow=0;
}
//---------------------------------------------------------------------


//---------------------------------------------------------------------
VSAddIdOp::~VSAddIdOp()
{
if(m_fArray!=NULL ) 
	if(m_fArray[0] != NULL) delete [] m_fArray[0];

if(m_fArray!=NULL) delete [] m_fArray;
}
//---------------------------------------------------------------------

//---------------------------------------------------------------------
void VSAddIdOp::printHelp()
{
	std::cout<<"Add a new column with a sequence of Ids in the input data table."<<std::endl<<std::endl;
	std::cout<<"Usage: VisIVOFilters --op addId [--outcol col_name] [--start start_number] [--help] [--file] inputFile.bin"<<std::endl;

	std::cout<<"Example: VisIVOFilters --op addId --outcol PtId  --file inputFile.bin"<<std::endl;

	std::cout<<"Note:"<<std::endl;
	std::cout<<"--start Starting Id. Default value is 1. Only an int value can be given."<<std::endl;
	std::cout<<"--outcol. Column name of the new id column. Default name is Id."<<std::endl;
	std::cout<<"--file  Input table filename."<<std::endl;
	std::cout<<"--help produce this output "<<std::endl;


	return;

}
//---------------------------------------------------------------------
//---------------------------------------------------------------------
bool VSAddIdOp::allocatefArray()
{
try
{
m_fArray=new  float*[1];
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
try
{
	m_fArray[0] = new float[m_nOfEle];
}
catch(std::bad_alloc &e)
{
	m_fArray[0]=NULL;
}
 		if(m_fArray[0]==NULL) 
 		{	
			goodAllocation=false;
			if(m_nOfEle==MIN_NUMBER_OF_ROW)
			{ 
				delete [] m_fArray;
				m_fArray=NULL;
				return false;
			}
			if(m_nOfEle<=MAX_NUMBER_TO_REDUCE_ROW)
				m_nOfEle=MIN_NUMBER_OF_ROW;
			else
				m_nOfEle=m_nOfEle-MAX_NUMBER_TO_REDUCE_ROW;
 		}
}
return true;
}
//---------------------------------------------------------------------
//---------------------------------------------------------------------
bool VSAddIdOp::execute()
{
std::string colNameOutput;
long long int offset=0;
if(getParameterAsString("outcol").empty() ||getParameterAsString("outcol")=="unknown")
  		colNameOutput="Id";
	 else
		colNameOutput=getParameterAsString("outcol");
if((m_tables[0]->getColId(colNameOutput))>=0)
{
	std::cerr<<"Invalid outcol name "<<colNameOutput<<". The column already exist. Operation Aborted."<<std::endl;
	return false;
}
if(isParameterPresent("start"))
{
	offset= (long long int) getParameterAsInt("start");
}

m_nOfRow=m_tables[0]->getNumberOfRows();
m_nOfCol=1;
if(m_nOfRow>getMaxNumberInt()) 
	m_nOfEle= getMaxNumberInt();
else 
	m_nOfEle= m_nOfRow;

if(!allocatefArray())
	return false;
VSTable tableId;
tableId.setLocator(m_tables[0]->getLocator());
m_realOutFilename.push_back(m_tables[0]->getLocator());

#ifdef VSBIGENDIAN
	std::string endianism="big";
	
#else	
	std::string endianism="little";
#endif

tableId.setEndiannes(endianism);

tableId.setType("float");
tableId.setNumberOfRows(m_tables[0]->getNumberOfRows());
tableId.setIsVolume(m_tables[0]->getIsVolume());
if(m_tables[0]->getIsVolume())
{
    tableId.setCellNumber(m_tables[0]->getCellNumber()[0],
                          m_tables[0]->getCellNumber()[1],
                          m_tables[0]->getCellNumber()[2]);
    tableId.setCellSize(m_tables[0]->getCellSize()[0],
                           m_tables[0]->getCellSize()[1],
                           m_tables[0]->getCellSize()[2]);
}
for(unsigned int i=0;i<m_tables[0]->getNumberOfColumns();i++)
	tableId.addCol(m_tables[0]->getColName(i));
tableId.addCol(colNameOutput); //append column
tableId.writeHeader();
 
unsigned int * colList;
colList=new unsigned int[1];
colList[0]=tableId.getNumberOfColumns()-1;

unsigned long long int index=0;
while(index< m_nOfRow)
{	
	for(int i=0;i<m_nOfEle;i++)
		m_fArray[0][i]=i+offset;
	offset=m_fArray[0][m_nOfEle-1]+1;
	tableId.putColumn(colList,1,index,index+m_nOfEle-1,m_fArray);
	index+= m_nOfEle;
	if(index+m_nOfEle>= m_nOfRow) 
		m_nOfEle= m_nOfRow-index;
}
delete [] colList;
return true;
}
//---------------------------------------------------------------------

