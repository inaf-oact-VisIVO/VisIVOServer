//
// C++ Implementation: vsahfhalolistop
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
#include <fstream>
#include <sstream>
#ifdef WIN32
	#include <time.h>
#endif
#include "vsahfhalolistop.h"
#include "vstable.h"

const unsigned int VSAhfHaloListOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSAhfHaloListOp::MIN_NUMBER_OF_ROW = 100;

//---------------------------------------------------------------------
VSAhfHaloListOp::VSAhfHaloListOp()
{ 
m_fArray=NULL;
m_nOfRow=0;
}
//---------------------------------------------------------------------

//---------------------------------------------------------------------
VSAhfHaloListOp::~VSAhfHaloListOp()
{
if(m_fArray!=NULL ) 
	if(m_fArray[0] != NULL) delete [] m_fArray[0];

if(m_fArray!=NULL) delete [] m_fArray;

}
//---------------------------------------------------------------------
//---------------------------------------------------------------------
void VSAhfHaloListOp::printHelp()
{
	std::cout<<"Specific filetr for Amiga Halo Finder."<<std::endl;
	std::cout<<"It extracts an ascii mutlist from a general input ascii multilist, following the Id in a given input table."<<std::endl<<std::endl;
	std::cout<<"Usage: VisIVOFilters --op ahfhalolist [--field col_name] --multilist in_mutlist_File [--out out_multilist_file] [--help] [--file] inputFile.bin"<<std::endl;

	std::cout<<"Example: VisIVOFilters --op ahfhalolist --field Id --multilist input.AHF_particles --out my_multilist.lst  --file sub_AHF_halos.bin"<<std::endl;

	std::cout<<"Note:"<<std::endl;
	std::cout<<"--field. Column name with contain id column. Default name is Id."<<std::endl;
	std::cout<<"--multilist. Ascii filename with a multilist."<<std::endl;
	std::cout<<"--out. Output ascii filename with a multilist."<<std::endl;
	std::cout<<"--file  Input table filename."<<std::endl;
	std::cout<<"--help produce this output "<<std::endl;


	return;

}
//---------------------------------------------------------------------
//---------------------------------------------------------------------
bool VSAhfHaloListOp::allocatefArray()
//---------------------------------------------------------------------
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
bool VSAhfHaloListOp::createNewList()
//---------------------------------------------------------------------
{
int dummy;
m_fileInput.seekg(0,std::ios::beg);
m_fileInput>> dummy;
for(int i=0;i<m_nOfEle;i++)
{
	if(m_fArray[0][i]>dummy || m_fArray[0][i]<=0)
	{
		std::cerr<<"SEVERE WARNING: Ignored list "<< m_fArray[0][i]<<std::endl;
		m_nOfRealList--;
		std::cerr<<"Please correct manualy the number of list in output file to "<< m_nOfRealList<<std::endl;
		continue;
	}
	int skip=(int) m_fArray[0][i]-1;
	for(int k=0;k<skip;k++)
	{
		int elements;
		m_fileInput>> elements;
		for(int j=0;j<elements;j++)
			m_fileInput>> dummy;
	}
	int totElements;
	int elements;
	m_fileInput>> totElements;
	m_fileOutput << totElements<<std::endl;
	for(int k=0;k<totElements;k++)
	{
		m_fileInput>> elements;
		m_fileOutput << elements<<std::endl;
	}
	m_fileInput.seekg(0,std::ios::beg);
	m_fileInput>> dummy;
}
return true;
}
//---------------------------------------------------------------------
bool VSAhfHaloListOp::execute()
{
if(!isParameterPresent("multilist"))
{
	std::cerr<<"No input option is given. Operation aborted."<<std::endl;
	return false;
}
std::string filename;
filename=getParameterAsString("multilist");

m_fileInput.open(filename.c_str());
if(!m_fileInput.is_open())
{
	std::cerr<<"Error opening input multilist file "<<filename<<" Operation aborted."<<std::endl;
	return false;
} 
std::stringstream fileNameOutputSStream;
std::string fileNameOutput;
if(!isParameterPresent("out") ||getParameterAsString("out").empty() || getParameterAsString("out")=="unknown")
{
  	std::string filenameInputTable=filename;
  	int len=filenameInputTable.length();
	time_t rawtime;
	struct tm * timeinfo;
	char buffer [80];
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	strftime (buffer,80,"%Y%m%d%H%M",timeinfo);
  	fileNameOutputSStream<<filenameInputTable.substr(0, len-4)<<"_AHFListExtraction_"<<buffer;  //QUI verificare
} else
	fileNameOutputSStream<<getParameterAsString("out");

fileNameOutput=fileNameOutputSStream.str();

m_fileOutput.open(fileNameOutput.c_str());
if(!m_fileOutput.is_open())
{
	std::cerr<<"Error opening output multilist file "<<fileNameOutput<<" Operation aborted."<<std::endl;
	return false;
} 


std::string colId;
if(getParameterAsString("field").empty() ||getParameterAsString("field")=="unknown")
  	colId="Id";
else
	colId=getParameterAsString("field");

if((m_tables[0]->getColId(colId))<0)
{
	std::cerr<<"Invalid field name "<<colId<<". The column NOT exist. Operation Aborted."<<std::endl;
	return false;
}

m_nOfRow=m_tables[0]->getNumberOfRows();
m_nOfCol=1;
if(m_nOfRow>getMaxNumberInt()) 
	m_nOfEle= getMaxNumberInt();
else 
	m_nOfEle= m_nOfRow;

if(!allocatefArray())
	return false;
unsigned int * colList;
colList=new unsigned int[1];
colList[0]=m_tables[0]->getColId(colId);

unsigned long long int index=0;
m_nOfRealList=m_nOfRow;
m_fileOutput<<m_nOfRow<<std::endl;
while(index< m_nOfRow)
{
	m_tables[0]->getColumn(colList,1,index,index+m_nOfEle-1,m_fArray);
	
	if(!createNewList())	
		return false;

	index+= m_nOfEle;
	if(index+m_nOfEle>= m_nOfRow) 
		m_nOfEle= m_nOfRow-index;

}
m_fileInput.close();
m_fileOutput.close();
return true;
}
//---------------------------------------------------------------------


