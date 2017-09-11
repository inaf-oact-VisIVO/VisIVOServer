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
#include "vsahfhalogalaxyextop.h"
#include "vstable.h"

const unsigned int VSAhfHaloGalaxyExtOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSAhfHaloGalaxyExtOp::MIN_NUMBER_OF_ROW = 100;

//---------------------------------------------------------------------
VSAhfHaloGalaxyExtOp::VSAhfHaloGalaxyExtOp()
{ 
m_fArray=NULL;
m_nOfRow=0;
m_threshold=1.0;
}
//---------------------------------------------------------------------

//---------------------------------------------------------------------
VSAhfHaloGalaxyExtOp::~VSAhfHaloGalaxyExtOp()
{
if(m_fArray!=NULL ) 
	if(m_fArray[0] != NULL) delete [] m_fArray[0];

if(m_fArray!=NULL) delete [] m_fArray;

}
//---------------------------------------------------------------------
//---------------------------------------------------------------------
void VSAhfHaloGalaxyExtOp::printHelp()
{
	std::cout<<"Specific filter for Amiga Halo Finder."<<std::endl;
	std::cout<<"It Extracts up to three ascii mutlist from a general input ascii multilist, following the Ids in a given input table, and using the threshold value given in select value."<<std::endl<<std::endl;
	std::cout<<"Usage: VisIVOFilters --op ahfhalogalaxyext --field col1 col2 col3 --multilist in_mutlist_File [--threshold value] [--out out_multilist_file] [--history] [--historyfile filename.xml] [--help] [--file] inputFile.bin"<<std::endl;

	std::cout<<"Example: VisIVOFilters --op ahfhalogalaxyext --field Id1 Id2 Id3 --multilist input.AHF_particles --out my_multilist.lst --file sub_AHF_halos.bin"<<std::endl;

	std::cout<<"Note:"<<std::endl;
	std::cout<<"The lists are extracted: my_multilist.lst_Id1 where Id1 >= 1.0, my_multilist.lst_Id2 where Id2 >= 1.0 and Id1 <1.0, my_multilist.lst_Id3 where Id3 >= 1.0 and Id1 <1.0, Id2<0,from input.AHF_particles"<<std::endl<<std::endl; 
	std::cout<<"--field. Up to three column names. Default name is Id."<<std::endl;
	std::cout<<"--multilist. Ascii filename with a multilist."<<std::endl;
	std::cout<<"--threshold. Selector value among the three columns. Default Value 1.0"<<std::endl;
	std::cout<<"--out. Output ascii root-filename with a multilist."<<std::endl;
    std::cout<<"--history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
    std::cout<<"--historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;
	std::cout<<"--file  Input table filename."<<std::endl;
	std::cout<<"--help produce this output "<<std::endl;


	return;

}
//---------------------------------------------------------------------
//---------------------------------------------------------------------
bool VSAhfHaloGalaxyExtOp::allocatefArray()
//---------------------------------------------------------------------
{
unsigned long long int tempLL=getMaxNumberInt()*2;
if(((unsigned long long int)m_nOfEle*m_countField)>tempLL) 
	m_nOfEle=(int)tempLL/m_countField;


try
{
m_fArray=new  float*[m_countField];
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
	for(unsigned int i=0;i<m_countField;i++)
	{
try
{
	m_fArray[i] = new float[m_nOfEle];
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
			if(m_nOfEle<=MAX_NUMBER_TO_REDUCE_ROW)
				m_nOfEle=MIN_NUMBER_OF_ROW;
			else
				m_nOfEle=m_nOfEle-MAX_NUMBER_TO_REDUCE_ROW;
 		}
        }
}
return true;

}
//---------------------------------------------------------------------
//---------------------------------------------------------------------
bool VSAhfHaloGalaxyExtOp::createNewList()
//---------------------------------------------------------------------
{
int dummy;
m_fileInput.seekg(0,std::ios::beg);
m_fileInput>> dummy;
for(int i=0;i<m_nOfEle;i++)
{
	int j=0;
	for(j=0;j<m_countField;j++)
	    if(m_fArray[j][i]>=m_threshold)
		break;   
	int elements;
	m_fileInput>> elements;
	if(j<m_countField)
	{
		m_fileOutput[j]<<elements<<std::endl;
		for(int k=0;k<elements;k++)
		{
			m_fileInput>> dummy;
			m_fileOutput[j]<<dummy<<std::endl;
		}
	}else
	{		
		std::clog<<"dummy list elements "<<elements<<std::endl;
		for(int k=0;k<elements;k++)
		{
			m_fileInput>> dummy;
			std::clog<<"dummy list! "<<dummy<<std::endl;
		}
	}
}	
return true;
}
//---------------------------------------------------------------------
bool VSAhfHaloGalaxyExtOp::execute()
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
std::string fileNameOutput[3];
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
  	fileNameOutputSStream<<filenameInputTable.substr(0, len-4)<<"_AHFhalogalaxyext_"<<buffer;  //QUI verificare
} else
	fileNameOutputSStream<<getParameterAsString("out");


if(isParameterPresent("threshold"))
	m_threshold=getParameterAsFloat("threshold");

std::string colId[3];
std::stringstream ssListparameters;
ssListparameters.str(getParameterAsString("field"));
m_countField=0;
while(!ssListparameters.eof() && m_countField<3)
{
	ssListparameters>>colId[m_countField];
	m_countField++;
}
if(m_countField==0)
{
	std::cerr<<"Operation Aborted. Valid fields must be given"<<std::endl;
	return false;
}
for(int i=0;i<m_countField;i++)
	if((m_tables[0]->getColId(colId[i]))<0)
	{
	  std::cerr<<"Invalid field name "<<colId[i]<<". The column does not exist. Operation Aborted."<<std::endl;
	 return false;
        }
for(int i=0;i<m_countField;i++)
	fileNameOutput[i]=fileNameOutputSStream.str()+"_"+colId[i]+"_particles";
for(int i=0;i<m_countField;i++)
{
	m_fileOutput[i].open(fileNameOutput[i].c_str());
	if(!m_fileOutput[i].is_open())
	{
		std::cerr<<"Error opening output multilist file "<<fileNameOutput[i]<<" Operation aborted."<<std::endl;
		return false;
	} 
}



m_nOfRow=m_tables[0]->getNumberOfRows();
m_nOfCol=3;
if(m_nOfRow>getMaxNumberInt()) 
	m_nOfEle= getMaxNumberInt();
else 
	m_nOfEle= m_nOfRow;

if(!allocatefArray())
	return false;
unsigned int * colList;
colList=new unsigned int[m_countField];
for(int i=0;i<m_countField;i++)
   colList[i]=m_tables[0]->getColId(colId[i]);

unsigned long long int index=0;

// count elelmnts for each output file
int *fileNOfEle;
try
{
	fileNOfEle = new int[m_countField];
}
catch(std::bad_alloc &e)
{
	fileNOfEle=NULL;
}

if (fileNOfEle==NULL)
	return false;

for(int i=0;i<m_countField;i++) fileNOfEle[i]=0;
while(index< m_nOfRow)
{
	m_tables[0]->getColumn(colList,m_countField,index,index+m_nOfEle-1,m_fArray);
	for(int i=0;i<m_nOfEle;i++)
		for(int j=0;j<m_countField;j++)
	    		if(m_fArray[j][i]>=m_threshold)
			{
				fileNOfEle[j]++;
				break;
			}


	index+= m_nOfEle;
	if(index+m_nOfEle>= m_nOfRow) 
		m_nOfEle= m_nOfRow-index;
}
index=0;
for(int i=0;i<m_countField;i++) m_fileOutput[i]<<fileNOfEle[i]<<std::endl;
if(m_nOfRow>getMaxNumberInt()) 
	m_nOfEle= getMaxNumberInt();
else 
	m_nOfEle= m_nOfRow;

while(index< m_nOfRow)
{
	m_tables[0]->getColumn(colList,m_countField,index,index+m_nOfEle-1,m_fArray);
	
	if(!createNewList())
        {
		delete [] fileNOfEle;
                delete [] colList;	
		return false;
	}
	index+= m_nOfEle;
	if(index+m_nOfEle>= m_nOfRow) 
		m_nOfEle= m_nOfRow-index;

}
m_fileInput.close();
for(int i=0;i<m_countField;i++)
  m_fileOutput[i].close();
delete [] fileNOfEle;
delete [] colList;	
return true;
}
//---------------------------------------------------------------------


