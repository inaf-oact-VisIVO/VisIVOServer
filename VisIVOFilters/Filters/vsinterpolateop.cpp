/***************************************************************************
 *   Copyright (C) 2009 by Ugo Becciani   *
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
#include <set>
#include <sstream>
#include <sstream>
#include <iostream>
#include <fstream>
#include "vsinterpolateop.h"
#include "vstable.h"
#include "vsstatisticop.h"
#include <cmath>
#ifdef WIN32
	#include <time.h>
#endif
#include "visivoutils.h"

const unsigned int VSInterpolateOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSInterpolateOp::MIN_NUMBER_OF_ROW = 100;

//---------------------------------------------------------------------
VSInterpolateOp::VSInterpolateOp()
//---------------------------------------------------------------------
{ 
 m_f1Array=NULL;
 m_f2Array=NULL;
 m_f3Array=NULL;
 m_nOfRows=0;
 m_nOfCols=0;	

}

//---------------------------------------------------------------------
VSInterpolateOp::~VSInterpolateOp()
//---------------------------------------------------------------------
{
if(m_f1Array!=NULL)	
	for(unsigned int i=0;i<m_nOfCols;i++)
	{
		if(m_f1Array[i]!=NULL) delete [] m_f1Array[i];
	}
if(m_f2Array!=NULL)	
	for(unsigned int i=0;i<m_nOfCols;i++)
	{
		if(m_f2Array[i]!=NULL) delete [] m_f2Array[i];
	}
if(m_f3Array!=NULL)	
	for(unsigned int i=0;i<m_nOfCols;i++)
	{
		if(m_f3Array[i]!=NULL) delete [] m_f3Array[i];
	}
if(m_f1Array!=NULL) delete [] m_f1Array;
if(m_f2Array!=NULL) delete [] m_f2Array;
if(m_f3Array!=NULL) delete [] m_f3Array;
}
//---------------------------------------------------------------------
void VSInterpolateOp::printHelp()
//---------------------------------------------------------------------
{
	std::cout<<"Interpolate two tables with a linear scale. The two table must have the same structure."<<std::endl<<std::endl;;

	std::cout<<"Usage: VisIVOFilters --op interpolate [--field columns_name]  [--numbin numberbin] [--periodic] [--interval from to] [--index column_name] [--out filename_out] [--history] [--historyfile filename.xml] [--help] --infiles file_start.bin file_end.bin"<<std::endl;

	std::cout<<"Example: VisIVOFilters --op interpolate --field X Y Z --numbin 20 --periodic --out mysequence  --infiles start.bin end.bin "<<std::endl<<std::endl;

	std::cout<<"Note: "<<std::endl;
	std::cout<<"--field a valid list of  columns name that must be exist on both the input tables. Default: all columns are considered"<<std::endl;
//	std::cout<<"--log generate frames with a log scale. Default is linear scale"<<std::endl;
	std::cout<<"--numbin is the number of bins between the starting and ending frames (input files)  or the interval given in the --interval option. Default value is 10. The number of created tables is equal to numberbin-1."<<std::endl;
	std::cout<<"--periodic applies a periodical boundary condition options."<<std::endl;
	std::cout<<"--interval. VisIVO assumes a distance of 1.0 between the starting frame and the ending tables. This option allow to produces the intermediate frames (table) in a subinterval between the two frames. 0.5 is the medium point of the interval. If the from value is lower than 0.0 it is considered 0.0. If the to value is lower than 1.0 it is considered 1.0. If from value is equal to to value the operation is not performed.  Default value from=0.0 to =1.0"<<std::endl;
	std::cout<<"--index is the column that contains points id for the interpolation. If this option is not given, the points order is  implicit in the table."<<std::endl;
	std::cout<<"--start is the suffix starting number used in --out option"<<std::endl;
	std::cout<<"--out is the root name of the new tables. Default name is given. New name is given by filename_out#.bin where # is the number of created table."<<std::endl;
    std::cout<<"--history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
    std::cout<<"--historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;

	std::cout<<"--infiles contains the names of starting and ending tables of the interpolation process."<<std::endl;
	std::cout<<"--help produce this output "<<std::endl;

	return;

}  

//---------------------------------------------------------------------
bool VSInterpolateOp::allocateArray()
//---------------------------------------------------------------------
{
	if(m_nOfRows>2*(getMaxNumberInt()/(m_nOfCols*3)))
		m_nOfRows=2*(getMaxNumberInt()/(m_nOfCols*3));
try
{
	m_f1Array=new  float*[m_nOfCols];
}
catch(std::bad_alloc & e)
{
	m_f1Array=NULL;
}

	if(m_f1Array == NULL)
		return false;

try
{
	m_f2Array=new  float*[m_nOfCols];
}
catch(std::bad_alloc & e)
{
	m_f2Array=NULL;
}

	if(m_f2Array == NULL)
	{
		delete [] m_f1Array;
		return false;
	}
try
{
	m_f3Array=new  float*[m_nOfCols];
}
catch(std::bad_alloc & e)
{
	m_f3Array=NULL;
}

	if(m_f3Array == NULL)
	{
		delete [] m_f1Array;
		delete [] m_f2Array;
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
			m_f1Array[i]=new  float[m_nOfRows];
}
catch(std::bad_alloc & e)
{
	m_f1Array[i]=NULL;
}
try
{
			m_f2Array[i]=new  float[m_nOfRows];
}
catch(std::bad_alloc & e)
{
	m_f2Array[i]=NULL;
}
try
{
			m_f3Array[i]=new  float[m_nOfRows];
}
catch(std::bad_alloc & e)
{
	m_f3Array[i]=NULL;
}
			if(m_f1Array[i]==NULL || m_f2Array[i]==NULL || m_f3Array[i]==NULL ) 
			{	
				goodAllocation=false;
				if(m_f1Array[i]!=NULL)
				{ 
					delete [] m_f1Array[i];
					m_f1Array[i]=NULL;
				}
				if(m_f2Array[i]!=NULL)
				{ 
					delete [] m_f2Array[i];
					m_f2Array[i]=NULL;
				}
				if(m_f3Array[i]!=NULL)
				{ 
					delete [] m_f3Array[i];
					m_f3Array[i]=NULL;
				}
				for(unsigned int j=0;j<i;j++)
				{ 
					delete [] m_f1Array[j];
					delete [] m_f2Array[j];
					delete [] m_f3Array[j];
					m_f1Array[j]=NULL;
					m_f2Array[j]=NULL;
					m_f3Array[j]=NULL;
				}
				if(m_nOfRows==MIN_NUMBER_OF_ROW)
				{ 
					delete [] m_f1Array;
					delete [] m_f2Array;
					delete [] m_f3Array;
					m_f1Array=NULL;
					m_f2Array=NULL;
					m_f3Array=NULL;
					return false;
				}
				m_nOfRows=m_nOfRows-MAX_NUMBER_TO_REDUCE_ROW;
				if(m_nOfRows<=MAX_NUMBER_TO_REDUCE_ROW) m_nOfRows=MIN_NUMBER_OF_ROW;
				break;
			}
		}
	}

	return true;

}
//---------------------------------------------------------------------
bool VSInterpolateOp::execute()
//---------------------------------------------------------------------
{

bool interval=false;
float intervalFrom=0.0;
float intervalTo=1.0;
int startNumber=0;
if(isParameterPresent("interval"))
{
	interval=true;
	std::stringstream tmp;
	tmp.str(getParameterAsString("interval"));
	tmp>>intervalFrom;
	tmp>>intervalTo;
	if(intervalFrom<0.0) intervalFrom=0.0;
	if(intervalTo>1.0) intervalTo=1.0;
	if(intervalFrom>=intervalTo)
	{
		std::cerr<<"Interpolate: Invalid range in interval option"<<std::endl;
		return false;
	}
}	

if(getParameterAsString("infiles").empty() || getParameterAsString("infiles")=="unknown" )
{
	std::cerr<<"vsinterpolateop: No input files are given"<<std::endl;
	return false;
}
std::stringstream ssInfileparameters;
ssInfileparameters.str(getParameterAsString("infiles"));
std::string file1="None";
std::string file2="None";
ssInfileparameters>>file1;
ssInfileparameters>>file2;
if(file1=="None" || file2=="None")
{
	std::cerr<<"Error: the interpolate operation requires two tables"<<std::endl; 
	return false;
}
if(file1.find(".bin") == std::string::npos)
    file1.append(".bin");
if(file2.find(".bin") == std::string::npos)
    file2.append(".bin");

int numbin=getParameterAsInt("numbin");
if(numbin<=0) numbin=10;


VSTable table1(file1);
VSTable table2(file2);
if(table1.getNumberOfRows()==0 || table2.getNumberOfRows()==0)
{
	std::cerr<<"Error: Tables do not exist or contains zero rows."<<std::endl;
	return false;

}
if(table1.getNumberOfRows()!=table2.getNumberOfRows())
{
	std::cerr<<"Error: Tables does not contain the same number of rows."<<std::endl;
	return false;

}
std::set<unsigned int> colNumberSet;
std::set<unsigned int>::iterator iter;
std::stringstream ssListparameters;
if(isParameterPresent("field"))
{
  
  if(!isParameterPresent("field"))
  {
	int tmpi=table1.getNumberOfColumns();
	for(unsigned int i=0;i<table1.getNumberOfColumns();i++)
	{
		if(table1.getColName(i)!=table2.getColName(i))
		{
			std::cerr<<"Error: Tables does not contain the same columns or the tables are not in the same order."<<std::endl;
 			colNumberSet.clear();
			return false;
		}		
		colNumberSet.insert(i);
		ssListparameters<<table1.getColName(i)<<" ";
	}

  } else
  {
    ssListparameters.str(getParameterAsString("field"));
    colNumberSet.clear();
    while (!ssListparameters.eof())
    {
	std::string paramField;
	ssListparameters>>paramField;
	if(table1.getColId(paramField)>=0)
		colNumberSet.insert(table1.getColId(paramField));
	if(table1.getColId(paramField)!=table2.getColId(paramField))
	{
		std::cerr<<"Error: Tables does not contain the same columns or the tables are not in the same order."<<std::endl;
 		colNumberSet.clear();
		return false;
	}
	
    }
}
  
  
} else  //obsolete
{
  if(!isParameterPresent("list"))
  {
	int tmpi=table1.getNumberOfColumns();
	for(unsigned int i=0;i<table1.getNumberOfColumns();i++)
	{
		if(table1.getColName(i)!=table2.getColName(i))
		{
			std::cerr<<"Error: Tables does not contain the same columns or the tables are not in the same order."<<std::endl;
 			colNumberSet.clear();
			return false;
		}		
		colNumberSet.insert(i);
		ssListparameters<<table1.getColName(i)<<" ";
	}

  } else
  {
    ssListparameters.str(getParameterAsString("list"));
    colNumberSet.clear();
    while (!ssListparameters.eof())
    {
	std::string paramField;
	ssListparameters>>paramField;
	if(table1.getColId(paramField)>=0)
		colNumberSet.insert(table1.getColId(paramField));
	if(table1.getColId(paramField)!=table2.getColId(paramField))
	{
		std::cerr<<"Error: Tables does not contain the same columns or the tables are not in the same order."<<std::endl;
 		colNumberSet.clear();
		return false;
	}
	
    }
}

}
bool periodic=false;
if(isParameterPresent("periodic")) periodic=true;

m_nOfCols=colNumberSet.size();
m_nOfRows=table1.getNumberOfRows();
float *maxMax,*minMin,*avgAvg;
try
{
	maxMax=new  float[m_nOfCols];
}
catch(std::bad_alloc &e)
{
	maxMax=NULL;
	std::cerr<<"vsinterpolateop: Not enough memory left"<<std::endl;
	return false;
}
try
{
	minMin=new  float[m_nOfCols];
}
catch(std::bad_alloc &e)
{
	minMin=NULL;
        delete [] maxMax;
	maxMax=NULL;
	std::cerr<<"vsinterpolateop: Not enough memory left"<<std::endl;
	return false;
}
try
{
	avgAvg=new  float[m_nOfCols];
}
catch(std::bad_alloc &e)
{
	avgAvg=NULL;
        delete [] maxMax;
        delete [] minMin;
	minMin=NULL;
	maxMax=NULL;
	std::cerr<<"vsinterpolateop: Not enough memory left"<<std::endl;
	return false;
}
bool index=false;
std::string indexColName;
if(!(getParameterAsString("index").empty() || getParameterAsString("index")=="unknown" ))
{
	std::cerr<<"Interpolate: index option is not available"<<std::endl;
	index=true;
	indexColName=getParameterAsString("index");
	int indexId= table1.getColId(indexColName);
	if(indexId==-1)
		std::cerr<<"Invalid index column name"<<std::endl;
	return false;
} 

/*** Execute the Statistic  OP **/
{
std::string listParameter;
if(isParameterPresent("list")) 
	listParameter=getParameterAsString("list");
else
	listParameter=ssListparameters.str();
VSStatisticOp op1;
VSStatisticOp op2;
op1.addParameter("list",listParameter);
op2.addParameter("list",listParameter);
op1.addParameter("silent","");
op2.addParameter("silent","");
op1.addInput(&table1);
op2.addInput(&table2);
op1.execute();
op2.execute();
float max1,min1,avg1,max2,min2,avg2;
for (unsigned int i=0;i<m_nOfCols;i++)
{
	float dummy;
	op1.getRange(i,max1,min1,avg1,dummy);
	op2.getRange(i,
max2,min2,avg2,dummy);
	if(min1<=min2) 
		minMin[i]=min1;
	else
		minMin[i]=min2;
	if(max1>=max2) 
		maxMax[i]=max1;
	else
		maxMax[i]=max2;
	avgAvg[i]=(maxMax[i]-minMin[i])/2.0;
}	
} //block to destroy statistic and op1 and op2
/*** END Statistic OP **/
unsigned int *colList=NULL;
int count=colNumberSet.size();
try
{
	colList=new unsigned int[count];
}
catch(std::bad_alloc &e)
{
	colList=NULL;
}

if(colList == NULL)
{
	std::cerr<<"Failed colList allocation. Interpolate terminated"<<std::endl;
	delete [] maxMax;
	delete [] minMin;
	delete [] avgAvg;
	return false;
}
count=0;
for (iter=colNumberSet.begin(); iter!=colNumberSet.end(); iter++)
{
	colList[count]=*iter;
	count++;
}
VSTable *tableInterpolate;
try
{
	tableInterpolate=new VSTable[numbin-1];
}
catch(std::bad_alloc &e)
{
	tableInterpolate=NULL;
}

if(tableInterpolate == NULL)
{
	std::cerr<<"Failed tableInterpolate allocation. Interpolate terminated"<<std::endl;
	delete [] maxMax;
	delete [] minMin;
	delete [] avgAvg;
	delete [] colList;
	return false;
}



unsigned int maxEle=m_nOfRows;
bool allocationArray=allocateArray();
if(!allocationArray)
{
	std::cerr<<"Failed allocateArray. Interpolate terminated"<<std::endl;
	delete [] maxMax;
	delete [] minMin;
	delete [] avgAvg;
	delete [] colList;
	delete [] tableInterpolate;
	return false;
}
maxEle=m_nOfRows;
unsigned long long int totRows=table1.getNumberOfRows();



//Set Output Table Filenames
std::stringstream fileNameOutputSStream;
fileNameOutputSStream<<getParameterAsString("out"); //QUI brutto
std::string fileNameOutputRoot;
if(fileNameOutputSStream.str()==""||fileNameOutputSStream.str()=="unknown")
{
 	std::string tmp=table1.getLocator();
	int len=tmp.length();
	fileNameOutputRoot=tmp.substr(0, len-4)+"_interpolate_";
}else
{	
	std::string tmp=fileNameOutputSStream.str();
	if(tmp.find(".bin") == std::string::npos)
   		fileNameOutputRoot=tmp;
	else
	{
		int len=tmp.length();
		fileNameOutputRoot=tmp.substr(0, len-4);
	}
}

std::string VSCycleFileList;
time_t rawtime;
struct tm * timeinfo;
char buffer [80];
time ( &rawtime );
timeinfo = localtime ( &rawtime );
strftime (buffer,80,"%Y%m%d%H%M%S",timeinfo);
std::stringstream fileCycleSStream;
fileCycleSStream<<getDir(fileNameOutputRoot)<<"VSInterpolatefilelist_"<<buffer<<".txt";  //QUI 
VSCycleFileList=fileCycleSStream.str();
std::ofstream ouFileList;
ouFileList.open(VSCycleFileList.c_str());


// creat numbin -1 headers

if(isParameterPresent("start"))
  startNumber=getParameterAsInt("start")-1;

for(int i=0; i<numbin-1; i++)
{	int isuffix=startNumber+i+1;
	std::string zeropad;
	if(isuffix<10) zeropad="000000";
	if(isuffix>=10 && isuffix<100 ) zeropad="00000";
	if(isuffix>=100 && isuffix<1000 ) zeropad="0000";
	if(isuffix>=1000 && isuffix<10000 ) zeropad="000";
	if(isuffix>=10000 && isuffix<100000 ) zeropad="00";
	if(isuffix>=100000 && isuffix<1000000 ) zeropad="0";
	if(isuffix>=1000000) zeropad="";
	std::stringstream tmp;
	tmp << fileNameOutputRoot<<"_"<<zeropad<<startNumber+i+1<<".bin";
	std::string fileNameOutput=tmp.str();
 	tableInterpolate[i].setLocator(fileNameOutput);
	ouFileList<<fileNameOutput<<std::endl;
	m_realOutFilename.push_back(fileNameOutput);
	
#ifdef VSBIGENDIAN
	std::string endianism="big";
	
#else	
	std::string endianism="little";
#endif

	tableInterpolate[i].setEndiannes(endianism);

 	tableInterpolate[i].setType("float");
 	tableInterpolate[i].setNumberOfRows(totRows);
	tableInterpolate[i].setIsVolume(table1.getIsVolume());
	if(table1.getIsVolume())
	{
	    tableInterpolate[i].setCellNumber(table1.getCellNumber()[0],
                           table1.getCellNumber()[1],
                           table1.getCellNumber()[2]);
            tableInterpolate[i].setCellSize(table1.getCellSize()[0],
                           table1.getCellSize()[1],
                           table1.getCellSize()[2]);
	}

	for (iter=colNumberSet.begin(); iter!=colNumberSet.end(); iter++)
		tableInterpolate[i].addCol(table1.getColName(*iter));	
 
	tableInterpolate[i].writeHeader();  //overwrite if table exist!
}
// write values in tables
ouFileList.close();
unsigned long long int totEle=totRows;
unsigned long long int fromRow, toRow, startCounter=0;

while(totEle!=0)
{		
	fromRow=startCounter;
	toRow=fromRow+maxEle-1;
	if(toRow>totRows-1)toRow=totRows-1;
	table1.getColumn(colList, m_nOfCols, fromRow, toRow, m_f1Array);
	table2.getColumn(colList, m_nOfCols, fromRow, toRow, m_f2Array);
// Interval section
	if(interval)
	{	
		for(int k=0;k<m_nOfCols;k++)
		{
			for(int j=0;j<=(toRow-fromRow);j++)
			{
				float step;
				float temp1=fabs(m_f2Array[k][j]-m_f1Array[k][j]);
				if(periodic && temp1>avgAvg[k])
				{
					if(m_f1Array[k][j]<=m_f2Array[k][j])
						step=-1*((m_f1Array[k][j]-minMin[k])+(maxMax[k]-m_f2Array[k][j]));
					else
						step=(m_f2Array[k][j]-minMin[k])+(maxMax[k]-m_f1Array[k][j]);
				}else
				{
					if(m_f1Array[k][j]<=m_f2Array[k][j])
						step=temp1;
					else
						step=-1*temp1;	
				}
					m_f1Array[k][j]=m_f1Array[k][j]+step*intervalFrom;	
					m_f2Array[k][j]=m_f2Array[k][j]-step*(1.0-intervalTo);
			}	
		}	
	}
//
	for(int i=0; i<numbin-1; i++)
	{
		for(int k=0;k<m_nOfCols;k++)
		{
			for(int j=0;j<=(toRow-fromRow);j++)
			{
				float step;
				float temp1=fabs(m_f2Array[k][j]-m_f1Array[k][j]);
				if(periodic && temp1>avgAvg[k])
				{
					if(m_f1Array[k][j]<=m_f2Array[k][j])
						step=-1*((m_f1Array[k][j]-minMin[k])+(maxMax[k]-m_f2Array[k][j]));
					else
						step=(m_f2Array[k][j]-minMin[k])+(maxMax[k]-m_f1Array[k][j]);
				}else
				{
					if(m_f1Array[k][j]<=m_f2Array[k][j])
						step=temp1;
					else
						step=-1*temp1;	
				}
		
				m_f3Array[k][j]=m_f1Array[k][j]+(i+1)*(step/(float)numbin);
				if(periodic && m_f3Array[k][j]>maxMax[k])
 					m_f3Array[k][j]=minMin[k]+(m_f3Array[k][j]-maxMax[k]);
				if(periodic && m_f3Array[k][j]<minMin[k])
 					m_f3Array[k][j]=maxMax[k]-(minMin[k]- m_f3Array[k][j]);
			}
		}
		tableInterpolate[i].putColumn(colList, m_nOfCols, fromRow, toRow, m_f3Array);
	}
	startCounter=toRow+1;
	totEle=totEle-(toRow-fromRow+1);
}

delete [] colList;
delete [] tableInterpolate;
delete [] maxMax;
delete [] minMin;
delete [] avgAvg;

return true;
}