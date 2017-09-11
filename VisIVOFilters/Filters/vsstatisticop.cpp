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
#include <set>
#include <sstream>
#include <fstream>
#include "vstable.h"
#include "vsstatisticop.h"

const unsigned int VSStatisticOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSStatisticOp::MIN_NUMBER_OF_ROW = 100;

VSStatisticOp::VSStatisticOp()
{
 m_fArray=NULL;
 m_nOfRow=0;
 m_nOfCol=0;	
 m_maxValue=NULL;
 m_minValue=NULL;
 m_sumValue=NULL;
 m_averageValue=NULL;
 m_histogram=NULL;
 m_numberOfHisto=0;

m_numberOfListParams=-100; // this is a "bad value"
}


VSStatisticOp::~VSStatisticOp()
{
	if(m_histogram!=NULL)
		for(unsigned int i=0;i<m_numberOfHisto;i++)
		{
			if(m_histogram[i]!=NULL) delete []m_histogram[i] ;
		}
	if(m_fArray!=NULL)
		for(unsigned int i=0;i<m_nOfCol;i++)
		{
			if(m_fArray[i]!=NULL) delete [] m_fArray[i];
		}
	if(m_fArray!=NULL) delete [] m_fArray;
	if(m_maxValue!=NULL) delete [] m_maxValue;
	if(m_minValue!=NULL) delete [] m_minValue;
	if(m_averageValue!=NULL) delete [] m_averageValue;
	if(m_sumValue!=NULL) delete [] m_sumValue;

}
//---------------------------------------------------------------------
bool VSStatisticOp::allocateArray()
//---------------------------------------------------------------------
{
unsigned long long int tempLL=getMaxNumberInt()*2;
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
		return false;

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
					delete [] m_fArray[j];
				if(m_nOfRow==MIN_NUMBER_OF_ROW)
				{ 
					delete [] m_fArray;
					m_fArray=NULL;
					return false;
				}
				m_nOfRow=m_nOfRow-MAX_NUMBER_TO_REDUCE_ROW;
				if(m_nOfRow<=MAX_NUMBER_TO_REDUCE_ROW) m_nOfRow=MIN_NUMBER_OF_ROW;
				break;
			}

		}
		
	}
	return true;
}
//--------------------------------------------------------------------
void VSStatisticOp::printHelp()
//---------------------------------------------------------------------
{
std::cout<<"Produce average, min and max value of field and creates an histogram of fields in the input table"<<std::endl<<std::endl;
std::cout<<"Usage: VisIVOFilters --op statistic [--list columns_name] [--histogram [bin]] [--range min max] [--out result.txt] [--history] [--historyfile filename.xml] [--help] [--file] inputFile.bin"<<std::endl;

std::cout<<"Example: VisIVOFilters --op statistic --list X Y --histogram 1000 --range 10.0 100.0  --out result.txt --file inputFile.bin"<<std::endl;

std::cout<<"Note:"<<std::endl;
std::cout<<"--list a valid list of  columns name. Default value all columns."<<std::endl;
std::cout<<"--histogram Produce histogram ascii file  with given number of bin. If bin number is not specified, the default value is fixed to 10% of the total rows of the input table."<<std::endl;
std::cout<<"--range produce the results only inside the specified interval."<<std::endl;
std::cout<<"--out Output ascii filename with histogram."<<std::endl;
std::cout<<"--historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;
std::cout<<"--filelist list of tables and valid column name to be merged. Wildcard * meabs all columns."<<std::endl;
std::cout<<"--file  Input table filename."<<std::endl;
std::cout<<"--help produce this output "<<std::endl;

return;
 }
//--------------------------------------------------------------------
bool VSStatisticOp::getRange(unsigned int parElement, float &max, float &min, float  &avg, float &sum)
//--------------------------------------------------------------------
{ 
	if(m_numberOfListParams<0)
		return false;
	if(parElement>(m_numberOfListParams-1))
		return false;
	max=m_maxValue[parElement]; 
	min= m_minValue[parElement]; 
	avg=m_averageValue[parElement]; 
	sum=m_sumValue[parElement]; 
	return true;
}
//---------------------------------------------------------------------
bool VSStatisticOp::execute()
//---------------------------------------------------------------------
// totRows is the total rows of the out tables
// maxRows is the bigger number of rows among input tables
// 
// maxEle maximum number allowed to upload/download data from a generic table: It is the lower value between maxRows, maxInt and fArray allocation fArray[0][maxEle]
// 
// nOfEle total number of element to upload for the specific input table.
// 
// totEle is a counter: the number of elements that I have still to  upload for the specific table. Starts from nOfEle, its value decrese of maxEle each cycle
{
	bool Api=false;
	if(isParameterPresent("Api")) Api=true;
	if(Api && !isParameterPresent("histogram"))
	  addParameter("histogram","unknown");
  
	bool silent=false;
	if(isParameterPresent("silent"))
		silent=true;
	float offset=0;
	float binInterval;
	int numberOfBin=0;
	std::string fileNameOutput;
	bool range=true;
	bool histoRequest=true;
	float minRange, maxRange;
	
	if(!isParameterPresent("histogram"))
			histoRequest=false;
	if(histoRequest)
	{
		numberOfBin=getParameterAsInt("histogram");
		if(numberOfBin<=0 )
		{	
			unsigned long long int numberOfDefaultBin=m_tables[0]->getNumberOfRows()/100;
			if(numberOfDefaultBin>10000) numberOfBin=10000;
			if(numberOfDefaultBin<10) numberOfBin=10;
			if(numberOfDefaultBin>=10 && numberOfDefaultBin<=10000) numberOfBin=numberOfDefaultBin;
		}
	}
	if(getParameterAsString("range").empty() || getParameterAsString("range")=="unknown" )
		range=false;
	else
	{
		range=true;
		std::stringstream ssListparameters;
		ssListparameters.str(getParameterAsString("range"));
		ssListparameters>>minRange;
		ssListparameters>>maxRange;
	}
	
	std::set<unsigned int> colNumberSet;
	std::stringstream ssListparameters;
	std::set<unsigned int>::iterator iter;
	colNumberSet.clear();
	if(isParameterPresent("field"))
	{
	  if(!isParameterPresent("field"))
	  {
		for(unsigned int i=0;i<m_tables[0]->getNumberOfColumns();i++)
			colNumberSet.insert(i);
	  }else{
		ssListparameters.str(getParameterAsString("field"));
		while (!ssListparameters.eof())
		{
			std::string paramField;
			ssListparameters>>paramField;
			if(m_tables[0] -> getColId(paramField)>=0)
				colNumberSet.insert(m_tables[0] -> getColId(paramField));
		}
	  }
	}else //obsolete
	{
	  if(!isParameterPresent("list"))
	  {
		for(unsigned int i=0;i<m_tables[0]->getNumberOfColumns();i++)
			colNumberSet.insert(i);
	  }else{
		ssListparameters.str(getParameterAsString("list"));
		while (!ssListparameters.eof())
		{
			std::string paramField;
			ssListparameters>>paramField;
			if(m_tables[0] -> getColId(paramField)>=0)
				colNumberSet.insert(m_tables[0] -> getColId(paramField));
		}
	  }
	}
	m_numberOfHisto=colNumberSet.size();
	m_numberOfListParams=colNumberSet.size();
	if(m_numberOfHisto==0)
	{
		std::cerr<<"vsstatisticop: Invalid columns name is given"<<std::endl;
 		return false;
	}

try
{
	m_maxValue= new float[m_numberOfHisto];
}
catch(std::bad_alloc &e)
{
	m_maxValue=NULL;
}

	if(m_maxValue==NULL)
	{
		std::cerr<<"vsstatisticop: Number of requested bin too high"<<std::endl;
 		return false;
	}
try
{
	m_minValue= new float[m_numberOfHisto];
}
catch(std::bad_alloc &e)
{
	m_minValue=NULL;
}

	if(m_minValue==NULL)
	{
		std::cerr<<"vsstatisticop: Number of requested bin too high"<<std::endl;
		delete [] m_maxValue;
 		return false;
	}
try
{
	m_averageValue= new float[m_numberOfHisto];
}
catch(std::bad_alloc &e)
{
	m_averageValue=NULL;
}

	if(m_averageValue==NULL)
	{
		std::cerr<<"vsstatisticop: Number of requested bin too high"<<std::endl;
 		delete [] m_maxValue;
		delete [] m_minValue;
		return false;
	}
try
{
	m_sumValue= new float[m_numberOfHisto];
}
catch(std::bad_alloc &e)
{
	m_sumValue=NULL;
}

	if(m_sumValue==NULL)
	{
		std::cerr<<"vsstatisticop: Number of requested bin too high"<<std::endl;
 		delete [] m_maxValue;
		delete [] m_minValue;
		delete [] m_averageValue;
		return false;
	}

	if(histoRequest)
	{
		std::stringstream fileNameOutputSStream;
		fileNameOutputSStream<<getParameterAsString("out");
		if(fileNameOutputSStream.str()==""||fileNameOutputSStream.str()=="unknown")
		{
  			std::string filenameInputTable=m_tables[0]->getLocator();
  			int len=filenameInputTable.length();
  			fileNameOutputSStream<<filenameInputTable.substr(0, len-4)<<"_statistic_gnuplotfile"<<getParameterAsString("field")<<".txt";
		}
		fileNameOutput=fileNameOutputSStream.str();
			
		m_realOutFilename.push_back(fileNameOutput);

		remove(fileNameOutput.c_str());
try
{
		m_histogram= new unsigned long long int*[m_numberOfHisto];
}
catch(std::bad_alloc &e)
{
	m_histogram=NULL;
}

		if(m_histogram==NULL)
		{
			std::cerr<<"vsstatisticop: Number of requested bin too high"<<std::endl;
  			delete [] m_maxValue;
			delete [] m_minValue;
			delete [] m_averageValue;
			delete [] m_sumValue;
			return false;
		}
		for(unsigned int i=0;i<m_numberOfHisto;i++)
		{
try
{
		  m_histogram[i]= new unsigned long long int[numberOfBin+1];
}
catch(std::bad_alloc &e)
{
	m_histogram[i]=NULL;
}

		  if(m_histogram[i]==NULL)
		  {
			std::cerr<<"vsstatisticop: Number of requested bin too high"<<std::endl;
		        for(unsigned int j=0;j<i;j++) 
			       delete [] m_histogram[j];
			delete [] m_histogram;
  			delete [] m_maxValue;
			delete [] m_minValue;
			delete [] m_averageValue;
			delete [] m_sumValue;
 			return false;
		  }
		  for(unsigned int j=0;j<=numberOfBin;j++) m_histogram[i][j]=0;
		}
	}

	unsigned long long int totRows=m_tables[0]->getNumberOfRows();
	unsigned long long int nOfEle=totRows;
	int maxInt=getMaxNumberInt();
	unsigned int maxEle;

	if(totRows>maxInt)
		maxEle=maxInt; 
	else
		maxEle=totRows;


	m_nOfCol = 1;
	m_nOfRow=maxEle;

 	bool allocationArray=allocateArray();
	maxEle=m_nOfRow;
	if(m_fArray==NULL ||  !allocationArray )	
	{
		std::cerr<<"Failed Array allocation. Showtable Operation terminated"<<std::endl;
		if(histoRequest)
		{
			for(unsigned int j=0;j<m_numberOfHisto;j++) 
			       delete [] m_histogram[j];
			delete [] m_histogram;
		}
  		delete [] m_maxValue;
		delete [] m_minValue;
		delete [] m_averageValue;
		delete [] m_sumValue;
		return false;
	}

	int counterCols=-1;
	std::set<unsigned int>::iterator end = colNumberSet.end();;
	for(iter = colNumberSet.begin(); iter != end; iter++) 	
	{
 	  unsigned long long int totEle=nOfEle;
	  unsigned long long int fromRow, toRow, startCounter=0;
	  unsigned int colList[1];	
	  colList[0]=*iter;
	  unsigned int nOfValidElement=0;
	  counterCols++;
	  m_sumValue[counterCols]=0.;
	  m_averageValue[counterCols]=0.;
	  while(totEle!=0)
	  { 
		fromRow=startCounter;
		toRow=fromRow+maxEle-1;
		if(toRow>totRows-1)toRow=totRows-1;
		
	  	m_tables[0]->getColumn(colList, m_nOfCol, fromRow, toRow, m_fArray);

		if(startCounter==0)
		{
			  m_maxValue[counterCols]=m_fArray[0][0];
			  m_minValue[counterCols]=m_fArray[0][0];
		}
		if(range)
		{
			for(unsigned int j=0;j<(toRow-fromRow+1);j++)
			{
			  if(range && (m_fArray[0][j]>maxRange || m_fArray[0][j]<minRange)) continue;
			  m_maxValue[counterCols]=m_fArray[0][j];
			  m_minValue[counterCols]=m_fArray[0][j];
			  break;
			}
		}

		for(unsigned int j=0;j<(toRow-fromRow+1);j++)
		{
			if(range && (m_fArray[0][j]>maxRange || m_fArray[0][j]<minRange)) continue;
			if(m_maxValue[counterCols]<m_fArray[0][j]) m_maxValue[counterCols]=m_fArray[0][j];
			if(m_minValue[counterCols]>m_fArray[0][j]) m_minValue[counterCols]=m_fArray[0][j];
			m_averageValue[counterCols]=m_averageValue[counterCols]+m_fArray[0][j];
			m_sumValue[counterCols]+=m_fArray[0][j];
			nOfValidElement++;
		}
		startCounter=toRow+1;
		totEle=totEle-(toRow-fromRow+1);
	  }
	  if(nOfValidElement>0)
		  m_averageValue[counterCols]=m_averageValue[counterCols]/nOfValidElement;
	  else
	    {
		std::cerr<<"vsstatisticop: Invalid range: no data in the specified range"<<std::endl;
		if(histoRequest)
		{
			for(unsigned int j=0;j<m_numberOfHisto;j++) 
			       delete [] m_histogram[j];
			delete [] m_histogram;
		}
  		delete [] m_maxValue;
		delete [] m_minValue;
		delete [] m_averageValue;
		delete [] m_sumValue;
		return false;
	    }
	}

	float maxMaxValue=m_maxValue[0];
	float minMinValue=m_minValue[0];
	for(unsigned int i=1;i<m_numberOfHisto;i++)
	{
		if(maxMaxValue<m_maxValue[i])maxMaxValue=m_maxValue[i];
		if(minMinValue>m_minValue[i])minMinValue=m_minValue[i];
	}

	if(histoRequest)
	{
		
		if(!range) 
			binInterval=((float)(maxMaxValue-minMinValue))/(float) numberOfBin;
		else
			binInterval=(maxRange-minRange)/numberOfBin;
		counterCols=-1;
		offset=minMinValue;
		if(range) offset=minRange;
		for(iter = colNumberSet.begin(); iter != end; iter++) 	
		{
 	  	  unsigned long long int totEle=nOfEle;
	  	  unsigned long long int fromRow, toRow, startCounter=0;
	  	  unsigned int colList[1];	
	  	  colList[0]=*iter;
		  counterCols++;
		  totEle=nOfEle;
		  startCounter=0;
		  int histoIndex;
		  while(totEle!=0)
		  {
			fromRow=startCounter;
			toRow=fromRow+maxEle-1;
			if(toRow>totRows-1)toRow=totRows-1;
	  		m_tables[0]->getColumn(colList, m_nOfCol, fromRow, toRow, m_fArray);
			for(unsigned int j=0;j<(toRow-fromRow+1);j++)
			{
				if(range && (m_fArray[0][j]>maxRange || m_fArray[0][j]<minRange)) continue;
				if(binInterval==0) 
				    histoIndex=0;
				else
				  histoIndex=(int) ((m_fArray[0][j]-offset)/binInterval);
				if(histoIndex>numberOfBin)
					histoIndex=numberOfBin;
				m_histogram[counterCols][histoIndex]++;
			}
			startCounter=toRow+1;
			totEle=totEle-(toRow-fromRow+1);
		  }
	          m_histogram[counterCols][numberOfBin-1]=m_histogram[counterCols][numberOfBin-1]+m_histogram[counterCols][numberOfBin];
		}
	}
	counterCols=-1;
	if(!silent)
	if(!Api)
	{
	  for(iter = colNumberSet.begin(); iter != end; iter++) 	
	  {	
	   counterCols++;
	   std::string colName=m_tables[0]->getColName(*iter);
	   std::cout<<" Column "<<colName<<" :";
	   std::cout<<"  max="<<m_maxValue[counterCols];
	   std::cout<<"  min="<<m_minValue[counterCols];
	   std::cout<<"  average="<<m_averageValue[counterCols];
	   std::cout<<"  sum="<<m_sumValue[counterCols]<<std::endl;
	  }
	}
  	if(histoRequest && !silent)
      	{
		std::ofstream fileOutput(fileNameOutput.c_str(), std::ios::out);
		for(unsigned int i=0;i<numberOfBin;i++)
		{
			fileOutput<<i*binInterval+offset;
			for(unsigned int j=0;j<m_numberOfHisto;j++)
			   fileOutput<<	"  "<<m_histogram[j][i];
			fileOutput<<std::endl;
		}
		fileOutput.close();
	}
return true;
}
