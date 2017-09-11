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
#include <set>
#include <fstream>
#include "vstable.h"
#include "vswrvotableop.h"
#include <time.h>

const unsigned int VSWriteVotableOp::MAXROWS = 1000000;
const unsigned int VSWriteVotableOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSWriteVotableOp::MIN_NUMBER_OF_ROW = 100;
const unsigned int VSWriteVotableOp::WIDTH = 10;
const unsigned int VSWriteVotableOp::PRECISION = 8;

VSWriteVotableOp::VSWriteVotableOp()
{
 m_fArray=NULL;
 m_nOfRows=0;
 m_nOfCols=0;	

}


VSWriteVotableOp::~VSWriteVotableOp()
{
 	if(m_fArray!=NULL)
		for(unsigned int i=0;i<m_nOfCols;i++)
		{
			if(m_fArray[i]!=NULL) delete [] m_fArray[i];
		}
	if(m_fArray!=NULL) delete [] m_fArray;

}
//---------------------------------------------------------------------
bool VSWriteVotableOp::allocateArray()
//---------------------------------------------------------------------
{
unsigned long long int tempLL=getMaxNumberInt()*2;
if(((unsigned long long int)m_nOfRows*m_nOfCols)>tempLL) 
	m_nOfRows=(int)tempLL/m_nOfCols;

try
{
	m_fArray=new  float*[m_nOfCols];
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
		for(unsigned int i=0;i<m_nOfCols;i++)
		{
try
{
			m_fArray[i]=new  float[m_nOfRows];
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
				if(m_nOfRows==MIN_NUMBER_OF_ROW)
				{ 
					delete [] m_fArray;
					m_fArray=NULL;
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

void VSWriteVotableOp::printHelp()
//---------------------------------------------------------------------
{
std::cout<<"Produce a VOTable 1.2 with selected field"<<std::endl<<std::endl;
std::cout<<"Usage: VisIVOFilters --op wrvotable [--field column_name] [--force] [--out filename_out.xml] [--history] [--historyfile filename.xml] [--help] [--file] inputFile.bin"<<std::endl<<std::endl;

std::cout<<"Example: VisIVOFilters --op wrvotable --field ra dec u   --out votable_out.xml --file inputFile.bin"<<std::endl;

std::cout<<"Note: "<<std::endl;
std::cout<<"--field Valid columns names. Default value ALL columns will be reported."<<std::endl;
std::cout<<"--force Must be given to force the VOTable creation when the input table has more than "<<MAXROWS<<" rows."<<std::endl;
std::cout<<"--out Output VOTable filename. Default name is given."<<std::endl;
std::cout<<"--history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
std::cout<<"--historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;
std::cout<<"--file  Input table filename."<<std::endl;
std::cout<<"--help produce this output "<<std::endl;

return;
 }
//---------------------------------------------------------------------
bool VSWriteVotableOp::writeVoHeader()
//---------------------------------------------------------------------
{		
time_t rawtime;
struct tm * timeinfo;
char buffer [80];
time ( &rawtime );
timeinfo = localtime ( &rawtime );
strftime (buffer,80,"%Y-%m-%dT%H:%M:%S",timeinfo);


m_fileOutput<<"<?xml version="<<'"'<<"1.0"<<'"'<<" encoding="<<'"'<<"UTF-8"<<'"'<<"?>"<<std::endl;

m_fileOutput<<"<VOTABLE version="<<'"'<<"1.2"<<'"'<< " xmlns:xsi="<<'"'<<"http://www.w3.org/2001/XMLSchema-instance"<<'"'<<std::endl;
m_fileOutput<<"  xmlns="<<'"'<<"http://www.ivoa.net/xml/VOTable/v1.2"<<'"'<<std::endl;

m_fileOutput<<"   xsi:schemaLocation="<<'"'<<"http://www.ivoa.net/xml/VOTable/v1.2 http://www.ivoa.net/xml/VOTable/v1.2"<<'"'<<">"<<std::endl;

m_fileOutput<<"   <DESCRIPTION>"<<std::endl<<"   VisIVO: visivo.oact.inaf.it"<<'\t'<<buffer<<std::endl;	
m_fileOutput<<"   In case of problem, please report to:	"<<'\t'<<" visivo-support@oact.inaf.it"<<std::endl; m_fileOutput<<"   </DESCRIPTION>"<<std::endl;

m_fileOutput<<"  <INFO name="<<'"'<<"rowcount, "<<m_tables[0]->getLocator()<<'"'<<" value="<<'"'<<m_tables[0]->getNumberOfRows()<<'"'<<"/>"<<std::endl;
m_fileOutput<<"  <RESOURCE>"<<std::endl<<"    <TABLE>"<<std::endl;
for(int i=0;i<m_colNumberSet.size();i++)
{
	m_fileOutput<<"      <FIELD datatype="<<'"'<<"float"<<'"'<<" name="<<'"'<<m_tables[0]->getColName(m_colNumberSet[i])<<'"'<<"/>"<<std::endl;
}
m_fileOutput<<"      <DATA>"<<std::endl<<"        <TABLEDATA>"<<std::endl;

return true;
}
//---------------------------------------------------------------------
bool VSWriteVotableOp::execute()
//---------------------------------------------------------------------
{	
	bool allColumns=false;

	
	unsigned long long int fromRow, toRow,startCounter=0;
	unsigned long long int nOfEle,totRows;

	if(getParameterAsString("field").empty() || getParameterAsString("field")=="unknown" )
		allColumns=true;
	

	std::stringstream ssListparameters;
	ssListparameters.str(getParameterAsString("field"));
	std::string paramField="";
	std::string paramFieldGlobal;
	std::vector<unsigned int>::iterator iter;
	if(allColumns)
		for(unsigned int i=0;i<m_tables[0] -> getNumberOfColumns();i++)
			m_colNumberSet.push_back(i);
	else
	{
		m_colNumberSet.clear();
		while (!ssListparameters.eof())
		{
			ssListparameters>>paramField;
			if(m_tables[0] -> getColId(paramField)>=0)
			{ 
				paramFieldGlobal.append(paramField);
				long int test=m_tables[0] -> getColId(paramField);
				m_colNumberSet.push_back(m_tables[0] -> getColId(paramField));
			}
		}

	}


	/*** Output Filename  **/
	std::stringstream fileNameOutputSStream;
	fileNameOutputSStream<<getParameterAsString("out");
	std::string fileNameOutput;

	if(fileNameOutputSStream.str()=="")
	{
		std::string filenameInputTable=m_tables[0]->getLocator();
		int len=filenameInputTable.length();
		fileNameOutputSStream<<filenameInputTable.substr(0, len-4)<<"_votable_"<< paramFieldGlobal<<".xml";
	}

	fileNameOutput=fileNameOutputSStream.str();
	
	m_realOutFilename.push_back(fileNameOutput);
	
//Clean existing tab

	totRows=m_tables[0]->getNumberOfRows();
        nOfEle=totRows;
	if(totRows>MAXROWS)
		if(!isParameterPresent("force"))
		{
			std::cerr<<"Operation aborted. The number of input table rows "<<totRows<<" exceed the maximum allowed "<<MAXROWS<<std::endl;
			std::cerr<<"Use the --force option to force the writing of the VOTable."<<std::endl;
			return false;
		}


	int maxInt=getMaxNumberInt();
	unsigned int maxEle;

	if(nOfEle>maxInt)
		maxEle=maxInt; 
	else
		maxEle=nOfEle;
	

	m_nOfCols = m_colNumberSet.size();
	m_nOfRows=maxEle;
	unsigned int *colList=NULL;	
try
{
	colList=new unsigned int[m_nOfCols];
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


 	bool allocationArray=allocateArray();
	maxEle=m_nOfRows;
	if(m_fArray==NULL ||  !allocationArray )	
	{
		std::cerr<<"Failed Array allocation. Showtable Operation terminated"<<std::endl;
		delete [] colList;
		return false;
	}
	remove(fileNameOutput.c_str());
// write VOTable header

        m_fileOutput.open(fileNameOutput.c_str(), std::ios::out);
	if(!writeVoHeader())
	{
		std::cerr<<"Failed Header Writer. WriteVOTableable operation aborted"<<std::endl;
		delete [] colList;
		return false;
	}
	

/*    	iter = m_colNumberSet.begin();
  	fileOutput.width(width);
 	fileOutput.precision(precision);
 	fileOutput.setf(std::ios::left);
	for(unsigned int i=0;i<m_nOfCols;i++)
	{ 
		colList[i]=*iter;
		fileOutput<<m_tables[0]->getColName(*iter).c_str()<<" "; 
  		fileOutput.width(width);  //QUI chiedere a Marco
  		fileOutput.precision(precision);
		iter++;
	}
	fileOutput<<std::endl;
 	fileOutput.width(width);
 	fileOutput.precision(precision);*/

// write VOTable data
	iter = m_colNumberSet.begin();
	for(unsigned int i=0;i<m_nOfCols;i++)
	{
		colList[i]=*iter;
		iter++;
	}
 	unsigned long long int totEle=nOfEle;

	while(totEle!=0)
	{
		fromRow=startCounter;
		toRow=fromRow+maxEle-1;
		if(toRow>totRows-1)toRow=totRows-1;
	  	m_tables[0]->getColumn(colList, m_nOfCols, fromRow, toRow, m_fArray);
		for(unsigned int j=0;j<(toRow-fromRow+1);j++)
		{
			m_fileOutput<<"<TR>"<<std::endl;
			for(unsigned int i=0;i<m_nOfCols;i++)
			{
				m_fileOutput<<"<TD>"<<m_fArray[i][j]<<"</TD>";
			}
			m_fileOutput<<std::endl<<"</TR>"<<std::endl;
//			fileOutput<<std::endl;
		}
		startCounter=toRow+1;
		totEle=totEle-(toRow-fromRow+1);
	}
	m_fileOutput<<"        </TABLEDATA>"<<std::endl;
	m_fileOutput<<"      </DATA>"<<std::endl;
	m_fileOutput<<"    </TABLE>"<<std::endl;
	m_fileOutput<<"  </RESOURCE>"<<std::endl;
	m_fileOutput<<"</VOTABLE>"<<std::endl;
	m_fileOutput.close();
delete [] colList;
return true;
}

