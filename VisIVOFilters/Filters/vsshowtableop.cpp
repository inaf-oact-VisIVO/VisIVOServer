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
#include "vsshowtableop.h"

const unsigned int VSShowTableOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSShowTableOp::MIN_NUMBER_OF_ROW = 100;
const unsigned int VSShowTableOp::WIDTH = 10;
const unsigned int VSShowTableOp::PRECISION = 8;

VSShowTableOp::VSShowTableOp()
{
 m_fArray=NULL;
 m_nOfRow=0;
 m_nOfCol=0;	

}


VSShowTableOp::~VSShowTableOp()
{
 	if(m_fArray!=NULL)
		for(unsigned int i=0;i<m_nOfCol;i++)
		{
			if(m_fArray[i]!=NULL) delete [] m_fArray[i];
		}
	if(m_fArray!=NULL) delete [] m_fArray;

}
//---------------------------------------------------------------------
bool VSShowTableOp::allocateArray()
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
//---------------------------------------------------------------------

void VSShowTableOp::printHelp()
//---------------------------------------------------------------------
{
std::cout<<"Produce an ASCII table with selected field of the first number of rows as specified in the --numrows parameter "<<std::endl<<std::endl;
std::cout<<"Usage: VisIVOFilters --op showtable [--field column_name] [--numrows num_of_rows] [--rangerows fromRow toRow] [--width format_width] [--precision format_precision] [--out filename_out.txt] [--help] [--file] inputFile.bin"<<std::endl<<std::endl;

std::cout<<"Example: VisIVOFilters --op showtable --field X Y Z --numrows 100  --out filename_out.txt --file inputFile.bin"<<std::endl;

std::cout<<"Note: "<<std::endl;
std::cout<<"--field Valid columns names. Default value ALL columns will be reported."<<std::endl;
std::cout<<"--numrows number of first Rows in the Ascii output file. Default value is equal to the number of Rows of the input table."<<std::endl;
std::cout<<"--rangerows Rows range of the inputFile that will be reported in the ascii output file. First row is 0 . Default range is equal to the all rows of the input table. It is ignored if numrows is specified."<<std::endl;
std::cout<<"--width field width in the Ascii output file. Default value is given."<<std::endl;
std::cout<<"--precision field precision in the Ascii output file. Default value is given."<<std::endl;
std::cout<<"--out Output ascii filename. Default name is given."<<std::endl;
std::cout<<"--file  Input table filename."<<std::endl;
std::cout<<"--help produce this output "<<std::endl;

return;
 }

//---------------------------------------------------------------------
bool VSShowTableOp::execute()
//---------------------------------------------------------------------
{	
	bool allColumns=false;

	
	int width=getParameterAsInt("width");
	int precision=getParameterAsInt("precision");
	unsigned long long int extrFromRow, extrToRow;
	unsigned long long int fromRow, toRow, startCounter=0;
	unsigned long long int nOfEle;

	if(width<=0) width=WIDTH;
	if(precision<=0) precision=PRECISION;

	if(getParameterAsString("field").empty() || getParameterAsString("field")=="unknown" )
		allColumns=true;
	

	std::vector<unsigned int> colNumberSet;
	std::stringstream ssListparameters;
	ssListparameters.str(getParameterAsString("field"));
	std::string paramField="";
	std::string paramFieldGlobal;
	std::vector<unsigned int>::iterator iter;
	if(allColumns)
		for(unsigned int i=0;i<m_tables[0] -> getNumberOfColumns();i++)
			colNumberSet.push_back(i);
	else
	{
		colNumberSet.clear();
		while (!ssListparameters.eof())
		{
			ssListparameters>>paramField;
//			std::clog<<paramField<<std::endl;

			if(m_tables[0] -> getColId(paramField)>=0)
			{ 
				paramFieldGlobal.append(paramField);
				colNumberSet.push_back(m_tables[0] -> getColId(paramField));
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
		fileNameOutputSStream<<filenameInputTable.substr(0, len-4)<<"_showtable_"<< paramFieldGlobal<<".ascii";
	}

	fileNameOutput=fileNameOutputSStream.str();
	m_realOutFilename.push_back(fileNameOutput);

//Clean existing tab

	unsigned long long int totRows=(unsigned long long int)getParameterAsInt("numrows");
	if(totRows==0)
	{ 
		totRows=m_tables[0]->getNumberOfRows();
		std::stringstream range;
		range<<getParameterAsString("rangerows");
		if(range.str()=="")  //QUI verificare
                  nOfEle=totRows;
		else
		{
			range >> extrFromRow;
			range >> extrToRow;
			if(extrToRow>m_tables[0]->getNumberOfRows()-1)
				extrToRow=m_tables[0]->getNumberOfRows()-1;
			if(extrToRow<extrFromRow)
			{
				std::cerr<<"vsshowtableop: Invalid rangerows parameters"<<std::endl;
				return false;
				
			}
			startCounter=extrFromRow;
			nOfEle=extrToRow-extrFromRow+1;
		}
	}else
	{	if(totRows>m_tables[0]->getNumberOfRows())
			totRows=m_tables[0]->getNumberOfRows();
		nOfEle=totRows;
	}
	int maxInt=getMaxNumberInt();
	unsigned int maxEle;

	if(nOfEle>maxInt)
		maxEle=maxInt; 
	else
		maxEle=nOfEle;
	

	m_nOfCol = colNumberSet.size();
	m_nOfRow=maxEle;
	unsigned int *colList=NULL;	
try
{
	colList=new unsigned int[m_nOfCol];
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
	maxEle=m_nOfRow;
	if(m_fArray==NULL ||  !allocationArray )	
	{
		std::cerr<<"Failed Array allocation. Showtable Operation terminated"<<std::endl;
		delete [] colList;
		return false;
	}
	remove(fileNameOutput.c_str());

        std::ofstream fileOutput(fileNameOutput.c_str(), std::ios::out);
  	iter = colNumberSet.begin();
  	fileOutput.width(width);
 	fileOutput.precision(precision);
 	fileOutput.setf(std::ios::left);
	for(unsigned int i=0;i<m_nOfCol;i++)
	{ 
		colList[i]=*iter;
		fileOutput<<m_tables[0]->getColName(*iter).c_str()<<" "; 
  		fileOutput.width(width);  //QUI chiedere a Marco
  		fileOutput.precision(precision);
		iter++;
	}
	fileOutput<<std::endl;
 	unsigned long long int totEle=nOfEle;
 	fileOutput.width(width);
 	fileOutput.precision(precision);
	while(totEle!=0)
	{
		fromRow=startCounter;
		toRow=fromRow+maxEle-1;
		if(toRow>totRows-1)toRow=totRows-1;
	  	m_tables[0]->getColumn(colList, m_nOfCol, fromRow, toRow, m_fArray);
		for(unsigned int j=0;j<(toRow-fromRow+1);j++)
		{
			for(unsigned int i=0;i<m_nOfCol;i++)
			{
				fileOutput<<m_fArray[i][j]<<" ";
  				fileOutput.width(width);  //QUI chiedere a Marco
  				fileOutput.precision(precision);
			}
			fileOutput<<std::endl;
		}
		startCounter=toRow+1;
		totEle=totEle-(toRow-fromRow+1);
	}
	fileOutput.close();
delete [] colList;
return true;
}

