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
#include <cstdio>
#include <cstring>
#include "vsextractsubvolumeop.h"
#include "vstable.h"
#include <iostream>
#include <sstream>
#ifdef WIN32
	#include <time.h>
#endif
#include "VisIVOFiltersConfigure.h"

const unsigned int VSExtractSubVolumeOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSExtractSubVolumeOp::MIN_NUMBER_OF_ROW = 100;
const unsigned int VSExtractSubVolumeOp::MAX_NUMBER_OF_BYTES = 536870912; // 512x512x512*sizeof(float)

//---------------------------------------------------------------------
VSExtractSubVolumeOp::VSExtractSubVolumeOp()
//---------------------------------------------------------------------
{
	m_fArray=NULL;
	m_cArray=NULL;
}

//---------------------------------------------------------------------
VSExtractSubVolumeOp::~VSExtractSubVolumeOp()
//---------------------------------------------------------------------
{
	if(m_fArray !=NULL)
	{
	 for(int i=0;i<m_nOfCol;i++)
		delete [] m_fArray[i];
	 delete [] m_fArray;
	}
	if(m_cArray !=NULL)
	{
	 for(int i=0;i<m_nOfCol;i++)
		delete [] m_cArray[i];
	 delete [] m_cArray;
	}

}
//---------------------------------------------------------------------
void VSExtractSubVolumeOp::printHelp()
//---------------------------------------------------------------------
{
std::cout<<"Produce a table which represent a subvolume from the original volume"<<std::endl<<std::endl;
std::cout<<"Usage: VisIVOFilters --op extractsubvolume   --startingcell X Y Z --resolution x_res y_res z_res [--field column_names] [--out filename_out.bin] [--help] [--file] inputFile.bin"<<std::endl<<std::endl;

std::cout<<"Example: VisIVOFilters --op extractsubvolume --startingcell 8 8 8 --field Mass Temperature --resolution 16 16 16  --out mysubvolume.bin --file inputFile.bin"<<std::endl<<std::endl;

std::cout<<"Note:"<<std::endl;
std::cout<<"--startingcell X Y Z number of the first cell to be extracted: 0 0 0 is the first cell of the original grid"<<std::endl;
std::cout<<"--resolution  grid size (3D) of the new subgrid"<<std::endl;
std::cout<<"--field valid columns name list to be reported in the new table."<<std::endl;
std::cout<<"--out Name of the new table. Default name is given."<<std::endl;
std::cout<<"--file Input table filename."<<std::endl;

std::cout<<"--help produce this output "<<std::endl;

return;


}
//---------------------------------------------------------------------
bool VSExtractSubVolumeOp::allocateArray()
//---------------------------------------------------------------------
{
unsigned long long int tempLL=getMaxNumberInt();
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
try
{
  m_cArray=new  float*[m_nOfCol];
}
catch(std::bad_alloc &e)
{
  m_cArray=NULL;
}

  if(m_fArray==NULL || m_cArray==NULL)
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
	for(unsigned int i=0;i<m_nOfCol;i++)
	{
try
{
		m_cArray[i]=new  float[m_nExtractOfRow];
}
catch(std::bad_alloc &e)
{
		m_cArray[i]=NULL;
}

		if(m_cArray[i]==NULL) 
		{	
			goodAllocation=false;
			for(unsigned int j=0;j<i;j++) 
				delete [] m_cArray[j];
			for(unsigned int j=0;j<m_nOfCol;j++) 
				delete [] m_cArray[j];
			if(m_nExtractOfRow==MIN_NUMBER_OF_ROW)
			{ 
				delete [] m_fArray;
				delete [] m_cArray;
				m_fArray=NULL;
				m_cArray=NULL;
				return false;
			}
			m_nOfRow=m_nOfRow-MAX_NUMBER_TO_REDUCE_ROW;
			if(m_nOfRow<=MAX_NUMBER_TO_REDUCE_ROW) 
				m_nOfRow=MIN_NUMBER_OF_ROW;
			m_nExtractOfRow=m_nExtractOfRow-MAX_NUMBER_TO_REDUCE_ROW;
			if(m_nExtractOfRow<=MAX_NUMBER_TO_REDUCE_ROW)
				m_nExtractOfRow=MIN_NUMBER_OF_ROW;
			break;
		}

	}

  }
	return true;

}
//---------------------------------------------------------------------
bool VSExtractSubVolumeOp::execute()
//---------------------------------------------------------------------
{
  if(!m_tables[0] -> getIsVolume())
  {
	std::cerr<<"vscoarsevolumeop: the input table is not a volume"<<std::endl;
	return false;
  }
  bool defaultperc= false;
  VSTable tableExtractGrid;
  int extractGrid[3];
  int resGrid[3];
  float dim[3];
  for(int i=0;i<3;i++)
	dim[i]=m_tables[0] -> getCellNumber()[i]; //QUI Mcomp come mai è float?

// starting cell coordinate
  std::stringstream startNameSStream;	
  startNameSStream<<getParameterAsString("startingcell");
  if(startNameSStream.str()==""||startNameSStream.str()=="unknown")
  {	
	std::cerr<<"vsextractsubvolumeop: the startingcell parameter has not valid values names"<<std::endl;
	return false;
   }
  else
  {
	int i=0;
 	while (!startNameSStream.eof())
  	{
		startNameSStream>>extractGrid[i];
	        i++;
		if(i==3) 
		  break;
	}
        if(i<3)
  	{
		std::cerr<<"vsextractsubvolumeop: the startingcell parameter has not valid values names"<<std::endl;
		return false;
  	}
  } 		
  for(int i=0;i<3;i++)
	if(extractGrid[i]<0 || extractGrid[i]>=dim[i])
	{
		std::cerr<<"vsextractsubvolumeop: the startingcell parameter has not valid values names"<<std::endl;
		return false;	
	}
// extractGrid resolution
  std::stringstream resNameSStream;	
  resNameSStream<<getParameterAsString("resolution");
  if(resNameSStream.str()==""||resNameSStream.str()=="unknown")
  {	
	std::cerr<<"vsextractsubvolumeop: the resolution parameter has not valid values names"<<std::endl;
	return false;
   }
  else
  {
	int i=0;
 	while (!resNameSStream.eof())
  	{
		resNameSStream>>resGrid[i];
	        i++;
		if(i==3) 
		  break;
	}  
	if(i<3)
  	{
		std::cerr<<"vsextractsubvolumeop: the resolution parameter has not valid values"<<std::endl;
		return false;
  	}

  } 		
  for(int i=0;i<3;i++)
	if(resGrid[i]<0 ||  resGrid[i]+extractGrid[i]>dim[i])
	{
		std::cerr<<"vsextractsubvolumeop: the resolution parameter has not valid values"<<std::endl;
		return false;	
	}


// ceck list of field
  std::stringstream fieldNameSStream;	
  fieldNameSStream<<getParameterAsString("field");
  if(fieldNameSStream.str()==""||fieldNameSStream.str()=="unknown")
	for(unsigned int i=0;i<m_tables[0] -> getNumberOfColumns();i++)
		m_fieldList.push_back((int) i);
  else
  {
 	while (!fieldNameSStream.eof())
  	{
		std::string paramField;
		fieldNameSStream>>paramField;
		if(m_tables[0] -> getColId(paramField)>=0)
			m_fieldList.push_back(m_tables[0] -> getColId(paramField));
	}
  } 		
  if(m_fieldList.size()==0)
  {
	std::cerr<<"vsextractsubvolumeop: the field parameter has not valid columns names"<<std::endl;
	return false;
  }
   
   m_nOfCol=m_fieldList.size();
   unsigned int *colList;
try
{
	colList=new unsigned int[m_nOfCol];	
}
catch(std::bad_alloc &e)
{
	std::cout<<"vsextractsubvolumeop: allocation failed"<<std::endl;
	return false;
}
   for(int i=0;i<m_nOfCol;i++)
	colList[i]=m_fieldList[i];
   m_nExtractOfRow=(unsigned long long int) resGrid[0];
   m_nExtractOfRow=m_nExtractOfRow*(unsigned long long int) resGrid[1];
   m_nExtractOfRow=m_nExtractOfRow*(unsigned long long int) resGrid[2];

//open file output	
  std::stringstream fileNameOutputSStream;
  fileNameOutputSStream<<getParameterAsString("out");
  std::string fileNameOutput;
  if(fileNameOutputSStream.str()==""||fileNameOutputSStream.str()=="unknown")
  {
	fileNameOutputSStream.str().erase(); //QUI verificare
 	std::string filenameInputTable=m_tables[0]->getLocator();
  	int len=filenameInputTable.length();
	time_t rawtime;
	struct tm * timeinfo;
	char buffer [80];
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	strftime (buffer,80,"%Y%m%d%H%M",timeinfo);
 	fileNameOutputSStream<<filenameInputTable.substr(0, len-4)<<"_extractsubvolume_"<<buffer<<".bin";  //QUI verificare
  }

  fileNameOutput=fileNameOutputSStream.str();
  if(fileNameOutput.find(".bin") == std::string::npos)
     fileNameOutput.append(".bin");
  m_realOutFilename.push_back(fileNameOutput);

//Clean existing tab
  remove(fileNameOutput.c_str());
  tableExtractGrid.setLocator(fileNameOutput);
#ifdef VSBIGENDIAN
  std::string endianism="big";	
#else	
  std::string endianism="little";
#endif

  tableExtractGrid.setEndiannes(endianism);
  tableExtractGrid.setType("float");
  for(int i=0;i<m_fieldList.size();i++)
	tableExtractGrid.addCol(m_tables[0] -> getColName(m_fieldList[i])); 	
  tableExtractGrid.setNumberOfRows(m_nExtractOfRow);
  tableExtractGrid.setIsVolume(true);
  tableExtractGrid.setCellNumber(resGrid[0],resGrid[1],resGrid[2]);
  tableExtractGrid.setCellSize(m_tables[0] -> getCellSize()[0],m_tables[0] -> getCellSize()[1],m_tables[0] -> getCellSize()[2]);
  tableExtractGrid.writeHeader();


   if(m_nExtractOfRow>getMaxNumberInt())
	m_nExtractOfRow=getMaxNumberInt();


   m_nOfRow=m_tables[0] -> getNumberOfRows();
   long long unsigned globalRows=m_nOfRow;

   if(m_nOfRow>getMaxNumberInt())
	m_nOfRow=getMaxNumberInt();

   m_fieldList.empty();
   if(!allocateArray())
   {
	std::cout<<"vsextractsubvolumeop: allocation failed"<<std::endl;
	delete [] colList;
	return false;
    }
   unsigned long long int totRows=m_tables[0]->getNumberOfRows();
   unsigned long long int startCounter=dim[0]*dim[1]*extractGrid[2];
   unsigned long long int fromRow,toRow;
   unsigned long long int fromExtractRow=0,toExtractRow;
   unsigned long long int totEle= totRows;
   int discard1Counter=extractGrid[0];
   int discard2Counter=(dim[0]-1)-(extractGrid[0]+resGrid[0]-1);
   int elementsForPlane=dim[0]*dim[1];
   int totLoadPlane=0;
   while(totEle!=0)
   {
        if(totLoadPlane==resGrid[2])
		break;
	int farrayCounter=0;
	int carrayCounter=0;
	fromRow=startCounter;
	toRow=fromRow+m_nOfRow-1;
	if(toRow>totRows-1)toRow=totRows-1;
// load a fix number of  z planes
	unsigned long long int nZplane=(toRow-fromRow+1)/(unsigned long long int) elementsForPlane;
	toRow=elementsForPlane*nZplane+fromRow-1;
	if(nZplane==0) 
	{
		std::cerr<<"vsextractsubvolumeop: Internal error or wrong table data"<<std::endl;
		delete colList;
		return false;
	}

	m_tables[0] -> getColumn(colList, m_nOfCol, fromRow,toRow,m_fArray);
	for(int l=0;l<nZplane;l++)
        {
 	  totLoadPlane++;
          if(totLoadPlane>resGrid[2])
		break;
/*          if(l+1>resGrid[2])
		break;*/
	   for(int j=0;j<dim[1];j++)
           {	
		if(j<extractGrid[1] || j>=resGrid[1]+extractGrid[1])
		{
			farrayCounter=farrayCounter+dim[0];
		 	continue;
		}
		farrayCounter=farrayCounter+discard1Counter;	
		for(int i=0;i<resGrid[0];i++)
		{
	   		for(int m=0;m<m_nOfCol;m++)
				m_cArray[m][carrayCounter]=m_fArray[m][farrayCounter];
			farrayCounter++;
			carrayCounter++;
		}
		farrayCounter=farrayCounter+discard2Counter;
	      
	    }
/*	    std::clog<<l<<" "<<totLoadPlane<<" "<<farrayCounter<<" "<<carrayCounter<<std::endl;*/
        }
        toExtractRow=carrayCounter-1+fromExtractRow;
       	tableExtractGrid.putColumn(colList,m_nOfCol, fromExtractRow,toExtractRow,m_cArray);  
	fromExtractRow=toExtractRow+1;
	totEle=totEle-(toRow-fromRow+1);
	startCounter=toRow+1;
	if(startCounter>=globalRows)totEle=0;
	if(totEle<0) totEle=0;
  }
  delete [] colList;

return true;
}
