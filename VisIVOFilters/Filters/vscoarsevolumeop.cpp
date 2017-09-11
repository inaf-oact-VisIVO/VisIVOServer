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
#include "vscoarsevolumeop.h"
#include "vstable.h"
#include "vsselectcolumnsop.h"
#include <map>
#include <iostream>
#include <sstream>
#ifdef WIN32
	#include <time.h>
#endif
#include "VisIVOFiltersConfigure.h"

const unsigned int VSCoarseVolumeOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSCoarseVolumeOp::MIN_NUMBER_OF_ROW = 100;
const unsigned int VSCoarseVolumeOp::MAX_NUMBER_OF_BYTES = 67108864; // 256*256*256*sizeof(float)

//---------------------------------------------------------------------
VSCoarseVolumeOp::VSCoarseVolumeOp()
//---------------------------------------------------------------------
{
	m_fArray=NULL;
	m_cArray=NULL;
}

//---------------------------------------------------------------------
VSCoarseVolumeOp::~VSCoarseVolumeOp()
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
void VSCoarseVolumeOp::printHelp()
//---------------------------------------------------------------------
{
std::cout<<"Produce a coarsed subvolume with plane extraction from the original volume "<<std::endl<<std::endl;
std::cout<<"Usage: VisIVOFilters --op coarsevolume  [--perc percentage] [--newres x_res y_res z_res]  [--field column_names] [--out filename_out.bin] [--history] [--historyfile filename.xml] [--help] [--file] inputFile.bin"<<std::endl<<std::endl;

std::cout<<"Example: VisIVOFilters --op coarsevolume  --perc 10.0 --field Mass Temperature   --out subvolume.bin --file inputFile.bin"<<std::endl;

std::cout<<"Note: --help produce this output "<<std::endl;
std::cout<<"      --file input table filename. A volume table must be given."<<std::endl;
std::cout<<"      --perc a percentage (from 0.0 to 100.0) subvolume  will be produced. Default value produce a  sub volume that  could be directly uploaded and visualised with VisIVOServer and VisIVODesktop applications."<<std::endl;
std::cout<<"      --newres a subvolume with new resolution will be produced. No default is given. This parameter is ignored if --perc option is given."<<std::endl;
std::cout<<"      --field list of columns contained in the original file. Default:  all columns will be extracted."<<std::endl;

std::cout<<"      --out filename of output result. Default name is given."<<std::endl;
std::cout<<"      --history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
std::cout<<"      --historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;

return;


}
//---------------------------------------------------------------------
bool VSCoarseVolumeOp::selectColumnsOnly()
//---------------------------------------------------------------------
{  
  std::map<std::string, std::string> appParameters;
  std::stringstream tableSStream;
  tableSStream<<getParameterAsString("file");
  if(tableSStream.str()==""||tableSStream.str()=="unknown")
  {
	std::cerr<<"vscoarsevolumeop: table filename is not given"<<std::endl;
		return false;
  }
  std::string fileNameBin;
  tableSStream>>fileNameBin;
  if(fileNameBin.find(".bin") == std::string::npos)
	    fileNameBin.append(".bin");
  VSTable table(fileNameBin);
  std::stringstream valueSStream;
  std::string temp=" ";
  for(int i=0;i<m_fieldList.size();i++)
	valueSStream <<" "<<m_tables[0] -> getColName(m_fieldList[i]);
  appParameters.insert(make_pair("--list", valueSStream.str()));
  appParameters.insert(make_pair("--extract", temp));
			
  VSSelectColumnsOp op;
  op.setParameters(appParameters);
  op.addInput(&table);
  op.execute();
  return true;
}
//---------------------------------------------------------------------
bool VSCoarseVolumeOp::allocateArray()
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
		m_cArray[i]=new  float[m_nCoarseOfRow];
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
			if(m_nCoarseOfRow==MIN_NUMBER_OF_ROW)
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
			m_nCoarseOfRow=m_nCoarseOfRow-MAX_NUMBER_TO_REDUCE_ROW;
			if(m_nCoarseOfRow<=MAX_NUMBER_TO_REDUCE_ROW)
				m_nCoarseOfRow=MIN_NUMBER_OF_ROW;
			break;
		}

	}

  }
	return true;

}
//---------------------------------------------------------------------
bool VSCoarseVolumeOp::execute()
//---------------------------------------------------------------------
{
		
  if(!m_tables[0] -> getIsVolume())
  {
	std::cerr<<"vscoarsevolumeop: the input table is not a volume"<<std::endl;
	return false;
  }
  bool defaultperc= false;
  bool newres= false;
  VSTable tableCoarseGrid;
  int coarseGrid[3];
  float perc=100.0;
  float dim[3];
  int **discardPlanes;

  for(int i=0;i<3;i++)
	dim[i]=m_tables[0] -> getCellNumber()[i]; //QUI Mcomp come mai è float?

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
	std::cerr<<"vscoarsevolumeop: the field parameter has not valid columns names"<<std::endl;
	return false;
  }
	
//ceck percentage to be reduced
  std::stringstream percSStream;
  percSStream<<getParameterAsString("perc");
  if(percSStream.str()==""||percSStream.str()=="unknown")
	defaultperc=true;
  else
  {
 	percSStream>> perc;
	if(perc>95.0) perc=95.0;
	if(perc<0.01) perc=0.01;
  }

// ceck for newres
  if(defaultperc)
  {	
  	std::stringstream resSStream;	
  	resSStream<<getParameterAsString("newres");
  	if(resSStream.str()!="" && resSStream.str()!="unknown")
        {
 		int i=0;
		while (!resSStream.eof())
  		{
			resSStream>>coarseGrid[i];
			if(coarseGrid[i]>dim[i])
			   coarseGrid[i]=dim[i];
			i++;
			if(i==3)
			   break;
		}
		if(i<3)
			std::cerr<<"vscoarsevolumeop: newres parameter contain invalid values. The parameter is ignored."<<std::endl;
		else
			newres=true;
		if(coarseGrid[0]>=dim[0]-1 || coarseGrid[1]>=dim[1]-1 || coarseGrid[2]>=dim[2]-1)
		{
			selectColumnsOnly();
			return true;	
		}
	}
  }
//  set for perc given or default perc
  if(defaultperc && !newres)	
  {
	unsigned long long int totalBytes=(unsigned long long int) (m_fieldList.size())*m_tables[0] ->getNumberOfRows()*sizeof(float);
	if(totalBytes>MAX_NUMBER_OF_BYTES*0.95)  // not operated if not lower than 0.95 of MAX_...
	{
		float yx=dim[1]/dim[0];		
		float zx=dim[2]/dim[0];		
		float temp;
		temp=(((float)(MAX_NUMBER_OF_BYTES)/sizeof(float))/(1+yx+zx));
		temp=temp/m_fieldList.size();
		coarseGrid[0]=(int) temp;
		coarseGrid[1]=(int) coarseGrid[0]*yx;
		coarseGrid[2]=(int) coarseGrid[0]*zx;
	} else
	{
		selectColumnsOnly();
		return true;
	}
  } 
  if(!defaultperc)
  {
	if(perc >= 95.0)
	{
		selectColumnsOnly();
		return true;
	}
	float yx=dim[1]/dim[0];		
	float zx=dim[2]/dim[0];		
	float temp;
	perc=perc/100;
	temp=(dim[0]*dim[1]*dim[2]*perc)/(1*yx*zx);
//search for cubic value fitting the number
	int xCoarse;
	for(xCoarse=dim[0];xCoarse>0;xCoarse--)
		if((xCoarse*xCoarse*xCoarse) < temp)
		   break;
	if(xCoarse <2)
  	{
		std::cerr<<"vscoarsevolumeop: the requested percentage value is too low"<<std::endl;
		return false;
  	}

	coarseGrid[0]=(int) xCoarse;
	coarseGrid[1]=(int) coarseGrid[0]*yx;
	coarseGrid[2]=(int) coarseGrid[0]*zx;
   }
	
// start plane Subtraction phase

   m_nOfCol=m_fieldList.size();
   for (int i=0;i<3;i++)
   	m_nOfDiscardPlanes[i]=m_tables[0] -> getCellNumber()[i] - coarseGrid[i];

// fill and check for discarded planes, and eventually correct m_nOfDiscardPlanes and coarseGrid values
  discardPlanes=new  int*[3];
  for(int i=0;i<3;i++) 
	discardPlanes[i]=new  int[m_nOfDiscardPlanes[i]];

   for (int i=0;i<3;i++)
   {
	float nDiv;
	nDiv= ((float)(dim[i]))/((float)(m_nOfDiscardPlanes[i]+1));
	int k=0;
	for(int j=0;j<m_nOfDiscardPlanes[i];j++)
	{	
		if(j>0)
		{
		   int temp=(int)(nDiv+nDiv*j);
		   if(temp != discardPlanes[i][k-1])
		   {
			discardPlanes[i][k]=temp;
			k++;
//			std::clog<<"discardPlanes["<<i<<"]["<<j<<"]="<<discardPlanes[i][j]<<std::endl;
		   }
		} else  //execeute the following only the first time
		{
			discardPlanes[i][k]=(int)(nDiv+nDiv*j);
			k++;
//			std::clog<<"discardPlanes["<<i<<"]["<<j<<"]="<<discardPlanes[i][j]<<std::endl;
		}
	}
        m_nOfDiscardPlanes[i]=k;
	coarseGrid[i]=dim[i]-k;
   }	

   unsigned int *colList;
try
{
	colList=new unsigned int[m_nOfCol];	
}
catch(std::bad_alloc &e)
{
	std::cout<<"vscoarsevolumeop: allocation failed"<<std::endl;
	return false;
}
   for(int i=0;i<m_nOfCol;i++)
	colList[i]=m_fieldList[i];
   m_nCoarseOfRow=(unsigned long long int) coarseGrid[0];
   m_nCoarseOfRow=m_nCoarseOfRow*(unsigned long long int) coarseGrid[1];
   m_nCoarseOfRow=m_nCoarseOfRow*(unsigned long long int) coarseGrid[2];
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
 	fileNameOutputSStream<<filenameInputTable.substr(0, len-4)<<"_coarsevolume_"<<buffer<<".bin";  //QUI verificare
  }

  fileNameOutput=fileNameOutputSStream.str();
  if(fileNameOutput.find(".bin") == std::string::npos)
	    fileNameOutput.append(".bin");
  m_realOutFilename.push_back(fileNameOutput);

//Clean existing tab
  remove(fileNameOutput.c_str());
  tableCoarseGrid.setLocator(fileNameOutput);
#ifdef VSBIGENDIAN
  std::string endianism="big";	
#else	
  std::string endianism="little";
#endif

  tableCoarseGrid.setEndiannes(endianism);
  tableCoarseGrid.setType("float");
  for(int i=0;i<m_fieldList.size();i++)
	tableCoarseGrid.addCol(m_tables[0] -> getColName(m_fieldList[i])); 	
  tableCoarseGrid.setNumberOfRows(m_nCoarseOfRow);
  tableCoarseGrid.setIsVolume(true);
  tableCoarseGrid.setCellNumber(coarseGrid[0],coarseGrid[1],coarseGrid[2]);
  tableCoarseGrid.setCellSize(m_tables[0] -> getCellSize()[0],m_tables[0] -> getCellSize()[1],m_tables[0] -> getCellSize()[2]);
  tableCoarseGrid.writeHeader();


   if(m_nCoarseOfRow>getMaxNumberInt())
	m_nCoarseOfRow=getMaxNumberInt();


   m_nOfRow=m_tables[0] -> getNumberOfRows();
   long long unsigned globalRows=m_nOfRow;
   if(m_nOfRow>getMaxNumberInt())
	m_nOfRow=getMaxNumberInt();

   m_fieldList.empty();
   if(!allocateArray())
   {
	std::cout<<"vscoarsevolumeop: allocation failed"<<std::endl;
	delete [] colList;
	return false;
   }
   
   unsigned long long int totRows=m_tables[0]->getNumberOfRows();
   unsigned long long int startCounter=0;
   unsigned long long int fromRow,toRow;
   unsigned long long int fromCoarseRow=0,toCoarseRow;
   unsigned long long int totEle= totRows;
   int zplaneCounter=0;
   int yrowCounter=0;
   int xeleCounter=0;
   int arrayCounter=0;
   int elementsForPlane=dim[0]*dim[1];
   unsigned int j;
   int loadedPlanes=0;
   int examinPlane=-1;
   while(totEle!=0)
   {
	j=0;
	arrayCounter=0;
	fromRow=startCounter;
	toRow=fromRow+m_nOfRow-1;
	if(toRow>totRows-1)toRow=totRows-1;
// load a fix number of  z planes
	unsigned long long int nZplane=(toRow-fromRow+1)/(unsigned long long int) elementsForPlane;
	toRow=elementsForPlane*nZplane+fromRow-1;
	int maxfArrayEle=toRow-fromRow;
	if(nZplane==0) 
	{
		std::cerr<<"vscoarsevolume: Internal error or wrong table data"<<std::endl;
		delete colList;
		return false;
	}
	
	m_tables[0]->getColumn(colList, m_nOfCol, fromRow,toRow,m_fArray);
	loadedPlanes=loadedPlanes+nZplane;
//	std::clog<<"LoadedPlanes= "<<loadedPlanes<<" totEle= "<<totEle<<std::endl;

	for(int zplane=0;zplane<nZplane;zplane++)
	{	
		examinPlane++;
		if(discardPlanes[2][zplaneCounter] == examinPlane)
		{
			arrayCounter=arrayCounter+elementsForPlane;
	 		zplaneCounter++;
			continue;
		}
		yrowCounter=0;
		for(int yrow=0;yrow<dim[1];yrow++)
		{
			if(discardPlanes[1][yrowCounter] == yrow)
			{
				arrayCounter=arrayCounter+dim[0];
	 			yrowCounter++;
				continue;
			}
			xeleCounter=0;
			for(int xele=0;xele<dim[0];xele++)
			{
				if(discardPlanes[0][xeleCounter] == xele)
				{
					arrayCounter++;
	 				xeleCounter++;
				} else
				{
					for(int i=0;i<m_nOfCol;i++)
						m_cArray[i][j]=m_fArray[i][arrayCounter];
					j++;
					arrayCounter++;
				}
			}
		}
		if(arrayCounter>maxfArrayEle)
			break;
	}
       	if(j>0) //there are new computed elements to write
	{
		toCoarseRow=fromCoarseRow+j-1;
       		tableCoarseGrid.putColumn(colList,m_nOfCol, fromCoarseRow,toCoarseRow,m_cArray);  
		fromCoarseRow=toCoarseRow+1;
	}
	totEle=totEle-(toRow-fromRow+1);
	startCounter=toRow+1;
	if(startCounter>=globalRows)totEle=0;
	if(totEle<0) totEle=0;
  }
  delete [] colList;
  return true;
}
