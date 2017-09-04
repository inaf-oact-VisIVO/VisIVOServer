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
#include <sstream>
#include <iostream>
#include <fstream>
//#include <ctime>
#include "time.h"
#include "vstable.h"
#include "vsvisualop.h"
#include "VisIVOFiltersConfigure.h"

const unsigned int VSVisualOp::DEFAULT_ROW_TO_VISUALIZE = 8000000;
const unsigned int VSVisualOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSVisualOp::MIN_NUMBER_OF_ROW = 100;

//---------------------------------------------------------------------
VSVisualOp::VSVisualOp()
//---------------------------------------------------------------------
{
   m_visualSize=DEFAULT_ROW_TO_VISUALIZE;
   m_colSize=0;  
   m_globalNumberOfRow=0;
   m_fVisualArray=NULL;
   m_fArray=NULL;
}
//---------------------------------------------------------------------
VSVisualOp::~VSVisualOp()
//---------------------------------------------------------------------
{
	if(m_fArray !=NULL)
	{
		if(m_fArray[0]!=NULL) delete [] m_fArray[0];
	 	delete [] m_fArray;
	}
	if(m_fVisualArray !=NULL)
	{
		for(int i=0;i<m_colSize;i++)
			if(m_fVisualArray[i]!=NULL) delete [] m_fVisualArray[i] ;
	 	delete [] m_fVisualArray;
	}

}
//---------------------------------------------------------------------
bool VSVisualOp::allocatefArray()
//---------------------------------------------------------------------
{

	m_fArray=new  float*[1];
try
{
	m_fVisualArray=new  float*[m_colSize];
}
catch(std::bad_alloc &e)
{
	m_fVisualArray=NULL;
}

	if(m_fVisualArray==NULL)
		return false;
	for(unsigned int i=0;i<m_colSize;i++)
	{
try
{
		m_fVisualArray[i]=new  float[m_visualSize];
}
catch(std::bad_alloc &e)
{
	m_fVisualArray[i]=NULL;
}

		if(m_fVisualArray[i]==NULL) 
		{	
			for(unsigned int j=0;j<i;j++) 
				delete [] m_fVisualArray[j];
			delete [] m_fVisualArray;
			return false;
		}
	}

	m_maxEle=getMaxNumberInt();
	if(m_maxEle>m_globalNumberOfRow)m_maxEle=m_globalNumberOfRow;

try
{
	m_fArray[0]=new float[m_maxEle];
}
catch(std::bad_alloc &e)
{
	m_fArray[0]=NULL;
}

	while (m_fArray[0]==NULL) 
	{			
		if(m_maxEle==MIN_NUMBER_OF_ROW)
		{
 			for(unsigned int i=0;i<m_colSize;i++)
				delete [] m_fVisualArray[i];
			delete [] m_fVisualArray;
			return false;
		}
		if(m_maxEle<=MAX_NUMBER_TO_REDUCE_ROW) 
			m_maxEle=MIN_NUMBER_OF_ROW;
		else
			m_maxEle=m_maxEle-MAX_NUMBER_TO_REDUCE_ROW;
try
{
		m_fArray[0]=new float[m_maxEle];
}
catch(std::bad_alloc &e)
{
	m_fArray[0]=NULL;
}

	}
	return true;
}
//---------------------------------------------------------------------
void VSVisualOp::printHelp()
//---------------------------------------------------------------------
{
	std::cout<<"Create an eventually randomized new table from  one or more input tables. All the input table must have the same number of rows. The new table can be used with VisIVOViewer. Cannot be applied to volume tables."<<std::endl<<std::endl;

	std::cout<<"Usage: VisIVOFilters --op visual  [--size number_of_elemnts] [--out filename_out.bin]  [--help] [--filelist] tab_selection_file.txt"<<std::endl;

	std::cout<<"Example: VisIVOFilters --op visual --out filename_out.bin --filelist tab_selection_file.txt"<<std::endl;

	std::cout<<"Note:"<<std::endl;
	std::cout<<"--size. Number of max rows in output table. Default is the minimum between 8000000 of rows and the rows of  input tables. Input table must have the same number of rows."<<std::endl;
	std::cout<<"--out Output table filename. Default name is given."<<std::endl;
	std::cout<<"--filelist  Input text file with a list of tables and columns."<<std::endl;
	std::cout<<"--help produce this output."<<std::endl;

	return;
}

//---------------------------------------------------------------------
bool VSVisualOp::execute()
//---------------------------------------------------------------------
// totRows is the minimum or the maximum number of rows among all tables
// 
// maxEle maximum number allowed to upload/download data from a generictable: It is the lower value between totRows, maxInt and fArray allocation fArray[0][maxEle]
// 
// nOfEle total number of ellement to upload for the specific input table. nOfEle is the lower value between the rows of the table and totRows
// 
// totEle is a counter: the numebr of elemnts that I have still to upload for the specific table. Starts from nOfEle, its value devrese of maxEle each cycle
// 
// nPadElements is the number of rows to pad (if any) for the specific input table

{
	std::string filename;
	int value=-1,tableValue;
	std::map<std::string, int>::iterator p;
	unsigned int totColumns=0;
	unsigned long long int totRows=0; 
	VSTable tableVisual;

	if(isParameterPresent("filelist"))
	{
	  if(getParameterAsString("filelist").empty() || getParameterAsString("filelist")=="unknown" )
	  {
		std::cerr<<"VSVisualOp: no file with table list is given"<<std::endl;
		return false;
	  } else
	  {
		filename=getParameterAsString("filelist");
	  }
	}else //obsolete
	{
	  if(getParameterAsString("file").empty() || getParameterAsString("file")=="unknown" )
	  {
		std::cerr<<"VSVisualOp: no file with table list is given"<<std::endl;
		return false;
	  } else
	  {
		filename=getParameterAsString("file");
	  }
	}
	std::ifstream fileInput(filename.c_str());
	if(!fileInput)
	{
		std::cerr<<"VSVisualOp: cannot open table list file"<<filename<<std::endl;
		return false;
	}

	if(!getParameterAsString("size").empty() && getParameterAsString("size")!="unknown" )
		m_visualSize=getParameterAsInt("size");
	if(m_visualSize<=0)	
	{
		std::cerr<<"VSVisualOp: invalid size "<<filename<<std::endl;
		return false;
	}
	if(m_visualSize>DEFAULT_ROW_TO_VISUALIZE) 
		std::cerr<<"Visual Operation. Visualization could fail. Lower the size of visualization to "<<DEFAULT_ROW_TO_VISUALIZE<<std::endl;

	int progressiveOutName=1;
	while(!fileInput.eof())
	{
		std::string  key,stringColumn;
		fileInput >> key;
		fileInput >> stringColumn;
		if(key!="")
		{
			VSTable testTable(key);
			if(testTable.tableExist())
			{
				p=m_listOfTables.find(key);
				if(p ==m_listOfTables.end()) //QUI verifica
				{
					value++;

					if(value==0) 
					{ 
						m_globalNumberOfRow=testTable.getNumberOfRows();
						if(m_globalNumberOfRow<m_visualSize)
						m_visualSize=m_globalNumberOfRow;		
					}
					if(value>0 && testTable.getNumberOfRows() != m_globalNumberOfRow)
					{ 
						std::cerr<<"Visual Operation cannot be executed: different value of rows in tables"<<std::endl;
						return false;
					}	

					if(testTable.getColId(stringColumn)>=0)
					{
						m_listOfTables.insert(make_pair(key,value));
						m_listColumns[value].push_back(testTable.getColId(stringColumn));
						std::stringstream sstmp;
						sstmp<<stringColumn<<"_visual_"<<progressiveOutName;
						progressiveOutName++;
						std::string test=sstmp.str();
						m_listColumnsOutName[value].push_back(sstmp.str());
					}
					if(stringColumn=="*")
					{
						m_listOfTables.insert(make_pair(key,value));
						unsigned int numberOfCols=testTable.getNumberOfColumns();
						for(int i=0;i<numberOfCols;i++)
						{
							m_listColumns[value].push_back(i);
							std::stringstream sstmp;
							sstmp<<testTable.getColName(i)<<"_visual_"<<progressiveOutName;
							progressiveOutName++;
							std::string test=sstmp.str();
							m_listColumnsOutName[value].push_back(sstmp.str());
						}
					}

					if(value==101) return false;

					continue;
				}else  //Table already added
				{
					int tableValue=p->second;
					int colId;
					colId=testTable.getColId(stringColumn);
					if(colId >= 0)
					{
						bool colExist=false;
						for(int i=0;i<m_listColumns[tableValue].size();i++)
						{
 						 	if(m_listColumns[tableValue][i]==colId)
							{
							  colExist=true;
							  break;
							}
						}
						if(!colExist)
						{
							m_listColumns[tableValue].push_back(colId);
							std::stringstream sstmp;
							sstmp<<stringColumn<<"_visual_"<<progressiveOutName;
							progressiveOutName++;
							std::string test=sstmp.str();
							m_listColumnsOutName[tableValue].push_back(sstmp.str());
						}
					}
					if(stringColumn=="*")
					{
						unsigned int numberOfCols=testTable.getNumberOfColumns();
						m_listColumns[tableValue].clear();
						for(unsigned int i=0;i<numberOfCols;i++)
						{
							m_listColumns[tableValue].push_back(i);
							std::stringstream sstmp;
							sstmp<<testTable.getColName(i)<<"_visual_"<<progressiveOutName;
							progressiveOutName++;
							std::string test=sstmp.str();
							m_listColumnsOutName[tableValue].push_back(sstmp.str());

						}
					}
				}
			}				
		}
	}
	fileInput.close();
	for(int i=0;i<=value;i++)
		totColumns=totColumns+m_listColumns[i].size();
	if(totColumns==0) return false;
	m_colSize=totColumns;
	if(m_colSize<3) m_colSize=3;
	bool allocationfArray=allocatefArray();
	if(m_fArray==NULL ||m_fVisualArray==NULL ||  !allocationfArray )	
	{
		std::cerr<<"Visual Op: Failed fArray allocation"<<std::endl;
		return false;
	}
//Set Null NONE columns
	if(m_colSize>totColumns)
	{	
		for( int j=totColumns+1;j<m_colSize;j++)
		   for(int i=0;i<m_visualSize;i++)
			m_fVisualArray[j][i]=0.0;
	}	
//Set Output Table Filename
	std::string fileNameOutput;

	if(getParameterAsString("out").empty()||getParameterAsString("out")=="unknown")
		fileNameOutput="visualTable.bin";
	else
		fileNameOutput=getParameterAsString("out");

	m_realOutFilename.push_back(fileNameOutput);
	
//Clean existing tab
	remove(fileNameOutput.c_str());

	tableVisual.setLocator(fileNameOutput);
#ifdef VSBIGENDIAN
	std::string endianism="big";
	
#else	
	std::string endianism="little";
#endif

	tableVisual.setEndiannes(endianism);
 	tableVisual.setType("float");

	unsigned long long int nOfEle;
	unsigned int colVisualize=0;
	unsigned long long int fromRow, toRow, visualEle, startCounter;
	int iran,rran;
	unsigned int nOfCol=1;
	unsigned int colList[1];
	unsigned int listElement;
	unsigned long long int nTotSkip=(m_globalNumberOfRow/m_visualSize)/m_maxEle;

	for(p=m_listOfTables.begin();p!=m_listOfTables.end();p++)
	{
		std::string locator=p->first;
		if(locator!="")
		{
			VSTable testTable(locator);
			tableValue=p->second;
			unsigned int size=m_listColumns[tableValue].size();
			for(unsigned int i=0;i<size;i++)
			{	
				colList[0]=m_listColumns[tableValue][i];
				startCounter=0;
				unsigned long long int totEle=m_globalNumberOfRow;
				visualEle=0;
				srand(0);
			        unsigned long long int nSkip=nTotSkip;
				while(totEle!=0)
				{
					fromRow=startCounter;
					toRow=fromRow+m_maxEle-1;
					if(toRow>m_globalNumberOfRow-1)toRow=m_globalNumberOfRow-1;
					if(nSkip==0)
					{
	  					testTable.getColumn(colList, nOfCol, fromRow, toRow, m_fArray);
						int numberOfSampleForRead=(int) ((toRow-fromRow+1)*m_visualSize/m_globalNumberOfRow);
						if(numberOfSampleForRead==0) 
							numberOfSampleForRead=1; //could be considered a BUG
						int numberOfValueChanc=(int) ((toRow-fromRow+1)/numberOfSampleForRead);
//fill m_fVisualArray
						for(unsigned int k=0; k<numberOfSampleForRead;k++)
						{
						  iran=rand();
        					  rran= (int)(numberOfValueChanc* (float)(iran) / RAND_MAX);
        					  listElement = rran+numberOfValueChanc*k;  // QUI verifica
						  if(listElement>=m_maxEle) 
							listElement=m_maxEle;
						  m_fVisualArray[colVisualize][visualEle]=m_fArray[0][listElement];
						  visualEle++;
						  if(visualEle==m_visualSize)
						  {
							totEle=toRow-fromRow+1; //put totEle=0 next rows 
							break;
						  }
						}
						nSkip=nTotSkip;
					}
					totEle=totEle-(toRow-fromRow+1);
					startCounter=toRow+1;
					if(nSkip>0) nSkip--;
				}
				colVisualize++;
			}
		}
		
	}
	tableVisual.setNumberOfRows(visualEle);
	std::string fileWebName=fileNameOutput+"_visualWeb";
	std::ofstream fileWeb(fileWebName.c_str());
	for(p=m_listOfTables.begin();p!=m_listOfTables.end();p++)
		fileWeb<<p->first<<" "<<p->second<<std::endl;
	fileWeb.close();
		
	for(p=m_listOfTables.begin();p!=m_listOfTables.end();p++)
	{
		std::string locator=p->first;
		if(locator!="")
		{
			VSTable testTable(locator);
			tableValue=p->second;
			unsigned int size=m_listColumnsOutName[tableValue].size();
			for(unsigned int i=0;i<size;i++)
				tableVisual.addCol(m_listColumnsOutName[tableValue][i]);
		}
	}
	unsigned int missCol=3;
	if(missCol <=tableVisual.getNumberOfColumns()) 
		missCol=0;
	else
		missCol=missCol-tableVisual.getNumberOfColumns();
	for(int i=0; i<missCol; i++)
		tableVisual.addCol("NONE");

	tableVisual.writeTable(m_fVisualArray);
	return true;
}


