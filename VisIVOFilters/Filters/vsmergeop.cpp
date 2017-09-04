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
#include "vsmergeop.h"
#include "VisIVOFiltersConfigure.h"

//---------------------------------------------------------------------
VSMergeOp::VSMergeOp()
//---------------------------------------------------------------------
{
 
}
//---------------------------------------------------------------------
VSMergeOp::~VSMergeOp()
//---------------------------------------------------------------------
{
 
}
//---------------------------------------------------------------------
void VSMergeOp::printHelp()
//---------------------------------------------------------------------
{
	std::cout<<"Merge  up to 100 tables"<<std::endl<<std::endl;;

	std::cout<<"Usage: VisIVOFilters --op merge  [--size HUGE/SMALLEST] [--pad value]  [--out filename_out.bin] [--help] [--filelist] table_param.txt"<<std::endl;

	std::cout<<"Example: VisIVOFilters --op merge --out out_table_file.bin --filelist table_param.txt"<<std::endl<<std::endl;

	std::cout<<"Note: "<<std::endl;
	std::cout<<"--size SMALLEST default option. Produce a new table having the size the smallest table"<<std::endl;
	std::cout<<"--pad 0 default option. Pad the table rows of smaller table with value if HUGE size is used"<<std::endl;
	std::cout<<"--out Name of the new table. Default name is given."<<std::endl;
	std::cout<<"--filelist list of tables and valid column name to be merged. Wildcard * meabs all columns."<<std::endl;
	std::cout<<"--help produce this output "<<std::endl;

	return;
}

//---------------------------------------------------------------------
bool VSMergeOp::execute()
//---------------------------------------------------------------------
// totRows is the minimum or the maximum number of rows among all tables
// 
// maxEle maximum number allowed to upload/download data from a generictable: It is the lower value between totRows, maxInt and fArray allocation fArray[0][maxEle]
// 
// nOfEle total number of ellement to upload for the specific input table. nOfEle is the lower value between the rows of the table and totRows
// 
// totEle is a counter: the numebr of elemnts that I have still to upload for the specific table. Starts from nOfEle, its value devrese of maxEle each cycle
// 

{
	bool sizeSmall=true;
	std::string filename;
	int value=-1,tableValue;
	std::map<std::string, int>::iterator p;
	unsigned int totColumns=0;
	unsigned long long int totRows=0; 
	VSTable tableMerged;
	float padValue=0.0;

	if(isParameterPresent("filelist"))
	{
	  if(getParameterAsString("filelist").empty() || getParameterAsString("filelist")=="unknown" )
	  {
		std::cerr<<"vsmergeop: No file with table list is given"<<std::endl;
		return false;
	  } else
	  {
		filename=getParameterAsString("filelist");
	  }
	
	} else
	{
	  std::clog<<"filelist:  isParameterPresent NOT given"<<std::endl;
	  if(getParameterAsString("file").empty() || getParameterAsString("file")=="unknown" )
	  {
		std::cerr<<"vsmergeop: No file with table list is given"<<std::endl;
		return false;
	  } else
	  {
		filename=getParameterAsString("file");
	  }
	}  
	if(getParameterAsString("size").empty() || getParameterAsString("size")=="unknown" )
	{
		sizeSmall=true;
	} else
	{
		std::string sizeString=getParameterAsString("size");
		if(sizeString=="HUGE") sizeSmall=false;
	}
	if(getParameterAsString("pad").empty())
		padValue=0.0;
	 else
		padValue=getParameterAsFloat("pad");

	std::ifstream fileInput(filename.c_str());
	if(!fileInput)
	{
		std::cerr<<"Cannot open table list file"<<filename<<std::endl;
		return false;
	}



	while(!fileInput.eof())
	{
		std::string key,stringColumn;
		fileInput >> key;
		fileInput >> stringColumn;
		if(key!="")
		{
                	if(key.find(".bin") == std::string::npos)
   				key.append(".bin");
			VSTable testTable(key);
			if(testTable.tableExist())
			{
				p=m_listOfTables.find(key);
				if(p ==m_listOfTables.end()) //QUI verifica
				{
					value++;
					if(testTable.getColId(stringColumn)>=0)
					{
						m_listOfTables.insert(make_pair(key,value));
						m_listColumns[value].push_back(testTable.getColId(stringColumn));
						if(value==0) totRows=testTable.getNumberOfRows();

						if(value>0 && sizeSmall && totRows>testTable.getNumberOfRows()) totRows=testTable.getNumberOfRows();

						if(value>0 && !sizeSmall && totRows<testTable.getNumberOfRows()) totRows=testTable.getNumberOfRows();

/*						std::stringstream sStringColumn;
						sStringColumn<<stringColumn<<"_tab_"<<value+1;
						tableMerged.addCol(sStringColumn.str()); */
						if(value==101) return false;
						
					}
					if(stringColumn=="*")
					{
						m_listOfTables.insert(make_pair(key,value));
						unsigned int numberOfCols=testTable.getNumberOfColumns();
						for(int i=0;i<numberOfCols;i++)
						{ 
							m_listColumns[value].push_back(i);
/*							std::stringstream sStringColumn;
							sStringColumn<<testTable.getColName(i)<<"_tab_"<<value+1;
							tableMerged.addCol(sStringColumn.str());*/
						}

						if(value==0) totRows=testTable.getNumberOfRows();
						if(value>0 && sizeSmall && totRows>testTable.getNumberOfRows()) totRows=testTable.getNumberOfRows();
						if(value>0 && !sizeSmall && totRows<testTable.getNumberOfRows()) totRows=testTable.getNumberOfRows();
						if(value==101) return false;
					}

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
/*							std::stringstream sStringColumn;
							sStringColumn<<stringColumn<<"_tab_"<<tableValue+1;
							tableMerged.addCol(sStringColumn.str());*/
						}
					}
					if(stringColumn=="*")
					{
						unsigned int numberOfCols=testTable.getNumberOfColumns();
						m_listColumns[tableValue].clear();
						for(unsigned int i=0;i<numberOfCols;i++)
						{
/*							std::stringstream sStringColumn;
							sStringColumn<<testTable.getColName(i)<<"_tab_"<<tableValue+1;*/
							m_listColumns[tableValue].push_back(i);
/*							tableMerged.addCol(sStringColumn.str());*/
						}
					}
				}
			}				
		}
	}
	fileInput.close();
	for(unsigned int i=0;i<=value;i++)
		totColumns=totColumns+m_listColumns[i].size();



	if(totColumns==0) return false;
//Set Output Table Filename
	std::stringstream fileNameOutputSStream;
	fileNameOutputSStream<<getParameterAsString("out"); //QUI brutto
	std::string fileNameOutput;

	if(fileNameOutputSStream.str()==""||fileNameOutputSStream.str()=="unknown")
	{
		p=m_listOfTables.begin();
  		std::string filenameInputTable=p->first;
  		int len=filenameInputTable.length();
		time_t rawtime;
		struct tm * timeinfo;
		char buffer [80];
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		strftime (buffer,80,"%Y%m%d%H%M",timeinfo);
  		fileNameOutputSStream<<filenameInputTable.substr(0, len-4)<<"_merge_"<<buffer<<".bin";  //QUI verificare
	}

	fileNameOutput=fileNameOutputSStream.str();
	if(fileNameOutput.find(".bin") == std::string::npos)
   		fileNameOutput.append(".bin");
//Clean existing tab
	remove(fileNameOutput.c_str());	
	m_realOutFilename.push_back(fileNameOutput);
	
	tableMerged.setLocator(fileNameOutput);
#ifdef VSBIGENDIAN
	std::string endianism="big";
	
#else	
	std::string endianism="little";
#endif

	tableMerged.setEndiannes(endianism);
 	tableMerged.setType("float");



	tableMerged.setNumberOfRows(totRows);
	int maxInt=getMaxNumberInt();
	unsigned long long int nOfEle;
	int maxEle;
	unsigned int mergedTableCol=-1;

	if(totRows>maxInt)
		maxEle=maxInt; 
	else
		maxEle=totRows;
	
	unsigned long long int fromRow, toRow, startCounter;
	int nOfCol=1;
	unsigned int colList[1];
	float **fArray=new  float*[1];
try
{
	fArray[0]=new  float[maxEle];
}
catch(std::bad_alloc &e)
{
	fArray[0]=NULL;
}

	while(fArray[0] == NULL)
	{
		maxEle=maxEle-10000;
		if(maxEle<=0)
		{
			std::cerr<<"Bad **fArray allocation. Select Column terminated"<<std::endl;
			return false;
		}
try
{
		fArray[0]=new  float[maxEle];
}
catch(std::bad_alloc &e)
{
	fArray[0]=NULL;
}

	}


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
				std::stringstream sStringColumn;
				sStringColumn<<testTable.getColName(m_listColumns[tableValue][i])<<"_tab_"<<tableValue+1;			
				tableMerged.addCol(sStringColumn.str());
			}
		}
	}
	tableMerged.writeTable();

	for(p=m_listOfTables.begin();p!=m_listOfTables.end();p++)
	{
		std::string locator=p->first;
		if(locator!="")
		{
			VSTable testTable(locator);
			tableValue=p->second;
			if(totRows>testTable.getNumberOfRows()) 
				nOfEle=testTable.getNumberOfRows();
			else
				nOfEle=totRows;			

			unsigned int size=m_listColumns[tableValue].size();
			for(unsigned int i=0;i<size;i++)
			{	
//				colList[0]=m_listColumns[tableValue][i];
				startCounter=0;
				unsigned long long int totEle=totRows;
				mergedTableCol++;

				while(totEle!=0)
				{
					fromRow=startCounter;
					toRow=fromRow+maxEle-1;
					if(toRow>totRows-1)toRow=totRows-1;
					if(fromRow<nOfEle)
					{
						if(toRow>nOfEle-1)toRow=nOfEle-1;
						colList[0]=m_listColumns[tableValue][i];
	  					testTable.getColumn(colList, nOfCol, fromRow, toRow, fArray);
					}else //pad
					{
						for(unsigned int j=0;j<(toRow-fromRow+1);j++) 
							fArray[0][j]=padValue;
					}
					totEle=totEle-(toRow-fromRow+1);
					colList[0]=mergedTableCol;
					tableMerged.putColumn(colList, nOfCol, fromRow, toRow, fArray);
					startCounter=toRow+1;
				}
				
			}
		}
		
	}
	delete [] fArray[0];
	delete [] fArray;

	return true;
}


