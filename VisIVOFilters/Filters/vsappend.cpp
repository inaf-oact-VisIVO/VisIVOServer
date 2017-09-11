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
//#include <string>
#include <sstream>
//#include <set>
#include <iostream>
#include <fstream>
#include <map>
#include <ctime>
#include "vstable.h"
//#include "vstableop.h"
#include "vsappend.h"
#include "VisIVOFiltersConfigure.h"

//---------------------------------------------------------------------
VSAppendOp::VSAppendOp()
//---------------------------------------------------------------------
{
 
}
//---------------------------------------------------------------------
VSAppendOp::~VSAppendOp()
//---------------------------------------------------------------------
{
 
}
//---------------------------------------------------------------------
void VSAppendOp::printHelp()
//---------------------------------------------------------------------
{
	std::cout<<"Create e new table appending  data from a list of existing tables"<<std::endl;
	std::cout<<"Append Filter can append  up to 100 tables with the same number of Columns"<<std::endl<<std::endl;

	std::cout<<"Usage: VisIVOFilters --op append [--out filename_out.bin] [--history] [--historyfile filename.xml] [--help] [--file] table_list"<<std::endl;

    std::cout<<"--history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
    std::cout<<"--historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;

	std::cout<<"Example: VisIVOFilters --op append  --out out_table.bin --file table_list.txt"<<std::endl<<std::endl;
	std::cout<<"tab_list.txt is a file that contain a list of a valid table name. The .bin extension is automatically added if the listed filename does not contain it"<<std::endl<<"The filter produce a new table (out_table.bin and out_table.bin.head). The column name are copied from the first table. An error is given if tables contain different number of columns."<<std::endl;
	std::cout<<"Note:"<<std::endl;

	std::cout<<"--out Name of the new table. Default name is given."<<std::endl;
	std::cout<<"--file an ascii file with a valid list of tables"<<std::endl;

	std::cout<<"--help produce this output "<<std::endl;


	return;
}

//---------------------------------------------------------------------
bool VSAppendOp::execute()
//---------------------------------------------------------------------
// totRows is the total rows of the out tables
// maxRows is the bigger number of rows among input tables
// 
// maxEle maximum number allowed to upload/download data from a generic table: It is the lower value between maxRows, maxInt and fArray allocation fArray[0][maxEle]
// 
// nOfEle total number of element to upload for the specific input table.
// 
// totEle is a counter: the number of elements that I have still to  upload for the specific table. Starts from nOfEle, its value decrese of maxEle each cycle
// 

{


	std::string filename;
	unsigned long long int totRows=0,nOfColsTab=-1;
	unsigned long long int maxRows=-1;
	VSTable tableAppend;
	if(isParameterPresent("filelist"))
	{  
	  std::cout<<getParameterAsString("filelist");
	  if(getParameterAsString("filelist").empty() || getParameterAsString("filelist")=="unknown" )
	  {
		std::cerr<<"vsappendop: No file with table list is given"<<std::endl;
		return false;
	  } else
	  {
		filename=getParameterAsString("filelist");
	  }
	} else
	{
	  if(getParameterAsString("file").empty() || getParameterAsString("file")=="unknown" )
	  {
		std::cerr<<"vsappendop: No file with table list is given"<<std::endl;
		return false;
	  } else
	  {
		filename=getParameterAsString("file");
	  }
	}
	std::ifstream fileInput(filename.c_str());
	if(!fileInput)
	{
		std::cerr<<"Cannot open table list file"<<filename<<std::endl;
		return false;
	}

// 	while(!fileInput.eof())
// 	{ 
// 		std::string tableName;
// 		fileInput>>tableName;
// 	}

	while(!fileInput.eof())
	{
		std::string tableName;
		fileInput >> tableName;
		if(tableName!="")
		{
			if(tableName.find(".bin") == std::string::npos)
	    			tableName.append(".bin");
			VSTable testTable(tableName);
			if(testTable.getIsVolume())
				std::cerr<<tableName<<" warning has a geometry"<< std::endl;
			if(testTable.tableExist())
			{
				m_listOfTables.push_back(tableName);
				
				if(nOfColsTab==-1)
				{
					nOfColsTab=testTable.getNumberOfColumns();
					for(unsigned int i=0;i<nOfColsTab;i++)
						tableAppend.addCol(testTable.getColName(i));
				}
				if(nOfColsTab>0 && nOfColsTab != testTable.getNumberOfColumns())
				{ 
					std::cerr<<"vsappendop: different number of columns in "<<tableName<<std::endl;
					return false;				
				}
				if(maxRows<testTable.getNumberOfRows())maxRows=testTable.getNumberOfRows();
				totRows=totRows+testTable.getNumberOfRows();
			}
		}
	}
	fileInput.close();
// Set Ouptut Table Filename

	std::stringstream fileNameOutputSStream;
	fileNameOutputSStream<<getParameterAsString("out");
	std::string fileNameOutput;

	if(fileNameOutputSStream.str()=="")
	{
  		std::string filenameInputTable=m_listOfTables[0];
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

	m_realOutFilename.push_back(fileNameOutput);
	
//Clean existing tab
	remove(fileNameOutput.c_str());

	tableAppend.setLocator(fileNameOutput);
#ifdef VSBIGENDIAN
	std::string endianism="big";
	
#else	
	std::string endianism="little";
#endif

	tableAppend.setEndiannes(endianism);
 	tableAppend.setType("float");
	tableAppend.setNumberOfRows(totRows);
	tableAppend.writeTable();


	int maxInt=getMaxNumberInt();
	unsigned long long int nOfEle;
	int maxEle;
	int mergedTableCol=0;

	if(totRows>maxInt)
		maxEle=maxInt; 
	else
		maxEle=totRows;
	
	unsigned long long int fromRow, toRow, startCounter;
	unsigned long long int fromDestRow, toDestRow;
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

	for(unsigned int j=0;j<tableAppend.getNumberOfColumns();j++)
	{	
		fromDestRow=0;
		colList[0]= j;
		for(unsigned int i=0; i<m_listOfTables.size();i++)
		{
			std::string tablename=m_listOfTables[i]; 
			if(tablename !="")
			{
				VSTable testTable(tablename); 
				startCounter=0;
				unsigned long long int totLocalRows=testTable.getNumberOfRows();
				unsigned long long int totEle=totLocalRows;
				while(totEle!=0)
				{
					fromRow=startCounter;
					toRow=totLocalRows-1;
					if(toRow-fromRow>maxEle)toRow=fromRow+maxEle-1;
	  				testTable.getColumn(colList, nOfCol, fromRow, toRow, fArray);
					toDestRow=fromDestRow+(toRow-fromRow);
					tableAppend.putColumn(colList, nOfCol, fromDestRow, toDestRow, fArray);
					fromDestRow=toDestRow+1;
					startCounter=toRow+1;
					totEle=totEle-(toRow-fromRow+1);

				}
			}
		}
	}
		
if(fArray!=NULL ) 
   if(fArray[0] != NULL) delete [] fArray[0];

if(fArray!=NULL) delete [] fArray;
	
return true;
}
