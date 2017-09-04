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
/*#include <string>*/
#include <sstream>
/*#include <set>*/
#include <fstream>
#include <exception>
#ifdef WIN32
	#include <time.h>
#endif
#include "vstable.h"
#include "vsselectfieldop.h"
#include "VisIVOFiltersConfigure.h"
#include "visivoutils.h"
#include <sse_utils.h>

const unsigned int VSSelectFieldOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSSelectFieldOp::MIN_NUMBER_OF_ROW = 100;


//---------------------------------------------------------------------
VSSelectFieldOp::VSSelectFieldOp()
//---------------------------------------------------------------------
{
 m_fArray=NULL;
 m_fArrayWrite=NULL;
 m_nOfRow=0;
 m_nOfCol=0;	
}
//---------------------------------------------------------------------
VSSelectFieldOp::~VSSelectFieldOp()
//---------------------------------------------------------------------
{
	if(m_fArray!=NULL)	
 	  for(unsigned int i=0;i<m_nOfCol;i++)
		{
			if(m_fArray[i]!=NULL) delete [] m_fArray[i];
			if(m_fArrayWrite[i]!=NULL) delete [] m_fArrayWrite[i];
		}
	if(m_fArray!=NULL) delete [] m_fArray;
	if(m_fArrayWrite!=NULL) delete [] m_fArrayWrite;
}
//---------------------------------------------------------------------
bool VSSelectFieldOp::allocatefArray()
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
	if(m_fArray==NULL)
	{
		return false;
	}
try
{
	m_fArrayWrite=new  float*[m_nOfCol];
}
catch(std::bad_alloc &e)
{
	m_fArrayWrite=NULL;
}

	if(m_fArrayWrite==NULL)
	{
		delete [] m_fArray;
		return false;
	}

	bool goodAllocation=false;
	while(!goodAllocation)
	{
		goodAllocation=true;
		for(unsigned int i=0;i<m_nOfCol;i++)
		{
try
{
			m_fArray[i] = new float[m_nOfRow];
}
catch(std::bad_alloc &e)
{
	m_fArray[i]=NULL;
}
 			if(m_fArray[i]==NULL) 
 			{	
				goodAllocation=false;
				for(unsigned int j=0;j<i;j++) 
				{	
					delete [] m_fArray[j];
					delete [] m_fArrayWrite[j];
				}
				if(m_nOfRow==MIN_NUMBER_OF_ROW)
				{ 
					delete [] m_fArray;
					delete [] m_fArrayWrite;
					m_fArray=NULL;
					m_fArrayWrite=NULL;
					return false;
				}
				if(m_nOfRow<=MAX_NUMBER_TO_REDUCE_ROW)
					m_nOfRow=MIN_NUMBER_OF_ROW;
				else
					m_nOfRow=m_nOfRow-MAX_NUMBER_TO_REDUCE_ROW;
				break;
 			}
try
{
			m_fArrayWrite[i]=new  float[m_nOfRow];
}
catch(std::bad_alloc &e)
{
	m_fArrayWrite[i]=NULL;
}

 			if(m_fArrayWrite[i]==NULL) 
 			{	
				goodAllocation=false;
				for(unsigned int j=0;j<i;j++)
				{ 	
					delete [] m_fArrayWrite[j];
					delete [] m_fArray[j];
				}
				delete [] m_fArray[i];
				if(m_nOfRow==MIN_NUMBER_OF_ROW)
				{ 
					delete [] m_fArray;
					delete [] m_fArrayWrite;
					m_fArray=NULL;
					m_fArrayWrite=NULL;
					return false;
				}
				if(m_nOfRow<=MAX_NUMBER_TO_REDUCE_ROW) 
					m_nOfRow=MIN_NUMBER_OF_ROW;
				else
					m_nOfRow=m_nOfRow-MAX_NUMBER_TO_REDUCE_ROW;
				break;
			}
		}
	}
	return true;
}
//---------------------------------------------------------------------
void VSSelectFieldOp::printHelp()
//---------------------------------------------------------------------
{
	std::cout<<"Create e new table using limits on one or more fields of a data table"<<std::endl<<std::endl;
	std::cout<<"Usage: VisIVOFilters --op selfield  --limits filename_limits [--operator AND/OR] [--outlist list_filename] [--format uns/int/ascii] [--out filename_out.bin] [--help] [--file] inputFile.bin"<<std::endl;

	std::cout<<"Example: VisIVOFilters --op selfield --limits limitsfile.txt --operator AND --out filename_out.bin --file inputFile.bin"<<std::endl;

	std::cout<<"Note:"<<std::endl;
	std::cout<<"-limits A file that has three columns: a valid column name and an interval that indicate the extraction limits."<<std::endl;
	std::cout<<"--operator Limits on all field listed in --limits option file are combined by default with logic AND operator. If this option is given with OR value the field limits are combined with logic OR operator "<<std::endl;
	std::cout<<"--out Output table filename. Default name is given."<<std::endl;
	std::cout<<"--outlist Output filename list of elements of selected points. Default name is given"<<std::endl;  
	std::cout<<"--format Output filename list data format: unsigned long long int, int, ascii. Default value unsigned long long int"<<std::endl;  
	std::cout<<"--file  Input table filename."<<std::endl;
	std::cout<<"--help produce this output "<<std::endl;


	return;
}

//---------------------------------------------------------------------
bool VSSelectFieldOp::execute()
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
std::stringstream sstmp1;
std::string  randat;
time_t rawtime;
struct tm * timeinfo;
char buffer [80];
time ( &rawtime );
timeinfo = localtime ( &rawtime );
strftime (buffer,80,"%Y%m%d%H%M%S",timeinfo);
sstmp1<<"_"<<rand()<<buffer<<"_";  
randat=sstmp1.str();

	std::string fileListName;
	bool list=false;
	std::ofstream fileList;
        bool andOp = true;
	std::string filename,tmpDir;
	unsigned long long int totRows;
	VSTable tableSelectField;
	if(m_tables[0]->getIsVolume())
	{
	 	std::cerr<<"vsselectfieldop: operation not allowed: a volume is selected"<<std::endl; 
		return false;
	}
	
	if(getParameterAsString("operator")=="OR") andOp = false;
	if(getParameterAsString("limits").empty() || getParameterAsString("limits")=="unknown" )
	{
		std::cerr<<"vsselectfieldop: No file with field limits is given"<<std::endl;
		return false;
	} else
	{
		filename=getParameterAsString("limits");
	}

	std::ifstream fileInput(filename.c_str());
	if(!fileInput)
	{
		std::cerr<<"Cannot open table list file"<<filename<<std::endl;
		return false;
	}
	limitField inputVar;
	while(!fileInput.eof()) 
	{
	        std::string colName,limInf,limSup;
		fileInput >> colName;
		fileInput >> limInf;
		fileInput >> limSup;

		if(m_tables[0]->getColId(colName)>=0 && limInf !="" && limSup!="")
		{
			inputVar.colId=m_tables[0]->getColId(colName);
			if(limInf=="unlimited") 
				inputVar.downUnlimited=true;
			else
			{
				inputVar.downUnlimited=false;
				inputVar.downLimit=atof(limInf.c_str());
			}
			if(limSup=="unlimited") 
				inputVar.upUnlimited=true;
			else
			{
				inputVar.upUnlimited=false;
				inputVar.upLimit=atof(limSup.c_str());
			}
			m_colVector.push_back(inputVar);	
		}
		
	}
	fileInput.close();	
	if(m_colVector.size()<=0) return false;
	int maxInt=getMaxNumberInt();
	unsigned int maxEle;
	unsigned int numberOfTempTables=0;
	unsigned int writeCounter=0;
	totRows=m_tables[0]->getNumberOfRows();
	unsigned long long int nOfEle=totRows;

	if(totRows>maxInt)
		maxEle=maxInt; 
	else
		maxEle=totRows;
	
	unsigned long long int fromRow, toRow, startCounter=0;

	int *IeleList=NULL;
	unsigned long long int *UeleList=NULL;
	if(isParameterPresent("outlist"))
	{
	  list=true;
//	  std::clog<<getParameterAsString("format")<<std::endl;
	  if(getParameterAsString("format")=="int") 
	    IeleList= new  int[1000000];
	  else 
	    UeleList= new  unsigned long long int[1000000];
	  
	  if(getParameterAsString("outlist")=="" | getParameterAsString("outlist")=="unknown")
	    fileListName="list_"+randat;
	  else
	    fileListName=getParameterAsString("outlist");
	  if(getParameterAsString("format")!="ascii")
	    fileList.open(fileListName.c_str(),std::ios::binary);
	  else
	    fileList.open(fileListName.c_str());
	    
	}
	unsigned int *colList=NULL;
	unsigned int nOfCol=m_tables[0]->getNumberOfColumns();
try
{
	colList=new unsigned int[nOfCol];
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

	m_nOfCol=(unsigned int) nOfCol;
	m_nOfRow=maxEle;
	bool allocationfArray=allocatefArray();
	maxEle=m_nOfRow;
	if(m_fArray==NULL ||m_fArrayWrite==NULL ||  !allocationfArray )	
	{
		std::cerr<<"Failed fArray allocation. Select Field terminated"<<std::endl;
		delete [] colList;
		return false;
	}

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
  		fileNameOutputSStream<<filenameInputTable.substr(0, len-4)<<"_selectfield_"<<buffer<<".bin";  //QUI verificare
	}

	fileNameOutput=fileNameOutputSStream.str();
	if(fileNameOutput.find(".bin") == std::string::npos)
   		fileNameOutput.append(".bin");
	
	m_realOutFilename.push_back(fileNameOutput);
	
	tmpDir=getDir(fileNameOutput);
//Clean existing tab
	remove(fileNameOutput.c_str());
	tableSelectField.setLocator(fileNameOutput);
#ifdef VSBIGENDIAN
	std::string endianism="big";
	
#else	
	std::string endianism="little";
#endif

	tableSelectField.setEndiannes(endianism);

 	tableSelectField.setType("float");
	for(unsigned int k=0;k<m_tables[0]->getNumberOfColumns();k++)
		tableSelectField.addCol(m_tables[0]->getColName(k));	
// 

// 	tableSelectField.


	for(int i=0;i<nOfCol;i++) colList[i]=i;

	unsigned long long int totEle=nOfEle;
	bool goodEle;
	unsigned long long int generalCounter=-1;
	int eleCounter=0;

	while(totEle!=0)
	{
		fromRow=startCounter;
		toRow=fromRow+maxEle-1;
		if(toRow>totRows-1)toRow=totRows-1;
	  	m_tables[0]->getColumn(colList, nOfCol, fromRow, toRow, m_fArray);
		for(unsigned int i=0;i<toRow-fromRow+1;i++)
		{
			generalCounter++;
			if(andOp)
			{
			  goodEle=true;
			  for(unsigned int j=0;j<m_colVector.size();j++)
			  {
				unsigned int colId=m_colVector[j].colId;
				if(!m_colVector[j].downUnlimited)
				{
					if(m_fArray[colId][i]<m_colVector[j].downLimit) 
					{	
						goodEle=false;
						break;
					}
				}
				if(!m_colVector[j].upUnlimited)
				{
					if(m_fArray[colId][i]>m_colVector[j].upLimit)
					{	
						goodEle=false;
						break;
					}
				}
			  }
			} else //OR limits
			{
			  goodEle=false;
			  for(unsigned int j=0;j<m_colVector.size();j++)
			  {
				unsigned int colId=m_colVector[j].colId;
				if(m_colVector[j].downUnlimited && m_colVector[j].upUnlimited)
				{
					goodEle=true;
					break; 
				}
				if(m_colVector[j].downUnlimited && m_fArray[colId][i]>m_colVector[j].upLimit)
					continue;

				if(m_colVector[j].upUnlimited && m_fArray[colId][i]<m_colVector[j].downLimit)
					continue;
				if(m_fArray[colId][i]>=m_colVector[j].downLimit && m_fArray[colId][i]<=m_colVector[j].upLimit)
				{	
						goodEle=true;
						break;
				}
			  }
			}


			if(goodEle)
			{
				if(list)
				{
				  if(getParameterAsString("format")=="int") 
				    IeleList[eleCounter]=generalCounter;
				  else 
				    UeleList[eleCounter]=generalCounter;
				  eleCounter++;
				  if(eleCounter==1000000)
				  {
				    if(getParameterAsString("format")=="ascii")
				      for(int k=0;k<1000000;k++) fileList<<UeleList[k]<<std::endl;
				      else
				      {
					if(getParameterAsString("format")=="int") 
					  fileList.write((char *) IeleList,sizeof(int)*eleCounter);
					else 
					  fileList.write((char *) UeleList,sizeof(unsigned long long int)*eleCounter);
				      }
				    eleCounter=0;
				  }
				  
				}
				// write temporary table
				if(writeCounter==maxEle ||writeCounter<0 )
				{
				        std::stringstream convertStream;
					convertStream<<tmpDir<<"SelFieldOp_temp__"<<randat<<numberOfTempTables<<".bin";
					std::string filenameTemp;
					convertStream>>filenameTemp;
					numberOfTempTables++;
					VSTable tableTemp;

					tableTemp.setLocator(filenameTemp);
					tableTemp.setEndiannes(endianism);
 					tableTemp.setType("float");
					unsigned long long int wrNOfRow=writeCounter;
 					tableTemp.setNumberOfRows(wrNOfRow);

					for(unsigned int k=0;k<m_tables[0]->getNumberOfColumns();k++)
						tableTemp.addCol(m_tables[0]->getColName(k));

					tableTemp.writeTable(m_fArrayWrite);
					writeCounter=0;
				}
				for(unsigned int j=0;j<nOfCol;j++) 
					m_fArrayWrite[j][writeCounter]=m_fArray[j][i];
				writeCounter++;
			}
		}		
		startCounter=toRow+1;
		totEle=totEle-(toRow-fromRow+1);
	}

	if(numberOfTempTables==0) // No temp tables are written
	{
		unsigned long long int wrNOfRow=writeCounter;
 		tableSelectField.setNumberOfRows(wrNOfRow);
		tableSelectField.writeTable(m_fArrayWrite);
	} else	
	{
		unsigned long long int fromWRow, toWRow=-1;
		unsigned long long int wrNOfRow=(unsigned long long int) (maxEle)*(numberOfTempTables)+writeCounter;

 		tableSelectField.setNumberOfRows(wrNOfRow);
		tableSelectField.writeTable();

		for(unsigned int i=0;i<numberOfTempTables;i++)
		{
		  
		  
			std::stringstream convertStream,convertStreamHead;
			convertStream<<tmpDir<<"SelFieldOp_temp__"<<randat<<i<<".bin";
			convertStreamHead<<tmpDir<<"SelFieldOp_temp__"<<randat<<i<<".bin.head";
			std::string filenameTemp;
			std::string filenameTempHead;
			convertStream>>filenameTemp;
			filenameTempHead=filenameTemp+".head";
			VSTable tableTemp(filenameTemp);
			fromWRow=toWRow+1;
			toWRow=fromWRow+maxEle-1;
			fromRow=0;
			toRow=maxEle-1;
			tableTemp.getColumn(colList,nOfCol, fromRow,  toRow, m_fArray);
			tableSelectField.putColumn(colList,nOfCol, fromWRow,  toWRow, m_fArray);
			remove(filenameTemp.c_str());
			remove(filenameTempHead.c_str());
		}
		if(writeCounter>0) //exist m_fArrayWrite partially filled
		{
			fromWRow=toWRow+1;
			toWRow=fromWRow+writeCounter-1;
			tableSelectField.putColumn(colList,nOfCol, fromWRow,  toWRow, m_fArrayWrite);
			
		}
	}
	if(list)
	{
	    if(getParameterAsString("format")=="ascii")
	      for(int k=0;k<eleCounter;k++) fileList<<UeleList[k]<<std::endl;
	    else
		{
		    if(getParameterAsString("format")=="int") 
		      fileList.write((char *) IeleList,sizeof(int)*eleCounter);
		    else 
		      fileList.write((char *) UeleList,sizeof(unsigned long long int)*eleCounter);
		}
	}
	
	fileList.close();
	
if(list)
{
  if(getParameterAsString("format")=="int") delete [] IeleList;
  else delete [] UeleList;
}
delete [] colList;
return true;		
}

