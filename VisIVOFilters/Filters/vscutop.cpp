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
#include "vscutop.h"
#include "VisIVOFiltersConfigure.h"

const unsigned int VSCutOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSCutOp::MIN_NUMBER_OF_ROW = 100;


//---------------------------------------------------------------------
VSCutOp::VSCutOp()
//---------------------------------------------------------------------
{
 m_fArray=NULL;
 m_fArrayWrite=NULL;
 m_nOfRow=0;
 m_nOfCol=0;	
}
//---------------------------------------------------------------------
VSCutOp::~VSCutOp()
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
bool VSCutOp::allocatefArray()
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
void VSCutOp::printHelp()
//---------------------------------------------------------------------
{
	std::cout<<"It fixes to a given value fields where a condition is satisfied."<<std::endl<<std::endl;
	std::cout<<"Usage: VisIVOFilters --op cut  [--field columns_list] --limits filename_limits [--threshold value] [--operator AND/OR] [--out filename_out.bin] [--history] [--historyfile filename.xml] [--help] [--file] inputFile.bin"<<std::endl; 

	std::cout<<"Example: VisIVOFilters --op cut --field A B C --limits limitsfile.txt --operator AND --out filename_out.bin --threshold 1.0 --file inputFile.bin"<<std::endl<<std::endl;
	std::cout<<"The command produces a new table (filename_out.bin and filename_out.bin.head) that contains all the data points and columns of the inputFile table. In any row  of the input table where limits are satisfied, the fields A B and C wil be changed with the threshold value 1.0. Other fields will be not changed."<<std::endl<<std::endl;

	std::cout<<"Note:"<<std::endl;
	std::cout<<"--limits A file that has three columns: a valid column name and an interval that indicate the  limits."<<std::endl;
	std::cout<<"--field valid columns name list to be reported in the new table that will be cut with the threshold. Default all columns will be cutted"<<std::endl;
	std::cout<<"--threshold. Value to be used to cut data. Default value is 0.0"<<std::endl;

	std::cout<<"--operator Limits on all field listed in --limits option file are combined by default with logic AND operator. If this option is given with OR value the field limits are combined with logic OR operator "<<std::endl;
	std::cout<<"--out Output table filename. Default name is given."<<std::endl;
    std::cout<<"--history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
    std::cout<<"--historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;
	std::cout<<"--file  Input table filename."<<std::endl;
	std::cout<<"--help produce this output "<<std::endl;


	return;
}

//---------------------------------------------------------------------
bool VSCutOp::execute()
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

        bool andOp = true;
	float thValue=0.0;
	std::string filename;
	unsigned long long int totRows;
	VSTable tableCut;
	

	bool allColumns=false;

	if(!isParameterPresent("field") || getParameterAsString("field").empty() || getParameterAsString("field")=="unknown")
		allColumns=true;

	if(isParameterPresent("threshold")) thValue=getParameterAsFloat("threshold");

	if(allColumns)
		for(unsigned int i=0;i<m_tables[0] -> getNumberOfColumns();i++)
			m_colNumberSet.push_back(i);
	else
	{
		std::stringstream ssListparameters;
		ssListparameters.str(getParameterAsString("field"));
		std::string paramField;
		while (!ssListparameters.eof())
		{
			ssListparameters>>paramField;
			if(m_tables[0] -> getColId(paramField)>=0)
				m_colNumberSet.push_back(m_tables[0] -> getColId(paramField));
		}
		if(m_colNumberSet.size()==0)
		{
			std::cerr<<"Cut. Invalid fields. Operation aborted"<<filename<<std::endl;
			return false;
		}

	}

	if(getParameterAsString("operator")=="OR") andOp = false;
	if(getParameterAsString("limits").empty() || getParameterAsString("limits")=="unknown" )
	{
		std::cerr<<"vsCutop: No file with field limits is given"<<std::endl;
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
	if(m_colVector.size()<=0) 
	{
		std::cerr<<"Cut. Invalid limits are given."<<std::endl;
		return false;
	}
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
	bool needTempTable=false;
	if(maxEle<m_tables[0]->getNumberOfRows()) needTempTable=true;
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
  		fileNameOutputSStream<<filenameInputTable.substr(0, len-4)<<"_Cut_"<<buffer<<".bin";  //QUI verificare
	}

	fileNameOutput=fileNameOutputSStream.str();
	if(fileNameOutput.find(".bin") == std::string::npos)
   		fileNameOutput.append(".bin");
	m_realOutFilename.push_back(fileNameOutput);
	
	
//Clean existing tab
	remove(fileNameOutput.c_str());
	tableCut.setLocator(fileNameOutput);
	tableCut.setNumberOfRows(m_tables[0]->getNumberOfRows());
#ifdef VSBIGENDIAN
	std::string endianism="big";
	
#else	
	std::string endianism="little";
#endif

	tableCut.setEndiannes(endianism);

	if(m_tables[0]->getIsVolume())
	{
		tableCut.setIsVolume(true);
		const unsigned int *temp;
		temp=new unsigned int[3];
		temp=m_tables[0]->getCellNumber();
		tableCut.setCellNumber(temp[0],temp[1],temp[2]);
		const float *tempf;
		tempf=new float[3];
		tempf=m_tables[0]->getCellSize();
		tableCut.setCellSize(tempf[0],tempf[1],tempf[2]);
	}
 	tableCut.setType("float");
	for(unsigned int k=0;k<m_tables[0]->getNumberOfColumns();k++)
		tableCut.addCol(m_tables[0]->getColName(k));	
// 

// 	tableCut.


	for(int i=0;i<nOfCol;i++) colList[i]=i;

	unsigned long long int totEle=nOfEle;
	bool goodEle;


	while(totEle!=0)
	{
		fromRow=startCounter;
		toRow=fromRow+maxEle-1;
		if(toRow>totRows-1)toRow=totRows-1;
	  	m_tables[0]->getColumn(colList, nOfCol, fromRow, toRow, m_fArray);
		for(unsigned int i=0;i<toRow-fromRow+1;i++)
		{
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

			for(unsigned int j=0;j<nOfCol;j++) m_fArrayWrite[j][i]=m_fArray[j][i];

			if(goodEle)
				for(unsigned int j=0;j<m_colNumberSet.size();j++) m_fArrayWrite[m_colNumberSet[j]][i]=thValue;

		}		
		startCounter=toRow+1;
		totEle=totEle-(toRow-fromRow+1);
		if(needTempTable)
		{
			// write temporary table
			std::stringstream convertStream;
			convertStream<<"__CutOp_temp__"<<randat<<numberOfTempTables<<".bin";
			std::string filenameTemp;
			std::string filenameTempHead;
			convertStream>>filenameTemp;
			convertStream>>filenameTempHead;
			remove(filenameTemp.c_str());
			filenameTempHead.append(".head");
			remove(filenameTempHead.c_str());
			numberOfTempTables++;
			VSTable tableTemp;
			tableTemp.setLocator(filenameTemp);
			tableTemp.setEndiannes(endianism);
 			tableTemp.setType("float");
			unsigned long long int numOfTempRows=toRow-fromRow+1;
 			tableTemp.setNumberOfRows(numOfTempRows);
			for(unsigned int k=0;k<m_tables[0]->getNumberOfColumns();k++)
						tableTemp.addCol(m_tables[0]->getColName(k));
			tableTemp.writeTable(m_fArrayWrite);
		}

	}

	if(!needTempTable) // No temp tables are written
		tableCut.writeTable(m_fArrayWrite);
	else	
	{
		unsigned long long int fromWRow=0, toWRow=-1;
		tableCut.writeTable();

		for(unsigned int i=0;i<numberOfTempTables;i++)
		{
			std::stringstream convertStream,convertStreamHead;
			convertStream<<"__CutOp_temp__"<<randat<<i<<".bin";
			convertStreamHead<<"__CutOp_temp__"<<randat<<i<<".bin.head";
			std::string filenameTemp;
			std::string filenameTempHead;
			convertStream>>filenameTemp;
			convertStreamHead>>filenameTempHead;
			VSTable tableTemp(filenameTemp);
			fromRow=0;
			toRow=tableTemp.getNumberOfRows()-1;
			fromWRow=toWRow+1;
			unsigned int tmp=tableTemp.getNumberOfRows();
			toWRow=fromWRow+tableTemp.getNumberOfRows()-1;
			tableTemp.getColumn(colList,nOfCol, fromRow,  toRow, m_fArray);
			tableCut.putColumn(colList,nOfCol, fromWRow,  toWRow, m_fArray);
			remove(filenameTemp.c_str());
			remove(filenameTempHead.c_str());
		}
	}
delete [] colList;
return true;		
}

