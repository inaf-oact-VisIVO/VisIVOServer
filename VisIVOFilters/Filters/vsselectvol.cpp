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
#include <fstream>
#include <iostream>
#include <sstream>

/*#include <string>*/
#include <sstream>
/*#include <set>*/
#include <fstream>
#include <exception>
#ifdef WIN32
	#include <time.h>
#endif
#include "vstable.h"
#include "vsselectvol.h"
#include "VisIVOFiltersConfigure.h"
#include "visivoutils.h"

const unsigned int VSSelectVolOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSSelectVolOp::MIN_NUMBER_OF_ROW = 100;
const unsigned int VSSelectVolOp::MAXOUT = 100000;


//---------------------------------------------------------------------
VSSelectVolOp::VSSelectVolOp()
//---------------------------------------------------------------------
{
 m_fArray=NULL;
 m_fArrayWrite=NULL; 
 m_nOfRow=0;
 m_nOfCol=0;	
 m_numCells=0;
}
//---------------------------------------------------------------------
VSSelectVolOp::~VSSelectVolOp()
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
bool VSSelectVolOp::allocatefArray()
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
				}
				if(m_nOfRow==MIN_NUMBER_OF_ROW)
				{ 
					delete [] m_fArray;
					m_fArray=NULL;
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
void VSSelectVolOp::printHelp()
//---------------------------------------------------------------------
{
	std::cout<<"This filter creates a new volume setting to zero all the volume cells that does not satisfy given limits."<<std::endl<<std::endl;
	std::cout<<"Optionally it writes an output ascii file containing the volume cells and their values that satisfy given limits."<<std::endl<<std::endl;
	std::cout<<"Usage: VisIVOFilters --op selvol [--field column_name] --limits filename_limits [--operator AND/OR]  [--numcells value] [--out filename_out.bin] [--outlist list_filename] [--help] [--file] inputFile.bin"<<std::endl; 

	std::cout<<"Example: VisIVOFilters --op selvol --field density --limits limitsfile.txt --operator AND --file inputFile.bin"<<std::endl<<std::endl;
	std::cout<<"The command produces a new volume with all cells where limits are satisfied, all the other cells are set to zero."<<std::endl<<std::endl;

	std::cout<<"Note:"<<std::endl;
	std::cout<<"--field valid column name  to be reported in the output table. Default value is all column"<<std::endl;
	std::cout<<"--limits A file that has three columns: a valid column name and an interval that indicate the  limits."<<std::endl;
	std::cout<<"--operator Limits on all field listed in --limits option file are combined by default with logic AND operator. If this option is given with OR value the field limits are combined with logic OR operator "<<std::endl;
	std::cout<<"--numcells Set the maximum number of cells that will be reported in the output ascii file."<<std::endl;
	std::cout<<"--out Output volume filename. Default name is given."<<std::endl;
	std::cout<<"--outlist Output list filename containing the volume cells satisfying the requested condition. Default name is given."<<std::endl;
	std::cout<<"--file  Input table filename."<<std::endl;
	std::cout<<"--help produce this output "<<std::endl;


	return;
}

//---------------------------------------------------------------------
bool VSSelectVolOp::execute()
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
	if(!m_tables[0]->getIsVolume())
	{
		std::cerr<<"VolLimit can be applied to volume only. Operation aborted"<<std::endl;
		return false;
		
	}
	bool andOp = true;
	std::string filename, tmpDir;
	unsigned long long int totRows;
	bool numCellSet=false;
	int numCell=0;
	int wrElements=0;
	VSTable tableSelectVol;
 
	if(isParameterPresent("numcells"))
	{
		numCell=getParameterAsInt("numcells");
		if(numCell<=0)
		{
			std::cerr<<"selvol. Invalid numcells. Operation aborted"<<std::endl;
			return false;

		}
		numCellSet=true;	
	}

	

	if(!isParameterPresent("field") || getParameterAsString("field").empty() || getParameterAsString("field")=="unknown")
	{
		for(unsigned int i=0;i<m_tables[0]->getNumberOfColumns();i++)
				m_colNumberSet.push_back(i);
	}else{
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
			std::cerr<<"selvol. Invalid field. Operation aborted"<<filename<<std::endl;
			return false;
		}

	}
	
	if(getParameterAsString("operator")=="OR") andOp = false;
	if(getParameterAsString("limits").empty() || getParameterAsString("limits")=="unknown" )
	{
		std::cerr<<"selvol: No file with field limits is given"<<std::endl;
		return false;
	} else
	{
		filename=getParameterAsString("limits");
	}

	std::ifstream fileInput(filename.c_str());
	if(!fileInput.is_open())
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
		std::cerr<<"VolLimit. Invalid limits are given."<<std::endl;
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
	if(m_fArray==NULL ||  !allocationfArray )	
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
                fileNameOutputSStream<<filenameInputTable.substr(0, len-4)<<"_selectvolume_"<<buffer<<".bin";  //QUI verificare
        }

        fileNameOutput=fileNameOutputSStream.str();
        if(fileNameOutput.find(".bin") == std::string::npos)
                fileNameOutput.append(".bin");
        m_realOutFilename.push_back(fileNameOutput);


        tmpDir=getDir(fileNameOutput);
//Clean existing tab
        remove(fileNameOutput.c_str());
        tableSelectVol.setLocator(fileNameOutput);
#ifdef VSBIGENDIAN
        std::string endianism="big";

#else   
        std::string endianism="little";
#endif

        tableSelectVol.setEndiannes(endianism);

        tableSelectVol.setType("float");
        for(unsigned int k=0;k<m_tables[0]->getNumberOfColumns();k++)
                tableSelectVol.addCol(m_tables[0]->getColName(k));    
	



	std::stringstream fileListOutputSStream;
	fileListOutputSStream<<getParameterAsString("outlist");
	std::string fileListOutput;

	if(fileListOutputSStream.str()==""||fileListOutputSStream.str()=="unknown")
	{
		fileListOutputSStream.str().erase(); //QUI verificare
  		std::string filenameInputTable=m_tables[0]->getLocator();
  		int len=filenameInputTable.length();
		time_t rawtime;
		struct tm * timeinfo;
		char buffer [80];
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		strftime (buffer,80,"%Y%m%d%H%M",timeinfo);
  		fileListOutputSStream<<filenameInputTable.substr(0, len-4)<<"_VolLimit_"<<buffer<<".txt";  //QUI verificare
	}

	fileListOutput=fileListOutputSStream.str();
	std::ofstream fileList(fileListOutput.c_str(), std::ios::out);
	fileList<<"X  Y  Z  ";

	for(int i=0;i<m_colNumberSet.size();i++)
		fileList<<m_tables[0]->getColName(m_colNumberSet[i])<<" ";
	fileList<<std::endl; 

//Clean existing tab
	for(int i=0;i<nOfCol;i++)
		colList[i]=i;

	unsigned long long int totEle=nOfEle;
	bool goodEle;
	int counterEle=0;
	int globalCounter=-1;
	const unsigned int *tableCells;
	tableCells=new unsigned int[3];
	tableCells=m_tables[0]->getCellNumber();
	
	while(totEle!=0)
	{
		fromRow=startCounter;
		toRow=fromRow+maxEle-1;
		if(toRow>totRows-1)toRow=totRows-1;
	  	m_tables[0]->getColumn(colList, nOfCol, fromRow, toRow, m_fArray);
		for(unsigned int i=0;i<toRow-fromRow+1;i++)
		{
			globalCounter++;
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
			}  // if (andop) else


			if(goodEle)
			{
				int z=globalCounter/(tableCells[0]*tableCells[1]);
				int y=globalCounter-z*(tableCells[0]*tableCells[1]);
				y=y/tableCells[0];
				int x=globalCounter-z*(tableCells[0]*tableCells[1])-y*tableCells[0];
				fileList<<x+1<<" ";
				fileList<<y+1<<" ";
				fileList<<z+1<<" ";
				for(int j=0;j<m_colNumberSet.size();j++){
					fileList<<m_fArray[m_colNumberSet[j]][i]<<" ";
					m_fArrayWrite[m_colNumberSet[j]][i]=m_fArray[m_colNumberSet[j]][i];
				}
				fileList<<std::endl;
				wrElements++;
				if(numCellSet && wrElements==numCell)
				{
					fileList.close();
					delete [] colList;
					return true;
				}
			}
			else{ // not good cell
                        	for(int j=0;j<m_colNumberSet.size();j++)
                                        m_fArrayWrite[m_colNumberSet[j]][i]=0.;
                            }
			
		}
		
		startCounter=toRow+1;
		totEle=totEle-(toRow-fromRow+1);
	}
	unsigned long long int wrNOfRow=m_nOfRow;
        tableSelectVol.setNumberOfRows(wrNOfRow);

	tableSelectVol.setIsVolume(true);
  	tableSelectVol.setCellNumber(m_tables[0] -> getCellNumber()[0],m_tables[0] -> getCellNumber()[1],m_tables[0] -> getCellNumber()[2]);
  	tableSelectVol.setCellSize(m_tables[0] -> getCellSize()[0],m_tables[0] -> getCellSize()[1],m_tables[0] -> getCellSize()[2]);

	tableSelectVol.writeTable(m_fArrayWrite);

fileList.close();
delete [] colList;
return true;		
}

