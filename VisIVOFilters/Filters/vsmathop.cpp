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
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cstring>
/*#include <set>*/
#include <fstream>
#ifdef WIN32
	#include <time.h>
#endif
#include "vstable.h"
#include "vsmathop.h"
#include "fparser.h"
#include "VisIVOFiltersConfigure.h"

const unsigned int VSMathOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSMathOp::MIN_NUMBER_OF_ROW = 100;

//---------------------------------------------------------------------
VSMathOp::VSMathOp()
//---------------------------------------------------------------------
{
 m_fArray=NULL;
 m_result=NULL;
 m_nOfRow=0;
 m_nOfCol=0;	
}
//---------------------------------------------------------------------
VSMathOp::~VSMathOp()
//---------------------------------------------------------------------
{	
	if(m_fArray!=NULL)
 		for(unsigned int i=0;i<m_nOfCol;i++)
		{
			if(m_fArray[i] != NULL) delete [] m_fArray[i];
		}
	if(m_result != NULL)
	{ 
		if(m_result[0] != NULL)
			delete [] m_result[0];
	}
	if(m_fArray!=NULL) delete [] m_fArray;
	if(m_result!=NULL) delete [] m_result;
}
//---------------------------------------------------------------------
bool VSMathOp::allocateArray()
//---------------------------------------------------------------------
{
unsigned long long int tempLL=getMaxNumberInt();
if(((unsigned long long int)m_nOfRow*m_nOfCol)>tempLL) 
	m_nOfRow=(int)tempLL/m_nOfCol;

	m_result=new float*[1];
try
{
	m_fArray=new  float*[m_nOfCol];
}
catch(std::bad_alloc &e)
{
	m_fArray=NULL;
}

	if(m_fArray == NULL)
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

			if(m_fArray[i] == NULL) 
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


			if(i==m_nOfCol-1)
			{
try
{
				m_result[0]=new  float[m_nOfRow];
}
catch(std::bad_alloc &e)
{
	m_result[0]=NULL;
}

				if(m_result[0] == NULL) 
				{	
					goodAllocation=false;
					for(unsigned int j=0;j<=i;j++)
						delete [] m_fArray[j];
					if(m_nOfRow==MIN_NUMBER_OF_ROW)
					{ 
						delete [] m_fArray;
						m_fArray=NULL;
						m_result=NULL;
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
		
	}
	return true;
}
//---------------------------------------------------------------------
void VSMathOp::printHelp()
//---------------------------------------------------------------------
{
	std::cout<<"Create a new field in a data table as the result of a mathematical operation between the existing fields"<<std::endl<<std::endl;
	std::cout<<"Usage: VisIVOFilters --op mathop [--expression math_expression.txt] [--compute <<expression>>] [--append] [--outcol col_name] [--out filename_out.bin] [--help] [--history] [--historyfile filename.xml] [--file] filename.bin"<<std::endl;
    std::cout<<"--history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
    std::cout<<"--historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;


	std::cout<<"Example: VisIVOFilters --op mathop --expression my_expression.txt --outcol SqrMod --out pos_math.bin --file pos.bin"<<std::endl;

	std::cout<<"Note: A New field is produced. Default options: a new table with only the new field is produced "<<std::endl;
	std::cout<<"      --expression a file with a valid expression with  Valid Column names. Ignored if the compute option is given."<<std::endl;
	std::cout<<"      --compute a  a valid expression with  Valid Column names. The expression must start with << and finish with >> characters. It has the priority on the  expression option. The expression  must contain the escape character control  for the << and >> symbols and the parentheses.  This option has the priority on the expression option. NOTE: the << , >> and escape characters MUST NOT BE GIVEN if  the parameter file is used."<<std::endl;
	std::cout<<"      --outcol. Column name of the new field"<<std::endl;
	std::cout<<"      --append. No new table will be cretaed. The original table will have the new field "<<std::endl;
	std::cout<<"      --out. Name of the new table. Ignored if -- append is specified."<<std::endl;
	std::cout<<"      --help produce this output "<<std::endl;


	return;
}
//---------------------------------------------------------------------
bool VSMathOp::nameSub()
//---------------------------------------------------------------------
{
	int newNameNumber=-1;
	bool still;
	 int found=0;
	bool newName;
	std::vector<int> positions;
for(unsigned int i=0;i<m_tables[0]->getNumberOfColumns();i++)
{
	still=true;
	newName=true;
	while(still)
	{
		bool foundName=false,beforeChar=false,afterChar=false;
		found=m_mathExpression.find(m_tables[0]->getColName(i));
//		std::clog<<m_tables[0]->getColName(i)<<" "<<m_mathExpression<<std::endl;
  		if (found!=std::string::npos) 
		{
			foundName=true;
			if(found==0)
				beforeChar=true;
			else
			{
				if(m_mathExpression.compare(found-1,1," ") ==0)
					beforeChar=true;
				else if(m_mathExpression.compare(found-1,1,"(") ==0)
					beforeChar=true;
				else if(m_mathExpression.compare(found-1,1,"-") ==0)
					beforeChar=true;
				else if(m_mathExpression.compare(found-1,1,"!") ==0)
					beforeChar=true;
				else if(m_mathExpression.compare(found-1,1,"*") ==0)
					beforeChar=true;
				else if(m_mathExpression.compare(found-1,1,"/") ==0)
					beforeChar=true;
				else if(m_mathExpression.compare(found-1,1,"%") ==0)
					beforeChar=true;
				else if(m_mathExpression.compare(found-1,1,">") ==0)
					beforeChar=true;
				else if(m_mathExpression.compare(found-1,1,"<") ==0)
					beforeChar=true;
				else if(m_mathExpression.compare(found-1,1,"&") ==0)
					beforeChar=true;
				else if(m_mathExpression.compare(found-1,1,"+") ==0)
					beforeChar=true;
				else if(m_mathExpression.compare(found-1,1,"=") ==0)
					beforeChar=true;
				else if(m_mathExpression.compare(found-1,1,"<") ==0)
					beforeChar=true;
				else if(m_mathExpression.compare(found-1,1,">") ==0)
					beforeChar=true;
				else if(m_mathExpression.compare(found-1,1,"&") ==0)
					beforeChar=true;
				else if(m_mathExpression.compare(found-1,1,"|") ==0)
					beforeChar=true;

			}
			//search for end string
			{
			  int lenString=m_tables[0]->getColName(i).size();
			  if((found+lenString)==m_mathExpression.size())
					afterChar=true;
			  else
			  {
				if(m_mathExpression.compare(found+lenString,1," ") ==0)
					afterChar=true;
				else if(m_mathExpression.compare(found+lenString,1,")") ==0)
					afterChar=true;
				else if(m_mathExpression.compare(found+lenString,1,"-") ==0)
					afterChar=true;
				else if(m_mathExpression.compare(found+lenString,1,"!") ==0)
					afterChar=true;
				else if(m_mathExpression.compare(found+lenString,1,"*") ==0)
					afterChar=true;
				else if(m_mathExpression.compare(found+lenString,1,"/") ==0)
					afterChar=true;
				else if(m_mathExpression.compare(found+lenString,1,"%") ==0)
					afterChar=true;
				else if(m_mathExpression.compare(found+lenString,1,">") ==0)
					afterChar=true;
				else if(m_mathExpression.compare(found+lenString,1,"<") ==0)
					afterChar=true;
				else if(m_mathExpression.compare(found+lenString,1,"&") ==0)
					afterChar=true;
				else if(m_mathExpression.compare(found+lenString,1,"+") ==0)
					afterChar=true;
				else if(m_mathExpression.compare(found+lenString,1,"=") ==0)
					afterChar=true;
				else if(m_mathExpression.compare(found+lenString,1,"<") ==0)
					afterChar=true;
				else if(m_mathExpression.compare(found+lenString,1,">") ==0)
					afterChar=true;
				else if(m_mathExpression.compare(found+lenString,1,"&") ==0)
					afterChar=true;
				else if(m_mathExpression.compare(found+lenString,1,"|") ==0)
					afterChar=true;
				else if(m_mathExpression.compare(found+lenString,1,"^") ==0)
					afterChar=true;
			  }
			}
		  if(foundName && beforeChar && afterChar)
		  {
			std::stringstream colNewName;
			if(newName)newNameNumber++;
			colNewName<<"VoOMaTh"<<newNameNumber; //QUI inserirre un controllo sulla NON presenza di tale nome nella stringa originaria: genera una strinnga casuale e verifica
			if(newName) m_colOldName.push_back(m_tables[0]->getColName(i));
			if(newName) m_colName.push_back(colNewName.str());
			newName=false;
			m_mathExpression.erase(found,m_tables[0]->getColName(i).size());
			m_mathExpression.insert(found,colNewName.str());
		  } else
		  {
			m_mathExpression.erase(found,m_tables[0]->getColName(i).size());
			positions.push_back(found);
		  }
		  
 		} else
			still=false;
	}
	for(int j=0;j<positions.size();j++)
		m_mathExpression.insert(positions[j],m_tables[0]->getColName(i));
	positions.clear();
}
if(m_colName.size()==0)
{
	std::cerr<<"vsmathop: Invalid expression is given"<<std::endl;
	return false;	
}   //String name replacement
return true;
}
//---------------------------------------------------------------------
bool VSMathOp::execute()
//---------------------------------------------------------------------
// totRows is the total rows of the out tables
// maxRows is the bigger number of rows among input tables
// 
// maxEle maximum number allowed to upload/download data from a generic table: It is the lower value between maxRows, maxInt and fArray allocation fArray[0][maxEle]
// 
// nOfEle total number of element to upload for the specific input table.
// 
// totEle is a counter: the number of elements that I have still to  upload for the specific table. Starts from nOfEle, its value decrese of maxEle each cycle
{
bool compute=false;
if(isParameterPresent("compute"))
	compute=true;
	
std::stringstream mathExpressionSS;
if(!compute)
{
	if(getParameterAsString("expression").empty() || getParameterAsString("expression")=="unknown" )
	{
		std::cerr<<"vsmathop: No File with mathematical expression is given"<<std::endl;
		return false;
	} else
	{
		std::string filename=getParameterAsString("expression");
		std::ifstream fileInput(filename.c_str());
		if(!fileInput)
		{
			std::cerr<<"Cannot open table list file "<<filename<<std::endl;
			return false;
		}

		while(!fileInput.eof())		//QUI si potrebbe fare un ck se ci sono 2 righe
		{ 
			std::string inputString;
			fileInput>>inputString;
			if(inputString !="")
				mathExpressionSS<<inputString<<" " ;
		}
		fileInput.close();
		m_mathExpression=mathExpressionSS.str(); 
	}
// 	bool append=true; // ** QUI strano ma non funziona!
// 	if(getParameterAsString("append").empty() || getParameterAsString("append")=="unknown" ) append=false;
} else
	m_mathExpression=getParameterAsString("compute");

std::cerr<<"Mathematical Expression: "<<m_mathExpression<<std::endl;
bool append=true; // ** QUI strano!
if(!isParameterPresent("append")) append=false;
//	 size_t found; // ** QUI strano! ma non funziona: lo prende come unsigned int!
	 int found=0;

	if(!nameSub())
		return false;

	m_nOfCol=m_colName.size(); 
	unsigned int *colList=NULL;	
try
{
	colList=new unsigned int[m_nOfCol];
}
catch(std::bad_alloc &e)
{
	colList=NULL;
}

	if(colList == NULL)
	{
		std::cerr<<"Failed colList allocation. Select Field terminated"<<std::endl;
		return false;
	}

	int maxInt=getMaxNumberInt();
	unsigned int maxEle;
	unsigned long long int totRows=m_tables[0]->getNumberOfRows();
	unsigned long long int nOfEle=totRows;

	if(totRows>maxInt)
		maxEle=maxInt; 
	else
		maxEle=totRows;
	
	unsigned long long int fromRow, toRow, startCounter=0;

	m_nOfRow=maxEle;
	bool allocationArray=allocateArray();
	maxEle=m_nOfRow;
	if(m_fArray==NULL ||m_result==NULL ||  !allocationArray )	
	{
		std::cerr<<"Failed Array allocation. Mathematical Operation terminated"<<std::endl;
		delete [] colList;
		return false;
	}
	std::string fileNameOutput;

	if(!append)
	{
		std::stringstream fileNameOutputSStream;
		fileNameOutputSStream<<getParameterAsString("out"); //QUI Brutto! modificare quie a ltri files come sotto per colNameOutput!

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
  			fileNameOutputSStream<<filenameInputTable.substr(0, len-4)<<"_mathop_"<<buffer<<".bin";  //QUI verificare
		}
		fileNameOutput=fileNameOutputSStream.str();
  		if(fileNameOutput.find(".bin") == std::string::npos)
	    		fileNameOutput.append(".bin");
		
	} else
		fileNameOutput=m_tables[0]->getLocator();
	m_realOutFilename.push_back(fileNameOutput);
	VSTable tableMath;
	tableMath.setLocator(fileNameOutput);

#ifdef VSBIGENDIAN
	std::string endianism="big";
	
#else	
	std::string endianism="little";
#endif

	tableMath.setEndiannes(endianism);

 	tableMath.setType("float");
 	tableMath.setNumberOfRows(totRows);
	tableMath.setIsVolume(m_tables[0]->getIsVolume());
	if(m_tables[0]->getIsVolume())
	{
	    tableMath.setCellNumber(m_tables[0]->getCellNumber()[0],
                           m_tables[0]->getCellNumber()[1],
                           m_tables[0]->getCellNumber()[2]);
            tableMath.setCellSize(m_tables[0]->getCellSize()[0],
                           m_tables[0]->getCellSize()[1],
                           m_tables[0]->getCellSize()[2]);
	}



	std::string colNameOutput;
	if(getParameterAsString("outcol").empty() ||getParameterAsString("outcol")=="unknown")
  		colNameOutput="__Mathop_";
	 else
		colNameOutput=getParameterAsString("outcol");
	for(unsigned int i=0;i<m_nOfCol;i++) 
		colList[i]=m_tables[0]->getColId(m_colOldName[i]);

	if(append)
		for(unsigned int i=0;i<m_tables[0]->getNumberOfColumns();i++) 
			tableMath.addCol(m_tables[0]->getColName(i));	
		
	if(!tableMath.addCol(colNameOutput))
	{
		std::cerr<<"Error: Invalid or duplicate column name in existing table: "<<colNameOutput<<std::endl;
		return false;
	}	
//	tableMath.writeHeader();  
	bool writeHeader=true;
 	unsigned long long int totEle=nOfEle;
	while(totEle!=0)
	{
		fromRow=startCounter;
		toRow=fromRow+maxEle-1;
		if(toRow>totRows-1)toRow=totRows-1;
	  	m_tables[0]->getColumn(colList, m_nOfCol, fromRow, toRow, m_fArray);
   	        FunctionParser fparser;
		std::stringstream varList;
		for(unsigned int i=0;i<m_nOfCol-1;i++)
		{
			varList<<m_colName[i]<<",";
		}
		varList<<m_colName[m_nOfCol-1];
    		int parseReuslt = fparser.Parse(m_mathExpression,varList.str());	
        	if(parseReuslt!=-1)
		{ 
			std::cerr<< m_mathExpression.substr(0,parseReuslt+1)<< "^"<<std::endl<< fparser.ErrorMsg()<<std::endl;
			delete [] colList;
			return false;
		} else
		{	if(writeHeader)
			{
				writeHeader=false;
				tableMath.writeHeader();//overwrite if table exist!
			}
		}
		fparser.setNumberOfRow(maxEle);
		fparser.setNumberOfCol(m_nOfCol);
		fparser.setArray(m_fArray,m_result[0]); //QUI VERIFIA
		fparser.Eval(NULL);
		unsigned int putColList[1];
		putColList[0]=tableMath.getColId(colNameOutput);
		tableMath.putColumn(putColList,1,fromRow, toRow, m_result);
		startCounter=toRow+1;
		totEle=totEle-(toRow-fromRow+1);
	}
delete [] colList;
return true;
}
