#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <sstream>
#include <set>
#include <fstream>
#include "vstable.h"
//#include "vstableop.h"
#include "vsmoduleop.h"
#include <cmath>

const unsigned int VSModuleOp::MIN_NUMBER_OF_ROW = 10000;
const unsigned int VSModuleOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;

//---------------------------------------------------------------------
VSModuleOp::VSModuleOp()
//---------------------------------------------------------------------
{
m_fArray=NULL; 
m_nOfCol=3;
m_nOfRows=0;
m_nOfEle=0;
}
//---------------------------------------------------------------------
VSModuleOp::~VSModuleOp()
//---------------------------------------------------------------------
{

if(m_fArray!=NULL)
{
	for(int i=0;i<m_nOfCol;i++)
	{
		if(m_fArray[i]!=NULL)
			delete [] m_fArray[i];
			
	}
	delete [] m_fArray; 
}
}

//---------------------------------------------------------------------
bool VSModuleOp::allocateArray()
//---------------------------------------------------------------------
{
unsigned long long int tempLL=getMaxNumberInt()*2;
if(((unsigned long long int)m_nOfEle*m_nOfCol)>tempLL) 
	m_nOfEle=(int)tempLL/m_nOfCol;

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
		for(int i=0;i<m_nOfCol;i++)
		{
try
{
			m_fArray[i]=new  float[m_nOfEle];
}
catch(std::bad_alloc &e)
{
	m_fArray[i]=NULL;
}

			if(m_fArray[i]==NULL) 
			{	
				goodAllocation=false;
				for(int j=0;j<i;j++) 
					delete [] m_fArray[j];
				if(m_nOfEle==MIN_NUMBER_OF_ROW)
				{ 
					delete [] m_fArray;
					m_fArray=NULL;
					return false;
				}
				m_nOfEle=m_nOfEle-MAX_NUMBER_TO_REDUCE_ROW;
				if(m_nOfEle<=MIN_NUMBER_OF_ROW) m_nOfEle=MIN_NUMBER_OF_ROW;
				break;
			}

		}
		
	}

	return true;
}

//---------------------------------------------------------------------
void VSModuleOp::printHelp()
//---------------------------------------------------------------------
{
	std::cout<<"Creates e new table (or a new field) computing the module of three fields of the input data table"<<std::endl<<std::endl;
	std::cout<<"Usage: VisIVOFilters --op module --field parameters [--append] [--outcol colname] [--out filename_out.bin][--help] [--file] inputFile.bin"<<std::endl;

	std::cout<<"Example: VisIVOFilters --op module --field Vx Vy Vz  --outcol Module --append --file inputFile.bin"<<std::endl;

	std::cout<<"Note: A New field is produced Vx, Vy and Vz sqrt(Vx^2+Vy^2+Vz^2) . Default options: a new table with only the new field is produced "<<std::endl;
	std::cout<<"      --field Three valid fields to compute the module"<<std::endl;
	std::cout<<"      --outcol. Column name of the new field"<<std::endl;
	std::cout<<"      --append. No new table will be cretaed. The original table will have the new field "<<std::endl;
	std::cout<<"      --out. Name of the new table. Ignored if --append is specified."<<std::endl;
    std::cout<<"      --history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
    std::cout<<"        --historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;
	std::cout<<"      --help produce this output "<<std::endl;


	return;
}

//---------------------------------------------------------------------
bool VSModuleOp::execute()
//---------------------------------------------------------------------
{
	if(getParameterAsString("field").empty() || getParameterAsString("field")=="unknown" )
	{
		std::cerr<<"VSModuleOp: No list of fileds is given"<<std::endl;
		return false;
	}

	std::stringstream ssField(getParameterAsString("field"));
	std::string *col;
	try
	{
		col = new std::string [m_nOfCol];
	}
	catch(std::bad_alloc &e)
	{
		col=NULL;
	}
	if (col == NULL)
		return false;
	unsigned int outColId[1];
	
	unsigned int * colId;
	try{
	colId = new unsigned int [m_nOfCol];
	}
	catch (std::bad_alloc &e) 
	{
	colId = NULL;
	}
	if (colId==NULL)
	{
		delete [] col;
		return false;
	}
	
	
	
	
	for(int i=0;i<m_nOfCol;i++)
	{
		if(ssField.eof())
		{
			if(col!=NULL) delete [] col;
			if(colId!=NULL) delete [] colId;
			return false;
		}
		ssField>>col[i];
		colId[i]=m_tables[0]->getColId(col[i]);
		if(colId[i]==-1)
		{
			std::cerr<<"Invalid filed name"<< col[i]<<std::endl;
			if(col!=NULL) delete [] col;
			if(colId!=NULL) delete [] colId;
			return false;
		}
	}
	std::string outColName;

	if(!isParameterPresent("outcol") || getParameterAsString("outcol")=="unknown")
		outColName="ModuleOp";
        else
		outColName=getParameterAsString("outcol");

	std::string outName;

	if(!isParameterPresent("out") || getParameterAsString("out")=="unknown")
		outName="VSModuleOp.bin";
        else
		outName=getParameterAsString("out");




        m_nOfRows=m_tables[0]->getNumberOfRows();

	VSTable newTable;
	if(!isParameterPresent("append"))
	{
                m_realOutFilename.push_back(outName);
		newTable.setType("float");
		newTable.setNumberOfRows(m_nOfRows);
		newTable.setIsVolume(m_tables[0]->getIsVolume());
		newTable.setCellSize(m_tables[0]->getCellSize()[0],m_tables[0]->getCellSize()[1],m_tables[0]->getCellSize()[2]);
		newTable.setCellNumber(m_tables[0]->getCellNumber()[0],m_tables[0]->getCellNumber()[1],m_tables[0]->getCellNumber()[2]);
		newTable.setEndiannes(m_tables[0]->getEndiannes());
		newTable.addCol(outColName);
		newTable.setLocator(outName);
		newTable.writeHeader();
	} else
	{
                m_realOutFilename.push_back(m_tables[0]->getLocator());
		m_tables[0]->addCol(outColName);
		m_tables[0]->writeHeader();
	}



	if(m_nOfRows>getMaxNumberInt())
		m_nOfEle=getMaxNumberInt();
	else
		m_nOfEle=m_nOfRows;

	bool goodAllocation=allocateArray();	
	if(!goodAllocation)
	{
		std::cerr<<"Invalid allocation"<<std::endl;
		if(col!=NULL) delete [] col;
		if(colId!=NULL) delete [] colId;
		return false;
	}
	unsigned long long int index=0;
	while(index<m_nOfRows)
	{

		m_tables[0]->getColumn(colId,m_nOfCol,index,index+m_nOfEle-1,m_fArray);

		for(int j=0;j<m_nOfCol;j++)
		{
			for(int i=0;i<m_nOfEle;i++)
			{
			     m_fArray[j][i]=m_fArray[j][i]*m_fArray[j][i];
			}
			
		}
		for(int j=1;j<m_nOfCol;j++)
		{
			for(int i=0;i<m_nOfEle;i++)
			{
			     m_fArray[0][i]+=m_fArray[j][i];
			}
		}
		for(int i=0;i<m_nOfEle;i++)
			m_fArray[0][i]=sqrt(m_fArray[0][i]);


		if(!isParameterPresent("append"))
 		{	
			outColId[0]=0;
			newTable.putColumn(outColId,1,index,index+m_nOfEle-1,m_fArray);		
		} else
		{
			outColId[0]=m_tables[0]->getNumberOfColumns()-1;
			m_tables[0]->putColumn(outColId,1,index,index+m_nOfEle-1,m_fArray);		
			
		}
		index+=m_nOfEle;
		if(index+m_nOfEle>=m_nOfRows)
			m_nOfEle=m_nOfRows-index;
		
	}
if(col!=NULL) delete [] col;
if(colId!=NULL) delete [] colId;
return true;
}

