#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <sstream>
#include <set>
#include <fstream>
#include "vstable.h"
//#include "vstableop.h"
#include "vsselectcolumnsop.h"

//---------------------------------------------------------------------
VSSelectColumnsOp::VSSelectColumnsOp()
//---------------------------------------------------------------------
{
 
}
//---------------------------------------------------------------------
VSSelectColumnsOp::~VSSelectColumnsOp()
//---------------------------------------------------------------------
{
 
}
//---------------------------------------------------------------------
void VSSelectColumnsOp::printHelp()
//---------------------------------------------------------------------
{
	std::cout<<"Create e new table using (or exluding) one or more fields of a data table. The default case produces output table including only listed fields."<<std::endl<<std::endl;
	std::cout<<"Usage: VisIVOFilters --op selcolumns --field parameters  [--delete][--out filename_out.bin] [--history] [--historyfile filename.xml] [--help] [--file] inputFile.bin"<<std::endl;

	std::cout<<"Example: VisIVOFilters --op selcolumns --field X Y --out filename_out.bin --file inputFile.bin"<<std::endl;

	std::cout<<"Note:"<<std::endl;
	std::cout<<"--field Valid columns names of the input table."<<std::endl;
	std::cout<<"--delete  Produce output table excluding only field listed in the --field option."<<std::endl;
	std::cout<<"--out Output table filename. Default name is given."<<std::endl;
    std::cout<<"--history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
    std::cout<<"--historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;
	std::cout<<"--file  Input table filename."<<std::endl;
	std::cout<<"--help produce this output "<<std::endl;

	return;
}

//---------------------------------------------------------------------
bool VSSelectColumnsOp::execute()
//---------------------------------------------------------------------
{
	if(isParameterPresent("field"))
	{
	  if(getParameterAsString("field").empty() || getParameterAsString("field")=="unknown" )
	  {
		std::cerr<<"vsselectcolumnsop: No list of fileds is given"<<std::endl;
		return false;
	  }
	}else
	{
	  if(getParameterAsString("list").empty() || getParameterAsString("list")=="unknown" )
	  {
		std::cerr<<"vsselectcolumnsop: No list of fileds is given"<<std::endl;
		return false;
	  }
	}
	if(m_tables[0] ->getType()!="float")
	{
		std::cerr<<"vsselectcolumnsop: This filter can be applied to float Table only"<<std::endl;
		return false;

	}
	std::set<unsigned int> colNumberSet;
	std::stringstream ssListparameters;
	
	if(isParameterPresent("field"))
	  ssListparameters.str(getParameterAsString("field"));
	else
	  ssListparameters.str(getParameterAsString("list"));
	  
	std::string paramField="";
	std::string paramFieldGlobal;
	std::set<unsigned int>::iterator iter;

	if(!isParameterPresent("delete")) 
	{
		colNumberSet.clear();
		while (!ssListparameters.eof())
		{
			ssListparameters>>paramField;

			if(m_tables[0] -> getColId(paramField)>=0)
			{ 
				paramFieldGlobal.append(paramField);
				colNumberSet.insert(m_tables[0] -> getColId(paramField));
			}
		}
	} else
	{
		paramFieldGlobal="NO";
		for(unsigned int i=0;i<m_tables[0] -> getNumberOfColumns();i++)colNumberSet.insert(i);
		while (!ssListparameters.eof())
		{
			ssListparameters>>paramField;
 
			if(m_tables[0] -> getColId(paramField)>=0)
			{
				paramFieldGlobal.append(paramField);
				iter=colNumberSet.find(m_tables[0] -> getColId(paramField));
				colNumberSet.erase(iter);
			}
		}   
   
	}
	int nCols = colNumberSet.size();

	/*** Output Filename  **/
	std::stringstream fileNameOutputSStream;
	fileNameOutputSStream<<getParameterAsString("out");
	std::string fileNameOutput;

	if(fileNameOutputSStream.str()=="")
	{
		std::string filenameInputTable=m_tables[0]->getLocator();
		int len=filenameInputTable.length();
		fileNameOutputSStream<<filenameInputTable.substr(0, len-4)<<"_selcols_"<< paramFieldGlobal<<".bin";
	}

	fileNameOutput=fileNameOutputSStream.str();
	if(fileNameOutput.find(".bin") == std::string::npos)
   		fileNameOutput.append(".bin");


	m_realOutFilename.push_back(fileNameOutput);
	//Clean existing tab

	remove(fileNameOutput.c_str());
	unsigned long long int totRows=m_tables[0]->getNumberOfRows();
	int maxInt=getMaxNumberInt();
	int maxEle;

	if(totRows>maxInt) 
		maxEle=maxInt; 
	else 
		maxEle=totRows;

	float **fArray=new  float*[1];
try
{
	fArray[0]=new  float[maxEle];
}
catch(std::bad_alloc &e)
{
	fArray[0]=NULL;
}

	while(fArray[0]==NULL)
	{
		maxEle=maxEle-10000;
		if(maxEle<=0)
		{
			std::cerr<<"Bad **fArray allocation. Select Column terminated"<<std::endl;
			return false;//QUI stranissimo il controllo ritorna QUI alla fine
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

	unsigned long long int nOfEle,fromRow,toRow;
	unsigned int colList[1];
	unsigned int nOfCol = 1;
	int startCounter;

	VSTable percTable;
	percTable.setLocator(fileNameOutput);
	percTable.setType("float");
	percTable.setNumberOfRows(m_tables[0]->getNumberOfRows());
	percTable.setIsVolume(m_tables[0]->getIsVolume());
	if(m_tables[0]->getIsVolume())
	{
  	 	percTable.setCellNumber(m_tables[0]->getCellNumber()[0],
                           m_tables[0]->getCellNumber()[1],
                           m_tables[0]->getCellNumber()[2]);
   		percTable.setCellSize(m_tables[0]->getCellSize()[0],
                           m_tables[0]->getCellSize()[1],
                           m_tables[0]->getCellSize()[2]);
	}

	percTable.setEndiannes(m_tables[0]->getEndiannes());

	std::set<unsigned int>::iterator end = colNumberSet.end();;
  
	for(iter = colNumberSet.begin(); iter != end; iter++) percTable.addCol(m_tables[0]->getColName(*iter));
	percTable.writeTable();

	unsigned long long int j=-1; 
	for(iter = colNumberSet.begin(); iter != end; iter++) 
	{	
		j++;
		startCounter=0;
		nOfEle=totRows;
		while(nOfEle!=0)
		{
			fromRow=startCounter;
			toRow=fromRow+maxEle-1;
			if(toRow>totRows-1)toRow=totRows-1;
			colList[0]=*iter;
	  		m_tables[0]->getColumn(colList, nOfCol, fromRow, toRow, fArray);
			colList[0]=j;
			percTable.putColumn(colList, nOfCol, fromRow, toRow, fArray);
			nOfEle=nOfEle-(toRow-fromRow+1);
			startCounter=toRow+1;
		}
	}
	delete [] fArray[0];
	delete [] fArray;
	return true;

}
