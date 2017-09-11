//
// C++ Implementation: vsextractlistrowsop
//
// Description: 
//
//
// Author:  <>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
// UNDEBUGGED for BINARY lists
//#include <cstdlib>
//#include <cstring>
#include <cstdlib>
#include <cstring>
#include <iostream>

/*#include <string>*/
#include <sstream>
/*#include <set>*/
#include <fstream>

#include "vsextractlistrowsop.h"
#include "vstable.h"
#ifdef WIN32
	#include <time.h>
#endif
const unsigned int VSExtractListRowsOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSExtractListRowsOp::MIN_NUMBER_OF_ROW = 100;
const long int VSExtractListRowsOp::MAX_SEEK = 500000000;
//---------------------------------------------------------------------
VSExtractListRowsOp::VSExtractListRowsOp()
{
m_multlistFormatBinary=true;
m_numberOfLists=0;
 m_fArray=NULL;
 m_wfArray=NULL;
 m_list=NULL;
 m_lllist=NULL;
 m_nOfRow=0;
 m_nOfCol=0;	
 m_nOfEle=0;
 m_totalListEle=0;
 m_oneList=false;
 m_numberLists=false;
 m_listElemnts=false;
m_elementsInLists=0;
 m_residualLoadList=0;
 m_startPoint=0;
m_listTotal=false;

}
//---------------------------------------------------------------------

//---------------------------------------------------------------------
VSExtractListRowsOp::~VSExtractListRowsOp()
{
if(m_fArray!=NULL)
   for(unsigned int i=0;i<m_nOfCol;i++)
	{
		if(m_fArray[i]!=NULL) delete [] m_fArray[i];
	}
if(m_fArray!=NULL) delete [] m_fArray;
if(m_wfArray!=NULL)
   for(unsigned int i=0;i<m_nOfCol;i++)
	{
		if(m_wfArray[i]!=NULL) delete [] m_wfArray[i];
	}
if(m_wfArray!=NULL) delete [] m_wfArray;
if(m_list!=NULL) delete [] m_list;
if(m_lllist!=NULL) delete [] m_lllist;
m_numberOfElements.clear();

}
//---------------------------------------------------------------------

//---------------------------------------------------------------------
void VSExtractListRowsOp::printHelp()
//---------------------------------------------------------------------
{
	std::cout<<"Create e new table extracting rows from a data table. A multi-list is given in ascii or binary format (unsigned long long int, or int), or, from yust a list (see --onelist option). The multi-list has the following structure:"<<std::endl<<std::endl; 
	std::cout<<"	Number  NL of lists;"<<std::endl; 
	std::cout<<"		NL sequences of:"<<std::endl;  
	std::cout<<"			1) Number N0 of  elements in the list,"<<std::endl; 
	std::cout<<"			2) N0 element ids (numbers of rows)."<<std::endl<<std::endl;
	std::cout<<"Option can be given to provide the NL number. In this case the multi-list file must not contain this information."<<std::endl;
	std::cout<<"Option can be given to provide the N0 number. In this case the multi-list file must not contain this information, and it is a multilist, each  list must contain N0 elements"<<std::endl;
	std::cout<<"Usage: VisIVOFilters --op extractlist  --multilist filename_list [--binaryint] [--asciilist] [--numberlists nl]  [--listelements n0] [--onelist] [--history] [--historyfile filename.xml] [--out filename_out.bin] [--help] [--file] inputFile.bin"<<std::endl<<std::endl;

	std::cout<<"Example: VisIVOFilters --op extractlist  --multilist mylist.txt --asciilist  --out list_extract.bin --file inputFile.bin"<<std::endl;

	std::cout<<"Note:"<<std::endl;
	std::cout<<"--multilist Is the multilist filename ."<<std::endl;
	std::cout<<"--binaryint The default multilist file format is binary unsigned long long int. If this parameter is specified the file is binary int."<<std::endl;
	std::cout<<"--asciilist The multilist file format is binary. If this parameter is specified the file is an ascii text."<<std::endl;
	std::cout<<"--numberlists The multilist file format is just a sequence of nl lists specified in this option.  Each list  starts with the number of elements in the list"<<std::endl;
	std::cout<<"--listelements The multilist file format is just a sequence of nl lists.  Each list has  the same number of n0 elements. This option requires that the --numberlists option is specified, otherwise it is ignored."<<std::endl;
	std::cout<<"--onelist If this option is given, the multlist file is considered as  one list only. Each element is the ID of the particle to be extracted. The --numberlists and --listelements options will be ignored."<<std::endl;
	std::cout<<"--out Output table filename. Default name is given."<<std::endl;
    std::cout<<"--history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
    std::cout<<"--historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;
    std::cout<<"--file  Input table filename."<<std::endl;
	std::cout<<"--help produce this output "<<std::endl;


	return;
}
//---------------------------------------------------------------------
bool VSExtractListRowsOp::allocatefArray()
//---------------------------------------------------------------------
{
unsigned long long int tempLL=getMaxNumberInt();
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
try
{
m_wfArray=new  float*[m_nOfCol];
}
catch(std::bad_alloc &e)
{
m_wfArray=NULL;
}
if(m_wfArray==NULL)
{
		delete [] m_fArray;
		m_fArray=NULL;	
		return false;
}
for(int i=0;i<m_nOfCol;i++)
{
	m_fArray[i]=NULL;
	m_wfArray[i]=NULL;
}
bool goodAllocation=false;
while(!goodAllocation)
{
	goodAllocation=true;
	for(unsigned int i=0;i<m_nOfCol;i++)
	{
try
{
		m_fArray[i] = new float[m_nOfEle];
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
			if(m_nOfEle==MIN_NUMBER_OF_ROW)
			{ 
				delete [] m_fArray;
				m_fArray=NULL;
				return false;
			}
			if(m_nOfEle<=MAX_NUMBER_TO_REDUCE_ROW)
				m_nOfEle=MIN_NUMBER_OF_ROW;
			else
				m_nOfEle=m_nOfEle-MAX_NUMBER_TO_REDUCE_ROW;
			break;
 		}
	}

	for(unsigned int i=0;i<m_nOfCol;i++)
	{
try
{
		m_wfArray[i] = new float[m_nOfEle];
}
catch(std::bad_alloc &e)
{
	m_wfArray[i]=NULL;
}
 		if(m_wfArray[i]==NULL) 
 		{	
			goodAllocation=false;
			for(unsigned int j=0;j<m_nOfCol;j++) 
				delete [] m_fArray[j];
			for(unsigned int j=0;j<i;j++) 
				delete [] m_wfArray[j];
			if(m_nOfEle==MIN_NUMBER_OF_ROW)
			{ 
				delete [] m_fArray;
				m_fArray=NULL;
				delete [] m_wfArray;
				m_wfArray=NULL;
				return false;
			}
			if(m_nOfEle<=MAX_NUMBER_TO_REDUCE_ROW)
				m_nOfEle=MIN_NUMBER_OF_ROW;
			else
				m_nOfEle=m_nOfEle-MAX_NUMBER_TO_REDUCE_ROW;
			break;
 		}
	}
	if(!goodAllocation)
		continue;


	if(isParameterPresent("binaryint"))
	{
try
{
	m_list = new unsigned int [m_nOfEle];
}
catch(std::bad_alloc &e)
{
	m_list=NULL;
}
 	if(m_list==NULL) 
 	{	
		goodAllocation=false;
		for(unsigned int j=0;j<m_nOfCol;j++) 
		{
			delete [] m_fArray[j];
			delete [] m_wfArray[j];
		}
		if(m_nOfEle==MIN_NUMBER_OF_ROW)
		{ 
			delete [] m_fArray;
			m_fArray=NULL;
			delete [] m_wfArray;
			m_wfArray=NULL;
			return false;
		}
		if(m_nOfEle<=MAX_NUMBER_TO_REDUCE_ROW)
			m_nOfEle=MIN_NUMBER_OF_ROW;
		else
			m_nOfEle=m_nOfEle-MAX_NUMBER_TO_REDUCE_ROW;
		break;
 	}
	}
	if(!isParameterPresent("binaryint"))
	{
try
{
	m_lllist = new unsigned long long int [m_nOfEle];
}
catch(std::bad_alloc &e)
{
	m_lllist=NULL;
}
 	if(m_lllist==NULL) 
 	{	
		goodAllocation=false;
		for(unsigned int j=0;j<m_nOfCol;j++) 
		{
			delete [] m_fArray[j];
			delete [] m_wfArray[j];
		}
		if(m_nOfEle==MIN_NUMBER_OF_ROW)
		{ 
			delete [] m_fArray;
			m_fArray=NULL;
			delete [] m_wfArray;
			m_wfArray=NULL;
			return false;
		}
		if(m_nOfEle<=MAX_NUMBER_TO_REDUCE_ROW)
			m_nOfEle=MIN_NUMBER_OF_ROW;
		else
			m_nOfEle=m_nOfEle-MAX_NUMBER_TO_REDUCE_ROW;
		break;
 	}
	}


}
return true;
}

//---------------------------------------------------------------------
bool VSExtractListRowsOp::readListsElements()
//---------------------------------------------------------------------
{
if(m_oneList)
{
	unsigned long long int tmp, tmpsave=0;
	m_totalListEle=0;
	while(!m_fileInput.eof()) 
	{   
		if(m_multlistFormatBinary)
		{
			m_fileInput.seekg (0, std::ios::end);
			unsigned long long int length=m_fileInput.tellg();
			if(isParameterPresent("binaryint"))
				m_totalListEle=length/sizeof(int);
			else
				m_totalListEle=length/sizeof(long long int);
		}
   		else
		{
			m_fileInput>> tmp;
//			std::clog<<tmp<<std::endl;
			if(tmp != tmpsave)
			{
			     m_totalListEle++;
			     tmpsave=tmp;
			}
		}		

	}
	m_numberOfElements.push_back(m_totalListEle);
m_fileInput.close();
return true;
}

if(m_numberOfLists==0)
{
   if(m_multlistFormatBinary)
	if(isParameterPresent("binaryint"))
	    m_fileInput.read((char *)&m_numberOfLists,sizeof(int));
	else
	{
	   unsigned long long int lltmp;
	    m_fileInput.read((char *)&lltmp,sizeof(unsigned long long int));
	    m_numberOfLists=lltmp;
	}
   else
	m_fileInput>> m_numberOfLists;	
}

for(int i=0;i<m_numberOfLists;i++)
{
	unsigned long long int numberOfElements;

	if(m_elementsInLists==0)
	{
	  if(m_multlistFormatBinary)
		if(isParameterPresent("binaryint"))
		{
			int itmp;
			m_fileInput.read((char *)&itmp,sizeof(int));
			numberOfElements=itmp;
		}else
			m_fileInput.read((char *)&numberOfElements,sizeof(unsigned long long int));
	  else
		m_fileInput>> numberOfElements;
	} else
		numberOfElements=m_elementsInLists;

	m_totalListEle+=numberOfElements;
	m_numberOfElements.push_back(numberOfElements);
	if(m_multlistFormatBinary)
	{
		int indexOffset;
		if(isParameterPresent("binaryint"))
			indexOffset=sizeof(int);
		else
			indexOffset=sizeof(long long int);
		
		std::streamoff seekShift=0;
		if(numberOfElements>MAX_SEEK)
			seekShift=MAX_SEEK;
		else
			seekShift=numberOfElements;
		while(numberOfElements>0)
		{
			m_fileInput.seekg(seekShift*indexOffset,std::ios::cur);
			numberOfElements-=seekShift;
			if(numberOfElements<seekShift)
				seekShift=(long int)numberOfElements;
		}
	}else
	{
		std::string null;
		for(unsigned long long int j=0;j<numberOfElements;j++)
			m_fileInput>> null;	
	}
}


m_fileInput.close();
return true;
}

//---------------------------------------------------------------------
bool VSExtractListRowsOp::execute()
//---------------------------------------------------------------------
{
if(m_tables[0]->getIsVolume())
{
 	std::cerr<<"vsextractlistrowsop: operation not allowed: a volume is selected"<<std::endl; 
	return false;
}

if(isParameterPresent("onelist"))
{
	m_oneList=true,
	m_numberOfLists=1;
} else
{
	if(isParameterPresent("numberlists"))
	{	
		m_numberLists=true;
		m_numberOfLists=getParameterAsInt("numberlists");
		if(m_numberOfLists<=0)
		{
			std::cerr<<"Invalid option numberlists. Operation aborted"<<std::endl;
			return false; 
		}
		if(isParameterPresent("listelements"))
		{
			m_listElemnts=true;
			m_elementsInLists=getParameterAsInt("listelements");
			if(m_elementsInLists<=0)
			{
				std::cerr<<"Invalid option listelements. Operation aborted"<<std::endl;
				return false; 
			}

			
		}
	}
}
std::string filename;

if(getParameterAsString("multilist").empty() || getParameterAsString("multilist")=="unknown" )
{
	std::cerr<<"vsextractlistrowsop: No multilist file is given"<<std::endl;
	return false;
} 
if(isParameterPresent("asciilist"))
	m_multlistFormatBinary=false;
		
filename=getParameterAsString("multilist");

if(m_multlistFormatBinary)
	m_fileInput.open(filename.c_str(),std::ios::binary);
else
	m_fileInput.open(filename.c_str());
if(!m_fileInput.is_open())
{
	std::cerr<<"vsextractlistrowsop: Error opening multilist file "<<filename<<" Operation aborted."<<std::endl;
	return false;
} 


std::stringstream fileNameOutputSStream;
std::string fileNameOutput;
if(!isParameterPresent("out") ||getParameterAsString("out").empty() || getParameterAsString("out")=="unknown")
{
  	std::string filenameInputTable=m_tables[0]->getLocator();
  	int len=filenameInputTable.length();
	time_t rawtime;
	struct tm * timeinfo;
	char buffer [80];
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	strftime (buffer,80,"%Y%m%d%H%M",timeinfo);
  	fileNameOutputSStream<<filenameInputTable.substr(0, len-4)<<"_multilist_"<<buffer<<".bin";  //QUI verificare
} else
	fileNameOutputSStream<<getParameterAsString("out");


fileNameOutput=fileNameOutputSStream.str();
if(fileNameOutput.find(".bin") == std::string::npos)
	fileNameOutput.append(".bin");

m_realOutFilename.push_back(fileNameOutput);

m_nOfRow=m_tables[0]->getNumberOfRows();
m_nOfCol= m_tables[0]->getNumberOfColumns();

if(m_nOfRow>getMaxNumberInt()) 
	m_nOfEle= getMaxNumberInt();
else 
	m_nOfEle= m_nOfRow;



unsigned int *colList=NULL;
try
{
	colList=new unsigned int[m_nOfCol];
}
catch(std::bad_alloc &e)
{
	colList=NULL;
}
if(colList==NULL)
{
	std::cerr<<"vsextractlistrowsop: Invalid colList allocation. Operation aborted"<<std::endl;
	return false;
} 
for(unsigned int i=0;i<m_nOfCol;i++)
	colList[i]=i;

if(!allocatefArray())
{
	delete [] colList;
	return false;
}
if(!readListsElements())
{
	delete [] colList;
	return false;
}
VSTable multiListTable;
multiListTable.setLocator(fileNameOutput);
multiListTable.setType("float");

#ifdef VSBIGENDIAN
	std::string endianism="big";
	
#else	
	std::string endianism="little";
#endif

multiListTable.setEndiannes(endianism);
multiListTable.setNumberOfRows(m_totalListEle);
for(unsigned int i=0;i<m_nOfCol;i++)
	multiListTable.addCol(m_tables[0]->getColName(i));	
multiListTable.writeHeader();


std::string strdummy;
unsigned long long int idummy;
unsigned long long int index=0;
unsigned long long int fromRow=0;
unsigned long long int toRow=0;
unsigned long long int fromwRow=0;
unsigned long long int towRow=0;
int wCount=0;
int lCount=0;
int  globalEle=m_nOfEle;
std::ifstream fileInp;
while(index< m_nOfRow)
{	
	if(m_multlistFormatBinary)
		fileInp.open(filename.c_str(),std::ios::binary);
	else
		fileInp.open(filename.c_str());
	if(!fileInp.is_open())
	{
		std::cerr<<"vsextractlistrowsop: Error opening multilist file "<<filename<<" Operation aborted."<<std::endl;
		delete [] colList;
		return false;
	} else
//	  std::clog<<"filename to reopen "<<filename<<std::endl;

	if(!m_oneList) 
 		if(!m_numberLists)
 		{ 
   			if(m_multlistFormatBinary)
				if(isParameterPresent("binaryint"))
					fileInp.read((char *)&idummy,sizeof(int));
				else
					fileInp.read((char *)&idummy,sizeof(unsigned long long int));
   			else
				fileInp>> strdummy;	
 		}
	m_tables[0]->getColumn(colList,m_nOfCol,index,index+globalEle-1,m_fArray);
	int nReadListEle=0;	
	for(int i=0;i<m_numberOfElements.size();i++)
	{
		if(!m_oneList)
		{ 
   		   if(m_multlistFormatBinary)
		   {
		       if(isParameterPresent("binaryint"))
			  fileInp.read((char *)&idummy,sizeof(int));
		       else
			  fileInp.read((char *)&idummy,sizeof(unsigned long long int));
		   } else
			fileInp>> strdummy;
		}
		bool stillRead=true;
		int residualEle=0;
		int alreadyRead=0;
		while(stillRead)
		{
		       if(residualEle==0)
				nReadListEle=m_numberOfElements[i];
		       else
			    	nReadListEle=residualEle;
		       if(nReadListEle+lCount > m_nOfEle)
		       {
			 stillRead=true;
			 nReadListEle=m_nOfEle-lCount;
			 residualEle=m_numberOfElements[i]-nReadListEle-alreadyRead;
		       } else
			  stillRead=false;
		       if(m_multlistFormatBinary)
		       {
		         if(isParameterPresent("binaryint"))
			    fileInp.read((char *)&m_list[lCount],sizeof(int)*nReadListEle);
		         else
			    fileInp.read((char *)&m_lllist[lCount],sizeof(unsigned long long int)*nReadListEle);
		       } else
			for(int j=0;j<nReadListEle;j++)
			{
				fileInp>> m_lllist[j+lCount];
//				std::clog<<j<<" "<< m_lllist[j+lCount]<<std::endl;
			}
		  alreadyRead+=nReadListEle;
		  lCount+=nReadListEle;
		  if(lCount==m_nOfEle)
		  {
		    for(int k=0;k<lCount;k++)
		       if(isParameterPresent("binaryint"))
		       {
		  		if(m_list[k]>=m_nOfRow || m_list[k]<0)
		  		{
					std::cerr<<"Invalid Id is given: "<<m_list[k]<<std::endl<<" Operation Aborted."<<std::endl;
					return false;			
		  		}

				if(m_list[k]<index+m_nOfEle && m_list[k]>=index)
				{
					for(int m=0;m<m_nOfCol;m++)
					{
					   m_wfArray[m][wCount]=m_fArray[m][m_list[k]-index];
//					   std::clog<<"m_wfArray["<<m<<"]["<<wCount<<"]"<< m_wfArray[m][wCount]<<std::endl;
					}
					wCount++;
				}			
			}else
			{
				if(m_lllist[k]>=m_nOfRow || m_lllist[k]<0)
		  		{
					std::cerr<<"Invalid Id is given: "<<m_lllist[k]<<std::endl<<" Operation Aborted."<<std::endl;
					return false;			
		  		}

				if(m_lllist[k]<index+m_nOfEle && m_lllist[k]>=index)
				{
					for(int m=0;m<m_nOfCol;m++)
					{
					   m_wfArray[m][wCount]=m_fArray[m][m_lllist[k]-index];
//					   std::clog<<"m_wfArray["<<m<<"]["<<wCount<<"]"<< m_wfArray[m][wCount]<<std::endl;
					}

					wCount++;
				}
			}
			if(wCount!=0)
			{
				towRow=fromwRow+wCount-1;
				multiListTable.putColumn(colList, m_nOfCol,fromwRow,towRow,m_wfArray);
				fromwRow=towRow+1;
			}
			wCount=0;
			lCount=0;
			
		    }
		}
	}
	for(int k=0;k<lCount;k++)
	{
		if(isParameterPresent("binaryint"))
		{
		  if(m_list[k]>=m_nOfRow || m_list[k]<0)
		  {
			std::cerr<<"Invalid Id is given: "<<m_list[k]<<std::endl<<" Operation Aborted."<<std::endl;
			return false;			
		  }
		  if(m_list[k]<index+globalEle && m_list[k]>=index)
		  {
			for(int m=0;m<m_nOfCol;m++)
			{
				m_wfArray[m][wCount]=m_fArray[m][m_list[k]-index];
//			       std::clog<<"m_wfArray["<<m<<"]["<<wCount<<"]"<< m_wfArray[m][wCount]<<std::endl;
			}
			wCount++;
		 }
		} else
		{
		  if(m_lllist[k]>=m_nOfRow || m_lllist[k]<0)
		  {
			std::cerr<<"Invalid Id is given: "<<m_lllist[k]<<std::endl<<" Operation Aborted."<<std::endl;
			return false;			
		  }

		  if(m_lllist[k]<index+globalEle && m_lllist[k]>=index)
		  {
			for(int m=0;m<m_nOfCol;m++)
			{
				m_wfArray[m][wCount]=m_fArray[m][m_lllist[k]-index];
//			       std::clog<<"m_wfArray["<<m<<"]["<<wCount<<"]"<< m_wfArray[m][wCount]<<std::endl;
			}
			wCount++;
		   }
		}
	}	
	if(wCount!=0)
	{
	   towRow=fromwRow+wCount-1;
	   multiListTable.putColumn(colList, m_nOfCol,fromwRow,towRow,m_wfArray);
	   fromwRow=towRow+1;
         }
	wCount=0;
	lCount=0;

	fileInp.close();
	index+= globalEle;
	if(index+globalEle>= m_nOfRow) 
		globalEle= m_nOfRow-index;
}

delete [] colList;
return true;
}