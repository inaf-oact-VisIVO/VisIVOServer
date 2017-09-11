#include <cstdlib>
#include <cstring>
#include "vspolarop.h"
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include "vstable.h"

 const unsigned int VSPolarOp::MIN_NUMBER_OF_ROW = 100; 
 const unsigned int VSPolarOp::MAX_NUMBER_TO_REDUCE_ROW = 10000;
//*******************************************************************
VSPolarOp::VSPolarOp() 
{
 m_fArray=NULL; 
 m_nOfCol=5; 
 m_nOfColIn=3; 
 m_nOfRow=0;
 m_nOfElements=0; 
}
//********************************************************************
VSPolarOp::~VSPolarOp() 
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
//******************************************************************
bool VSPolarOp::allocateArray()
//******************************************************************
{
unsigned long long int tempLL=getMaxNumberInt()*2;
if(((unsigned long long int)m_nOfElements*m_nOfCol)>tempLL) 
	m_nOfElements=(int)tempLL/m_nOfCol;

 try
 {
  m_fArray=new  float*[m_nOfCol]; //crea l'array facendo attenzione al possibile errore di allocazione
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
	for(unsigned int i=0;i<m_nOfCol;i++) 
	{
           try
           {
		m_fArray[i]=new  float[m_nOfElements]; 
           }
           catch(std::bad_alloc &e)
           {
		m_fArray[i]=NULL; //in caso di errore
           }

	   if(m_fArray[i]==NULL) 
		{	
			goodAllocation=false;
			for(unsigned int j=0;j<i;j++) 
				delete [] m_fArray[j]; 
			if(m_nOfElements==MIN_NUMBER_OF_ROW) 
			{ 
				delete [] m_fArray;
				m_fArray=NULL;
				return false;
			}
			m_nOfElements=m_nOfElements-MAX_NUMBER_TO_REDUCE_ROW; 
			if(m_nOfElements<=MAX_NUMBER_TO_REDUCE_ROW) 
				m_nOfElements=MIN_NUMBER_OF_ROW;
				break;
		}

	}

  }
	return true;

}
//********************************************************************
bool VSPolarOp::execute()
{
 if( !isParameterPresent("field") ) 
  {
    std::cerr <<"No field given"<<std::endl;
    return false;
   }
 std::stringstream ssField(getParameterAsString("field")); 
 std::string Col[3],outName; 
 unsigned int i, colId[3],colPutId[3]; 
 for (i=0;i<3;i++)
  {
  if (ssField.eof()) 
    {
     std::cerr <<"Invalid field parameters"<<std::endl;
     return false;
    }
  ssField>>Col[i]; 
  colId[i]=m_tables[0]->getColId(Col[i]); 
  if (colId[i]==-1) 
    {
     std::cerr <<"Wrong field "<<Col[i]<<std::endl;
     return false;
    }
  }
  m_nOfRow=m_tables[0]->getNumberOfRows(); 
  if(m_nOfRow>getMaxNumberInt()) 
    {
     m_nOfElements=getMaxNumberInt(); 
    }
   else
    {
     m_nOfElements=m_nOfRow;
    }
  if(!allocateArray()) 
    {
     std::cerr <<"Invalid allocation"<<std::endl;
     return false;
    }

   if( !isParameterPresent("outcol") ) 
     {
      std::cerr <<"No output columns specified, default names = rho, theta, phi"<<std::endl;
      Col[0]="rho"; // se non ci sono metto valori default
      Col[1]="theta";
      Col[2]="phi";
     }
    else
     {
      std::stringstream ssCols(getParameterAsString("outcol")); 
      for (i=0;i<3;i++)
       {
       if (ssCols.eof())
        {
         std::cerr <<"Invalid number of output column names parameters"<<std::endl;
         return false;
        }
       ssCols>>Col[i];
      if (Col[i]=="unknown")
        {
         std::cerr <<"Invalid columns name, default name are rho, theta, phi"<<std::endl;
         switch(i)
          {
           case 0: { Col[i]="rho"; break; }
           case 1: { Col[i]="theta"; break; }
           case 2: Col[i]="phi";
          }
        }
       }
     }
    VSTable newTable; 

    if(isParameterPresent("append"))
      {
      	outName=m_tables[0]->getLocator(); 
        std::clog<<m_tables[0]->getNumberOfColumns()<<std::endl;
      	for (i=0;i<3;i++)
      	{
         	m_tables[0]->addCol(Col[i]);
        std::clog<<m_tables[0]->getNumberOfColumns()<<std::endl;
         	colPutId[i]=m_tables[0]->getColId(Col[i]); 
      	}
        m_tables[0]->writeHeader();
      }
     else
      {
        if( !isParameterPresent("out") ) 
        {
           std::cerr <<"No output table specified, default name = VSPolarOp.bin"<<std::endl;
           outName="VSPolarOp.bin"; 
        }
        else
        {
           outName=getParameterAsString("out");
           if (outName=="unknown")
           {
              std::cerr <<"Invalid table name, default name = VSPolarOp.bin"<<std::endl;
              outName="VSPolarOp.bin";
           }
        }
	m_realOutFilename.push_back(outName);        
        newTable.setType("float"); 
        newTable.setLocator(outName);
        newTable.setNumberOfRows(m_nOfRow);
        newTable.setIsVolume(m_tables[0]->getIsVolume());
        newTable.setCellNumber(m_tables[0]->getCellNumber()[0], m_tables[0]->getCellNumber()[1], m_tables[0]->getCellNumber()[2]);
        newTable.setCellSize(m_tables[0]->getCellSize()[0], m_tables[0]->getCellSize()[1], m_tables[0]->getCellSize()[2]);
        newTable.setEndiannes(m_tables[0]->getEndiannes());
        for (i=0;i<3;i++)
        {
           newTable.addCol(Col[i]); 
           colPutId[i]=newTable.getColId(Col[i]); 
        }
        newTable.writeHeader();
      }

  unsigned long long int index=0;
  
  while(index<m_nOfRow) 
   {
   m_tables[0]->getColumn(colId,m_nOfColIn,index,index+m_nOfElements-1,m_fArray); 
   // theta=zeta copio la colonna 3 nella colonna 4 
   for (unsigned long long int j=0;j<m_nOfElements;j++){
     m_fArray[3][j]=m_fArray[2][j];
     }
   // phi=ypsilon coulmn copy 1 on 4
   for (unsigned long long int j=0;j<m_nOfElements;j++){
     m_fArray[4][j]=m_fArray[1][j];
     }
  // phi=arctan (phi/x) 
   for (unsigned long long int j=0;j<m_nOfElements;j++){
     m_fArray[4][j]=atan(m_fArray[4][j]/m_fArray[0][j]);
     }
   // rho calculus 
   for (unsigned long long int k=0;k<m_nOfColIn;k++){
    for (unsigned long long int j=0;j<m_nOfElements;j++){
     m_fArray[k][j]=m_fArray[k][j]*m_fArray[k][j]; 
     }
    }
   for (unsigned long long int k=1;k<m_nOfColIn;k++){
    for (unsigned long long int j=0;j<m_nOfElements;j++){
     m_fArray[0][j]=m_fArray[0][j]+m_fArray[k][j]; 
     }
    }
   //rho on column 0
   for (unsigned long long int j=0;j<m_nOfElements;j++){
     m_fArray[0][j]=sqrt(m_fArray[0][j]);
     }   

   // theta=arccos (theta/x)    
   for (unsigned long long int j=0;j<m_nOfElements;j++){
     m_fArray[3][j]=acos(m_fArray[3][j]/m_fArray[0][j]);
     }
  // theta on column 1
   for (unsigned long long int j=0;j<m_nOfElements;j++){
     m_fArray[1][j]=m_fArray[3][j];
     }

  // phi on column 2

   for (unsigned long long int j=0;j<m_nOfElements;j++){
     m_fArray[2][j]=m_fArray[4][j];
     }

   if(isParameterPresent("append"))
     {
     	m_tables[0]->putColumn(colPutId,3,index,index+m_nOfElements-1,m_fArray); 
     }
    else 
      {
       newTable.putColumn(colPutId,3,index,index+m_nOfElements-1,m_fArray); 
      }
  index+=m_nOfElements; 
  if(index+m_nOfElements>=m_nOfRow) m_nOfElements=m_nOfRow-index; 
   }   
 return true;
}
//********************************************************************
void VSPolarOp::printHelp()
{
 std::cout<<"It creates three new fields in a data table as the result of the shperical polar transformation of three existing fields."<<std::endl<<std::endl;
 std::cout<<"Usage: VisIVOFilters --op cartesian2polar --field X Y Z [--append] [--outcol rho theta phi] [--out filename_out.bin] [--history] [--historyfile filename.xml] --file input"<<std::endl;

 std::cout<<"Example: VisIVOFilters --op cartesian2polar --field X Y Z --append --file inputFile.bin"<<std::endl;

 std::cout<<"Note: Three new fields are produced."<<std::endl;
 std::cout<<"      --field. Three cartesian coordinates  to compute the polar coordinates vector."<<std::endl;
 std::cout<<"      --outcol. Three column names of the new fields. Default names are given."<<std::endl;
 std::cout<<"      --append. No new table will be created. The original table will have  new fields "<<std::endl;
 std::cout<<"      --out. Name of the new table. Ignored if -- append is specified. Default nam is given"<<std::endl;
 std::cout<<"      --history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
 std::cout<<"      --historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;
 std::cout<<"      --file. Input file name"<<std::endl;
 std::cout<<"      --help produce this output "<<std::endl;
return;
}
//********************************************************************

