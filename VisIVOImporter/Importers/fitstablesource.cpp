/***************************************************************************
 *   Copyright (C) 2008 by Gabriella Caniglia *
 *  UBE 03/10/11 *
 *  gabriella.caniglia@oact.inaf.it *
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

#include "fitstablesource.h"

#include "stdio.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "visivoutils.h"

extern "C" {
#include "fitsio.h"
}

fitsfile *pFile;

//----------------------------------------------------------------------------
int FitsTableSource::readHeader()
//----------------------------------------------------------------------------

{
  return 0;
}


//----------------------------------------------------------------------------
int FitsTableSource::readData()
//----------------------------------------------------------------------------
//!CHECK ENDINAISM PLEASE!!!!
//fits_read_colnull(pFile, typecode, i+1 , firstrow,firstelem, nLoad, readColData, nullArray, &anynul, &status);
// pFile Fits table pointer; typecode Dtata Type in the Column; i+1: number of Column (1= first column), 
//firstrow= row to be read in the column, firstelem First element in the array where data are stored
// (readColData), nLoad number of lement to be read, readColData, array where data are stored, nullArry: missing 
//element=1
{ 
  std::ofstream outfile(m_pointsBinaryName.c_str(),std::ofstream::binary );  //!open out binary file
  
 
  int i = 0;
  int j = 0;
 
  long nRows;
  char *val = 0;
  char value[1000];
  value[0] = '\0';
  char nullstr[1] = {'\0'};
  int anynul = 1;
  val = value;
  int hdunum = 0;
  int hdunum1 = 0;
  int hdutype = 0;
  int typecode;
  unsigned long long int res=0;
      
  long firstrow = 1;
  long firstelem = 1;
  long nelements;
  float *DataArrayScal = 0;
  char *nullArray = 0;
  int err_ret=0;  
  status=0;
  err_ret=fits_open_file(&pFile, m_pointsFileName.c_str(), READONLY, &status); //!open FITS input file
  if(err_ret!=0)
  {
	std::cerr<<"Invalid open of fitstable file. Importer aborted"<<std::endl;
	return -1;
  } 
  fits_get_num_hdus(pFile, &hdunum, &status);
  fits_get_hdu_num(pFile,  &hdunum1); 
  fits_get_hdu_type(pFile, &hdutype, &status);
  fits_movabs_hdu(pFile, 2, &hdutype, &status);
  fits_get_hdu_num(pFile,  &hdunum1); 
  fits_get_num_rows(pFile, &nRows, &status);
  fits_get_num_cols(pFile, &m_nCols, &status);
  status = 0;
    
  m_nRows=int(nRows);
 int nLoad=MAX_LOAD;
  
     if(m_nRows<=nLoad)
     nLoad=m_nRows;

        //!get col names into vecotr m_fields
  for(int is = 1; is <= m_nCols; is++) 
  {
    char colname[70];
    char templ[2] = {'*', '\0'};
    fits_get_colname(pFile, casesen, templ, colname, &is, &status);
    std::string nomeColonna = colname;
    m_fieldNames.push_back(nomeColonna);
  }
  try
  {
    DataArrayScal = new float[nLoad];
  }
  catch (std::bad_alloc & e)
  {
    return 1;
  }
  try
  {  
    nullArray = new char[nLoad];
  }
  catch (std::bad_alloc & e )
  {
    return 1;
  }
  
  for(i = 0; i < nLoad; i++)
    nullArray[i] = 0;
        
  nelements = (long)m_nRows;
  status = 0;

//  int nSelScalFields = scalFields.size();
  int nSelScalFields = 0;
  int numeroScal = -1;
  std::string unit = "";
  std::string type = "";
  std::string format = "";

  char **sData = 0;
  long pos=0;
    
    
  for(i = 0; i < m_nCols; i++)
  {
     firstrow=1;  
     int nLoad=MAX_LOAD;

     if(m_nRows<=nLoad)
        nLoad=m_nRows;
    
     res=m_nRows;

//    int nLoad/*=(int)(10000000/m_nCols)*/;
     while(res!=0)
     { 
         if(nLoad>res)
     	 	nLoad=res; 


           type = GetColType(i,&typecode);
           unit = GetColUnit(i);

           if(!(iCompare(type, "Type: String"))) //!string file reading
           {
             int stringSize = 0;

             format = GetColFormat(i);

             int localHDUType = GetHDUType();
             if(localHDUType == 2)
             {
       
               format.erase(format.size() - 1, 1);
             }
             else
             {
               format.erase(0, 1);
             }

             stringSize = atoi(format.c_str()) + 1;

             if(stringSize < 2)
             {
               delete [] sData;
               sData = 0;
          

               break;
             }

             sData = new char*[nLoad];

             for(j = 0; j < nLoad; j++)
               sData[j] = new char[stringSize];

             fits_read_col(pFile, TSTRING, i , firstrow,
                           firstelem, nelements, &nullstr, sData, 
                           &anynul, &status);
				firstelem=firstelem+nLoad;
           }
    
           else
           {
  		switch(typecode)
  		{
    		  case 11:
			{ 
			unsigned char *readColData;
try
  {
    			readColData = new unsigned char[nLoad];
  }
catch (std::bad_alloc & e)
  {
    return 1;
  }
             		fits_read_colnull(pFile, typecode, i+1 , firstrow,firstelem, nLoad, readColData, nullArray, &anynul, &status);
             		for(j = 0; j < nLoad; j++)
				DataArrayScal[j]=readColData[j];
			firstelem=firstelem+nLoad;
			delete [] readColData;
			}
			break;
    		  case 21:
			{
			short *readColData;
try
  {
    			readColData = new short[nLoad];
  }
catch (std::bad_alloc & e)
  {
    return 1;
  }
             		fits_read_colnull(pFile, typecode, i+1 , firstrow,firstelem, nLoad, readColData, nullArray, &anynul, &status);
             		for(j = 0; j < nLoad; j++)
				DataArrayScal[j]=readColData[j];
			firstelem=firstelem+nLoad;
			delete [] readColData;
			}
			break;
    		  case 30:
			{
			unsigned int *readColData;
try
  {
    			readColData = new unsigned int[nLoad];
  }
catch (std::bad_alloc & e)
  {
    return 1;
  }
             		fits_read_colnull(pFile, typecode, i+1 , firstrow,firstelem, nLoad, readColData, nullArray, &anynul, &status);
             		for(j = 0; j < nLoad; j++)
				DataArrayScal[j]=readColData[j];
			firstelem=firstelem+nLoad;
			delete [] readColData;
			}
			break;
    		  case 31:
			{
			int *readColData;
try
  {
    			readColData = new int[nLoad];
  }
catch (std::bad_alloc & e)
  {
    return 1;
  }
             		fits_read_colnull(pFile, typecode, i+1 , firstrow,firstelem, nLoad, readColData, nullArray, &anynul, &status);
             		for(j = 0; j < nLoad; j++)
				DataArrayScal[j]=readColData[j];
			firstelem=firstelem+nLoad;
			delete [] readColData;
			}
			break;
    		  case 40:
			{
			unsigned long int *readColData;
try
  {
    			readColData = new unsigned long int[nLoad];
  }
catch (std::bad_alloc & e)
  {
    return 1;
  }
             		fits_read_colnull(pFile, typecode, i+1 , firstrow,firstelem, nLoad, readColData, nullArray, &anynul, &status);
             		for(j = 0; j < nLoad; j++)
				DataArrayScal[j]=readColData[j];
			firstelem=firstelem+nLoad;
			delete [] readColData;
			}
			break;
    		  case 41:
			{
			long int *readColData;
try
  {
    			readColData = new long int[nLoad];
  }
catch (std::bad_alloc & e)
  {
    return 1;
  }
             		fits_read_colnull(pFile, typecode, i+1 , firstrow,firstelem, nLoad, readColData, nullArray, &anynul, &status);
             		for(j = 0; j < nLoad; j++)
				DataArrayScal[j]=readColData[j];
			firstelem=firstelem+nLoad;
			delete [] readColData;
			}
			break;
    		  case 42:
             		fits_read_colnull(pFile, typecode, i+1 , firstrow,firstelem, nLoad, DataArrayScal, nullArray, &anynul, &status);
			firstelem=firstelem+nLoad;
			break;
    		  case 81:
			{
			long long int *readColData;
try
  {
    			readColData = new long long int[nLoad];
  }
catch (std::bad_alloc & e)
  {
    return 1;
  }
             		fits_read_colnull(pFile, typecode, i+1 , firstrow,firstelem, nLoad, readColData, nullArray, &anynul, &status);
             		for(j = 0; j < nLoad; j++)
				DataArrayScal[j]=readColData[j];
			firstrow=firstrow+nLoad;
			delete [] readColData;
			}
			break;
    		  case 82:
			{
			double *readColData;
try
  {
    			readColData = new double[nLoad];
  }
catch (std::bad_alloc & e)
  {
    return 1;
  }
             		fits_read_colnull(pFile, typecode, i+1 , firstrow,firstelem, nLoad, readColData, nullArray, &anynul, &status);
             		for(j = 0; j < nLoad; j++)
				DataArrayScal[j]=readColData[j];
			firstrow=firstrow+nLoad;
			delete [] readColData;
			}
			break;
    		  default:  
			std::cerr<<"Invalid "<<type<<" in column "<<numeroScal<<std::endl;
             		for(j = 0; j < nLoad; j++)
				DataArrayScal[j]=TEXT_VALUE;
			break;
  		}//switch
          }//else

    //!to substitute null values with a fake value equals to
    //!the mean of the valid values
//            double sum = 0;
//            double nValidValues = 0;
//            double fakeValue = 0;
// 
//            for(j = 0; j < nLoad; j++)
//              if(nullArray[j] == 0)
//            {
//              sum += DataArrayScal[j];
//              nValidValues++;
//            }
// 
//            fakeValue = sum/nValidValues;

//            for(j = 0; j <nLoad; j++)
//              if(nullArray[j] == 1)
//                DataArrayScal[j] = fakeValue;  //!here there is the substitution

           if(!(iCompare(unit, "mas")))
             masToRad(DataArrayScal, nLoad);
           else if(!(iCompare(unit, "\"h:m:s\"")))
             hmsToRad(sData, DataArrayScal, nLoad);
           else if(!(iCompare(unit, "\"d:m:s\"")))
             dmsToRad(sData, DataArrayScal, nLoad);
           else if(!(iCompare(unit, "deg")))
             degToRad(DataArrayScal, nLoad);

     
           for(j = 0; j <nLoad; j++)
             if(nullArray[j] == 1)
               DataArrayScal[j] = MISSING_VALUE;  //!here there is the       

      
           outfile.write((char*)(DataArrayScal), sizeof(float)*nLoad); 
           res=res-nLoad;

           if(sData)
           {
             for(j = 0; j < nLoad; j++)
               delete [] sData[j];

             delete [] sData;
             sData = 0;
           }
    
     } //while(res!=0)
  }



  delete [] DataArrayScal;
  delete [] nullArray;
        //!close fits file
  fits_close_file(pFile, &status);
  outfile.close();
  makeHeader(m_nRows,m_pointsBinaryName,m_fieldNames,m_cellSize,m_cellComp,m_volumeOrTable);
 

  return 0;
}

//----------------------------------------------------------------------------
  long FitsTableSource::GetNumRows(int ntab)
//----------------------------------------------------------------------------
{
  Ntable = ntab;
  fits_movabs_hdu(pFile, Ntable, &hdutype, &status);
	
  long nm_nRows;
  fits_get_num_rows(pFile, &nm_nRows, &status);

  return nm_nRows;
}
//----------------------------------------------------------------------------
  int FitsTableSource::GetNumColumns(int ntab)
//----------------------------------------------------------------------------
{
  Ntable = ntab;
	
	//! Move to the selected table (starting from 1, which is dummy)
	
  fits_movabs_hdu(pFile, Ntable, &hdutype, &status);
	
  int nm_nCols;
  fits_get_num_cols(pFile, &nm_nCols, &status);
	
  return nm_nCols;
}

//----------------------------------------------------------------------------
  int FitsTableSource::GetHDUType()
//----------------------------------------------------------------------------
{
  int type;
  fits_get_hdu_type(pFile, &type, &status);
	
  return type;
}
//----------------------------------------------------------------------------
  std::string FitsTableSource::GetColUnit(int ncol)
//----------------------------------------------------------------------------
{
  std::string fieldunit = "";
  char value[81]; //! temp
  status = 0;

  std::stringstream skey;
  std::string key;
  skey << "TUNIT" << (ncol + 1);
  key = skey.str();
  char *keyword = const_cast<char *>(key.c_str());
	
  fits_read_key(pFile, TSTRING, keyword, value, 0, &status);

  if(status == 0)
  {
    fieldunit = value;
    return fieldunit;
  }

  char err[31];
  fits_get_errstatus(status, err);
  status = 0;
	
  return "";
}

//----------------------------------------------------------------------------
  std::string FitsTableSource::GetColFormat(int ncol)
//----------------------------------------------------------------------------
{
  std::string fieldFormat = "";
  char value[81]; //! temp
  status = 0;

  std::stringstream skey;
  std::string key;
  skey << "TFORM" << (ncol + 1);
  key = skey.str();
  char *keyword = const_cast<char *>(key.c_str());

  fits_read_key(pFile, TSTRING, keyword, value, 0, &status);

  if(status == 0)
  {
    fieldFormat = value;
    return fieldFormat;
  }

  return "";
}

//----------------------------------------------------------------------------
  std::string FitsTableSource::GetColName(int ncol)
//----------------------------------------------------------------------------
{
  if(m_fields.size() > 0)
    return m_fields[ncol];

  char *templ = "*";
  char colname[70];
  fits_get_colname(pFile, casesen, templ, colname, &ncol, &status);
  std::string fieldname = colname;

  return fieldname;
}
//----------------------------------------------------------------------------
  std::string FitsTableSource::GetColType(int ncol, int *typecode)
//----------------------------------------------------------------------------
{
  
  long repeat;
  long width;
  std::string fieldtype = "";
  status = 0;
  fits_get_coltype(pFile, ncol+1, typecode, &repeat, &width, &status);

  switch(*typecode)
  {
    case 1: fieldtype = "Type: Bit"; break;
    case 11: fieldtype = "Type: Byte"; break;
    case 12: fieldtype = "Type: SByte"; break;
    case 14: fieldtype = "Type: Logical"; break;
    case 16: fieldtype = "Type: String"; break;
    case 20: fieldtype = "Type: UShort"; break;
    case 21: fieldtype = "Type: Short"; break;
    case 30: fieldtype = "Type: UInt"; break;
    case 31: fieldtype = "Type: Int"; break;
    case 40: fieldtype = "Type: ULong"; break;
    case 41: fieldtype = "Type: Long"; break;
    case 42: fieldtype = "Type: Float";  break;
    case 81:  fieldtype = "Type: LongLong"; break;
    case 82:  fieldtype = "Type: Double"; break;
    case 83: fieldtype = "Type: Complex"; break;
    case 163: fieldtype = "Type: Double Complex"; break;
    default:  fieldtype = "Type: Unknown";
  }


  return fieldtype;
}

//----------------------------------------------------------------------------
  void FitsTableSource::SetScalFields(int col)
//----------------------------------------------------------------------------
{
  scalFields.clear();

  for(int i = 0; i < col; i++)
    scalFields.push_back(i);
 
  return;
}

