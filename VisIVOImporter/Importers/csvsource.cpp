/***************************************************************************
 *   Copyright (C) 2008 by Gabriella Caniglia *
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

#include "csvsource.h"
#include "visivoutils.h"

#include <string>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>



//---------------------------------------------------------------------
int CSVSource::readData()
//---------------------------------------------------------------------
{
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //      in case of missing data we put 0 as default value
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  int i = 0;
  int j = 0;

  std::string::size_type index = std::string::npos;
  std::string fields;

  std::vector<std::string> lineData;

  std::ifstream inFile;

  unsigned long long int sum=0;
  std::ios::off_type pos=0;

                
  std::ofstream outfile(m_pointsBinaryName.c_str(),std::ofstream::binary ); 
        //std::clog<<m_pointsBinaryName.c_str()<<endl;
        // Read Data

  inFile.open(m_pointsFileName.c_str());

  if(!inFile)
    return 1;

  std::string tmp = "";

  getline(inFile, tmp);

  findAndReplace(tmp, '\t', ' ');
  tmp = trim(tmp);

  int indexp = tmp.find('#');

  if(indexp != std::string::npos)
    tmp.erase(indexp,indexp+1);

  if(tmp.compare(""))
    fields=tmp;

  std::string data;
  std::stringstream ss;
  ss << fields;

  while(!ss.eof())
  {
    getline(ss, data, ',');  //!instruction for comma separated data format 

    if(data.compare(""))
    {
     data=trim(data);
      m_fieldNames.push_back(data.c_str());
      data = "";
    }
  }

  m_nCols=m_fieldNames.size();
//   std::clog<<"m_nRows="<< m_nRows<<endl;
//   std::clog<<"m_nCols="<< m_nCols<<endl;
  int nLoad=(int)(10000000/m_nCols);
 
  float **matrix=NULL;
  try
  {
    matrix = new float*[m_nCols];
  }
  catch (std::bad_alloc e)
  {
    return 1;
  }
 
  for(i = 0; i < m_nCols; i++)
    
  {
    try
    {
      matrix[i] = new float[m_nRows];
    }
    catch (std::bad_alloc e)
    {
      return 1;
    }
  } 

  while(!inFile.eof())
  {
    tmp = "";

    getline(inFile, tmp);

    findAndReplace(tmp, '\t', ' ');
    tmp = trim(tmp);

    indexp = tmp.find('#');

    if(indexp != std::string::npos)
      tmp.erase(indexp,indexp+1);

    if(tmp.compare(""))
      lineData.push_back(tmp);
    
   
    if((lineData.size()==nLoad || inFile.eof()))   
    {
//       std::clog<<lineData.size()<<endl;
            
      for(i = 0; i <lineData.size(); i++)
      {
        std::stringstream ss;
        ss << lineData[i];


        for(j = 0; j < m_nCols; j++)
        {
          getline(ss, data, ','); //!fill data from ss up to the comma
	
	  bool numeric=false;	
	  char *cstr;
	  cstr = new char [data.size()+1];
	  strcpy (cstr, data.c_str());
	  for(unsigned int k = 0; k < data.length(); k++)
          {
		numeric=true;
		if(!(isdigit(cstr[k])))
                {
		   if ((k==0 && data.compare(k,1,"-")!=0)  && (k==0 && data.compare(k,1,"+")!=0 ))
		 	if(data.compare(k,1,".")!=0 )
		 	{ 
				numeric=false;
		 		break;
		 	}
		 if(k>0 && data.compare(k,1,".")!=0 )
		 { 
			numeric=false;
		 	break;
		 }
		}
	   }
	  
          if(numeric) 
		matrix[j][i]=atof(data.c_str());
          else
		if(data.length()==0)
			matrix[j][i]=MISSING_VALUE;
		else
			matrix[j][i]=TEXT_VALUE;
        }

                   //std::cerr << std::endl;
      }

      for (j=0;j<m_nCols;j++) 
      {
        pos=(sum*(sizeof(float)))+((sizeof(float))*((unsigned long long int)j*m_nRows));
        outfile.seekp(pos);
//         std::clog<<"pos="<<pos<<endl;
        outfile.write((char*)(matrix[j]), sizeof(float)*lineData.size()); 
      }

      sum=sum+lineData.size();
//       std::clog <<"sum="<<sum<<endl;
      lineData.clear();

      pos=outfile.tellp();
     // std::clog<<"posfinal="<<pos<<endl;
    }
  }
  
  if(matrix)
  {
    for(int i = 0; i < m_nCols; i++)
    {
      if(matrix[i])
      {
        delete [] matrix[i];
      }
    }
     
    delete [] matrix;

  }
  
  inFile.close();
  outfile.close();
  
  makeHeader(sum,m_pointsBinaryName,m_fieldNames,m_cellSize,m_cellComp,m_volumeOrTable);

  return 1;
}
//---------------------------------------------------------------------
int CSVSource::readHeader()
//---------------------------------------------------------------------
{
  m_fieldNames.clear();

  int i = 0;
  int rows=-1;

  std::string::size_type index = std::string::npos;

  std::ifstream inFile;
  inFile.open(m_pointsFileName.c_str());  if(!inFile) return 1;

  std::string tmp = "";

  while(!inFile.eof())
  {
    getline(inFile, tmp);


    if(tmp.compare(""))
      ++rows;
  }
  m_nRows=rows;

	//std::clog<<"rowsheader="<<m_nRows<<endl;

  inFile.close();
  if(rows==0)
	return 1;

  return 0;
}

