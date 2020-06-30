/***************************************************************************
 *   Copyright (C) 2014 by Fabio Vitello *
 *  fabio.vitello@oact.inaf.it *
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


#include "vialacteasource.h"

#include "visivoutils.h"
#include "base64.h"


#include <iostream>
#include <fstream>
#include <sstream>

//---------------------------------------------------------------------
int VialacteaSource::readData()
//---------------------------------------------------------------------
{
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //	in case of missing data we put 0 as default value
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    int i = 0;
    int j = 0;
    
    
    std::string::size_type index = std::string::npos;
    std::string fields;
    
    std::vector<std::string> lineData; //!will contain each row of file
    
    std::ifstream inFile;
    
    unsigned long long int  sum=0;
    std::ios::off_type pos=0;
    
    
    //std::ofstream outfile(m_pointsBinaryName.c_str(),std::ofstream::binary );
    std::ofstream outfile(m_pointsBinaryName.c_str() );
    //std::clog<<m_pointsBinaryName.c_str()<<endl;
    // Read Data
    
    if(!outfile)
        return 1;
    
    inFile.open(m_pointsFileName.c_str());
    
    if(!inFile)
        return 1;
    
    std::string tmp = "";
    bool discardLine=true;
    
    getline(inFile, tmp); //!read first line
    tmp = trim(tmp);
    
    int indexp = tmp.find('#');
    
    while(indexp == 0)
    {
        getline(inFile, tmp);
        indexp = tmp.find('#');
    }
    
    findAndReplace(tmp, '\t', ' ');
    findAndReplace(tmp, '#', ' ');
    tmp = trim(tmp);
    
    
    if(tmp.compare(""))
        fields=tmp;
    
    std::string data;
    std::stringstream ss;
    ss << fields;
    int suffix=1;
    while(!ss.eof()) //!fill m_fieldNames: name of columns
    {
        ss >> data;
        
        if(data.compare(""))
        {
            for(int k=0;k< m_fieldNames.size();k++)
            {
                if(data==m_fieldNames[k])
                {
                    std::cerr<<"Warning. Duplicate column name "<<data;
                    std::stringstream ksst;
                    ksst<<suffix;
                    data=data+"_"+ksst.str();
                    suffix++;
                    std::cerr<<" is changed with "<<data<<std::endl;
                }
            }
            m_fieldNames.push_back(data.c_str());
            data = "";
        }
    }
    
    m_nCols=m_fieldNames.size();
    
    //   std::clog<<"m_nRows="<< m_nRows<<endl;
    //   std::clog<<"m_nCols="<< m_nCols<<endl;
    
    int nLoad=(int)(10000000/m_nCols);  //!nLoad is the number of rows that I will read
    
    std::string **matrix=NULL;
    
    try
    {
        matrix = new std::string*[m_nCols];
    }
    catch (std::bad_alloc e)
    {
        return 1;
    }
    
    for(i = 0; i < m_nCols; i++)
        
    {
        try
        {
            matrix[i] = new std::string[m_nRows*2];
        }
        catch (std::bad_alloc e)
        {
            return 1;
        }
    }
    
    while(!inFile.eof())  //!read file content
    {
        tmp = "";
        
        getline(inFile, tmp); //!read row
        
        //    std::clog<<tmp<<std::endl;
        
        findAndReplace(tmp, '\t', ' ');
        tmp = trim(tmp);
        
        indexp = tmp.find('#');
        
        if(indexp == 0)
            continue;
        
        if(tmp.compare(""))
            lineData.push_back(tmp);  //!fill lineData with the row
        
        
        
        if( inFile.eof())  //!arrived to maximum allowed
        {
            //std::clog<<lineData.size()<<endl;
            
            for(i = 0; i <lineData.size(); i++)  //!for each vector element
            {
                std::stringstream ss;
                ss << lineData[i];  //!line extraction
                
                
                for(j = 0; j < m_nCols; j++)
                {
                    if(ss.eof())
                    {
                        matrix[j][i]=MISSING_VALUE;
                        continue;
                    }
                    std::string data;
                    ss>>data;
                    //	  std::clog<<data<<std::endl;

                    char *cstr;
                    cstr = new char [data.size()+1];
                    strcpy (cstr, data.c_str());

                    
                    
                    //matrix[j][i]=base64_encode(reinterpret_cast<const unsigned char*>(data.c_str()), data.length());
                   // matrix[j][i]=ss1.str();
                    matrix[j][i]= data.c_str();

                   // std::cout<< "i: "<<i<<" - "<<matrix[j][i].c_str()<<std::endl;

                
                    if(data.length()==0)
                        {
                            std::cerr<<"WARNIG. Missing value "<<data<<" is found. Replaced with VALUE="<<MISSING_VALUE <<std::endl;
                            matrix[j][i]=MISSING_VALUE;
                        }
                }
                
            }
            
            //std::cout<<"linedatasize: "<<lineData.size()<<" col: "<<m_nCols<<" size: "<<sizeof(std::string)<<" m_nRows: "<<m_nRows;
            
           for (j=0;j<m_nCols;j++) //!write matrix on file
            {
                for (int i=0; i<m_nRows; i++)
                {
                    
                   // outfile.write(reinterpret_cast<const char *>(&len), sizeof(size_t));
                    
                    outfile<<matrix[j][i].c_str()<<std::endl;
                    
                   // outfile.write((matrix[j][i+1].c_str()), matrix[j][i+1].length());
                   //                        std::cout<< "i: "<<i<<" - "<<matrix[j][i].c_str()<<std::endl;

                }
            }
    
            sum=sum+lineData.size();
            //std::clog <<"sum="<<sum<<endl;
            lineData.clear();
            // pos=outfile.tellp();
            //std::clog<<"posfinal="<<pos<<endl;
        }
    }
    
    
    
    
    if(matrix)  //! delete matrix
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
    
    makeHeader(sum,m_pointsBinaryName,m_fieldNames,m_cellSize,m_cellComp,m_volumeOrTable); //!create header file (in utility)
    
    return 1;
}

//---------------------------------------------------------------------
int VialacteaSource::readHeader()
//---------------------------------------------------------------------
{
    
    m_fieldNames.clear();
    
    int i = 0;
    int rows=-1;
    
    std::vector<std::string> counting;
    std::string::size_type index = std::string::npos;
    
    std::ifstream inFile;
    inFile.open(m_pointsFileName.c_str());  if(!inFile) return 1;
    
    std::string tmp = "";
    
    while(!inFile.eof())  //! count the number of row: rows variable is the toltal numer (excluding the first one
    {
        getline(inFile, tmp);
        
        index = tmp.find('#');
        
        if(index != std::string::npos)
        {
            tmp.erase(index);
            tmp = trimRight(tmp);
        }
        
        if(tmp.compare(""))
        {
            
            counting.push_back(tmp);
            ++rows;
            counting.clear();
        }
    }
    m_nRows=rows;
    //std::clog<<"rowsheader="<<m_nRows<<endl;
    //!tmp now is the line containing the field names
    
    inFile.close();
    
    if(rows==0)
        return 1;
    
    return 0;
}

