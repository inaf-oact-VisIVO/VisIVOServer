/***************************************************************************
 *   Copyright (C) 2011 U. Becciani
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

#ifndef PARFILE_H
#define PARFILE_H
#include <cstdlib>
#include <ctype.h>
#include <map>
#include <sstream>
#include <string>
#include <vector>

class Parfile
{ 
      //store downloaded (temporary) file names: <remote filename , local filename> 
      std::map<std::string, std::string> m_pars; 
      std::map<std::string, std::string>::iterator m_iterpars; 
      bool m_validFile;
      
      void cleanEnd(std::string &str)
      {
 	   std::string noTab="";
	   size_t sizeStr=str.size();
	  std::string::iterator it=str.end()-1;
	  size_t indStr;
	  for (indStr=sizeStr-1 ; indStr >=0; indStr-- )
	  {
//	    if ( !isblank ( *it ) )  break;
	    std::stringstream sstmp;
	    sstmp<<*it;
	    if (!(sstmp.str()==" " || sstmp.str()=="\t")) break;
	    it--;
	  }
	  
	  if(indStr>=0) 
	    noTab=str.substr(0,indStr+1);
	  
         str=noTab;
	 return;
      };
      
  public:

    Parfile(std::string filename)
    { 
      std::ifstream inputParFile;
      inputParFile.open(filename.c_str());
      if(inputParFile.is_open())
      {
	m_validFile=true;
	while(!inputParFile.eof())
	{
	  std::string line,key,value;
	  std::size_t found;
  	  getline(inputParFile, line);
	  found=line.find("#");
	  if(found==0)
	    continue;
	    
	  found=line.find("=");
	  if(found!=std::string::npos)
	  {
	    key=line.substr(0,found);
	    value=line.substr(found+1,line.size());
//	    cleanEnd(value);
	    m_pars.insert(make_pair(key,value));
	  }
	  
	}
      } else
	m_validFile=false;
      return;};
      
      
    template<typename T> T find
      (const std::string &key, const T &deflt)
      {
	T result;
	m_iterpars=m_pars.find(key);
	if(m_iterpars!=m_pars.end())
	{
	  std::stringstream sstmp;
	  sstmp<<m_iterpars->second;
	  sstmp>>result;
	}else
	   result=deflt;
	return result;
      };
    template<typename T> T find
      (const std::string &key, const std::string &deflt)
      {
	T result;
	m_iterpars=m_pars.find(key);
	if(m_iterpars!=m_pars.end()){
	  result=m_iterpars->second;
	  cleanEnd(result);
	 }
	else
	   result=deflt;
	return result;
      };
    
    bool param_present(std::string par)
    {
      if(m_pars.find(par)!=m_pars.end()) return true;
      return false;
    };
    bool valid(){return m_validFile;};
};

#endif 
