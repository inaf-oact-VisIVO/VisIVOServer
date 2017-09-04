/***************************************************************************
 *   Copyright (C) 2008 by Marco Comparato   *
 *   marco.comparato@oact.inaf.it   *
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
#include "parametersparser.h"
#include <cstdlib>
#include <cstring>

#include <sstream>
#include <iostream>
#include <fstream>



ParametersParser::ParametersParser(std::string parameters,int method)
{
if(method==0)
{
  parameters = trim(parameters);
  
  std::string dummy = "";
  std::string key   = "";
  std::string value = "";
  
  std::stringstream ss(parameters);
  std::streampos pos;

  while(!ss.eof())
  {
    ss >> dummy;
//    std::clog<<dummy<<std::endl;

    if(dummy[0] != '-')
      break; //return false;

    dummy.erase(0, 1);

    if(dummy[0] == '-')
      dummy.erase(0, 1);

    key = dummy;
//    std::clog<<key<<std::endl;

    if(ss.eof())
    {
      value = "unknown";

      m_parameters.insert(make_pair(key, value));
//     std::clog<<value<<std::endl;
     
      break;
    }

    pos = ss.tellg();
    
    ss >> dummy;
//    std::clog<<dummy<<std::endl;

    if(dummy.find("<<")!=std::string::npos)
    {
	 bool endExpression=false;
	std::stringstream sstmp;
	dummy.erase(dummy.find("<<"), 2);
//      std::clog<<dummy<<std::endl;

	while(!endExpression)
	{
	  if (dummy.find(">>")!=std::string::npos)
	  {
		dummy.erase(dummy.find(">>"),2);
//    		std::clog<<dummy<<std::endl;

		endExpression=true;
	  }
	  sstmp<<dummy;
	  if(!endExpression) ss>> dummy;
//	    std::clog<<dummy<<std::endl;

	}
	value=sstmp.str();      
//        std::clog<<key<<" "<<value<<std::endl;

	m_parameters.insert(make_pair(key, value));
	continue; // next while cycle.	
    }

    if(dummy[0] == '-')
    {
      if(!(isdigit(dummy[1])))
      {
        if(ss.eof())
          ss.clear();

        ss.seekg(pos);

        value = "unknown";

        m_parameters.insert(make_pair(key, value));

        continue;
      }
    }

    value = dummy;
//    std::clog<<value<<std::endl;
    
    while(!ss.eof())
    {
      pos = ss.tellg();
       
      ss >> dummy;
//    std::clog<<dummy<<std::endl;
       
      if(dummy[0] == '-')
      {
        if(!(isdigit(dummy[1])))
        {
          if(ss.eof())
            ss.clear();

          ss.seekg(pos);

          break;
        }
      }

      value += " ";
      value += dummy;
    }
//    std::clog<<key<<" "<<value<<std::endl;
    m_parameters.insert(make_pair(key, value));
  }
} else if(method==1){


       std::ifstream parameterFile(parameters.c_str());
       if(!parameterFile.is_open())
       {
	   std::cerr<<"Invalid Filter parameters filename"<<std::endl;
	   exit(1);
       } 
     
      while(!parameterFile.eof())  //!read file content
      {
	std::string key,value;
        std::string tmp = "";
        getline(parameterFile, tmp); //!read row
	if(tmp.compare(0,1,"#")==0) continue; //!comment line
 	tmp = trim(tmp);
	if(tmp.size()==0) continue; //!ignore blank line
	if(tmp.find("=")==std::string::npos)
	{
		std::cerr<<"Invalid row in parameter file: "<<tmp<<std::endl;
		exit(1);
	}
	key=tmp.substr(0,tmp.find("="));
	key=trim(key);
	value=tmp.substr(tmp.find("=")+1,tmp.size()-1);
	value=trim(value);
	m_parameters.insert(make_pair(key, value));

      }
      parameterFile.close();	  
} else 
	std::cerr<<"Invalid parameters method "<<method<<std::endl;


}


ParametersParser::~ParametersParser()
{
}


void ParametersParser::printSelf()
{
  std::map<std::string, std::string>::iterator p;
  std::map<std::string, std::string>::iterator end = m_parameters.end();

  std::cout << "Printing parameters..." << std::endl;

  for(p = m_parameters.begin(); p != end; p++)
    std::cout << p->first << ": " << p->second << std::endl;

  return;
}

//---------------------------------------------------------------------
std::string ParametersParser::trimRight(const std::string & source, const std::string & t /*= " "*/)
//---------------------------------------------------------------------
{
  std::string str = source;
  return str.erase(str.find_last_not_of(t) + 1);
}

//---------------------------------------------------------------------
std::string ParametersParser::trimLeft(const std::string & source, const std::string & t /*= " "*/)
//---------------------------------------------------------------------
{
  std::string str = source;
  return str.erase(0, str.find_first_not_of(t));
}

//---------------------------------------------------------------------
std::string ParametersParser::trim(const std::string & source, const std::string & t /*= " "*/)
//---------------------------------------------------------------------
{
  std::string str = source;
  return trimLeft(trimRight(str, t), t);
}
