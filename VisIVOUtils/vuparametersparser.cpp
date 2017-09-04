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
#include "vuparametersparser.h"

#include <sstream>
#include <iostream>



VUParametersParser::VUParametersParser(std::string parameters)
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

    if(dummy[0] != '-')
      break; //return false;

    dummy.erase(0, 1);

    if(dummy[0] == '-')
      dummy.erase(0, 1);

    key = dummy;

    if(ss.eof())
    {
      value = "unknown";

      m_parameters.insert(make_pair(key, value));
      
      break;
    }

    pos = ss.tellg();
    
    ss >> dummy;

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
    
    while(!ss.eof())
    {
      pos = ss.tellg();
       
      ss >> dummy;
       
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

    m_parameters.insert(make_pair(key, value));
  }
}


VUParametersParser::~VUParametersParser()
{
}


void VUParametersParser::printSelf()
{
  std::map<std::string, std::string>::iterator p;
  std::map<std::string, std::string>::iterator end = m_parameters.end();

  std::cout << "Printing parameters..." << std::endl;

  for(p = m_parameters.begin(); p != end; p++)
    std::cout << p->first << ": " << p->second << std::endl;

  return;
}

//---------------------------------------------------------------------
std::string VUParametersParser::trimRight(const std::string & source, const std::string & t /*= " "*/)
//---------------------------------------------------------------------
{
  std::string str = source;
  return str.erase(str.find_last_not_of(t) + 1);
}

//---------------------------------------------------------------------
std::string VUParametersParser::trimLeft(const std::string & source, const std::string & t /*= " "*/)
//---------------------------------------------------------------------
{
  std::string str = source;
  return str.erase(0, str.find_first_not_of(t));
}

//---------------------------------------------------------------------
std::string VUParametersParser::trim(const std::string & source, const std::string & t /*= " "*/)
//---------------------------------------------------------------------
{
  std::string str = source;
  return trimLeft(trimRight(str, t), t);
}
