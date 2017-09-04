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
#ifndef PARAMETERPARSER_H
#define PARAMETERPARSER_H

/**
        @author Marco Comparato <marco.comparato@oact.inaf.it>
 */

#include <map>
#include <string>
//! The constructor read input command line options (method=0) or options from a file (method = 1)
//! All options are "keyvalue" in the m_parameters map with the value given by the user.
class ParametersParser
{
  std::string trimRight(const std::string & source, const std::string & t = " ");
  std::string trimLeft(const std::string & source, const std::string & t = " ");
  std::string trim(const std::string & source, const std::string & t = " ");

  std::map<std::string, std::string> m_parameters;
  std::string m_paramFile;


  public:
    ParametersParser(std::string parameters,int method=0); //! method=0 commandline, method=1 read from file 
    ~ParametersParser();

    std::map<std::string, std::string> getParameters() {return m_parameters;};

    void printSelf();
};

#endif
