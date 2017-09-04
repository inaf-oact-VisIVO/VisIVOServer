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
#ifndef VSTABLEOP_H
#define VSTABLEOP_H

/**
	@author Marco Comparato <marco.comparato@oact.inaf.it>
*/

#include "vsobject.h"

#include <string>
#include <vector>
#include <map>

class VSTable;

class VSTableOp : public VSObject{

  static const unsigned int MAX_NUMBER_INT; //! maximum number of integer

  std::map<std::string, std::string> m_parameters;//! parameters for the operations: key-value
  int m_maxNumberInt; //! maximum number of integer
  protected:
    std::vector<VSTable *> m_tables; //!table pointers used in the operation
    std::vector<std::string> m_realOutFilename; 
    
  public:
    VSTableOp();
    ~VSTableOp();

    //bool setParameters(std::string parameters);
    bool setParameters(std::map<std::string, std::string> parameters) {m_parameters = parameters; return true;}; //! set parameters for the operation 
    bool addParameter(std::string key,std::string value); //! add an element to m_parametrs but  clean "-" signs if any.
    bool addInput(VSTable *table); //! add a table pointer to m_tables

     int isParameterPresent(std::string parameter) {return m_parameters.count(parameter);};//! return the number of the parameter or zero if not present
    std::vector<std::string> realOutFilename() {return m_realOutFilename;};//! return the number of the parameter or zero if not present
    
    std::string getParameterAsString(std::string parameter); //!return parameter value as a string from a given key
    int getParameterAsInt(std::string parameter);//!return parameter value as an int from a given key
    int getMaxNumberInt(); 
    float getParameterAsFloat(std::string parameter);//!return parameter value as a float from a given key

    std::map<std::string, std::string> getParameters() {return m_parameters;};//!return the map (all options)
    
    virtual bool execute() = 0; 
    virtual void printHelp() = 0;
};

#endif
