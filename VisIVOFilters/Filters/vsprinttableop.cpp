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
#include <cstdlib>
#include <cstring>
#include "vsprinttableop.h"
#include "vstable.h"

#include <iostream>

VSPrintTableOp::VSPrintTableOp()
{
}


VSPrintTableOp::~VSPrintTableOp()
{
}


bool VSPrintTableOp::execute()
{
  std::map<std::string, std::string> parameters = getParameters();
  
  std::map<std::string, std::string>::iterator p;
  std::map<std::string, std::string>::iterator end = parameters.end();

  std::cout << "Printing operation parameters..." << std::endl;

  for(p = parameters.begin(); p != end; p++)
    std::cout << p->first << ": " << p->second << std::endl;

//   unsigned int size = m_tables.size();
// 
//   for(int i = 0; i < size; ++i)
//     m_tables[i]->printSelf();
// 
//   VSTable t3;
//   t3.setName("lillo_test");
//   t3.setLocator("/home/marlin/Scrivania/prova2.bin");
//   t3.setDescription("desc_prova");
//   t3.setNumberOfRows(4);
//   t3.setType("double");
//   t3.setEndiannes("big");
//   t3.addCol("lillo1");
//   t3.addCol("lillo2");
//   t3.writeHeader();

  return true;
}

