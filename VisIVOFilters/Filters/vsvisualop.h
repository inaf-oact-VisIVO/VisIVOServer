/***************************************************************************
 *   Copyright (C) 2008 by Ugo Becciani   *
 *   ugo.becciani@oact.inaf.it   *
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
#ifndef VSVISUAL_H
#define VSISUAL_H

/**
	@author Ugo Becciani <ugo.becciani@oact.inaf.it>
*/
//#include <map>
//#include <vector>

#include "vstableop.h"

class VSVisualOp : public VSTableOp
{

  std::map<std::string, int> m_listOfTables;
  std::vector<unsigned int> m_listColumns[101];
  std::vector<std::string> m_listColumnsOutName[101];
  static const unsigned int DEFAULT_ROW_TO_VISUALIZE;
  static const unsigned int MAX_NUMBER_TO_REDUCE_ROW;
  static const unsigned int MIN_NUMBER_OF_ROW;

  unsigned int m_visualSize;
  unsigned int m_colSize;
  unsigned long long int m_globalNumberOfRow;
  unsigned int m_maxEle;

  float **m_fVisualArray;
  float **m_fArray;
  bool allocatefArray();

public:
    VSVisualOp();
    ~VSVisualOp();
    void printHelp();
    bool execute();

};

#endif
