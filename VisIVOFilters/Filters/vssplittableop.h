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
#ifndef __VSSPLITTABLE_H
#define __VSSPLITTABLE_H

#include "vstableop.h"


class VSSplitTableOp : public VSTableOp
{

static const unsigned int MAX_NUMBER_TO_REDUCE_ROW;
static const unsigned int MIN_NUMBER_OF_ROW;


  float **m_fArray;

  unsigned int m_nOfCols;
  unsigned long long int m_nOfRows, m_nOfLocalRow;
  int m_nOfEle;


  bool allocatefArray();

  bool splitVolume();
  bool splitPoints();
  bool fastTableSplit();
  std::string m_field;
  int m_vDir;
  unsigned int m_fieldId;
   int m_numOfTables;
 float m_lowLimit;
float m_upLimit;
 int m_tableNumber;
std::string m_fileNameOutput;

  bool writeTable();

  public:
    VSSplitTableOp();
    ~VSSplitTableOp();
    void printHelp();
    bool execute();
};

#endif
