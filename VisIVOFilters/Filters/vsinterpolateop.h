/***************************************************************************
 *   Copyright (C) 2009 by Ugo Becciani   *
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
*/

#ifndef VSINTERPOLATEOP_H
#define VSINTERPOLATEOP_H

#include "vstableop.h"


class VSInterpolateOp: public VSTableOp
{  
static const unsigned int MAX_NUMBER_TO_REDUCE_ROW;
static const unsigned int MIN_NUMBER_OF_ROW;

  float **m_f1Array;
  float **m_f2Array;
  float **m_f3Array;
  unsigned int m_nOfCols;
  unsigned int m_nOfRows;

  bool allocateArray();

public:
    VSInterpolateOp();
    ~VSInterpolateOp();
    void printHelp();
    bool execute();

};

#endif
