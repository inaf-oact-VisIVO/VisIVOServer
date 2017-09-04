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
#ifndef VSPOINTPROPERTYOP_H
#define VSPOINTPROPERTYOP_H
#include "vstableop.h"

/**
	@author Ugo Becciani <ugo.becciani@oact.inaf.it>
*/  
class VSPointPropertyOp: public VSTableOp
{
  static const unsigned int MAX_NUMBER_TO_REDUCE_ROW;
  static const unsigned int MIN_NUMBER_OF_ROW;
  unsigned long long int m_numNewPts; // numer of points grid input value and adjusted after allocation

  unsigned int m_nOfCol;
  unsigned int m_nOfRow; // numer of rows grid input value and adjusted  after allocation

  float **m_fArray;
  float **m_fresult;
  float **m_grid;
  bool allocateArray();
  float m_origin[3];
  float m_spacing[3];
public:
    VSPointPropertyOp();

    ~VSPointPropertyOp();
    void printHelp();
    bool execute();

};

#endif
