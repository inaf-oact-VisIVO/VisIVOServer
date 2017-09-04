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
#ifndef VSEXTRACTSUBVOLUMEOP_H
#define VSEXTRACTSUBVOLUMEOP_H

/**
	@author Ugo Becciani <ugo.becciani@oact.inaf.it>
*/
#include "vstableop.h"



class VSExtractSubVolumeOp: public VSTableOp
{
  static const unsigned int MAX_NUMBER_TO_REDUCE_ROW;
  static const unsigned int MIN_NUMBER_OF_ROW;
  static const unsigned int MAX_NUMBER_OF_BYTES;


  	float **m_fArray;  //input data
	float **m_cArray;  //extract array output data
	unsigned int m_nOfCol;
	unsigned long long int m_nOfRow;
	unsigned long long int  m_nExtractOfRow;

	
   std::vector<unsigned int> m_fieldList;

public:
    VSExtractSubVolumeOp();
    ~VSExtractSubVolumeOp();
    bool allocateArray();

    void printHelp();
    bool execute();

};


#endif
