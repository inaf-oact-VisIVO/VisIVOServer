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
#ifndef __VSCUT_H
#define __VSCUT_H

#include "vstableop.h"


class VSCutOp : public VSTableOp
{

static const unsigned int MAX_NUMBER_TO_REDUCE_ROW;
static const unsigned int MIN_NUMBER_OF_ROW;

	struct limitField{
		unsigned int colId;
		bool downUnlimited, upUnlimited;
		float 	downLimit, upLimit;		
	};

  std::vector<limitField> m_colVector;	

  float **m_fArray;
  float **m_fArrayWrite;

  unsigned int m_nOfCol;
  unsigned int m_nOfRow;
  std::vector<unsigned int> m_colNumberSet;

  bool allocatefArray();

  public:
    VSCutOp();
    ~VSCutOp();
    void printHelp();
    bool execute();
};

#endif
