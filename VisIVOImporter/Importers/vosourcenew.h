/***************************************************************************
 *   Copyright (C) 2008 by Gabriella Caniglia *
 *  gabriella.caniglia@oact.inaf.it *
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

#ifndef VOSOURCENEW_h
#define VOSOURCENEW_h

#include "abstractsource.h"

#include <string>
#include <vector>
#include <fstream>

class VOSourcenew : public AbstractSource
{
  static const unsigned int MAX_NUMBER_INT; //! maximum number of integer
  static const unsigned int MAX_NUMBER_TO_REDUCE_ROW;
  static const unsigned int MIN_NUMBER_OF_ROW;
	std::ifstream m_inputFile;
	std::ofstream m_outfile;
	std::vector<std::string> m_fieldName;
	float **m_fArray;
	char *m_buffer;
	bool allocateArray();
	bool writeData();
	unsigned long long int m_nOfRow;
	unsigned long long int m_alreadyWritten;
	int m_nOfEle;
	int m_toWrite;
	int m_nOfCol;
	int m_length;
	bool m_notgood;

  public:
    int readHeader();
    int readData();
    VOSourcenew();
    ~VOSourcenew();


};

#endif
