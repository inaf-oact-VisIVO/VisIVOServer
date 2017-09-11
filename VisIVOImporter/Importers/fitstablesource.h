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

#include  "abstractsource.h"

#include <vector>
extern "C" {
#include "fitsio.h"
}


class FitsTableSource : public AbstractSource 
{
  public:
    int readHeader();
    int readData();
    
  private:
    long GetNumRows(int ntab);
    int GetNumColumns(int ntab);
    int GetHDUType();
    void SetScalFields(int col);
  
    std::string GetColName(int ncol);
    std::string GetColType(int ncol, int *typecode, long *repeat);
    std::string GetColUnit(int ncol);
    std::string GetColFormat(int ncol);
    fitsfile *pFile;

  protected:
    int Ntable;
    int hdunum;
    int hdutype;
    int casesen;
    int status; 
  
    std::vector<std::string> m_fields;
    std::vector <int> scalFields;
	
};


