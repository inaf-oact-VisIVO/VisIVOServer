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

#ifndef FITSIMAGE5OURCE_H
#define FITIMAGESOURCE_H


#include "abstractsource.h"


extern "C" {
#include "fitsio.h"
}


class FitsImageSource : public AbstractSource 
{
  public:
 
  // Read the header file and set the basic table parameters
    int readHeader();
    int readData();
private:
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
#endif
