/***************************************************************************
 *   Copyright (C) 2011 U. Becciani
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

#ifndef GLITEGW_H
#define GLITEGW_H
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

class gLiteGw
{ 
      //store downloaded (temporary) file names: <remote filename , local filename> 
      std::map<std::string, std::string> m_lfnDownloadedFile; 

  public:
//    gLiteGw ( ) {};
//    ~gLiteGw ( ){};
    std::string readVF (std::map<std::string, std::string> arguments,int idOp=-1);
    bool writeVF (std::map<std::string, std::string> arguments,int idOp=-1);
    void rmv();

};

#endif // GLITEGW_H
