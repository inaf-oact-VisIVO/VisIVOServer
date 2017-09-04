/***************************************************************************
 *   Copyright (C) 2008 by Ugo Becciani *
 *  ugo.becciani@oact.inaf.it *
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
 
#ifndef XMLSOURCE_H
#define XMLSOURCE_H


#include "abstractsource.h"
#include "voFieldParam.h"
#include "VOTableParser.h"
#include <vector>
#include <map>

class voFieldParam;


enum itype {fl, dou, ldou, i, li,lli}; 
class XmlSource : public AbstractSource
{    

static const unsigned int MAX_NUMBER_TO_REDUCE_ROW;
static const unsigned int MIN_NUMBER_OF_ROW;
static const unsigned int MAXINT;

	std::vector<voFieldParam> m_tvoXmlList;  //!descrizione completa campo <FILED />
	std::vector<int> m_customRows;
        VOTableParser m_sPP;
	std::map<std::string, unsigned long long int> m_speciesMapele;
	std::map<std::string, std::string> m_speciesMapgeo;
  	std::map<std::string, unsigned long long int>::iterator m_it;
	std::map<std::string, std::string> m_outfileMap;
	std::map<std::string, std::string> m_listFileMap;
	int m_numIdSpecies;
	int m_numIdOffset;
	int m_numIdHref;
	int m_numIdField;
	int m_numIdPrecision;
	int m_numIdType;
	int m_numIdStride;
	int m_numIdFormat;
	int m_numIdEndianism;
	int m_numIdRank;
	int m_numIdArraysize;
	int m_maxEle;
	float *m_fArray;
	char *m_fReadArray;
	enum itype m_inputType;

	bool downloadFiles(); 
	bool allocateArray();
	bool readSpecifiedFormat();
	void writeListFiles();

  public:
 
  //! Read the header file and set the basic table parameters
    int readHeader();
    int readData();
     ~XmlSource();

};
#endif

