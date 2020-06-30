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

#ifndef RAMSESSOURCE_h
#define RAMSESSOURCE_h

#include "abstractsource.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

#define QWORD_PRECISION 16
#define DOUBLE_PRECISION 8
#define SINGLE_PRECISION 4


// File header
struct ramsesHeader{
	int 	ncpu;
	int 	ndim;
	int 	npart;
	// If local seed is fortran real(16) store in two parts
	double 	localseed1;
	double 	localseed2;
	int 	nstar_tot;
	// Assume real(8) for mstars.
	double 	mstar_tot;
	double 	mstar_lost;
	int 	nsink;
};

// Temp particle holder
struct particle{
	float x,y,z,vx,vy,vz,m;
	int id,lvl;
};

class RamsesSource : public AbstractSource
{
  public:
    int readHeader();
    int readData();

private:

	std::vector<particle> Read(std::string);
	int computeTotalParticles();
	std::string getCurrentFilename();
	std::vector<particle> ReadIndivid(std::string filename);
	ramsesHeader ReadHeader(std::ifstream&);
	void ReadFloatData(ramsesHeader&, std::vector<particle>&, std::ifstream&);
	void ReadDoubleData(ramsesHeader&, std::vector<particle>&, std::ifstream&);

	std::string infile;
	int numCpus;
	int nFields;
	int fidx;
	int currentFile;

};

// Read a particular type's worth of data, convert to said type and return.
// Ramses outputs fortran unformatted file, so do not use file>>container it wont work.

template<typename T>
T ByteReader(std::ifstream& file)
{
	int dataSize = sizeof(T);

	char* c = new char(dataSize);

	file.read(c,dataSize);

	T data = *((T*)c);

	delete(c);

	return(data);
}

// Typedef specialisations of this function
typedef int (*ByteReadI_type)(std::ifstream&);
ByteReadI_type const IByteRead = &ByteReader<int>;

typedef float (*ByteReadF_type)(std::ifstream&);
ByteReadF_type const FByteRead = &ByteReader<float>;

typedef double (*ByteReadD_type)(std::ifstream&);
ByteReadD_type const DByteRead = &ByteReader<double>;

// Skip count bytes of unformatted data
void ByteSkipper(std::ifstream& file, int count);

#endif
