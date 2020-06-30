#ifndef RAMSES_READER_H
#define RAMSES_READER_H

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

#include "cxxsupport/paramfile.h"
#include "splotch/splotchutils.h"


#define QWORD_PRECISION 16
#define DOUBLE_PRECISION 8
#define SINGLE_PRECISION 4
#define INT_SIZE 4

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

class RamsesReader
{
public:
	RamsesReader();
	~RamsesReader();

	void Read(paramfile& params,std::vector<particle_sim>&);

private:

	int calculateTotalParticles();
	std::string getCurrentFilename();
	std::vector<particle_sim> ReadIndividual(std::string filename);
	ramsesHeader ReadHeader(std::ifstream&);
	void ReadSingleData(ramsesHeader&, std::vector<particle_sim>&, std::ifstream&);
	void ReadDoubleData(ramsesHeader&, std::vector<particle_sim>&, std::ifstream&);

	std::string infile;
	int numCpus;
	int fidx;
	int currentFile;

	// splotch extra params
	float scale3d;
	float scaleRGB;
	unsigned short type;
	float radius;
	float intensity;
	float active;
};

// Read a particular type's worth of data, convert to said type and return.
// Ramses outputs fortran unformatted file, so do not use file>>container it wont work.

template<typename T>
T ByteRead(std::ifstream& file)
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
ByteReadI_type const IByteRead = &ByteRead<int>;

typedef float (*ByteReadF_type)(std::ifstream&);
ByteReadF_type const FByteRead = &ByteRead<float>;

typedef double (*ByteReadD_type)(std::ifstream&);
ByteReadD_type const DByteRead = &ByteRead<double>;

// Skip count bytes of unformatted data
void ByteSkip(std::ifstream& file, int count);

#endif
