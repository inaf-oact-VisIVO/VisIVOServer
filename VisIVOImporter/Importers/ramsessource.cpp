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
#include <cstdlib>
#include <cstring>


#include "ramsessource.h"
#include "visivoutils.h"
#include <iostream>
#include <fstream>
#include <sstream>

//---------------------------------------------------------------------
int RamsesSource::readData()
//---------------------------------------------------------------------
{
  	// Vector for particle data
	std::vector<particle> particles;
	
	std::string infile=m_pointsFileName.c_str();
	std::ofstream outfile(m_pointsBinaryName.c_str(),std::ofstream::binary ); 
    particles=Read(infile);
	if( particles.size()==0 )
    {
        return -1;
    }
	// std::cout << "N Fields = " << nFields << std::endl;
    if(nFields < 1){
		return -1;
	}	
	m_fieldNames.push_back("X");
	m_fieldNames.push_back("Y");
	m_fieldNames.push_back("Z");
	m_nCols=m_fieldNames.size();
  	
	int nPart=particles.size();
	makeHeader(nPart,m_pointsBinaryName,m_fieldNames,m_cellSize,m_cellComp,m_volumeOrTable); //!crea$
	
	for(int i = 0; i < nPart; i++)
	{
	  //std::cout<<particles[i].x<<std::endl;
	  outfile.write((char*)&(particles[i].x), sizeof(float)); 
	}

	for(int i = 0; i < nPart; i++)
	{
	  // std::cout<<particles[i].y<<std::endl;
		if(nFields < 2){
	  		particles[i].y=MISSING_VALUE; 
		}
	  		outfile.write((char*)&(particles[i].y), sizeof(float)); 
	}

	for(int i = 0; i < nPart; i++)
	{
	  // std::cout<<particles[i].z<<std::endl;
		if(nFields < 3){
	  		particles[i].z=MISSING_VALUE; 
		}
	  outfile.write((char*)&(particles[i].z), sizeof(float)); 
	}
	
	outfile.close();
	return 0;
}	

//---------------------------------------------------------------------
  int RamsesSource::readHeader()
//---------------------------------------------------------------------
{
  	
	return 0;  
}

//---------------------------------------------------------------------
std::vector<particle> RamsesSource::Read(std::string inputName)
//---------------------------------------------------------------------
{
	// Initialise
	std::vector<particle> points;
	infile = inputName; //params get infile (default part_)
	fidx = 11; //params get fidx (default 0? or throw error?)
	currentFile = 1;

	// Check how many particles in total
	int total = calculateTotalParticles();
    if(total <=0) return points;
	// std::cout << "Total particles in Ramses files: " << total << std::endl;

	// Read all files
	while(currentFile <= numCpus)
	{	
		// Read current file
		std::vector<particle> p = ReadIndividual(getCurrentFilename());

		// Add these particles to full list
		// (do this in read to allow indicidual reads to be mpi parallelised?)
		points.insert(points.end(),p.begin(),p.end());	

		currentFile++;
	}	

	// // Test output some particles
	// for(int i = 1463; i < 1473; i++)
	// {
	// 	std::cout << "Particle[" << i << "].x: " << points[i].x << std::endl;
	// 	std::cout << "Particle[" << i << "].y: " << points[i].y << std::endl;
	// 	std::cout << "Particle[" << i << "].z: " << points[i].z << std::endl;
	// 	std::cout << "Particle[" << i << "].vx: " << points[i].vx << std::endl;
	// 	std::cout << "Particle[" << i << "].vy: " << points[i].vy << std::endl;
	// 	std::cout << "Particle[" << i << "].vz: " << points[i].vz << std::endl;
	// 	std::cout << "Particle[" << i << "].m: " << points[i].m << std::endl;
	// 	std::cout << "Particle[" << i << "].id: " << points[i].id << std::endl;
	// 	std::cout << "Particle[" << i << "].lvl: " << points[i].lvl << std::endl;
	// }

//	std::cout << "Vector size: " << points.size() << std::endl;
	return points;
}

//---------------------------------------------------------------------
int RamsesSource::calculateTotalParticles()
//---------------------------------------------------------------------
{
	// Work out total number of particles used in Ramses simulation
	int totalParticles = 0;
	currentFile = 1;

	// Open first file, check number of processors used
	// This is how many files need reading
	// Add nparts from first file to total while we're at it
	std::string fn = getCurrentFilename();


	std::ifstream file(fn.c_str(), std::ios::in | std::ios::binary);
    if(!file.is_open())
    {
        std::cerr << "Invalid filename: " << fn << std::endl;
        return -1;
    }
	ramsesHeader h = ReadHeader(file);
	file.close();
	numCpus = h.ncpu;
    if(numCpus<=0)
    {
        std::cerr << "Invalid cpu number: " << numCpus << std::endl;
        return -1;
    }
	totalParticles += h.npart;
    if(totalParticles<=0)
    {
        std::cerr << "Invalid total particles number: " << totalParticles << std::endl;
        return -1;
    }
	currentFile++;

	// Check headers of rest of files and add nparts to total
	while(currentFile <= numCpus)
	{
			std::string f = getCurrentFilename();
			std::ifstream cFile(f.c_str(), std::ios::in | std::ios::binary);
 			h = ReadHeader(cFile);
 			file.close();
			totalParticles += h.npart;
			currentFile++;
	}
	
	//Reset currentFile
	currentFile = 1;

	return totalParticles;
}
//---------------------------------------------------------------------
std::string RamsesSource::getCurrentFilename()
//---------------------------------------------------------------------
{
	std::string filename = infile;
	std::stringstream ss;
	//ss.width(5);
	//ss.fill('0');
	//ss << std::right << fidx; //fidx from paramfile
	filename += (ss.str() + ".out"); 
	ss.str("");
	ss.width(5);
	ss.fill('0');
	ss << std::right << currentFile; //currentFile is ramses file
	filename += ss.str(); 	
	return filename;
}

//---------------------------------------------------------------------
std::vector<particle> RamsesSource::ReadIndividual(std::string filename)
//---------------------------------------------------------------------
{
	// Vector for particle data
	std::vector<particle> particles;

	// Open file
	std::ifstream file(filename.c_str(), std::ios::in | std::ios::binary);

	// Read header
	ramsesHeader h = ReadHeader(file);

	// // Check read
	//std::cout << "ncpu = " << h.ncpu << std::endl;
	//std::cout << "ndim = " << h.ndim << std::endl;
	//std::cout << "npart = " << h.npart << std::endl;
	//std::cout << "localseed1 = " << h.localseed1 << std::endl;
	//std::cout << "localseed2 = " << h.localseed2 << std::endl;
	//std::cout << "nstar_tot = " << h.nstar_tot << std::endl;
	//std::cout << "mstar_tot = " << h.mstar_tot << std::endl;
	//std::cout << "mstar_lost = " << h.mstar_lost << std::endl;
	//std::cout << "nsink = " << h.nsink << std::endl;

	// Allocate space for data
	particles.resize(h.npart);

	// Read xPos delimiter
	double xPosData = IByteRead(file);

	// Check precision (check for other precisions?)
	float precision = (xPosData/(double)h.npart == DOUBLE_PRECISION) ? DOUBLE_PRECISION : SINGLE_PRECISION;

	// Read data
	if(precision==SINGLE_PRECISION)
	{	
		//std::cout << "Data is single precision" << std::endl;
		ReadFloatData(h, particles, file);
	}
	else
	{	
		//std::cout << "Data is double precision" << std::endl;
		ReadDoubleData(h, particles, file);
	}
	
	nFields = h.ndim;	

	return particles;
}

//---------------------------------------------------------------------
ramsesHeader RamsesSource::ReadHeader(std::ifstream& file)
//---------------------------------------------------------------------
{
	// Each element of header data is flanked by two sets of 4 bytes, indicating the size of the data field
	// Eg, for a 4byte integer field storing the value 80, data would be arranged as so
	// Delimiter  |  Data     |  Delimiter
	// 04 00 00 00 50 00 00 00 04 00 00 00
	// Currently little endian, is this likely to be different depending on source data?  
	// Assumption is currently made that integers are all 4 bytes. Checks are made for fortran real(4) vs real(8)
	// And in the case of localseed, real(8) vs real(16). Currently returns nan for real(16), not needed atm though.
	// Replace sizeofint with constant defined at runtime when making class!
	// Error checking...

	ramsesHeader h;

	int dSize;

	// For each integer value, read/discard delimiting 4bytes and read data
	ByteSkip(file,sizeof(int));
	h.ncpu = IByteRead(file);

	ByteSkip(file,sizeof(int)*2);
	h.ndim = IByteRead(file);;

	ByteSkip(file,sizeof(int)*2);
	h.npart = IByteRead(file);

	ByteSkip(file,sizeof(int));

	// For localseed, check for quadwords -real(16) - and store in second double if necessary
	// Returns nan for real(16) at the moment (probably cant store in doubles...)
	dSize = IByteRead(file);
	h.localseed1 = DByteRead(file);
	h.localseed2 = (dSize == QWORD_PRECISION) ?  DByteRead(file) : 0;
	
	ByteSkip(file,sizeof(int)*2);
	h.nstar_tot = IByteRead(file);

	ByteSkip(file,sizeof(int));

	// For doubles, check for real(8) vs real(4)
	dSize = IByteRead(file);
	h.mstar_tot = (dSize == DOUBLE_PRECISION) ?  DByteRead(file) : FByteRead(file);
	ByteSkip(file,sizeof(int));

	dSize = IByteRead(file);
	h.mstar_lost = (dSize == DOUBLE_PRECISION) ?  DByteRead(file) : FByteRead(file);

	ByteSkip(file,sizeof(int)*2);
	h.nsink = IByteRead(file);

	ByteSkip(file,sizeof(int));

	return h;
}

//---------------------------------------------------------------------
void RamsesSource::ReadFloatData(ramsesHeader& h, std::vector<particle>& particles, std::ifstream& file)
//---------------------------------------------------------------------
{
	// Read x data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].x = FByteRead(file);
	}

	// Skip delimiters
	ByteSkip(file,sizeof(int)*2);

	// Read y data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].y = FByteRead(file);
	}

	// Skip delimiters
	ByteSkip(file,sizeof(int)*2);

	// Read z data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].z = FByteRead(file);
	}			

	// Skip delimiters
	ByteSkip(file,sizeof(int)*2);

	// Read vx data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].vx = FByteRead(file);
	}	

	// Skip delimiters
	ByteSkip(file,sizeof(int)*2);

	// Read vy data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].vy = FByteRead(file);
	}	

	// Skip delimiters
	ByteSkip(file,sizeof(int)*2);

	// Read vz data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].vz = FByteRead(file);
	}	

	// Skip delimiters
	ByteSkip(file,sizeof(int)*2);

	// Read mass data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].m = FByteRead(file);
	}		

	// Skip delimiters
	ByteSkip(file,sizeof(int)*2);

	// Read identity data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].id = IByteRead(file);
	}	

	// Skip delimiters
	ByteSkip(file,sizeof(int)*2);

	// Read level data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].lvl = IByteRead(file);
	}		
}

//---------------------------------------------------------------------
void RamsesSource::ReadDoubleData(ramsesHeader& h, std::vector<particle>& particles, std::ifstream& file)
//---------------------------------------------------------------------
{
	// Read x data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].x = DByteRead(file);
	}

	// Skip delimiters
	ByteSkip(file,sizeof(int)*2);

	// Read y data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].y = DByteRead(file);
	}

	// Skip delimiters
	ByteSkip(file,sizeof(int)*2);

	// Read z data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].z = DByteRead(file);
	}			

	// Skip delimiters
	ByteSkip(file,sizeof(int)*2);

	// Read vx data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].vx = DByteRead(file);
	}	

	// Skip delimiters
	ByteSkip(file,sizeof(int)*2);

	// Read vy data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].vy = DByteRead(file);
	}	

	// Skip delimiters
	ByteSkip(file,sizeof(int)*2);

	// Read vz data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].vz = DByteRead(file);
	}	

	// Skip delimiters
	ByteSkip(file,sizeof(int)*2);

	// Read mass data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].m = DByteRead(file);
	}

	// Skip delimiters
	ByteSkip(file,sizeof(int)*2);

	// Read identity data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].id = IByteRead(file);
	}	

	// Skip delimiters
	ByteSkip(file,sizeof(int)*2);

	// Read level data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].lvl = IByteRead(file);
	}			
}

//---------------------------------------------------------------------
void ByteSkip(std::ifstream& file, int count)
//---------------------------------------------------------------------
{
	char* c = new char(count);
	file.read(c,count);
	delete(c);
}
