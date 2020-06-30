#include "reader.h"
#include "ramses_reader.h"

// Read function to fit in with current splotch reading conventions

void ramses_reader(paramfile& params, std::vector<particle_sim>& points)
{
	RamsesReader rr;
	rr.Read(params,points);
}

//
// Ramses reader class
//
// This reader currently only reads *3* dimensional particle data from Ramses file formats
// Checks are made for single vs double precision locations/velocities and read appropriately 
//  
// This Splotch implementation reads positions, and velocities (stored as colour data)
//
// Optional readable data is commented out in ReadSingleData and ReadDoubleData functions
// This is currently unnecessary for Splotch 
//
// Ramses data that currently isnt read (and isnt in the optional readable data sections) 
// birth_epoch, metallicity, sink properties
//
//


RamsesReader::RamsesReader() {}

RamsesReader::~RamsesReader() {}

void RamsesReader::Read(paramfile& params, std::vector<particle_sim>& points)
{
	std::cout << "Reading Ramses file" << std::endl;

	// Initialise from splotch params
	infile = params.find<std::string>("infile","part_"); 
	fidx = params.find<int>("fidx",1); 
	scale3d = params.find<float>("ramses_scale3d",10000);
	scaleRGB  = params.find<float>("ramses_scaleRGB",500);
	type = params.find<float>("ramses_type",0);
	radius = params.find<float>("ramses_radius",100);
	intensity = params.find<float>("ramses_intensity",1);
	active = params.find<float>("ramses_active",1);


	currentFile = 1;

	// Check how many particles in total
	int total = calculateTotalParticles();
	std::cout << "Total particles in Ramses files: " << total << std::endl;

	// Read all files
	while(currentFile <= numCpus)
	{	
		// Read current file (this could be parallelised...)
		std::vector<particle_sim> p = ReadIndividual(getCurrentFilename());

		// Add this file's particles to full list
		points.insert(points.end(),p.begin(),p.end());	

		currentFile++;
	}	

	std::cout << "Total particles read in: " << points.size() << std::endl;

	std::cout << "Read complete" << std::endl;
}

int RamsesReader::calculateTotalParticles()
{	
	// Work out total number of particles used in Ramses simulation
	int totalParticles = 0;
	currentFile = 1;

	// Open first file, check number of processors used
	// This is how many files need reading
	// Add nparts from first file to total while we're at it
	std::string fn = getCurrentFilename();

	std::cout << "filename: " << fn << std::endl;

	std::ifstream file(fn.c_str(), std::ios::in | std::ios::binary);

	// Check file
	if(!file)
		planck_fail("File not found: "+fn);

	ramsesHeader h = ReadHeader(file);
	file.close();
	numCpus = h.ncpu;
	totalParticles += h.npart;
	currentFile++;

	// Check headers of rest of files and add nparts to total
	while(currentFile <= numCpus)
	{
			std::string f = getCurrentFilename();
			std::ifstream cFile(f.c_str(), std::ios::in | std::ios::binary);

			// Check file
			if(!cFile)
				planck_fail("File not found: "+f);

 			h = ReadHeader(cFile);
 			file.close();
			totalParticles += h.npart;
			currentFile++;
	}
	
	//Reset currentFile
	currentFile = 1;

	return totalParticles;
}

std::string RamsesReader::getCurrentFilename()
{
	// Work out filename according to ramses convention and infile/fidx
	// Eg. for infile = part_ and fidx = 11 and currentFile = 1, filename = part_00011.out.00001 

	// Get infile
	std::string filename = infile;

	// Append fidx with 5 digit precision
	std::stringstream ss;
	ss.width(5);
	ss.fill('0');
	ss << std::right << fidx; 

	// Append .out
	filename += (ss.str() + ".out"); 
	
	//Append id of current file to be read
	ss.str("");
	ss.width(5);
	ss.fill('0');
	ss << std::right << currentFile;
	filename += ss.str(); 

	return filename;
}

std::vector<particle_sim> RamsesReader::ReadIndividual(std::string filename)
{
	// Vector for particle data
	std::vector<particle_sim> particles;

	// Open file
	std::ifstream file(filename.c_str(), std::ios::in | std::ios::binary);

	// Check file
	if(!file)
		planck_fail("File not found: "+filename);

	// Read header
	ramsesHeader h = ReadHeader(file);

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
		ReadSingleData(h, particles, file);
	}
	else
	{	
		//std::cout << "Data is double precision" << std::endl;
		ReadDoubleData(h, particles, file);
	}

	return particles;
}

ramsesHeader RamsesReader::ReadHeader(std::ifstream& file)
{
	// Each element of header data is flanked by two sets of 4 bytes, indicating the size of the data field
	// Eg, for a 4byte integer field storing the value 80, data would be arranged as so
	// Delimiter  |  Data     |  Delimiter
	// 04 00 00 00 50 00 00 00 04 00 00 00
	//
	// Currently little endian only. Will we potentially need to swap endianness? 
	// Assumption is currently made that integers are all 4 bytes. Checks are made for fortran real(4) vs real(8)
	// And in the case of localseed, real(8) vs real(16). Currently returns nan for real(16), not needed atm though.

	ramsesHeader h;

	int dSize;

	// For each integer value, read/discard delimiting 4bytes and read data
	ByteSkip(file,INT_SIZE);
	h.ncpu = IByteRead(file);

	ByteSkip(file,INT_SIZE*2);
	h.ndim = IByteRead(file);;

	ByteSkip(file,INT_SIZE*2);
	h.npart = IByteRead(file);

	ByteSkip(file,INT_SIZE);

	// For localseed, check for quadwords -real(16) - and store in second double if necessary
	// Returns nan for real(16) at the moment (probably cant store in doubles...)
	dSize = IByteRead(file);
	h.localseed1 = DByteRead(file);
	h.localseed2 = (dSize == QWORD_PRECISION) ?  DByteRead(file) : 0;
	
	ByteSkip(file,INT_SIZE*2);
	h.nstar_tot = IByteRead(file);

	ByteSkip(file,INT_SIZE);

	// For doubles, check for real(8) vs real(4), read appropriately
	dSize = IByteRead(file);
	h.mstar_tot = (dSize == DOUBLE_PRECISION) ?  DByteRead(file) : FByteRead(file);
	ByteSkip(file,INT_SIZE);

	dSize = IByteRead(file);
	h.mstar_lost = (dSize == DOUBLE_PRECISION) ?  DByteRead(file) : FByteRead(file);

	ByteSkip(file,INT_SIZE*2);
	h.nsink = IByteRead(file);

	ByteSkip(file,INT_SIZE);

	return h;
}

// Read data as single precision (4bytes per field)
void RamsesReader::ReadSingleData(ramsesHeader& h, std::vector<particle_sim>& particles, std::ifstream& file)
{
	// Read x data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].x = FByteRead(file) * scale3d;
	}

	// Skip delimiters
	ByteSkip(file,INT_SIZE*2);

	// Read y data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].y = FByteRead(file) * scale3d;
	}

	// Skip delimiters
	ByteSkip(file,INT_SIZE*2);

	// Read z data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].z = FByteRead(file) * scale3d;
	}			

	// Skip delimiters
	ByteSkip(file,INT_SIZE*2);

	// ----------------------- Read velocities into colour slots for splotch ------------------------

	// Read vx data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].e.r = FByteRead(file) * scaleRGB;
	}	

	// Skip delimiters
	ByteSkip(file,INT_SIZE*2);

	// Read vy data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].e.g = FByteRead(file) * scaleRGB;
	}	

	// Skip delimiters
	ByteSkip(file,INT_SIZE*2);

	// Read vz data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].e.b = FByteRead(file) * scaleRGB;
	}

	// -----------------------------Additional data for splotch particles -----------------------------

	// Fillers for splotch particle data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].r = radius;
		particles[i].I = intensity;
		particles[i].type = type;
		particles[i].active = active;
	}



	// ----------------------------- Other ramses data that could be read -----------------------------

	// // Skip delimiters
	// ByteSkip(file,INT_SIZE*2);

	// // Read mass data
	// for(int i = 0; i < h.npart; i++)
	// {
	// 	particles[i].m = FByteRead(file);
	// }		

	// // Skip delimiters
	// ByteSkip(file,INT_SIZE*2);

	// // Read identity data
	// for(int i = 0; i < h.npart; i++)
	// {
	// 	particles[i].id = IByteRead(file);
	// }	

	// // Skip delimiters
	// ByteSkip(file,INT_SIZE*2);

	// // Read level data
	// for(int i = 0; i < h.npart; i++)
	// {
	// 	particles[i].lvl = IByteRead(file);
	// }		
}

// Read data as double precision (8bytes per field)
void RamsesReader::ReadDoubleData(ramsesHeader& h, std::vector<particle_sim>& particles, std::ifstream& file)
{
	// Read x data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].x = DByteRead(file) * scale3d;
	}

	// Skip delimiters
	ByteSkip(file,INT_SIZE*2);

	// Read y data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].y = DByteRead(file) * scale3d;
	}

	// Skip delimiters
	ByteSkip(file,INT_SIZE*2);

	// Read z data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].z = DByteRead(file) * scale3d;
	}			

	// Skip delimiters
	ByteSkip(file,INT_SIZE*2);

	// ----------------------- Read velocities into colour slots for splotch ------------------------

	// Read vx data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].e.r = DByteRead(file) * scaleRGB;
	}	

	// Skip delimiters
	ByteSkip(file,INT_SIZE*2);

	// Read vy data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].e.g = DByteRead(file) * scaleRGB;
	}	

	// Skip delimiters
	ByteSkip(file,INT_SIZE*2);

	// Read vz data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].e.b = DByteRead(file) * scaleRGB;
	}	

	// -----------------------------Additional data for splotch particles -----------------------------

	// Fillers for splotch particle data
	for(int i = 0; i < h.npart; i++)
	{
		particles[i].r = radius;
		particles[i].I = intensity;
		particles[i].type = type;
		particles[i].active = active;
	}


	// ----------------------------- Other ramses data that could be read -----------------------------

	// // Skip delimiters
	// ByteSkip(file,INT_SIZE*2);

	// // Read mass data
	// for(int i = 0; i < h.npart; i++)
	// {
	// 	particles[i].m = DByteRead(file);
	// }

	// // Skip delimiters
	// ByteSkip(file,INT_SIZE*2);

	// // Read identity data
	// for(int i = 0; i < h.npart; i++)
	// {
	// 	particles[i].id = IByteRead(file);
	// }	

	// // Skip delimiters
	// ByteSkip(file,INT_SIZE*2);

	// // Read level data
	// for(int i = 0; i < h.npart; i++)
	// {
	// 	particles[i].lvl = IByteRead(file);
	// }			
}

void ByteSkip(std::ifstream& file, int count)
{
	char* c = new char(count);
	file.read(c,count);
	delete(c);
}
