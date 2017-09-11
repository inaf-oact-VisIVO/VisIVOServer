//
// C++ Implementation: vsvbt2ahf
//
// Description: 
//
//
// Author:  <>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <iostream>
#include <string>
#include <sstream>
#include <set>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <fstream>

#include "vsvbt2ahf.h"
#include "vstable.h"

#include "AHFstep/param.h"
#include "AHFstep/tdef.h"
#include "AHFstep/common.c"

const unsigned int VSvbt2ahf::MAX_NUMBER_TO_REDUCE_ROW = 10000;
const unsigned int VSvbt2ahf::MIN_NUMBER_OF_ROW = 100;
const unsigned int VSvbt2ahf::NUM_AMIGA_OUT = 1000*6; //IMPORTANT multple of six!

//---------------------------------------------------------------------
VSvbt2ahf::VSvbt2ahf()
//---------------------------------------------------------------------
{
 m_fArray=NULL;
 m_Amiga=NULL;
 m_nOfRow=0;
 m_nOfCol=0;	
 m_vMean=0.0;
}

//---------------------------------------------------------------------
VSvbt2ahf::~VSvbt2ahf()
//---------------------------------------------------------------------
{	if(m_fArray!=NULL)
		for(unsigned int i=0;i<m_nOfCol;i++)
		{
			if(m_fArray[i]!=NULL) delete [] m_fArray[i];
		}
	if(m_fArray!=NULL) delete [] m_fArray;
	if(m_Amiga!=NULL) delete [] m_Amiga;

}
//---------------------------------------------------------------------
void VSvbt2ahf::printHelp()
//---------------------------------------------------------------------
{	
	std::cout<<"It creates an Amiga halofinder file (float) from six columns of a table, typically a snapshot of a cosmological simulation"<<std::endl<<std::endl;;

	std::cout<<"Usage: VisIVOFilters --op vbt2ahf [--field columns_name]  --par parfile  [--out filename_out] [--history] [--historyfile filename.xml] [--help] --file inputFile.bin"<<std::endl;

	std::cout<<"Example: VisIVOFilters --op vbt2ahf --field X Y Z Vx Vy Vz  --par myFile.par  [--out filename_out] [--help] --file input.bin"<<std::endl<<std::endl;

	std::cout<<"Note: "<<std::endl;
	std::cout<<"--field a valid list of  columns names. Default names   X Y Z Vx Vy Vz are considered"<<std::endl;
	std::cout<<"--par  A valid input parameter file. It must contain the following data (space separated):"<<std::endl;
	std::cout<<"	np (int, number of elements), box (double,  box size), mass_frac (double, mass of each element), vel_fact (float, factor to transform velocity field in Km/sec), omega0, lambda0, omegab, gamma, H_frac, T_initial, z_initial, z_current (double, cosmological parameters), numberOfTimesteps (int), header (string, header identifier)"<<std::endl;
	std::cout<<"--out is the output amiga data format file. Default name is given."<<std::endl;
    std::cout<<"--history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
    std::cout<<"--historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;
	std::cout<<"--file. Input table. Tipically a snapshot of a cosmological simulation."<<std::endl;
	std::cout<<"--help produce this output "<<std::endl;
	std::cout<<"WARNING for AMIGA HF users. No DONT_KIK_ME, EQ_MOT2,  ART,  ISOLATED, ISOLATED_P, XDOT are included in this filter."<<std::endl;

}
//---------------------------------------------------------------------
bool VSvbt2ahf::allocateArray()
//---------------------------------------------------------------------
{	
unsigned long long int tempLL=getMaxNumberInt()*2;
if(((unsigned long long int)m_nOfEle*m_nOfCol)>tempLL) 
	m_nOfEle=(int)tempLL/m_nOfCol;

try 
{
		m_Amiga= new float[NUM_AMIGA_OUT];
}
	
catch(std::bad_alloc &e) 
{
		m_Amiga=NULL;
}

if(m_Amiga==NULL) return false;


try 
{
		m_fArray= new float*[m_nOfCol];
}
	
catch(std::bad_alloc &e) 
{
		m_fArray=NULL;
}

if(m_fArray==NULL) return false;
bool goodAllocation=false;
while(!goodAllocation)
{
	goodAllocation=true;
	for(unsigned int i=0;i<m_nOfCol;i++)
	{
try
{
	  m_fArray[i]=new  float[m_nOfEle];
}
catch(std::bad_alloc &e)
{
	  m_fArray[i]=NULL;
}

	  if(m_fArray[i]==NULL) 
	  {	
		goodAllocation=false;
		for(unsigned int j=0;j<i;j++) 
			delete [] m_fArray[j];
		if(m_nOfEle==MIN_NUMBER_OF_ROW)
		{ 
			delete [] m_fArray;
			m_fArray=NULL;
			return false;
		}
		m_nOfEle=m_nOfEle-MAX_NUMBER_TO_REDUCE_ROW;
		if(m_nOfEle<=MAX_NUMBER_TO_REDUCE_ROW) m_nOfEle=MIN_NUMBER_OF_ROW;
				break;
	  }
	}
		
}
return true;   
}
//---------------------------------------------------------------------
bool VSvbt2ahf::execute()
//---------------------------------------------------------------------
{   
std::string filename;
if(getParameterAsString("par").empty() || getParameterAsString("par")=="unknown" )
{
	std::cerr<<"VSvbt2ahf: no file with parameters is given"<<std::endl;
	return false;
} else
{
	filename=getParameterAsString("par");
}
std::ifstream fileInput(filename.c_str());
if(!fileInput)
{
	std::cerr<<"VSvbt2ahf: cannot open parameter file "<<filename<<std::endl;
	return false;
}
std::set<unsigned int> colNumberSet;


if (getParameterAsString("field").empty() || getParameterAsString("field")=="unknown" || !isParameterPresent("field"))
{
	if(m_tables[0] -> getColId("X")>=0)
		colNumberSet.insert(m_tables[0] -> getColId("X"));
	if(m_tables[0] -> getColId("Y")>=0)
		colNumberSet.insert(m_tables[0] -> getColId("Y"));
	if(m_tables[0] -> getColId("Z")>=0)
		colNumberSet.insert(m_tables[0] -> getColId("Z"));
	if(m_tables[0] -> getColId("Vx")>=0)
		colNumberSet.insert(m_tables[0] -> getColId("Vx"));
	if(m_tables[0] -> getColId("Vy")>=0)
		colNumberSet.insert(m_tables[0] -> getColId("Vy"));
	if(m_tables[0] -> getColId("Vz")>=0)
		colNumberSet.insert(m_tables[0] -> getColId("Vz"));
}else
{	
	std::stringstream ssListparameters;
	ssListparameters.str(getParameterAsString("field"));
	while (!ssListparameters.eof())
	{	
		std::string paramField;
		ssListparameters>>paramField;
		
		if(m_tables[0] -> getColId(paramField)>=0)
			colNumberSet.insert(m_tables[0] -> getColId(paramField));
	}
}

m_nOfCol = colNumberSet.size();
if(m_nOfCol!=6)
{
	std::cerr<<"VSvbt2amiga. Invalid field names. Operation aborted"<<std::endl;
	return false;
}
m_nOfRow=m_tables[0] -> getNumberOfRows();
if(m_nOfRow>getMaxNumberInt()) m_nOfEle= getMaxNumberInt();
else m_nOfEle= m_nOfRow;

bool goodAllocation= allocateArray();
if(!goodAllocation){
	std::cerr<<"VSvbt2amiga. Invalid allocation"<<std::endl;
	return false;	
} 	

long long unsigned int np;
double box,m_fac,omega0,lambda0,omegab,gamma,H_frac,T_init,z_initial,z_current;
int  no_timestep;
char header[HEADERSTRING];
fileInput>> np;
fileInput>> box;
fileInput>>m_fac;
fileInput>>m_velFactor;
fileInput>>omega0;
fileInput>>lambda0;
fileInput>>omegab;
fileInput>>gamma;
fileInput>>H_frac;
fileInput>>T_init;
fileInput>>z_initial;
fileInput>>z_current;
fileInput>>no_timestep;
fileInput>>header;

io.header.no_vpart = np;

// riga 236 fly2amiga
std::string outName;
if(!isParameterPresent("out")||getParameterAsString("out")=="unknown") 	
	outName="VSvbt2amiga.amg";
else 
	outName= getParameterAsString("out");

std::ofstream fileOutput(outName.c_str(), std::ios::out | std::ios::binary);
if(!fileOutput)
{
      std::cerr<<"Cannot open binary file "<< outName <<std::endl;
      return false;
}

int machine_sizeof_long = sizeof(long);
fileOutput.write((char *) &machine_sizeof_long,1*sizeof(int));

double   a_initial = 1.0/(1.0+z_initial);
double   a_current = 1.0/(1.0+z_current);
strcpy(io.header.header,header);
io.header.multi_mass       = 0;
io.header.double_precision = 0;
io.header.no_part          = np;
io.header.boxsize          = box;
io.header.omega0           = omega0;
io.header.lambda0          = lambda0;
io.header.omegab           = omegab;
io.header.gamma            = gamma;
io.header.H_frac           = H_frac;
io.header.T_init           = T_init;
io.header.a_initial        = a_initial;
io.header.a_current        = a_current;
io.header.no_timestep      = no_timestep;
io.header.K_initial        = 0.0;
io.header.U_initial        = 0.0;
io.header.K_current        = 0.0;
io.header.U_current        = 0.0;
io.header.Eintegral        = 0.0;
io.header.Econst           = 0.0;
io.header.pmass            = m_fac;
io.header.version          = VERSION;
io.header.built            = BUILT;
fileOutput.write((char *) &(io.header), 1*sizeof(io.header));
unsigned long long int index=0;
unsigned long long int index_sel=0;
unsigned int * colId;
try {
	colId = new unsigned int[m_nOfCol];
}
catch(std::bad_alloc  &e){
	colId = NULL;
	}
if (colId == NULL) return false;
std::set<unsigned int>::iterator iter;
iter = colNumberSet.begin();
float  x_fac=1.0/ box;
float  v_fac= a_current / (box*H0); // not Isolated, not Isolated_P, not XDOT
int counterWriteAmiga=0;
for(int i=0;i<m_nOfCol;i++)
{
	colId[i]=*iter;
	iter++;
}
while(index< m_nOfRow)
{	
	m_tables[0]->getColumn(colId,m_nOfCol,index,index+m_nOfEle-1,m_fArray);
	for(int i=0;i< m_nOfEle;i++)
	{

// 		if(i==0) std::clog<<"Start 1. "<<m_fArray[0][i]<<" "<<m_fArray[1][i]<<" "<<m_fArray[2][i]<<" "<<m_fArray[3][i]<<" "<<m_fArray[4][i]<<" "<<m_fArray[5][i]<<std::endl;
// 		if(i==m_nOfEle-1) std::clog<<"Start 2. "<<m_fArray[0][i]<<" "<<m_fArray[1][i]<<" "<<m_fArray[2][i]<<" "<<m_fArray[3][i]<<" "<<m_fArray[4][i]<<" "<<m_fArray[5][i]<<std::endl;


		m_fArray[0][i] =fmod(m_fArray[0][i]*x_fac+1.0,1.0);
		m_fArray[1][i] =fmod(m_fArray[1][i]*x_fac+1.0,1.0);
		m_fArray[2][i] =fmod(m_fArray[2][i]*x_fac+1.0,1.0);
		m_fArray[3][i] *=m_velFactor;
		m_fArray[4][i] *=m_velFactor;
		m_fArray[5][i] *=m_velFactor;
		m_vMean += sqrt(m_fArray[3][i]*m_fArray[3][i]+m_fArray[4][i]*m_fArray[4][i]+m_fArray[5][i]*m_fArray[5][i]);
		m_fArray[3][i] *=v_fac;
		m_fArray[4][i] *=v_fac;
		m_fArray[5][i] *=v_fac;
		m_Amiga[0+counterWriteAmiga*6]=m_fArray[0][i];
		m_Amiga[1+counterWriteAmiga*6]=m_fArray[3][i];
		m_Amiga[2+counterWriteAmiga*6]=m_fArray[1][i];
		m_Amiga[3+counterWriteAmiga*6]=m_fArray[4][i];
		m_Amiga[4+counterWriteAmiga*6]=m_fArray[2][i];
		m_Amiga[5+counterWriteAmiga*6]=m_fArray[5][i];
/* scale masses to min. particle mass: here the two coincide, because there is only one mass, so: dw[6] = dw[6]/m_fac = mass/m_fac = 1 */
/*		if(i==0) std::clog<<"FINALE1. "<<m_fArray[0][i]<<" "<<m_fArray[3][i]<<" "<<m_fArray[1][i]<<" "<<m_fArray[4][i]<<" "<<m_fArray[2][i]<<" "<<m_fArray[5][i]<<std::endl;
		if(i==m_nOfEle-1) std::clog<<"FINALE2. "<<m_fArray[0][i]<<" "<<m_fArray[3][i]<<" "<<m_fArray[1][i]<<" "<<m_fArray[4][i]<<" "<<m_fArray[2][i]<<" "<<m_fArray[5][i]<<std::endl;*/
		counterWriteAmiga++;
		if(counterWriteAmiga==NUM_AMIGA_OUT/6 - 1 || i==m_nOfEle-1)
		{
			fileOutput.write((char *) m_Amiga, counterWriteAmiga*6*sizeof(float));
			counterWriteAmiga=0;						
		}
	}
	index+= m_nOfEle;
	if(index+m_nOfEle>= m_nOfRow) m_nOfEle= m_nOfRow-index;
}
	std::cout<<"Done. Mean velocity in Km/sec is "<<m_vMean/np<<" "<<m_vMean<<std::endl;
delete [] colId;
return true;
}


