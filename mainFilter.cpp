/***************************************************************************
 *   Copyright (C) 2008 by Ugo Becciani   *
 *   ugo.becciani@oact.inaf.it   *
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

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include "parametersparser.h"
#include "startFilter.h"

int main(int argc, char *argv[])
{

std::string paramFilename;

std::stringstream commandParametersSStream;
std::map<std::string, std::string> appParameters;
std::map<std::string, std::string>::iterator iter;

bool paramFileGiven=false;
if(argc==2)  // The parameter File is given!
{
	paramFileGiven=true;
	std::string argStr;
	argStr.assign(argv[1]);
	if (argStr=="--op"||argStr=="-op")
		paramFileGiven=false;
	if (argStr=="-help"||argStr=="--help")
		paramFileGiven=false;
	if(paramFileGiven)
	  	paramFilename=argStr;
}
 
if(paramFileGiven) //fill  appParameters with parameter file content
{
  ParametersParser myparser(paramFilename,1); 
  appParameters=myparser.getParameters(); //! local copy of the options given
} else
{
  bool fileIsAtTheEnd=true;
  for (int i=1;i<argc;i++) 
  {	
	std::string argStr;
	argStr.assign(argv[i]);
	if (argStr=="--file"||argStr=="-file")
	{	
		fileIsAtTheEnd=false;
		break;
	}
	if (argStr=="--op" || argStr=="-op")
	{	
		std::string argStrop;
		if(i+1>=argc)
		{
			std::cerr<<"Invalid options"<<std::endl;
			exit(1);
		}	
		argStrop.assign(argv[i+1]);
		if(argStrop=="interpolate" ||argStrop=="AHFstep" || 
		  argStrop=="merge" ||argStrop=="visualop" ||argStrop=="append")
		{
			fileIsAtTheEnd=false;
			break;
		}
	}
  }

  if (fileIsAtTheEnd)
	{
		for (int i=1;i<argc-1;i++) commandParametersSStream<< argv[i]<<" ";
		commandParametersSStream<<"--file "<< argv[argc-1]<<" ";
	}
  else
	for (int i=1;i<argc;i++) commandParametersSStream<< argv[i]<<" ";

// TEST ZONE 
//commandParametersSStream<<"--op extractlist --multilist ugo_oneList --asciilist --onelist --out bodyInField_newone.bin --file VSout_fl_8Ml_0.300.bin"<<" "; 
//commandParametersSStream<<"--op changecolname --field X Y  --newnames A B C --file VSout_fl_8Ml_0.300.bin"<<" "; 
//commandParametersSStream<<"--op mathop --expression math_halo_volume.txt --append --outcol HaloVolume1 --file VS_out_fl_0.3000_halos.bin"<<" "; 
//commandParametersSStream<<"--op mathop --compute <<(Rvir_10*Rvir_10*Rvir_10)/(1000000000)>> --append --outcol HaloVolume --file VS_out_fl_0.3000_halos.bin"<<" "; 
//commandParametersSStream<<"--op showtable --out ugo.txt --file VS_out_fl_0.3000_halos.bin"<<" "; 
//commandParametersSStream<<"--op pointdistribute --resolution 70  70 70  --points X Y Z --out VSgrid.bin --tsc --periodic --file VSout_fl_8Ml_0.300.bin"<<" "; 
//commandParametersSStream<<"--op splittable --field X --numoftables 4 --volumesplit 1 --out VSGsplit  --file VSgrid.bin"<<" "; 

//std::cerr << "TEST0: Command str -> " << commandParametersSStream.str() << std::endl;
//std::cerr << "TEST1: " << iter2->first << " " << iter2->second << std::endl;
//std::cerr << "TEST2: valOutFilename[i] = " << valOutFilename[i] << std::endl;
// END TEST ZONE
//

  ParametersParser myparser(commandParametersSStream.str());
  appParameters=myparser.getParameters();
} //if paramFileGiven ... else

startFilter startFilter(appParameters);
return EXIT_SUCCESS;
}

//#ifdef HAVE_CONFIG_H
//#include <config.h>
//#endif


