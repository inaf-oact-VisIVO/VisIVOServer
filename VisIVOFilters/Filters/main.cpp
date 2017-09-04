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
#include<sstream>
#include<string>

#include "vstable.h"

#include "vsselectcolumnsop.h"
#include "vsrandomizertableop.h"
#include "parametersparser.h"
#include "vsmergeop.h"
#include "vsappend.h"
#include "vsselectfieldop.h"
#include "vsmathop.h"
#include "vsdecimatorop.h"
#include "vsvisualop.h"
#include "vssphereop.h"
#include "vsstatisticop.h"
#include "vsshowtableop.h"
#include "vspointdistributeop.h"
#include "vscoarsevolumeop.h"
#include "vsextractsubvolumeop.h"
#include "vsexampleop.h"
#include "vspointpropertyop.h"
#include "vsinterpolateop.h"
#include "vsmoduleop.h"
#include "vssigmacontoursop.h"
#include "vspolarop.h"
#include "AHFstep/main_AHFstep.h"
#include "vsvbt2ahf.h"
#include "vsgrid2pointdistr.h"
#include "vsswapop.h"
#include "vsextractlistrowsop.h"
#include "vsaddidop.h"
#include "vsahfhalolistop.h"
#include "vschangecolnameop.h"
#include "vsahfhalogalaxyextop.h"
#include "vssplittableop.h"
#include "vswrvotableop.h"
#include "vscutop.h"
#include "vsvollimit.h"
#include "vsincludeop.h"


int main(int argc, char *argv[])
{

std::string filename,paramFilename;
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
		if(argStrop=="interpolate" ||argStrop=="AHFstep")
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
// END TEST ZONE
//

  ParametersParser myparser(commandParametersSStream.str());
  appParameters=myparser.getParameters();
} //if paramFileGiven ... else

iter =appParameters.find("op");
  if( iter == appParameters.end())
  {
    iter =appParameters.find("help");
      if( iter == appParameters.end())
	 std::cerr <<"No operation is requested"<<std::endl;
    std::cerr<<"VisIVOFilters version 1.2  January 24th 2011"<<std::endl<<std::endl;
    std::cerr <<"Syntax1: VisIVOFilters --op operation  [PARAMETERS] [--help]"<<std::endl;
    std::cerr <<"Syntax2: VisIVOFilters parameterFile"<<std::endl;
    std::cerr <<"valid operations: randomizer selcolumns merge append selfield mathop decimator extraction visualop showtable statistic pointdistribute pointproperty coarsevolume extractsubvolume interpolate module sigmacontours cartesian2polar AHFstep vbt2ahf grid2point extractlist addId ahfhalolist ahfhalogalaxyext changecolname splittable wrvotable include"<<std::endl;
    return EXIT_SUCCESS;
   }
std::stringstream sstreamOp(iter->second);
appParameters.erase(iter);
int idOp=-1;

if(sstreamOp.str()=="randomizer") idOp=1;
if(sstreamOp.str()=="selcolumns") idOp=2;
if(sstreamOp.str()=="merge") idOp=3;
if(sstreamOp.str()=="append") idOp=4;
if(sstreamOp.str()=="selfield") idOp=5;
if(sstreamOp.str()=="mathop") idOp=6;
if(sstreamOp.str()=="decimator") idOp=7;
if(sstreamOp.str()=="visualop") idOp=8;
if(sstreamOp.str()=="extraction") idOp=9;
if(sstreamOp.str()=="showtable") idOp=10;
if(sstreamOp.str()=="statistic") idOp=11;
if(sstreamOp.str()=="pointdistribute") idOp=12;
if(sstreamOp.str()=="coarsevolume") idOp=13;
if(sstreamOp.str()=="extractsubvolume") idOp=14;
if(sstreamOp.str()=="example") idOp=15;
if(sstreamOp.str()=="pointproperty") idOp=16;
if(sstreamOp.str()=="interpolate") idOp=17;
if(sstreamOp.str()=="module") idOp=18;
if(sstreamOp.str()=="sigmacontours") idOp=19;
if(sstreamOp.str()=="cartesian2polar") idOp=20;
if(sstreamOp.str()=="AHFstep") idOp=21;
if(sstreamOp.str()=="vbt2ahf") idOp=22;
if(sstreamOp.str()=="grid2point") idOp=23;
if(sstreamOp.str()=="swap") idOp=24;
if(sstreamOp.str()=="extractlist") idOp=25;
if(sstreamOp.str()=="addId") idOp=26;
if(sstreamOp.str()=="ahfhalolist") idOp=27;
if(sstreamOp.str()=="changecolname") idOp=28;
if(sstreamOp.str()=="ahfhalogalaxyext") idOp=29;
if(sstreamOp.str()=="splittable") idOp=30;
if(sstreamOp.str()=="wrvotable") idOp=31;
if(sstreamOp.str()=="cut") idOp=32;
if(sstreamOp.str()=="showvol") idOp=33;
if(sstreamOp.str()=="include") idOp=34;

switch(idOp)
{
/*** Randomizer OP **/
case 1:
{
iter =appParameters.find("help");
  if( iter != appParameters.end())
  {
  VSRandomizerTableOp op;
  op.printHelp();
  return 1;
  }

iter =appParameters.find("file");
  if( iter == appParameters.end())
  {
    std::cerr <<"No input file table is provided"<<std::endl;
    return 1;
   }
std::stringstream sFilename(iter->second);
appParameters.erase(iter);

sFilename>>filename;
if(filename.find(".bin") == std::string::npos)
	    filename.append(".bin");
VSTable table(filename);

if(!table.tableExist())
  {
    std::cerr <<"No valid input file table is provided"<<std::endl;
    return 1;
   }


VSRandomizerTableOp op;
op.setParameters(appParameters);
op.addInput(&table);
op.execute();
break;
}
/***END Randomizer OP **/
/*** Select Columns OP **/
case 2:
{
iter =appParameters.find("help");
  if( iter != appParameters.end())
  {
  VSSelectColumnsOp op;
  op.printHelp();
  return 1;
  }
iter =appParameters.find("file");
  if( iter == appParameters.end())
  {
    std::cerr <<"No input file table is provided"<<std::endl;
    return 1;
   }
std::stringstream sFilename(iter->second);
appParameters.erase(iter);
sFilename>>filename;
if(filename.find(".bin") == std::string::npos)
	    filename.append(".bin");
VSTable table(filename);
if(!table.tableExist())
  {
    std::cerr <<"No valid input file table is provided"<<std::endl;
    return 1;
   }


VSSelectColumnsOp op;
op.setParameters(appParameters);
op.addInput(&table);
op.execute();

break;
}
/*** END Select Columns OP **/
/*** Merge Tables OP **/
case 3:
{
iter =appParameters.find("help");
  if( iter != appParameters.end())
  {
  VSMergeOp op;
  op.printHelp();
  return 1;
  }
VSMergeOp op;
op.setParameters(appParameters);
op.execute();

break;
}
/*** END Merge Tables OP **/
/*** Append Tables OP **/
case 4:
{
iter =appParameters.find("help");
  if( iter != appParameters.end())
  {
  VSAppendOp op;
  op.printHelp();
  return 1;
  }
VSAppendOp op;
op.setParameters(appParameters);
op.execute();

break;
}
/*** END Append Tables OP **/
/*** SelectField  OP **/
case 5:
{
iter =appParameters.find("help");
  if( iter != appParameters.end())
  {
  VSSelectFieldOp op;
  op.printHelp();
  return 1;
  }

iter =appParameters.find("file");
  if( iter == appParameters.end())
  {
    std::cerr <<"No input file table is provided"<<std::endl;
    return 1;
   }
std::stringstream sFilename(iter->second);
appParameters.erase(iter);
sFilename>>filename;
if(filename.find(".bin") == std::string::npos)
	    filename.append(".bin");
VSTable table(filename);
if(!table.tableExist())
  {
    std::cerr <<"No valid input file table is provided"<<std::endl;
    return 1;
   }

VSSelectFieldOp op;
op.setParameters(appParameters);
op.addInput(&table);
op.execute();

break;
}
/*** END SelectField  OP **/
/*** Math  OP **/
case 6:
{
iter =appParameters.find("help");
  if( iter != appParameters.end())
  {
  VSMathOp op;
  op.printHelp();
  return 1;
  }

iter =appParameters.find("file");
  if( iter == appParameters.end())
  {
    std::cerr <<"No input file table is provided"<<std::endl;
    return 1;
   }
std::stringstream sFilename(iter->second);
appParameters.erase(iter);
sFilename>>filename;
if(filename.find(".bin") == std::string::npos)
	    filename.append(".bin");
VSTable table(filename);
if(!table.tableExist())
  {
    std::cerr <<"No valid input file table is provided"<<std::endl;
    return 1;
   }

VSMathOp op;
op.setParameters(appParameters);
op.addInput(&table);
op.execute();

break;
}
/*** END Math  OP **/
/*** Decimator OP **/
case 7:
{
iter =appParameters.find("help");
  if( iter != appParameters.end())
  {
  VSDecimatorTableOp op;
  op.printHelp();
  return 1;
  }

iter =appParameters.find("file");
  if( iter == appParameters.end())
  {
    std::cerr <<"No input file table is provided"<<std::endl;
    return 1;
   }


std::stringstream sFilename(iter->second);
appParameters.erase(iter);

sFilename>>filename;
if(filename.find(".bin") == std::string::npos)
	    filename.append(".bin");
VSTable table(filename);
if(!table.tableExist())
  {
    std::cerr <<"No valid input file table is provided"<<std::endl;
    return 1;
   }

VSDecimatorTableOp op;
op.setParameters(appParameters);
op.addInput(&table);
op.execute();
break;
}
/***END Randomizer OP **/
/*** Visual OP **/
case 8:
{
iter =appParameters.find("help");
  if( iter != appParameters.end())
  {
  VSVisualOp op;
  op.printHelp();
  return 1;
  }
VSVisualOp op;
op.setParameters(appParameters);
op.execute();

break;
}
/*** END Visual OP **/
/*** Sphere  OP **/
case 9:
{
iter =appParameters.find("help");
  if( iter != appParameters.end())
  {
  VSSphereOp op;
  op.printHelp();
  return 1;
  }

iter =appParameters.find("file");
  if( iter == appParameters.end())
  {
    std::cerr <<"No input file table is provided"<<std::endl;
    return 1;
   }
std::stringstream sFilename(iter->second);
appParameters.erase(iter);
sFilename>>filename;
if(filename.find(".bin") == std::string::npos)
	    filename.append(".bin");
VSTable table(filename);
if(!table.tableExist())
  {
    std::cerr <<"No valid input file table is provided"<<std::endl;
    return 1;
   }

VSSphereOp op;
op.setParameters(appParameters);
op.addInput(&table);
op.execute();

break;
}
/*** END SelectField  OP **/
/*** ShowTable  OP **/
case 10:
{
iter =appParameters.find("help");
  if( iter != appParameters.end())
  {
  VSShowTableOp op;
  op.printHelp();
  return 1;
  }

iter =appParameters.find("file");
  if( iter == appParameters.end())
  {
    std::cerr <<"No input file table is provided"<<std::endl;
    return 1;
   }
std::stringstream sFilename(iter->second);
appParameters.erase(iter);
sFilename>>filename;
if(filename.find(".bin") == std::string::npos)
	    filename.append(".bin");
VSTable table(filename);
if(!table.tableExist())
  {
    std::cerr <<"No valid input file table is provided"<<std::endl;
    return 1;
   }

VSShowTableOp op;
op.setParameters(appParameters);
op.addInput(&table);
op.execute();

break;
}
/*** END ShowTable  OP **/
/*** Statistic  OP **/
case 11:
{
iter =appParameters.find("help");
  if( iter != appParameters.end())
  {
  VSStatisticOp op;
  op.printHelp();
  return 1;
  }

iter =appParameters.find("file");
  if( iter == appParameters.end())
  {
    std::cerr <<"No input file table is provided"<<std::endl;
    return 1;
   }
std::stringstream sFilename(iter->second);
appParameters.erase(iter);
sFilename>>filename;
if(filename.find(".bin") == std::string::npos)
	    filename.append(".bin");
VSTable table(filename);
if(!table.tableExist())
  {
    std::cerr <<"No valid input file table is provided"<<std::endl;
    return 1;
   }

VSStatisticOp op;
op.setParameters(appParameters);
op.addInput(&table);
op.execute();

break;
}
/*** END Statistic OP **/
/*** PointDistribute  OP **/
case 12:
{
iter =appParameters.find("help");
  if( iter != appParameters.end())
  {
  VSPointDistributeOp op;
  op.printHelp();
  return 1;
  }

iter =appParameters.find("file");
  if( iter == appParameters.end())
  {
    std::cerr <<"No input file table is provided"<<std::endl;
    return 1;
   }
std::stringstream sFilename(iter->second);
appParameters.erase(iter);
sFilename>>filename;
if(filename.find(".bin") == std::string::npos)
	    filename.append(".bin");
VSTable table(filename);
if(!table.tableExist())
  {
    std::cerr <<"No valid input file table is provided"<<std::endl;
    return 1;
   }

VSPointDistributeOp op;
op.setParameters(appParameters);
op.addInput(&table);
op.execute();

break;
}
/*** END PointDistribute OP **/
/*** CoarseVolume  OP **/
case 13:
{
iter =appParameters.find("help");
  if( iter != appParameters.end())
  {
  VSCoarseVolumeOp op;
  op.printHelp();
  return 1;
  }

iter =appParameters.find("file");
  if( iter == appParameters.end())
  {
    std::cerr <<"No input file table is provided"<<std::endl;
    return 1;
   }
std::stringstream sFilename(iter->second);
appParameters.erase(iter);
sFilename>>filename;
if(filename.find(".bin") == std::string::npos)
	    filename.append(".bin");
VSTable table(filename);
if(!table.tableExist())
  {
    std::cerr <<"No valid input file table is provided"<<std::endl;
    return 1;
   }

VSCoarseVolumeOp op;
op.setParameters(appParameters);
op.addInput(&table);
op.execute();

break;
}
/*** END CoarseVolume OP **/
/*** ExtractSubVolume  OP **/
case 14:
{
iter =appParameters.find("help");
  if( iter != appParameters.end())
  {
  VSExtractSubVolumeOp op;
  op.printHelp();
  return 1;
  }

iter =appParameters.find("file");
  if( iter == appParameters.end())
  {
    std::cerr <<"No input file table is provided"<<std::endl;
    return 1;
   }
std::stringstream sFilename(iter->second);
appParameters.erase(iter);
sFilename>>filename;
if(filename.find(".bin") == std::string::npos)
	    filename.append(".bin");
VSTable table(filename);
if(!table.tableExist())
  {
    std::cerr <<"No valid input file table is provided"<<std::endl;
    return 1;
   }

VSExtractSubVolumeOp op;
op.setParameters(appParameters);
op.addInput(&table);
op.execute();

break;
}
/*** END ExtracSubVolume OP **/
/*** Example  OP Used for developers only **/
case 15:
{
iter =appParameters.find("help");  //! if --help is given, print help
  if( iter != appParameters.end())
  {
  VSExampleOp op;
  op.printHelp();
  return 1;
  }

iter =appParameters.find("file");  //!check for --file parameter
  if( iter == appParameters.end())
  {
    std::cerr <<"No input file table is provided"<<std::endl;
    return 1;
   }
std::stringstream sFilename(iter->second);
sFilename>>filename; //! name of the file binary data format
if(filename.find(".bin") == std::string::npos)  //! table name must contain .bin extension
	    filename.append(".bin");
VSTable table(filename); //! create an object VSTable
if(!table.tableExist()) //! check that the file exist
  {
    std::cerr <<"No valid input file table is provided"<<std::endl;
    return 1;
   }

VSExampleOp op; //! specific operation to be executed (defined in a cpp class. All classes must be derived from VSTableOp)

op.setParameters(appParameters); //!parser of input line command: fill a map(key=string option, value=string value) defined in VSTableOP. The map is a private member  of the VSTableOP. Contained values can be obtained with the method getParameterAsString, getParameterAsFloat, getParameterAsInt public functions of the class. Parameters can be also set with addParameter method. 

op.addInput(&table); //! add the created table to the operation (defined in VSTableOp)    			std::vector<VSTable *> m_tables; m_tables[i] is a pointer to a 			VSTable object.

op.execute(); //! execute the action specific of the operation

break;
}
/*** END Example OP **/
/*** PointProperty  OP **/
case 16:
{
iter =appParameters.find("help");
  if( iter != appParameters.end())
  {
  VSPointPropertyOp op;
  op.printHelp();
  return 1;
  }

iter =appParameters.find("file");
  if( iter == appParameters.end())
  {
    std::cerr <<"No input file table is provided"<<std::endl;
    return 1;
   }
std::stringstream sFilename(iter->second);
appParameters.erase(iter);
sFilename>>filename;
if(filename.find(".bin") == std::string::npos)
	    filename.append(".bin");
VSTable table(filename);
if(!table.tableExist())
  {
    std::cerr <<"No valid input file table is provided"<<std::endl;
    return 1;
   }

VSPointPropertyOp op;
op.setParameters(appParameters);
op.addInput(&table);
op.execute();

break;
}
/*** END PointProperty OP **/
/*** Interpolate OP **/
case 17:
{
iter =appParameters.find("help");
  if( iter != appParameters.end())
  {
  VSInterpolateOp op;
  op.printHelp();
  return 1;
  }
VSInterpolateOp op;
op.setParameters(appParameters);
op.execute();

break;
}
/*** END Merge Tables OP **/
case 18:
{
iter =appParameters.find("help");
if( iter != appParameters.end())
{
  VSModuleOp op;
  op.printHelp();
  return 1;
}
iter =appParameters.find("file");
if( iter == appParameters.end())
{
    std::cerr <<"No input file table is provided"<<std::endl;
    return 1;
}
std::stringstream sFilename(iter->second);
appParameters.erase(iter);
sFilename>>filename;
if(filename.find(".bin") == std::string::npos)
	    filename.append(".bin");
VSTable table(filename);
if(!table.tableExist())
  {
    std::cerr <<"No valid input file table is provided"<<std::endl;
    return 1;
   }
VSModuleOp op;
op.setParameters(appParameters);
op.addInput(&table);
op.execute();

break;
}
case 19:
{
	iter =appParameters.find("help");
  if( iter != appParameters.end()){
  	VSSigmaContoursOp op;
  	op.printHelp();
  	return EXIT_SUCCESS;
  }
  
    
	VSSigmaContoursOp op;
	op.setParameters(appParameters);
	iter= appParameters.find("file");
	if(iter==appParameters.end()){
		std::cerr<<"No input file table is provided!"<<std::endl;
		return 1;
	}//close if
	std::stringstream sFilename(iter->second);
	appParameters.erase(iter);
	sFilename>>filename;
	if(filename.find(".bin") == std::string::npos) filename.append(".bin");
	VSTable table(filename);
	if(!table.tableExist()){
		std::cerr<<"No valid input file table is provided"<<std::endl;
		return 1;
	}
	
	op.addInput(&table);
	op.execute();
break;
}
/*** END SIGMA CONTOURS OP **/
/*** END SIGMA CONTOURS OP **/
/*** END SIGMA CONTOURS OP **/
case 20:
{
	iter =appParameters.find("help");
  if( iter != appParameters.end()){
  	VSPolarOp op;
  	op.printHelp();
  	return EXIT_SUCCESS;
  }
  
    
	VSPolarOp op;
	op.setParameters(appParameters);
	iter= appParameters.find("file");
	if(iter==appParameters.end()){
		std::cerr<<"No input file table is provided!"<<std::endl;
		return 1;
	}//close if
	std::stringstream sFilename(iter->second);
	appParameters.erase(iter);
	sFilename>>filename;
	if(filename.find(".bin") == std::string::npos) filename.append(".bin");
	VSTable table(filename);
	if(!table.tableExist()){
		std::cerr<<"No valid input file table is provided"<<std::endl;
		return 1;
	}
	
	op.addInput(&table);
	op.execute();
break;
}
case 21:
{
#ifdef WIN32
	std::cerr<<"Operation not allowed on this platform";
	return EXIT_SUCCESS;
#else
	iter =appParameters.find("help");
	if( iter != appParameters.end())
	{
		std::cout<<"VisIVOFilters --op AHFstep --input AmigaParameterFile"<<std::endl;
		return EXIT_SUCCESS;
	}
	iter= appParameters.find("input");
	if(iter==appParameters.end()){
		std::cerr<<"No input Amiga Parameter File is provided!"<<std::endl;
		return 1;
	}//close if
	std::stringstream sFilename(iter->second);

	std::string fileInput;
	sFilename>>fileInput;
	char *parFile=new char[fileInput.size()+1];
	parFile[fileInput.size()]=0;
	memcpy(parFile,fileInput.c_str(),fileInput.size());
	main_AHFstep(parFile,&argc,&argv);
	delete [] parFile;
#endif
break;
}
case 22:
{
	iter =appParameters.find("help");
  if( iter != appParameters.end()){
  	VSvbt2ahf op;
  	op.printHelp();
  	return EXIT_SUCCESS;
  }
  
    
	VSvbt2ahf op;
	op.setParameters(appParameters);
	iter= appParameters.find("file");
	if(iter==appParameters.end()){
		std::cerr<<"No input file table is provided!"<<std::endl;
		return 1;
	}//close if
	std::stringstream sFilename(iter->second);
	appParameters.erase(iter);
	sFilename>>filename;
	if(filename.find(".bin") == std::string::npos) filename.append(".bin");
	VSTable table(filename);
	if(!table.tableExist()){
		std::cerr<<"No valid input file table is provided"<<std::endl;
		return 1;
	}
	
	op.addInput(&table);
	op.execute();
break;
}
case 23:
{
	iter =appParameters.find("help");
  if( iter != appParameters.end()){
  	VSGrid2PointDistr op;
  	op.printHelp();
  	return EXIT_SUCCESS;
  }
  
    
	VSGrid2PointDistr op;
	op.setParameters(appParameters);
	iter= appParameters.find("file");
	if(iter==appParameters.end()){
		std::cerr<<"No input file table is provided!"<<std::endl;
		return 1;
	}//close if
	std::stringstream sFilename(iter->second);
	appParameters.erase(iter);
	sFilename>>filename;
	if(filename.find(".bin") == std::string::npos) filename.append(".bin");
	VSTable table(filename);
	if(!table.tableExist()){
		std::cerr<<"No valid input file table is provided"<<std::endl;
		return 1;
	}
	
	op.addInput(&table);
	op.execute();
break;
}
/******* Swap */
case 24:
{
	iter =appParameters.find("help");
  if( iter != appParameters.end()){
  	VSSwapOp op;
  	op.printHelp();
  	return EXIT_SUCCESS;
  }
  
    
	VSSwapOp op;
	op.setParameters(appParameters);
	iter= appParameters.find("file");
	if(iter==appParameters.end()){
		std::cerr<<"No input file table is provided!"<<std::endl;
		return 1;
	}//close if
	std::stringstream sFilename(iter->second);
	appParameters.erase(iter);
	sFilename>>filename;
	if(filename.find(".bin") == std::string::npos) filename.append(".bin");
	VSTable table(filename);
	if(!table.tableExist()){
		std::cerr<<"No valid input file table is provided"<<std::endl;
		return 1;
	}
	
	op.addInput(&table);
	op.execute();
break;
}
/******* ExtractList */
case 25:
{
	iter =appParameters.find("help");
  if( iter != appParameters.end()){
  	VSExtractListRowsOp op;
  	op.printHelp();
  	return EXIT_SUCCESS;
  }
  
    
	VSExtractListRowsOp op;
	op.setParameters(appParameters);
	iter= appParameters.find("file");
	if(iter==appParameters.end()){
		std::cerr<<"No input file table is provided!"<<std::endl;
		return 1;
	}//close if
	std::stringstream sFilename(iter->second);
	appParameters.erase(iter);
	sFilename>>filename;
	if(filename.find(".bin") == std::string::npos) filename.append(".bin");
	VSTable table(filename);
	if(!table.tableExist()){
		std::cerr<<"No valid input file table is provided"<<std::endl;
		return 1;
	}
	
	op.addInput(&table);
	op.execute();
break;
}
/******* AddId */
case 26:
{
	iter =appParameters.find("help");
  if( iter != appParameters.end()){
  	VSAddIdOp op;
  	op.printHelp();
  	return EXIT_SUCCESS;
  }
  
    
	VSAddIdOp op;
	op.setParameters(appParameters);
	iter= appParameters.find("file");
	if(iter==appParameters.end()){
		std::cerr<<"No input file table is provided!"<<std::endl;
		return 1;
	}//close if
	std::stringstream sFilename(iter->second);
	appParameters.erase(iter);
	sFilename>>filename;
	if(filename.find(".bin") == std::string::npos) filename.append(".bin");
	VSTable table(filename);
	if(!table.tableExist()){
		std::cerr<<"No valid input file table is provided"<<std::endl;
		return 1;
	}
	
	op.addInput(&table);
	op.execute();
break;
}
/******* AhfHaloId */
case 27:
{
	iter =appParameters.find("help");
  if( iter != appParameters.end()){
  	VSAhfHaloListOp op;
  	op.printHelp();
  	return EXIT_SUCCESS;
  }
  
    
	VSAhfHaloListOp op;
	op.setParameters(appParameters);
	iter= appParameters.find("file");
	if(iter==appParameters.end()){
		std::cerr<<"No input file table is provided!"<<std::endl;
		return 1;
	}//close if
	std::stringstream sFilename(iter->second);
	appParameters.erase(iter);
	sFilename>>filename;
	if(filename.find(".bin") == std::string::npos) filename.append(".bin");
	VSTable table(filename);
	if(!table.tableExist()){
		std::cerr<<"No valid input file table is provided"<<std::endl;
		return 1;
	}
	
	op.addInput(&table);
	op.execute();
break;
}
/******* ChangeColName*/
case 28:
{
	iter =appParameters.find("help");
  if( iter != appParameters.end()){
  	VSChangeColNameop op;
  	op.printHelp();
  	return EXIT_SUCCESS;
  }
  
    
	VSChangeColNameop op;
	op.setParameters(appParameters);
	iter= appParameters.find("file");
	if(iter==appParameters.end()){
		std::cerr<<"No input file table is provided!"<<std::endl;
		return 1;
	}//close if
	std::stringstream sFilename(iter->second);
	appParameters.erase(iter);
	sFilename>>filename;
	if(filename.find(".bin") == std::string::npos) filename.append(".bin");
	VSTable table(filename);
	if(!table.tableExist()){
		std::cerr<<"No valid input file table is provided"<<std::endl;
		return 1;
	}
	
	op.addInput(&table);
	op.execute();
break;
}
/******* AhfHaloGalaxyExt */
case 29:
{
	iter =appParameters.find("help");
  if( iter != appParameters.end()){
  	VSAhfHaloGalaxyExtOp op;
  	op.printHelp();
  	return EXIT_SUCCESS;
  }
  
    
	VSAhfHaloGalaxyExtOp op;
	op.setParameters(appParameters);
	iter= appParameters.find("file");
	if(iter==appParameters.end()){
		std::cerr<<"No input file table is provided!"<<std::endl;
		return 1;
	}//close if
	std::stringstream sFilename(iter->second);
	appParameters.erase(iter);
	sFilename>>filename;
	if(filename.find(".bin") == std::string::npos) filename.append(".bin");
	VSTable table(filename);
	if(!table.tableExist()){
		std::cerr<<"No valid input file table is provided"<<std::endl;
		return 1;
	}
	
	op.addInput(&table);
	op.execute();
break;
}
/******* SplitTable */
case 30:
{
	iter =appParameters.find("help");
  if( iter != appParameters.end()){
  	VSSplitTableOp op;
  	op.printHelp();
  	return EXIT_SUCCESS;
  }
  
    
	VSSplitTableOp op;
	op.setParameters(appParameters);
	iter= appParameters.find("file");
	if(iter==appParameters.end()){
		std::cerr<<"No input file table is provided!"<<std::endl;
		return 1;
	}//close if
	std::stringstream sFilename(iter->second);
	appParameters.erase(iter);
	sFilename>>filename;
	if(filename.find(".bin") == std::string::npos) filename.append(".bin");
	VSTable table(filename);
	if(!table.tableExist()){
		std::cerr<<"No valid input file table is provided"<<std::endl;
		return 1;
	}
	
	op.addInput(&table);
	op.execute();
break;
}
/******* wrvotrable */
case 31:
{
	iter =appParameters.find("help");
  if( iter != appParameters.end()){
  	VSWriteVotableOp op;
  	op.printHelp();
  	return EXIT_SUCCESS;
  }
  
    
	VSWriteVotableOp op;
	op.setParameters(appParameters);
	iter= appParameters.find("file");
	if(iter==appParameters.end()){
		std::cerr<<"No input file table is provided!"<<std::endl;
		return 1;
	}//close if
	std::stringstream sFilename(iter->second);
	appParameters.erase(iter);
	sFilename>>filename;
	if(filename.find(".bin") == std::string::npos) filename.append(".bin");
	VSTable table(filename);
	if(!table.tableExist()){
		std::cerr<<"No valid input file table is provided"<<std::endl;
		return 1;
	}
	
	op.addInput(&table);
	op.execute();
break;
}
/******* cut */
case 32:
{
	iter =appParameters.find("help");
  if( iter != appParameters.end()){
  	VSCutOp op;
  	op.printHelp();
  	return EXIT_SUCCESS;
  }
  
    
	VSCutOp op;
	op.setParameters(appParameters);
	iter= appParameters.find("file");
	if(iter==appParameters.end()){
		std::cerr<<"No input file table is provided!"<<std::endl;
		return 1;
	}//close if
	std::stringstream sFilename(iter->second);
	appParameters.erase(iter);
	sFilename>>filename;
	if(filename.find(".bin") == std::string::npos) filename.append(".bin");
	VSTable table(filename);
	if(!table.tableExist()){
		std::cerr<<"No valid input file table is provided"<<std::endl;
		return 1;
	}
	
	op.addInput(&table);
	op.execute();
break;
}

/******* vollimit */
case 33:
{
	iter =appParameters.find("help");
  if( iter != appParameters.end()){
  	VSVolLimitOp op;
  	op.printHelp();
  	return EXIT_SUCCESS;
  }
  
    
	VSVolLimitOp op;
	op.setParameters(appParameters);
	iter= appParameters.find("file");
	if(iter==appParameters.end()){
		std::cerr<<"No input file table is provided!"<<std::endl;
		return 1;
	}//close if
	std::stringstream sFilename(iter->second);
	appParameters.erase(iter);
	sFilename>>filename;
	if(filename.find(".bin") == std::string::npos) filename.append(".bin");
	VSTable table(filename);
	if(!table.tableExist()){
		std::cerr<<"No valid input file table is provided"<<std::endl;
		return 1;
	}
	
	op.addInput(&table);
	op.execute();
break;
}
/******* include */
case 34:
{
	iter =appParameters.find("help");
  if( iter != appParameters.end()){
  	VSIncludeOp op;
  	op.printHelp();
  	return EXIT_SUCCESS;
  }
VSIncludeOp op;
op.setParameters(appParameters);
iter=appParameters.find("file");
if(iter==appParameters.end())
{
	std::cerr<<"No input file table is provided!"<<std::endl;
	return 1;
}
std::stringstream sFilename(iter->second);
sFilename>>filename;
if(filename.find(".bin") == std::string::npos) filename.append(".bin");
VSTable table(filename);
if(!table.tableExist())
{
	std::cerr<<"No valid input file table is provided"<<std::endl;
	return 1;
}
op.addInput(&table);
op.execute();
break;
}


/*** Default **/
default:
{

    std::cerr <<"No valid operation was given"<<std::endl;
    std::cerr<<"VisIVOFilters version 1.2  January 24th 2011 "<<std::endl<<std::endl;
    std::cerr <<"Syntax1: VisIVOFilters --op operation  [PARAMETERS] [--help]"<<std::endl;
    std::cerr <<"Syntax2: VisIVOFilters parameterFile"<<std::endl;
    std::cerr <<"An operation code is expected: randomizer selcolumns merge append selfield mathop decimator extraction visualop showtable statistic pointdistribute pointproperty coarsevolume extractsubvolume interpolate module sigmacontours cartesian2polar AHFstep vbt2ahf grid2point extractlist addId ahfhalolist ahfhalogalaxyext changecolname splittable wrvotable include"<<std::endl;

    return 1;


}
/*** END Default  OP **/
}
return EXIT_SUCCESS;
}

//#ifdef HAVE_CONFIG_H
//#include <config.h>
//#endif

