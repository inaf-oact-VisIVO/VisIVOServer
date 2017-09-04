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
#include <iostream>
#include <fstream>
#include <sstream>
#include "vstable.h"
#include "vscreateslices.h"
#include "parametersparser.h"

//---------------------------------------------------------------------
VSCreateSlicesUT::VSCreateSlicesUT()
//---------------------------------------------------------------------
{
 m_sliceFrom=0;
 m_sliceTo=0;
 m_step=1;
 m_outFile="cycle.par";
}
//---------------------------------------------------------------------
void VSCreateSlicesUT::printHelp()
//---------------------------------------------------------------------
{
	std::cout<<"Append or create a file with 1 columns with the slice poisition in the volume table used in VisIVOViewer --cycle option for slice visualization"<<std::endl<<std::endl;;

	std::cout<<"Usage: VisIVOUtils --op slices  --pos from to [--xplane] [--yplane] [--zplane] [--step stepvalue] [--out filename] [--help] --file inputFile.bin"<<std::endl;

	std::cout<<"Example: VisIVOUtils --op slices  --pos  0 64 --step 1 --out my_cycle.par --file inputFile.bin"<<std::endl<<std::endl;

	std::cout<<"Note: "<<std::endl;
	std::cout<<"--pos Slice position (integer) from-to in the volume. Values outsides  the volume size are ignored."<<std::endl;
	std::cout<<"--xplane Set the direction x to be considered. Default is x"<<std::endl;
	std::cout<<"--yplane Set the direction y to be considered. Ignored if --xplane is given"<<std::endl;
	std::cout<<"--zplane Set the direction z to be considered. Ignored if --xplane or --yplane is given"<<std::endl;
	std::cout<<"--step Step increment (integer)  for slice position in the volume. Default value 1"<<std::endl;
	std::cout<<"--out output filename. Default filename cycle.par. The file is opened in append mode"<<std::endl;
	std::cout<<"--file Input volume table "<<std::endl;
	std::cout<<"--help produce this output "<<std::endl;

	return;

}

//---------------------------------------------------------------------
bool VSCreateSlicesUT::execute()
//---------------------------------------------------------------------
{
  bool fileFlag=true;
  if(isParameterPresent("file"))
  {
	m_inputFile=getParameterAsString("file");
  } else 
    fileFlag=false;
  if(m_inputFile.find(".bin") == std::string::npos && fileFlag)
	    m_inputFile.append(".bin");
  VSTable sliceTable;
  if(fileFlag)
  {  
    sliceTable.setLocator(m_inputFile);
    sliceTable.readHeader();
  }   
  
  
  if(!sliceTable.tableExist() && fileFlag)
  {
	std::cerr<<"Invalid filename "<<m_inputFile<<" Table does not exist. Operation aborted."<<std::endl;
	return false;
  }
  if(!sliceTable.getIsVolume() && fileFlag)
  {
	std::cerr<<"Table "<<m_inputFile<<" in not a volume. Operation aborted."<<std::endl;
	return false;
  }
  const unsigned int *tableCells;
  tableCells=new unsigned int[3];
  tableCells[0]==0;
  tableCells[1]==0;
  tableCells[2]==0;
  if(fileFlag) tableCells=sliceTable.getCellNumber();

  if(isParameterPresent("pos"))
  {
	std::stringstream sstmp;
	sstmp.str(getParameterAsString("pos"));
	sstmp>>m_sliceFrom;
	sstmp>>m_sliceTo;
  }
  if(m_sliceFrom<0 || m_sliceTo<0 )
  {
	std::cerr<<"Negative pos values: From "<<m_sliceFrom<<" To "<<m_sliceTo<<" Operation aborted."<<std::endl;
	return false;
  }
  if(isParameterPresent("step"))
	m_step=getParameterAsInt("step");
  if(m_step<=0)
  {
	std::cerr<<"Invalid step value: "<<m_step<<" Operation aborted."<<std::endl;
	return false;
  }
  int planeSelector=0;

  if(isParameterPresent("yplane") && !isParameterPresent("xplane")) planeSelector=1;
  if(isParameterPresent("zplane") && !isParameterPresent("yplane") && !isParameterPresent("xplane")) planeSelector=2;
  
  if(m_sliceTo>tableCells[planeSelector] && fileFlag)
  {
	 m_sliceTo=tableCells[planeSelector];
	 std::cerr<<"Invalid pos TO parameter. Lowered to the maximum size in the volume "<<m_sliceTo<<std::endl; 
  }
  if(m_sliceFrom>=m_sliceTo && fileFlag)
  {
	std::cerr<<"Invalid pos: From "<<m_sliceFrom<<" To "<<m_sliceTo<<" Operation aborted."<<std::endl;
	return false;
  }

  if(!(getParameterAsString("out").empty() || getParameterAsString("out")=="unknown" ))
  {
	m_outFile=getParameterAsString("out");
	
  }
  int outValue=m_sliceFrom;
  bool writeLast=false;
  std::ofstream outFile;
  outFile.open(m_outFile.c_str(), std::ios::app);
  while(true)
  {
        outFile<<outValue<<std::endl;
	if(outValue==m_sliceTo) writeLast=true;
	outValue+=m_step;
	if(outValue>m_sliceTo)
		break;
  }
  if(!writeLast)
	outFile<<m_sliceTo<<std::endl;
  outFile.close();
  return true;
}
