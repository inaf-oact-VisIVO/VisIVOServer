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
#include "vscreategenericslices.h"
#include "parametersparser.h"

//---------------------------------------------------------------------
VSCreateGenericSlicesUT::VSCreateGenericSlicesUT()
//---------------------------------------------------------------------
{
  m_movedown=false;
  m_size=1; 
  m_outFile="cycle.par";
}
//---------------------------------------------------------------------
void VSCreateGenericSlicesUT::printHelp()
//---------------------------------------------------------------------
{
	std::cout<<"Append or create a file with  six columns. The point position (plane point) is increased (decreased) of step_size for n steps."; 
	std::cout<<"The plane point is moved along the normal axis. The product step*size determines";
	std::cout<<" the movement of the plain point. If step*size is equal to 1, at the end";
	std::cout<<" the plane point will be at the same point of the normal point";
	std::cout<<std::endl<<std::endl;;

	std::cout<<"Usage: VisIVOUtils --op genericslices  --point x y z --normal x y z --step n [--size step_size]  [--movedown]  [--out filename] [--help]"<<std::endl;

	std::cout<<"Example: VisIVOUtils --op genericslices  --point 1 1 1  --normal 10 10 10 --size 0.05 --step 20  --out cyclefile"<<std::endl<<std::endl;

	std::cout<<"Note: "<<std::endl;
	std::cout<<"--point The three coordinates of a point in the plane."<<std::endl;
	std::cout<<"--normal The three coordinate fixing the normal axis to the plane"<<std::endl;
	std::cout<<"--step Number of generated new point positions along the normal axis"<<std::endl;
	std::cout<<"--size Value of increased (decreased) point coordinate. Default value 1"<<std::endl;
	std::cout<<"--movedown The plane point is moving in the opposte side of the normal point"<<std::endl;
	std::cout<<"--out output filename. Default filename cycle.par. The file is opened in append mode"<<std::endl;
	std::cout<<"--help produce this output "<<std::endl;

	return;

}

//---------------------------------------------------------------------
bool VSCreateGenericSlicesUT::execute()
//---------------------------------------------------------------------
{

  if(isParameterPresent("point"))
  {
	std::stringstream sstmp;
	sstmp.str(getParameterAsString("point"));
	sstmp>>m_point[0];
	sstmp>>m_point[1];
	sstmp>>m_point[2];
  } else
  {
    std::cerr<<"Error. The point option is not given."<<std::endl;
    return false;
  }
  if(isParameterPresent("normal"))
  {
	std::stringstream sstmp;
	sstmp.str(getParameterAsString("normal"));
	sstmp>>m_plane[0];
	sstmp>>m_plane[1];
	sstmp>>m_plane[2];
  } else
  {
    std::cerr<<"Error. The normal option is not given."<<std::endl;
    return false;
  }
  if(m_point[0]==m_plane[0] && m_point[1]==m_plane[1] && m_point[2]==m_plane[2])
  {
      std::cerr<<"Error. The normal and the point have the same vaules."<<std::endl;
       return false;
  
  }
  if(isParameterPresent("step"))
	m_step=getParameterAsInt("step");
  else
  {
    std::cerr<<"Error. The step option is not given."<<std::endl;
    return false;
  }
  if(m_step<=0)
  {
	std::cerr<<"Invalid step value: "<<m_step<<" Operation aborted."<<std::endl;
	return false;
  }
  if(isParameterPresent("size"))
	m_size=getParameterAsFloat("size");
  if(isParameterPresent("movedown"))
    m_movedown=true;
  
  

  if(!(getParameterAsString("out").empty() || getParameterAsString("out")=="unknown" ))
  {
	m_outFile=getParameterAsString("out");
	
  }
  float par[3];
  for(int i=0;i<3;i++) 
    par[i]=m_plane[i]-m_point[i];

  std::ofstream outFile;
  outFile.open(m_outFile.c_str(), std::ios::app);
  for(int i=0;i<=m_step;i++)
  {  
    outFile<<m_point[0]<<" "<<m_point[1]<<" "<<m_point[2]<<" ";
    outFile<<m_plane[0]<<" "<<m_plane[1]<<" "<<m_plane[2]<<std::endl;
    for(int j=0;j<3;j++)
    {
      if(m_movedown) m_point[j]=m_point[j]-par[j]*m_size;
      else m_point[j]=m_point[j]+par[j]*m_size;
    }  
  } 
  outFile.close();
  return true;
}
