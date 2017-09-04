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
#include "vscreatepath.h"
#include "parametersparser.h"
const double VSCreatePathUT::INVALID_CAM = -123456789.31;

//---------------------------------------------------------------------
VSCreatePathUT::VSCreatePathUT()
//---------------------------------------------------------------------
{
 m_azimuthFrom=0.0;
 m_azimuthTo=0.0;
 m_elevationFrom=0.0;
 m_elevationTo=0.0;
 m_zoomFrom=1.0;
 m_zoomTo=1.0;
 m_length=10;
 m_framePerSecond=10;
 m_zStepFrame=0.2;
 m_zoomend=false;
 m_outFile="cycle.par";
 m_cameraPoisitionStart[0]=INVALID_CAM;
 m_cameraPoisitionStart[1]=INVALID_CAM;
 m_cameraPoisitionStart[2]=INVALID_CAM;
 m_cameraPoisitionEnd[0]=INVALID_CAM;
 m_cameraPoisitionEnd[1]=INVALID_CAM;
 m_cameraPoisitionEnd[2]=INVALID_CAM;
 m_cameraFocalPointStart[0]=INVALID_CAM;
 m_cameraFocalPointStart[1]=INVALID_CAM;
 m_cameraFocalPointStart[2]=INVALID_CAM;
 m_cameraFocalPointEnd[0]=INVALID_CAM;
 m_cameraFocalPointEnd[1]=INVALID_CAM;
 m_cameraFocalPointEnd[2]=INVALID_CAM;
 m_cameraRollStart=INVALID_CAM;
 m_cameraRollEnd=INVALID_CAM;
}
//---------------------------------------------------------------------
void VSCreatePathUT::printHelp()
//---------------------------------------------------------------------
{
	std::cout<<"Create a file with 4,7 or 8 columns used in VisIVOViewer --cycle option"<<std::endl<<std::endl;;

	std::cout<<"Usage: VisIVOUtils --op createpath --type value [--azimuth from to] [--elevation from to] [--zoom from to]  [--zoomend [stepframe]] [--campos from to] [--camfp from to] [--camroll from to] [--framesec value] [--length value]  [--out filename] [--help]"<<std::endl;

	std::cout<<"Example: VisIVOUtils --op createpath --type 1 --azimuth 0.0 60. --elevation 0.0 10.0 --zoom 1.0 1.5  --zoomend --camfp 35 35 35  0 0 0 --length 20 --out my_cycle.par"<<std::endl<<std::endl;

	std::cout<<"Note: "<<std::endl;
	std::cout<<"--type . 0 Create path for azimuth, elevation, zoom and roll. Default value."<<std::endl; 
	std::cout<<"         1 Create path for azimuth, elevation, zoom, focal point and roll"<<std::endl;
	std::cout<<"         2 Create path for zoom, camera position, focal point and roll"<<std::endl;
	std::cout<<"         3 Create path for azimuth, elevation, zoom, camera position, focal point and roll"<<std::endl;
	std::cout<<"--azimuth Movement from to. Default values 0.0 and 0.0."<<std::endl;
	std::cout<<"--elevation Movement from to. Default values 0.0 and 0.0."<<std::endl;
	std::cout<<"--zoom Zoom from to. Default values 1.0 and 1.0."<<std::endl;
	std::cout<<"--zoomend The zoom is given at the end. The value stepframe represent the step for zooming . Default stepframe is 0.2 If this option is given priority with zoom will be ignored. The final zooming is added to global the length."<<std::endl;
	std::cout<<"--campos Movement from to. Three vale for starting point and three value for ending point are expected."<<std::endl;
	std::cout<<"--camfp Movement from to. Three vale for starting point and three value for ending point are expected."<<std::endl;
	std::cout<<"--camroll movement from to. To values starting and ending degrees are extpected"<<std::endl;
	std::cout<<"--framesec Number of frame values for each second. Default value is 10 sec."<<std::endl;
	std::cout<<"--length value in seconds. Default value is 10 sec."<<std::endl;
	std::cout<<"--out output filename. Default filename cycle.par. The file is opened in append mode"<<std::endl;
	std::cout<<"--help produce this output "<<std::endl;

	return;

}

//---------------------------------------------------------------------
bool VSCreatePathUT::execute()
//---------------------------------------------------------------------
{	
  bool path_enabled=false;
  bool camposIs=false;
  bool camfpIs=false;
  bool camrollIs=false;
  double camposStep[3];
  double camfpStep[3];
  double camrollStep;
  int type;
  type=getParameterAsInt("type");
  if(type!=0 && type!=1 && type!=2 && type!=3) type=0;
  
  if(!(getParameterAsString("azimuth").empty() || getParameterAsString("azimuth")=="unknown" ))
  {
	std::stringstream sstmp;
	sstmp.str(getParameterAsString("azimuth"));
	sstmp>>m_azimuthFrom;
	sstmp>>m_azimuthTo;
  }
  if(!(getParameterAsString("elevation").empty() || getParameterAsString("elevation")=="unknown" ))
  {
	std::stringstream sstmp;
	sstmp.str(getParameterAsString("elevation"));
	sstmp>>m_elevationFrom;
	sstmp>>m_elevationTo;
	if(m_elevationFrom>90.0 || m_elevationFrom<-90.0 || m_elevationTo>90 || m_elevationTo<-90.0)
		std::cerr<<"Elevation invalid range "<< m_elevationFrom<<" "<< m_elevationTo<<" Values are set to the interval extreme [-90,90]."<<std::endl; 
	if(m_elevationFrom>90.0)
		m_elevationFrom=90.0;
	if(m_elevationFrom<-90.0)
		m_elevationFrom=-90.0;
	if(m_elevationTo>90.0)
		m_elevationTo=90.0;
	if(m_elevationTo<-90.0)
		m_elevationTo=-90.0;
  }
  if(!(getParameterAsString("zoom").empty() || getParameterAsString("zoom")=="unknown" ))
  {
	std::stringstream sstmp;
	sstmp.str(getParameterAsString("zoom"));
	sstmp>>m_zoomFrom;
	sstmp>>m_zoomTo;
  }
  if(!(getParameterAsString("zoomend").empty())) 
	m_zoomend=true;
  if(!(getParameterAsString("zoomend").empty() || getParameterAsString("zoomend")=="unknown" ))
  {
	std::stringstream sstmp;
	sstmp.str(getParameterAsString("zoomend"));
	sstmp>>m_zStepFrame;
//	std::clog<<m_zStepFrame<<std::endl;
  }
  if(!(getParameterAsString("campos").empty() || getParameterAsString("campos")=="unknown" ))
  {
	std::stringstream sstmp;
	sstmp.str(getParameterAsString("campos"));
	sstmp>>m_cameraPoisitionStart[0];
	sstmp>>m_cameraPoisitionStart[1];
	sstmp>>m_cameraPoisitionStart[2];
	sstmp>>m_cameraPoisitionEnd[0];
	sstmp>>m_cameraPoisitionEnd[1];
	sstmp>>m_cameraPoisitionEnd[2];
	path_enabled=true;
	camposIs=true;
	for(int i=0;i<3;i++)
	{
	  if(m_cameraPoisitionStart[i]==INVALID_CAM) camposIs=false;
	  if(m_cameraPoisitionEnd[i]==INVALID_CAM) camposIs=false;
	}
  }
  if(!(getParameterAsString("camfp").empty() || getParameterAsString("camfp")=="unknown" ))
  {
	std::stringstream sstmp;
	sstmp.str(getParameterAsString("camfp"));
	sstmp>>m_cameraFocalPointStart[0];
	sstmp>>m_cameraFocalPointStart[1];
	sstmp>>m_cameraFocalPointStart[2];
	sstmp>>m_cameraFocalPointEnd[0];
	sstmp>>m_cameraFocalPointEnd[1];
	sstmp>>m_cameraFocalPointEnd[2];
	path_enabled=true;
	camfpIs=true;
	for(int i=0;i<3;i++)
	{
	  if(m_cameraFocalPointStart[i]==INVALID_CAM) camfpIs=false;
	  if(m_cameraFocalPointEnd[i]==INVALID_CAM) camfpIs=false;
	}
  }
  if(!(getParameterAsString("camroll").empty() || getParameterAsString("camfp")=="unknown" ))
  {
	std::stringstream sstmp;
	sstmp.str(getParameterAsString("camroll"));
	sstmp>>m_cameraRollStart;
	sstmp>>m_cameraRollEnd;
	path_enabled=true;
	camrollIs=true;
	if(m_cameraRollStart==INVALID_CAM) camrollIs=false;
	if(m_cameraRollEnd==INVALID_CAM) camrollIs=false;

  }
  if(!(getParameterAsString("framesec").empty() || getParameterAsString("framesec")=="unknown" ))
  {
	std::stringstream sstmp;
	sstmp.str(getParameterAsString("framesec"));
	sstmp>>m_framePerSecond;
  }
  if(!(getParameterAsString("length").empty() || getParameterAsString("length")=="unknown" ))
  {
	std::stringstream sstmp;
	sstmp.str(getParameterAsString("length"));
	sstmp>>m_length;
  }
  if(!(getParameterAsString("out").empty() || getParameterAsString("out")=="unknown" ))
  {
	m_outFile=getParameterAsString("out");
	
  }
  int numberOfFrames=m_length*m_framePerSecond;
  if(numberOfFrames<=0)
  {
  	std::cerr<<"Wrong parameters: length= "<<m_length<<" Frame for seconds= "<< m_framePerSecond<<std::endl;
	return false;
  }
  float azimuthStep=(m_azimuthTo-m_azimuthFrom)/numberOfFrames;
  float elevationStep=(m_elevationTo-m_elevationFrom)/numberOfFrames;
  if(path_enabled)
  {
    camposStep[0]=(m_cameraPoisitionEnd[0]-m_cameraPoisitionStart[0])/numberOfFrames;
    camposStep[1]=(m_cameraPoisitionEnd[1]-m_cameraPoisitionStart[1])/numberOfFrames;
    camposStep[2]=(m_cameraPoisitionEnd[2]-m_cameraPoisitionStart[2])/numberOfFrames;
    camfpStep[0]=(m_cameraFocalPointEnd[0]-m_cameraFocalPointStart[0])/numberOfFrames;
    camfpStep[1]=(m_cameraFocalPointEnd[1]-m_cameraFocalPointStart[1])/numberOfFrames;
    camfpStep[2]=(m_cameraFocalPointEnd[2]-m_cameraFocalPointStart[2])/numberOfFrames;
    camrollStep=(m_cameraRollEnd-m_cameraRollStart)/numberOfFrames;
  }
  float zoomStep=0;
  if(!m_zoomend) 
	zoomStep=(m_zoomTo-m_zoomFrom)/numberOfFrames;
  std::ofstream outFile;
  outFile.open(m_outFile.c_str(), std::ios::app);
  float azimuthValue=m_azimuthFrom;
  float elevationValue=m_elevationFrom;
  float zoomValue=m_zoomFrom;
  
  double camposValue[3];
  double camfpValue[3];
  double camrollValue;
  
  camposValue[0]=m_cameraPoisitionStart[0];
  camposValue[1]=m_cameraPoisitionStart[1];
  camposValue[2]=m_cameraPoisitionStart[2];
  camfpValue[0]=m_cameraFocalPointStart[0];
  camfpValue[1]=m_cameraFocalPointStart[1];
  camfpValue[2]=m_cameraFocalPointStart[2];
  camrollValue=m_cameraRollStart;
  
switch(type)
{
/*** Randomizer OP **/
case 0:
{
   for(int i=0;i<numberOfFrames;i++)
   {
	outFile<<azimuthValue<<"  "<<elevationValue<<"  "<<zoomValue;
	  if(camrollIs)
	    outFile<<" "<<camrollValue;
	  else
	    outFile<<" NULL";
	
	outFile<<std::endl;
	azimuthValue=azimuthValue+azimuthStep; 
	elevationValue=elevationValue+elevationStep; 
	zoomValue=zoomValue+zoomStep; 
	camrollValue=camrollValue+camrollStep;
   }
   if(!m_zoomend)
   {  
       outFile<<m_azimuthTo<<"  "<<m_elevationTo<<"  "<<m_zoomTo;
	  if(camrollIs)
	    outFile<<" "<<m_cameraRollEnd;
	  else
	    outFile<<" NULL";
       outFile<<std::endl;
   } else {  
       outFile<<m_azimuthTo<<"  "<<m_elevationTo<<"  "<<m_zoomFrom;
	  
	  if(camrollIs)
	    outFile<<" "<<m_cameraRollEnd;
	  else
	    outFile<<" NULL";
       outFile<<std::endl;
  
   float numbeOfZoomSteps=(m_zoomTo-m_zoomFrom)/m_zStepFrame;
    if(numbeOfZoomSteps<0)
	numbeOfZoomSteps=-1*numbeOfZoomSteps;
    if(m_zoomFrom>m_zoomTo)
	m_zStepFrame=-1.0*m_zStepFrame;
    for(int i=0;i<numbeOfZoomSteps;i++)
    {
	zoomValue=zoomValue+m_zStepFrame;
        outFile<<m_azimuthTo<<"  "<<m_elevationTo<<"  "<<zoomValue;
	  if(camrollIs)
	    outFile<<" "<<m_cameraRollEnd;
	  else
	    outFile<<" NULL";
       outFile<<std::endl;
    }
    }
   break;
} // case=0
case 1:
{

   for(int i=0;i<numberOfFrames;i++)
   {
	outFile<<azimuthValue<<"  "<<elevationValue<<"  "<<zoomValue;
	  if(camfpIs) 
	    outFile<<" "<<camfpValue[0]<<" "<<camfpValue[1]<<" "<<camfpValue[2];
	  else
	    outFile<<" NULL NULL NULL";
	  
	  if(camrollIs)
	    outFile<<" "<<camrollValue;
	  else
	    outFile<<" NULL";
	  
	  camfpValue[0]=camfpValue[0]+camfpStep[0];
	  camfpValue[1]=camfpValue[1]+camfpStep[1];
	  camfpValue[2]=camfpValue[2]+camfpStep[2];
	  camrollValue=camrollValue+camrollStep;
	
	outFile<<std::endl;
	azimuthValue=azimuthValue+azimuthStep; 
	elevationValue=elevationValue+elevationStep; 
	zoomValue=zoomValue+zoomStep; 
  }
   if(!m_zoomend)
   {  
       outFile<<m_azimuthTo<<"  "<<m_elevationTo<<"  "<<m_zoomTo;
	  if(camfpIs) 
	    outFile<<" "<<m_cameraFocalPointEnd[0]<<" "<<m_cameraFocalPointEnd[1]<<" "<<m_cameraFocalPointEnd[2];
	  else
	    outFile<<" NULL NULL NULL";
	  
	  if(camrollIs)
	    outFile<<" "<<m_cameraRollEnd;
	  else
	    outFile<<" NULL";

      
       
       outFile<<std::endl;
  }
   else
   {  
       outFile<<m_azimuthTo<<"  "<<m_elevationTo<<"  "<<m_zoomFrom;
	  if(camfpIs) 
	    outFile<<" "<<m_cameraFocalPointEnd[0]<<" "<<m_cameraFocalPointEnd[1]<<" "<<m_cameraFocalPointEnd[2];
	  else
	    outFile<<" NULL NULL NULL";
	  
	  if(camrollIs)
	    outFile<<" "<<m_cameraRollEnd;
	  else
	    outFile<<" NULL";
       outFile<<std::endl;
  }
   if(m_zoomend)
   {
    float numbeOfZoomSteps=(m_zoomTo-m_zoomFrom)/m_zStepFrame;
    if(numbeOfZoomSteps<0)
	numbeOfZoomSteps=-1*numbeOfZoomSteps;
    if(m_zoomFrom>m_zoomTo)
	m_zStepFrame=-1.0*m_zStepFrame;
    for(int i=0;i<numbeOfZoomSteps;i++)
    {
	zoomValue=zoomValue+m_zStepFrame;
        outFile<<m_azimuthTo<<"  "<<m_elevationTo<<"  "<<zoomValue;
	    
	  if(camfpIs) 
	    outFile<<" "<<m_cameraFocalPointEnd[0]<<" "<<m_cameraFocalPointEnd[1]<<" "<<m_cameraFocalPointEnd[2];
	  else
	    outFile<<" NULL NULL NULL";
	  
	  if(camrollIs)
	    outFile<<" "<<m_cameraRollEnd;
	  else
	    outFile<<" NULL";
       outFile<<std::endl;
    }
   }
   break;
} // case=1
case 2:
{

   for(int i=0;i<numberOfFrames;i++)
   {
	outFile<<zoomValue;
	  if(camposIs) 
	    outFile<<" "<<camposValue[0]<<" "<<camposValue[1]<<" "<<camposValue[2];
	  else
	    outFile<<" NULL NULL NULL";
	    
	  if(camfpIs) 
	    outFile<<" "<<camfpValue[0]<<" "<<camfpValue[1]<<" "<<camfpValue[2];
	  else
	    outFile<<" NULL NULL NULL";
	  
	  if(camrollIs)
	    outFile<<" "<<camrollValue;
	  else
	    outFile<<" NULL";
	  
	  camposValue[0]=camposValue[0]+camposStep[0];
	  camposValue[1]=camposValue[1]+camposStep[1];
	  camposValue[2]=camposValue[2]+camposStep[2];
	  camfpValue[0]=camfpValue[0]+camfpStep[0];
	  camfpValue[1]=camfpValue[1]+camfpStep[1];
	  camfpValue[2]=camfpValue[2]+camfpStep[2];
	  camrollValue=camrollValue+camrollStep;
	
	outFile<<std::endl;
	azimuthValue=azimuthValue+azimuthStep; 
	elevationValue=elevationValue+elevationStep; 
	zoomValue=zoomValue+zoomStep; 
  }
   if(!m_zoomend)
   {  
       outFile<<m_zoomTo;
       	  if(camposIs) 
	    outFile<<" "<<m_cameraPoisitionEnd[0]<<" "<<m_cameraPoisitionEnd[1]<<" "<<m_cameraPoisitionEnd[2];
	  else
	    outFile<<" NULL NULL NULL";
	    
	  if(camfpIs) 
	    outFile<<" "<<m_cameraFocalPointEnd[0]<<" "<<m_cameraFocalPointEnd[1]<<" "<<m_cameraFocalPointEnd[2];
	  else
	    outFile<<" NULL NULL NULL";
	  
	  if(camrollIs)
	    outFile<<" "<<m_cameraRollEnd;
	  else
	    outFile<<" NULL";

        
       outFile<<std::endl;
  }
   else
   {  
       outFile<<m_zoomFrom;
          if(camposIs) 
	    outFile<<" "<<m_cameraPoisitionEnd[0]<<" "<<m_cameraPoisitionEnd[1]<<" "<<m_cameraPoisitionEnd[2];
	  else
	    outFile<<" NULL NULL NULL";
	    
	  if(camfpIs) 
	    outFile<<" "<<m_cameraFocalPointEnd[0]<<" "<<m_cameraFocalPointEnd[1]<<" "<<m_cameraFocalPointEnd[2];
	  else
	    outFile<<" NULL NULL NULL";
	  
	  if(camrollIs)
	    outFile<<" "<<m_cameraRollEnd;
	  else
	    outFile<<" NULL";
       outFile<<std::endl;
  }
   if(m_zoomend)
   {
    float numbeOfZoomSteps=(m_zoomTo-m_zoomFrom)/m_zStepFrame;
    if(numbeOfZoomSteps<0)
	numbeOfZoomSteps=-1*numbeOfZoomSteps;
    if(m_zoomFrom>m_zoomTo)
	m_zStepFrame=-1.0*m_zStepFrame;
    for(int i=0;i<numbeOfZoomSteps;i++)
    {
	zoomValue=zoomValue+m_zStepFrame;
        outFile<<zoomValue;
          if(camposIs) 
	    outFile<<" "<<m_cameraPoisitionEnd[0]<<" "<<m_cameraPoisitionEnd[1]<<" "<<m_cameraPoisitionEnd[2];
	  else
	    outFile<<" NULL NULL NULL";
	    
	  if(camfpIs) 
	    outFile<<" "<<m_cameraFocalPointEnd[0]<<" "<<m_cameraFocalPointEnd[1]<<" "<<m_cameraFocalPointEnd[2];
	  else
	    outFile<<" NULL NULL NULL";
	  
	  if(camrollIs)
	    outFile<<" "<<m_cameraRollEnd;
	  else
	    outFile<<" NULL";
       
       outFile<<std::endl;
    }
   }
    break;
} // case=2
case 3:
{

   for(int i=0;i<numberOfFrames;i++)
   {
       outFile<<azimuthValue<<"  "<<elevationValue<<"  "<<zoomValue;
	  if(camposIs) 
	    outFile<<" "<<camposValue[0]<<" "<<camposValue[1]<<" "<<camposValue[2];
	  else
	    outFile<<" NULL NULL NULL";
	    
	  if(camfpIs) 
	    outFile<<" "<<camfpValue[0]<<" "<<camfpValue[1]<<" "<<camfpValue[2];
	  else
	    outFile<<" NULL NULL NULL";
	  
	  if(camrollIs)
	    outFile<<" "<<camrollValue;
	  else
	    outFile<<" NULL";
	  
	  camposValue[0]=camposValue[0]+camposStep[0];
	  camposValue[1]=camposValue[1]+camposStep[1];
	  camposValue[2]=camposValue[2]+camposStep[2];
	  camfpValue[0]=camfpValue[0]+camfpStep[0];
	  camfpValue[1]=camfpValue[1]+camfpStep[1];
	  camfpValue[2]=camfpValue[2]+camfpStep[2];
	  camrollValue=camrollValue+camrollStep;
	
	outFile<<std::endl;
	azimuthValue=azimuthValue+azimuthStep; 
	elevationValue=elevationValue+elevationStep; 
	zoomValue=zoomValue+zoomStep; 
  }
   if(!m_zoomend)
   {  
       outFile<<m_azimuthTo<<"  "<<m_elevationTo<<"  "<<m_zoomTo;
       	  if(camposIs) 
	    outFile<<" "<<m_cameraPoisitionEnd[0]<<" "<<m_cameraPoisitionEnd[1]<<" "<<m_cameraPoisitionEnd[2];
	  else
	    outFile<<" NULL NULL NULL";
	    
	  if(camfpIs) 
	    outFile<<" "<<m_cameraFocalPointEnd[0]<<" "<<m_cameraFocalPointEnd[1]<<" "<<m_cameraFocalPointEnd[2];
	  else
	    outFile<<" NULL NULL NULL";
	  
	  if(camrollIs)
	    outFile<<" "<<m_cameraRollEnd;
	  else
	    outFile<<" NULL";

        
       outFile<<std::endl;
  }
   else
   {  
       outFile<<m_azimuthTo<<"  "<<m_elevationTo<<"  "<<m_zoomFrom;
          if(camposIs) 
	    outFile<<" "<<m_cameraPoisitionEnd[0]<<" "<<m_cameraPoisitionEnd[1]<<" "<<m_cameraPoisitionEnd[2];
	  else
	    outFile<<" NULL NULL NULL";
	    
	  if(camfpIs) 
	    outFile<<" "<<m_cameraFocalPointEnd[0]<<" "<<m_cameraFocalPointEnd[1]<<" "<<m_cameraFocalPointEnd[2];
	  else
	    outFile<<" NULL NULL NULL";
	  
	  if(camrollIs)
	    outFile<<" "<<m_cameraRollEnd;
	  else
	    outFile<<" NULL";
       outFile<<std::endl;
  }
   if(m_zoomend)
   {
    float numbeOfZoomSteps=(m_zoomTo-m_zoomFrom)/m_zStepFrame;
    if(numbeOfZoomSteps<0)
	numbeOfZoomSteps=-1*numbeOfZoomSteps;
    if(m_zoomFrom>m_zoomTo)
	m_zStepFrame=-1.0*m_zStepFrame;
    for(int i=0;i<numbeOfZoomSteps;i++)
    {
	zoomValue=zoomValue+m_zStepFrame;
        outFile<<m_azimuthTo<<"  "<<m_elevationTo<<"  "<<zoomValue;
          if(camposIs) 
	    outFile<<" "<<m_cameraPoisitionEnd[0]<<" "<<m_cameraPoisitionEnd[1]<<" "<<m_cameraPoisitionEnd[2];
	  else
	    outFile<<" NULL NULL NULL";
	    
	  if(camfpIs) 
	    outFile<<" "<<m_cameraFocalPointEnd[0]<<" "<<m_cameraFocalPointEnd[1]<<" "<<m_cameraFocalPointEnd[2];
	  else
	    outFile<<" NULL NULL NULL";
	  
	  if(camrollIs)
	    outFile<<" "<<m_cameraRollEnd;
	  else
	    outFile<<" NULL";
       
       outFile<<std::endl;
    }
   }
    break;
} // case=3

  
  
}


  outFile.close();
  return true;
}
