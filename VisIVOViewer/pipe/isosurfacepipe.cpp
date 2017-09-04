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
#include <sstream>

#include "isosurfacepipe.h"
#include "visivoutils.h"

#include <vtkGenericRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkImageData.h>
#include <vtkImageMathematics.h>
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkColorTransferFunction.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include <vtkVolumeRayCastIsosurfaceFunction.h>
#include <vtkVolumeRayCastMapper.h>
#include <vtkVolumeProperty.h>
#include <vtkPiecewiseFunction.h>
#include "vtkImageCast.h"


#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkContourFilter.h"

#include "vtkImageGaussianSmooth.h"

#include "vtkPolyDataNormals.h"
#include "vtkProperty.h"
// #include "vtkStructuredPointsWriter.h"

//---------------------------------------------------------------------
IsosurfacePipe::IsosurfacePipe ( VisIVOServerOptions options)
//---------------------------------------------------------------------
{
  m_visOpt=options;
  constructVTK();
  m_imageData               = vtkImageData::New();
  
}
//---------------------------------
IsosurfacePipe::~IsosurfacePipe()
//---------------------------------
{
  destroyVTK();

  m_imageData->Delete();
}
//---------------------------------
void IsosurfacePipe::destroyAll()
//---------------------------------
{  
  m_imageData->Delete();

}
//-----------------------------------------------------------------------------------
int IsosurfacePipe::createPipe ()
//------------------------------------------------------------------------------------
{
   
 
  m_range[0]= 0;
  m_range[1] = 0;
   
  float tmp[1];
  std::ifstream inFile;
  vtkFloatArray *volumeField =vtkFloatArray::New();
  volumeField->SetNumberOfTuples(m_visOpt.nRows);
  
  int iVol;

  if(m_visOpt.dataRead && m_visOpt.goodAllocation)
  {
     std::map<std::string, int>::iterator p;
     for(p=m_visOpt.columns.begin();p!=m_visOpt.columns.end();p++)
	if(p->first==m_visOpt.isosurfaceField) iVol=p->second;

     for(int i=0; i<m_visOpt.nRows;i++)
       volumeField->  SetValue(i,m_visOpt.tableData[iVol][i]); //! this is the row that fill xAxis array
  } else
  {

  inFile.open(m_visOpt.path.c_str(), ios::binary);
 // std::clog<<m_visOpt.path.c_str()<<std::endl;
      
  if(!inFile)
    return -1;
     inFile.seekg((m_visOpt.nIsosurfaceField * m_visOpt.nRows)* sizeof(float));

  for (int i=0;i<m_visOpt.nRows;i++)
  {
     inFile.read((char *)(tmp ),  sizeof(float));
         
    if(m_visOpt.needSwap)
      tmp[0]=floatSwap((char *)(&tmp[0]));
          
    volumeField->  SetValue(i,tmp[0]);
  }
  inFile.close();
  }//else
  volumeField->SetName(m_visOpt.vRenderingField.c_str());
 
  m_imageData->SetScalarTypeToFloat();
  m_imageData->SetExtent(1,(int)m_visOpt.comp[0],1,(int)m_visOpt.comp[1],1,(int)m_visOpt.comp[2]);
  m_imageData->SetSpacing(m_visOpt.size[0],m_visOpt.size[1],m_visOpt.size[2]);
  m_imageData->GetPointData()->SetScalars(volumeField);
  m_imageData->GetScalarRange(m_range);
// vtkStructuredPointsWriter* spw = vtkStructuredPointsWriter::New();
//    spw->SetInput(m_imageData);
//    spw->SetFileName("/home/ube/bovolo.vtk");
//    spw->SetFileTypeToASCII();
//    spw->Write();

vtkImageGaussianSmooth *igs = vtkImageGaussianSmooth::New();
   igs->SetInput(m_imageData);
   igs->SetDimensionality(3);
   if(m_visOpt.isoSmooth=="medium")
   {
   igs->SetRadiusFactors(0.8 , 0.8 , 0.8); // 0 0 0 = escluso     0.8=medio    1.1=alto
   igs->SetStandardDeviations(1.6 , 1.6 , 1.6); //  0 0 0 = escluso  1.2=medio  2=alto
   }   
   else if(m_visOpt.isoSmooth=="high")
   {
   igs->SetRadiusFactors(1.1 , 1.1 , 1.1); // 0 0 0 = escluso     0.8=medio    1.1=alto
   igs->SetStandardDeviations(2 , 2 , 2); //  0 0 0 = escluso  1.2=medio  2=alto
   }
   else
   {
    igs->SetRadiusFactors(0 , 0 , 0); // 0 0 0 = escluso     0.8=medio    1.1=alto
    igs->SetStandardDeviations(0 , 0 , 0); //  0 0 0 = escluso  1.2=medio  2=alto
   }
// isosurface value  is normalized between 0 and 255 (it should be not normalized...only to avoid user mistakes)
if(m_visOpt.isosurfaceValue>255.0)
	std::cerr<<"Warning: Isosurface value out of range (0-255)"<<std::endl;
m_visOpt.isosurfaceValue=(m_visOpt.isosurfaceValue/255)*(m_range[1]-m_range[0]);
 
vtkContourFilter *cf = vtkContourFilter::New();
//   cf->SetInput(m_imageData );
   cf->SetInput((vtkDataSet *) igs->GetOutput() );
   cf->SetValue(0 , m_visOpt.isosurfaceValue);  // valore della prima (ed unica) isosuperficie
    cf->SetComputeNormals(1);
    cf->SetComputeScalars(1);

vtkPolyDataNormals *pdn = vtkPolyDataNormals::New();
   pdn->SetInput((vtkPolyData *) cf->GetOutput());
    pdn->SetFeatureAngle(90);
    pdn->SetComputeCellNormals(0);
    pdn->SetComputePointNormals(1);
    pdn->SetFlipNormals(0);
    pdn->SetSplitting(0);

vtkPolyDataMapper *pdm = vtkPolyDataMapper::New();
   pdm->SetInput((vtkPolyData *) pdn->GetOutput());
   pdm->SetNumberOfPieces(1);
   pdm->SetScalarRange(m_range[0] , m_range[1]);  
   pdm->SetColorMode(0);
   pdm->SetResolveCoincidentTopology(0);
   pdm->SetScalarMode(0);
   pdm->SetImmediateModeRendering(1);
   pdm->SetScalarVisibility(0);
   pdm->SetUseLookupTableScalarRange(0);

vtkProperty *P = vtkProperty::New();
     P->SetColor(1, 1 ,1);  //default is white
     if(m_visOpt.oneColor=="yellow") P->SetColor(1, 1 ,0);  
     if(m_visOpt.oneColor=="red") P->SetColor(1, 0 ,0);  
     if(m_visOpt.oneColor=="green") P->SetColor(0, 1 ,0);  
     if(m_visOpt.oneColor=="blue") P->SetColor(0, 0 ,1);  
     if(m_visOpt.oneColor=="cyan") P->SetColor(0, 1 ,1);  
     if(m_visOpt.oneColor=="violet") P->SetColor(1, 0 ,1);  
     if(m_visOpt.oneColor=="black") P->SetColor(0, 0 ,0);  
   //P->SetRepresentation(0);  // points
//   P->SetRepresentation(1);  // wireframe
     P->SetRepresentation(2);  // solid
     if(m_visOpt.wireframe) P->SetRepresentation(1); 

vtkActor *a = vtkActor::New();
   a->SetMapper(pdm);
   a->SetProperty(P);
   a->SetOrigin(0 , 0 , 0);
   a->SetPosition(0 , 0 , 0);
   a->SetScale(1 , 1 , 1);
/*   a->SetPickable(1);*/
   a->SetVisibility(1);

   P->Delete();
   pdm->Delete();
   pdn->Delete();
   cf->Delete();
   igs->Delete();

  
  

  setBoundingBox(m_imageData);
  m_pRenderer->AddActor(a);
  m_pRenderWindow->AddRenderer(m_pRenderer);

  m_pRenderer->SetBackground ( 0.0,0.0,0.0 );
     if(m_visOpt.backColor=="yellow") m_pRenderer->SetBackground (1, 1 ,0);  
     if(m_visOpt.backColor=="red")m_pRenderer->SetBackground (1, 0 ,0);  
     if(m_visOpt.backColor=="green") m_pRenderer->SetBackground (0, 1 ,0);  
     if(m_visOpt.backColor=="blue") m_pRenderer->SetBackground (0, 0 ,1);  
     if(m_visOpt.backColor=="cyan") m_pRenderer->SetBackground (0, 1 ,1);  
     if(m_visOpt.backColor=="violet") m_pRenderer->SetBackground (1, 0 ,1);  
     if(m_visOpt.backColor=="white") m_pRenderer->SetBackground (1, 1 ,1);  
  
  if(m_visOpt.imageSize=="small") 
	m_pRenderWindow->SetSize ( 512,365 );
  else if(m_visOpt.imageSize=="large") 
	m_pRenderWindow->SetSize ( 1024,731 );
  else 
	m_pRenderWindow->SetSize ( 792,566 );

  m_pRenderWindow->SetWindowName ("VisIVOServer View");

  if(m_visOpt.stereo)
  {  
      m_pRenderWindow->StereoRenderOn();
     if(m_visOpt.stereoMode=="RedBlue") 
	m_pRenderWindow->SetStereoTypeToRedBlue();

     else if(m_visOpt.stereoMode=="Anaglyph")
     {  
	m_pRenderWindow->SetStereoTypeToAnaglyph();
	m_pRenderWindow->SetAnaglyphColorSaturation(m_visOpt.anaglyphsat);
	m_visOpt.anaglyphmask=trim(m_visOpt.anaglyphmask);
	std::stringstream anatmp(m_visOpt.anaglyphmask);
	int rightColorMask=4, leftColorMask=3;
	anatmp>>rightColorMask;
	anatmp>>leftColorMask;
	m_pRenderWindow->SetAnaglyphColorMask(rightColorMask,leftColorMask);
     }	
     
     else
     {
       m_pRenderWindow->SetStereoTypeToCrystalEyes();

       if(m_visOpt.stereoImg==0)
        { 
	  m_pRenderWindow->SetStereoTypeToRight();
        }else 
        {     
	  m_pRenderWindow->SetStereoTypeToLeft();
        }
     }   
      m_pRenderWindow->StereoUpdate();
  }

       
  
//open view
//------------------------------------------------------------------
  vtkGenericRenderWindowInteractor* inter=vtkGenericRenderWindowInteractor::New();
  m_pRenderWindow->SetInteractor(inter);
  m_pRenderWindow->Render();
//--------------------------------------------------------------------

  setCamera ();
  double *bounds;
  bounds=new double[6];
  bounds[0]=0;
  bounds[1]=m_visOpt.comp[0];
  bounds[2]=0;
  bounds[3]=m_visOpt.comp[1];
  bounds[4]=0;
  bounds[5]=m_visOpt.comp[2];
  m_visOpt.xField="X";
  m_visOpt.yField="Y";
  m_visOpt.zField="Z";
  
  if(m_visOpt.showAxes)  setAxes(m_imageData,bounds);
  delete [] bounds;   
 //open view
 //-----------------------------------
  inter->Start();
  inter->ExitEvent();
 //-----------------------------------

  if(inter!=0)
    inter->Delete();
    volumeField->Delete();
    a->Delete();
}
