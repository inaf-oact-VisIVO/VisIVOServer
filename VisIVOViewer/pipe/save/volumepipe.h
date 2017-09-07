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

#ifndef VolumePipe_H
#define VolumePipe_H

#include "pipe.h"
#include "optionssetter.h"
#include "visivoutils.h"

   class vtkImageData;
   class vtkVolumeRayCastMapper;
   class vtkVolume;
   class vtkColorTransferFunction;
   class vtkImageMathematics;
   class vtkPiecewiseFunction;
   class vtkImageCast;
   class vtkVolumeProperty;
   class vtkVolumeRayCastCompositeFunction;
        
        
   class VolumePipe: public Pipe
{

  public:
    VolumePipe( VisIVOServerOptions options);
    ~VolumePipe();
   
  protected:
   
    int createPipe();
    bool setLookupTable();
    double m_range[2]; 
   void destroyAll();
    
    vtkColorTransferFunction *m_colorTransferFunction;
    vtkImageData *m_imageData;
    vtkImageMathematics *m_math;
    vtkImageMathematics *m_math2;
    vtkImageCast *m_charData;
    vtkPiecewiseFunction *m_opacityTransferFunction;
    vtkVolumeProperty *m_volumeProperty;
    vtkVolumeRayCastCompositeFunction * m_rayCastCompositFunction ;
    vtkVolumeRayCastMapper *m_rayCastMapper ;
    vtkVolume *m_volume ;
};
#endif
