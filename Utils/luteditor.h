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

#ifndef LUTEDITOR_H
#define LUTEDITOR_H
#include "optionssetter.h"

class  optionsSetter;
class vtkLookupTable;

void SelectLookTable(VisIVOServerOptions *opt, vtkLookupTable *lut=NULL);

void lutVolRenGreen(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutVolRenGlow(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutTenStep(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutTemperature(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutSar(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutPhysicsContour(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutGlow(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutEField(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutDefault(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutDefaultStep(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutGray(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutRun1(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutRun2(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutPureRed(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutPureGreen(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutPureBlue(VisIVOServerOptions *opt,vtkLookupTable *lut);

void lutAllYellow(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutAllCyane(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutAllViolet(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutAllWhite(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutAllBlack(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutAllRed(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutAllGreen(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutAllBlu(VisIVOServerOptions *opt,vtkLookupTable *lut);

void lutMinMax(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutVolRenRGB(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutVolRenTwoLev(VisIVOServerOptions *opt,vtkLookupTable *lut);
void lutFile(VisIVOServerOptions *opt,vtkLookupTable *lut);

#endif


