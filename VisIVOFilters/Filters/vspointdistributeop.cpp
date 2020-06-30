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
#include <fstream>
#include <cmath>
#ifdef WIN32
#include <time.h>
#endif
#include "vspointdistributeop.h"
#include "vstable.h"
#include "VisIVOFiltersConfigure.h"

const unsigned int VSPointDistributeOp::MAX_NUMBER_TO_REDUCE_ROW = 100000;
const unsigned int VSPointDistributeOp::MIN_NUMBER_OF_ROW = 100;

//---------------------------------------------------------------------
VSPointDistributeOp::VSPointDistributeOp()
//---------------------------------------------------------------------
{
    m_modelBounds[0] = 0.0;
    m_modelBounds[1] = 0.0;
    m_modelBounds[2] = 0.0;
    m_modelBounds[3] = 0.0;
    m_modelBounds[4] = 0.0;
    m_modelBounds[5] = 0.0;
    
    m_sampleDimensions[0] = 50;
    m_sampleDimensions[1] = 50;
    m_sampleDimensions[2] = 50;
    
    m_nullValue = 0.0;
    
    m_splattedScalar = -1;
    
    m_useConstant = false;
    m_constValue=1.0;
    m_fArray= NULL;
    m_grid=NULL;
    m_executeDone= false;
    m_tsc=false;
    m_cic=true;
    m_ngp=false;
    m_OriginSet=false;
    m_SpacingSet=false;
    m_gridSpacing=false;
    m_avg=false;
}

//---------------------------------------------------------------------
VSPointDistributeOp::~VSPointDistributeOp()
//---------------------------------------------------------------------
{
	if(m_fArray !=NULL)
	{
        for(int i=0;i<3;i++)
            delete [] m_fArray[i];
        delete [] m_fArray;
	}
	if(m_grid != NULL)
		delete [] m_grid;
}
//---------------------------------------------------------------------
void VSPointDistributeOp::printHelp()
//---------------------------------------------------------------------
{
    std::cout<<"It produces a table which represent a volume from selected fields of the input table that are distributed using NGP , CIC (default) or TSC algorithm"<<std::endl<<std::endl;
    std::cout<<"Usage: VisIVOFilters --op pointdistribute  --resolution x_res y_res z_res --points x_col y_col z_col [--field column_names] [--nodensity] [--avg] [--out filename_out.bin] [--tsc] [--ngp] [--gridOrigin xg0 xg1 xg2] [--gridSpacing sg0 sg1 sg2]  [--box length] [--periodic] [--history] [--historyfile filename.xml] [--help] [--file] inputFile.bin"<<std::endl<<std::endl;
    
    std::cout<<"Example: VisIVOFilters --op pointdistribute --resolution 16 16 16 --points X Y Z --field Mass Temperature   --out filename_out.bin --file inputFile.bin"<<std::endl;
    
    std::cout<<"Note:  "<<std::endl;
    std::cout<<"--resolution  3D mesh size."<<std::endl;
    std::cout<<"--points Columns to be assumed for points coordinates."<<std::endl;
    std::cout<<"--field Valid columns name list to be distributed in the grid."<<std::endl;
    std::cout<<"--constant Assign a constant to all points to be distributed in the grid Ignored if field option is given. Default value is a 1.0 for all points."<<std::endl;
    std::cout<<"--nodensity Overrides the default behavior. The field distribution is not divided for the cell volume."<<std::endl;
    std::cout<<"--avg Distributes the first field on the volume grid and compute the aritmethic average value on each cell, of the first field. The output volume table will have three field. For each cell the number of total elements in the cell NumberOfElements, the sum of total field value fieldSum, and the aritmetic average value fieldAvg. Only the ngp algorithm will be applied."<<std::endl;
    std::cout<<"--out Name of the new table. Default name is given."<<std::endl;
    std::cout<<"--tsc. The TSC algorithm is adopted."<<std::endl;
    std::cout<<"--ngp. The NGP algorithm is adopted."<<std::endl;
    std::cout<<"--gridOrigin. It specifies the coordinates of the lower left corner of the grid. Default values are assumed from the box of inputFile.bin"<<std::endl;
    std::cout<<"--gridSpacing. It specifies the length of each cell dimension in arbitray unit. This parameter is ignored if the box option is given. Default vaules are assumed from the box of inputFile.bin"<<std::endl;
    std::cout<<"--box. It specifies the length of a box. Default value is assumed from the box of inputFile.bin if the gridSpacing option is not given"<<std::endl;
    std::cout<<"--periodic. It specifies the box is periodic. Particles outside the box limits are considered inside on the other side."<<std::endl;
    std::cout<<"--history (optional) create an XML file which contains the history of operations performed (default create hist.xml file)"<<std::endl;
    std::cout<<"--historyfile [filename.xml]   (optional) Change default history file name  and or directory "<<std::endl;

    std::cout<<"--file Input table filename."<<std::endl;
    
    std::cout<<"--help produce this output."<<std::endl;
    
    return;
    
}
//---------------------------------------------------------------------
bool VSPointDistributeOp::getOrigin(float *origin)
//---------------------------------------------------------------------
{
	if(!m_executeDone)
		return false;
	origin[0]=m_origin[0];
	origin[1]=m_origin[1];
	origin[2]=m_origin[2];
	return true;
}
//---------------------------------------------------------------------
bool VSPointDistributeOp::getSpacing(float *spacing)
//---------------------------------------------------------------------
{
	if(!m_executeDone)
		return false;
	spacing[0]=m_spacing[0];
	spacing[1]=m_spacing[1];
	spacing[2]=m_spacing[2];
	return true;
    
}
//---------------------------------------------------------------------
bool VSPointDistributeOp::allocateArray(int nField)
//---------------------------------------------------------------------
{
    unsigned long long int tempLL=getMaxNumberInt();
    if(((unsigned long long int)m_nOfRow*m_nOfCol)>tempLL)
        m_nOfRow=(int)tempLL/m_nOfCol;
    
    try
    {
        m_fArray=new  float*[m_nOfCol];
    }
    catch(std::bad_alloc &e)
    {
        m_fArray=NULL;
    }
    
	if(m_fArray==NULL)
		return false;
    if(m_avg) nField=3;
    try
    {
        m_grid=new  float*[nField];
    }
    catch(std::bad_alloc &e)
    {
        m_grid=NULL;
    }
    
	if(m_grid==NULL)
		return false;
    
	bool goodAllocation=false;
	while(!goodAllocation)
	{
		goodAllocation=true;
		for(unsigned int i=0;i<m_nOfCol;i++)
		{
            try
            {
                m_fArray[i]=new  float[m_nOfRow];
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
				if(m_nOfRow==MIN_NUMBER_OF_ROW)
				{
					delete [] m_fArray;
					m_fArray=NULL;
					return false;
				}
				m_nOfRow=m_nOfRow-MAX_NUMBER_TO_REDUCE_ROW;
				if(m_nOfRow<=MAX_NUMBER_TO_REDUCE_ROW) m_nOfRow=MIN_NUMBER_OF_ROW;
				break;
			}
            //		std::clog<<i<<" " <<m_nOfRow<<std::endl;
		}
		if(!goodAllocation)
			continue;
		for(unsigned int i=0;i<nField;i++)
		{
            try
            {
                m_grid[i]=new  float[m_numNewPts];
            }
            catch(std::bad_alloc &e)
            {
                m_grid[i]=NULL;
            }
            
			if(m_grid[i]==NULL)
			{
				goodAllocation=false;
				for(unsigned int j=0;j<i;j++)
					delete [] m_grid[j];
				for(unsigned int j=0;j<m_nOfCol;j++)
					delete [] m_fArray[j];
				if(m_numNewPts==MIN_NUMBER_OF_ROW)
				{
					delete [] m_fArray;
					delete [] m_grid;
					m_fArray=NULL;
					m_grid=NULL;
					return false;
				}
				m_nOfRow=m_nOfRow-MAX_NUMBER_TO_REDUCE_ROW;
				if(m_nOfRow<=MAX_NUMBER_TO_REDUCE_ROW)
					m_nOfRow=MIN_NUMBER_OF_ROW;
				m_numNewPts=m_numNewPts-MAX_NUMBER_TO_REDUCE_ROW;
				if(m_numNewPts<=MAX_NUMBER_TO_REDUCE_ROW)
					m_numNewPts=MIN_NUMBER_OF_ROW;
				break;
			}
            
		}
        
	}
	return true;
}

//---------------------------------------------------------------------
bool VSPointDistributeOp::computeModelBounds()
//---------------------------------------------------------------------
// Compute ModelBounds from input geometry.
// bounds are the lower and upper coordinates of the particles
// in code units
// serach for min and max coordinates in the table
{
    
    unsigned int counterCols=3;
    unsigned long long int totRows=m_tables[0]->getNumberOfRows();
    
    unsigned long long int totEle=totRows;
    unsigned long long int fromRow, toRow, startCounter=0;
    unsigned int nOfValidElement=0;
    float maxValue[3],minValue[3];
    while(totEle!=0)
    {
        fromRow=startCounter;
        toRow=fromRow+m_nOfRow-1;
        if(toRow>totRows-1)
            toRow=totRows-1;
        m_tables[0]->getColumn(m_colList, counterCols, fromRow, toRow, m_fArray);
        if(startCounter==0)
        {
            for(int k=0;k<3;k++)
            {
                maxValue[k]=m_fArray[k][0];
                minValue[k]=m_fArray[k][0];
            }
        }
        for(unsigned int j=0;j<(toRow-fromRow+1);j++)
        {
            for(int k=0;k<3;k++)
            {
                if(maxValue[k]<m_fArray[k][j]) maxValue[k]=m_fArray[k][j];
                if(minValue[k]>m_fArray[k][j]) minValue[k]=m_fArray[k][j];
            }
        }
		startCounter=toRow+1;
		totEle=totEle-(toRow-fromRow+1);
		if(totEle<0) totEle=0;
    }
    //  for(int i=0; i<3; i++) std::clog<<"i="<<i<<" min="<<minValue[i]<<std::endl;//AA
    //  for(int i=0; i<3; i++) std::clog<<"i="<<i<<" max="<<maxValue[i]<<std::endl;//AA
    
    //// END search
    
    for(int i=0; i<3; i++) m_modelBounds[2*i] = minValue[i];
    for(int i=0; i<3; i++) m_modelBounds[2*i+1] = maxValue[i];
    
    // Set volume origin and data spacing
    
    for (int i=0; i<3; i++)
    {
        if(!m_OriginSet)
            m_origin[i] = m_modelBounds[2*i];
        if(!m_SpacingSet)
            m_spacing[i] = (m_modelBounds[2*i+1] - m_origin[i]) / m_sampleDimensions[i];
    }
    return true;
}

//---------------------------------------------------------------------
bool VSPointDistributeOp::setOrigin()
{
    
    if(getParameterAsString("gridOrigin").empty()||getParameterAsString("gridOrigin")=="unknown")
    {
        std::cerr<<"No valid input gridOrigin parameter."<<std::endl;
        return false;
    }
    std::stringstream ssOrigin;
    ssOrigin.str(getParameterAsString("gridOrigin"));
    int count=0;
    while (!ssOrigin.eof())
    {
        ssOrigin>>m_origin[count];
        if(count==2)
            break;
        count++;
    }
    
    return true;
}
//---------------------------------------------------------------------
bool VSPointDistributeOp::setSpacing()
{
    if(isParameterPresent("box"))
    {
        float box=getParameterAsFloat("box");
        if(box<=0.0)
        {
            std::cerr<<"No valid input box parameter."<<std::endl;
            return false;
        }
        for(int i=0;i<3;i++)
        {
            if(m_sampleDimensions[i]<=0)
            {
                std::cerr<<"No valid grid cell number."<<std::endl;
                return false;
            }
            m_spacing[i]=box/m_sampleDimensions[i];
        }
        return true;
    }
    if(isParameterPresent("gridSpacing"))
    {
        if(getParameterAsString("gridSpacing").empty()||getParameterAsString("gridSpacing")=="unknown")
        {
            std::cerr<<"No valid input gridSpacing parameter. Operation aborted"<<std::endl;
            return false;
        }
        std::stringstream ssSpacing;
        ssSpacing.str(getParameterAsString("gridSpacing"));
        int count=0;
        while (!ssSpacing.eof())
        {
            ssSpacing>>m_spacing[count];
            if(m_spacing[count] <=0.)
            {
                std::cerr<<"No valid input gridSpacing parameter."<<std::endl;
                return false;
            }
            if(count==2)
                break;
            count++;
        }
    }
    return true;
}
//---------------------------------------------------------------------

//---------------------------------------------------------------------
bool VSPointDistributeOp::execute()
//---------------------------------------------------------------------
{
    if(isParameterPresent("avg"))
        m_avg=true;
    bool periodic=false;
    if(isParameterPresent("periodic"))
        periodic=true;
    VSTable tableGrid;
    std::vector<int> fieldList;
    // check for points  coordinate columns
    int counterCols=0;
    std::stringstream ssListparameters;
    ssListparameters.str(getParameterAsString("points"));
    while (!ssListparameters.eof())
    {
        std::string paramField;
        ssListparameters>>paramField;
        if(m_tables[0] -> getColId(paramField)>=0)
        {
            m_colList[counterCols]=m_tables[0] -> getColId(paramField);
            counterCols++;
            if(counterCols==3)
                break;
        }
    }
    if(counterCols !=3)
    {
        std::cerr<<"VSPointDistributeOp: Invalid columns in --points argument is given"<<std::endl;
        return false;
    }
    // check for TSC
    if(isParameterPresent("tsc"))
    {
        m_tsc=true;
        m_ngp=false;
        m_cic=false;
    }
    if(isParameterPresent("ngp"))
    {
        if(m_tsc)
            std::cerr<<"Ignored --tsc parameter"<<std::endl;
        m_tsc=false;
        m_cic=false;
        m_ngp=true;
    }
    if(m_avg)	// force ngp
    {
        m_tsc=false;
        m_cic=false;
        m_ngp=true;
    }
    
    //check grid resolution
    
    counterCols =0;
    ssListparameters.clear();
    ssListparameters.str(getParameterAsString("resolution"));
    while (!ssListparameters.eof())
    {
        std::string paramField;
        ssListparameters>>paramField;
        m_sampleDimensions[counterCols]=atoi(paramField.c_str()); //set resolution
        if(m_sampleDimensions[counterCols]<=0)
        {
            std::cerr<<"VSPointDistributeOp: Invalid resolution is given"<<std::endl;
            return false;
        }
        counterCols++;
        if(counterCols==3)
            break;
    }
    
    if(counterCols<3)
    {
        std::cerr<<"VSPointDistributeOp: Invalid resolution is given"<<std::endl;
        return false;
    }
    std::stringstream fieldNameSStream;
    fieldNameSStream<<getParameterAsString("field");
    if(fieldNameSStream.str()==""||fieldNameSStream.str()=="unknown")
    {
        m_useConstant=true;
        fieldList.push_back(-1);
        if(isParameterPresent("constant"))
            m_constValue=getParameterAsFloat("constant");
    }else
    {
        while (!fieldNameSStream.eof())
        {
            std::string paramField;
            fieldNameSStream>>paramField;
            if(m_tables[0] -> getColId(paramField)>=0)
                fieldList.push_back(m_tables[0] -> getColId(paramField));
            if(m_avg) break;
        }
    }
    if(fieldList.size()<=0)
    {
        std::cerr<<"Pointdistribute. Invalid field is given"<<std::endl;
        return false;
    }
    //prepare colList: list of columns to be read
    unsigned int *colList;
    unsigned int nOfField= fieldList.size();
    m_nOfCol=3+nOfField;
    if(m_useConstant) m_nOfCol=3;
    try
    {
        colList= new unsigned int[m_nOfCol];
    }
    catch(std::bad_alloc &e)
    {
        colList=NULL;
    }
    
    if(colList==NULL)
    {
        std::cerr<<"Failed Array allocation. vspointdistribute Operation terminated"<<std::endl;
        return false;
    }
    
    // allocate m_arrays
    
    unsigned long long int totRows=m_tables[0]->getNumberOfRows();
    int maxInt=getMaxNumberInt();
    
    if(totRows>maxInt)
        m_nOfRow=maxInt;
    else
        m_nOfRow=totRows;
    
    unsigned long long int gridPts = m_sampleDimensions[0] * m_sampleDimensions[1] * m_sampleDimensions[2];
    if(gridPts>maxInt)
        m_numNewPts=maxInt;
    else
        m_numNewPts=gridPts;
    
    bool allocationArray=allocateArray((int) fieldList.size());
    
    
    if(m_fArray==NULL || m_grid==NULL || !allocationArray )
    {
        std::cerr<<"Failed Array allocation. vspointdistribute Operation terminated"<<std::endl;
        delete [] colList;
        return false;
    }
    
    
    
    if(isParameterPresent("gridOrigin"))
        m_OriginSet=setOrigin();
    if(isParameterPresent("box"))
        m_SpacingSet=setSpacing();
    if(isParameterPresent("gridSpacing") && !m_SpacingSet)
    {
        m_gridSpacing=true;
        m_SpacingSet=true;
        std::stringstream ssgridSpacing;
        ssgridSpacing.str(getParameterAsString("gridSpacing"));
        counterCols=0;
        while (!ssgridSpacing.eof())
        {
            ssgridSpacing>>m_spacing[counterCols];
            if(m_spacing[counterCols]<=0.)
            {
                std::cerr<<"Invalid gridSpacing values .Pointdistribute Operation terminated"<<std::endl;
                delete [] colList;
                return false;
            }
            counterCols++;
            if(counterCols==3)
                break;
        }
    }
    if(!m_SpacingSet || !m_OriginSet)
        if(!computeModelBounds())
        {
            delete [] colList;
            return false;
        }
    // initialize the grid
    //open file output
    std::stringstream fileNameOutputSStream;
    fileNameOutputSStream<<getParameterAsString("out");
    std::string fileNameOutput;
    if(fileNameOutputSStream.str()==""||fileNameOutputSStream.str()=="unknown")
    {
        fileNameOutputSStream.str().erase(); //QUI verificare
        std::string filenameInputTable=m_tables[0]->getLocator();
        int len=filenameInputTable.length();
        time_t rawtime;
        struct tm * timeinfo;
        char buffer [80];
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        strftime (buffer,80,"%Y%m%d%H%M",timeinfo);
        fileNameOutputSStream<<filenameInputTable.substr(0, len-4)<<"_pointdistribute_"<<buffer<<".bin";  //QUI verificare
    }
    
    fileNameOutput=fileNameOutputSStream.str();
    if(fileNameOutput.find(".bin") == std::string::npos)
   		fileNameOutput.append(".bin");
    m_realOutFilename.push_back(fileNameOutput);
    
    //Clean existing tab
    remove(fileNameOutput.c_str());
    tableGrid.setLocator(fileNameOutput);
#ifdef VSBIGENDIAN
    std::string endianism="big";
#else
    std::string endianism="little";
#endif
    
    tableGrid.setEndiannes(endianism);
    tableGrid.setType("float");
    if(m_avg)
    {
        std::stringstream fileNameColSStream;
        if(m_useConstant)
            fileNameColSStream<<"Constant";
        else
            fileNameColSStream<<m_tables[0] -> getColName(fieldList[0]);
        tableGrid.addCol("NumberOfElements");
        std::string avgColName=fileNameColSStream.str()+"Avg";
        std::string sumColName=fileNameColSStream.str()+"Sum";
        tableGrid.addCol(sumColName);
        tableGrid.addCol(avgColName);
    }
    else if(m_useConstant)
        tableGrid.addCol("Constant");
    else
        for(int i=0;i<fieldList.size();i++)
            tableGrid.addCol(m_tables[0] -> getColName(fieldList[i]));
    
    tableGrid.setNumberOfRows(gridPts);
    tableGrid.setIsVolume(true);
    tableGrid.setCellNumber(m_sampleDimensions[0],m_sampleDimensions[1],m_sampleDimensions[2]);
    float spacing[3];
    if(m_gridSpacing)
    {
        spacing[0]=m_spacing[0];
        spacing[1]=m_spacing[1];
        spacing[2]=m_spacing[2];
        
    } else{
        spacing[0]=1.0;
        spacing[1]=m_spacing[1]/m_spacing[0];
        spacing[2]=m_spacing[2]/m_spacing[0];
    }
    tableGrid.setCellSize(spacing[0],spacing[1],spacing[2]); //fixed cellSize QUI could be modified
    tableGrid.writeHeader();
    
    //  Create Empty Binary File
    // use only the first col of m_grid for comodity
    int tmp=(int) fieldList.size();
    unsigned int *gridList ;
    try {
        gridList = new unsigned int[tmp];
    }
    catch(std::bad_alloc &e)
    {
        gridList=NULL;
    }
    if (gridList==NULL)
    {
        delete [] colList;
        return false;
    }
    int numOfField=fieldList.size();
    if(m_avg) numOfField=3;
	
    for(int k=0;k<numOfField;k++)
    {
        gridList[k]=k;
        for(unsigned int i=0;i<(unsigned int) m_numNewPts;i++)
            m_grid[k][i]=0.0;
    }
    // fill the table with zero
    
    unsigned long long int startCounter=0;
    unsigned long long int fromRow,toRow;
    unsigned long long int totEle= gridPts;
    while(totEle!=0)
    {
        fromRow=startCounter;
        toRow=fromRow+m_numNewPts-1;
        if(toRow>gridPts-1)toRow=gridPts-1;
        tableGrid.putColumn(gridList, nOfField, fromRow,toRow,m_grid);
        totEle=totEle-(toRow-fromRow+1);
        startCounter=toRow+1;
        if(totEle<0) totEle=0;
    }
    
    //////////////////////
    /////
    // Start Point Distribution
    ////
    
    // set cashed grid on memory
    unsigned long long int gridIndex[2];  //limits of cashed grid
    gridIndex[0]=0;
    gridIndex[1]=m_numNewPts-1;
    totEle=totRows;
    startCounter=0;
    colList[0]=m_colList[0];
    colList[1]=m_colList[1];
    colList[2]=m_colList[2];
    if(!m_useConstant)
        for(int i=0;i<nOfField;i++)
            colList[3+i]=fieldList[i];
    
    int jkFactor = m_sampleDimensions[0]*m_sampleDimensions[1];
    int jFactor = m_sampleDimensions[0];
    float norm;
    
    
    float cellVolume=m_spacing[0]*m_spacing[1]*m_spacing[2];
    if(isParameterPresent("nodensity") || m_avg)
        cellVolume=1.0;
    
    while(totEle!=0)
    {
        // Table douwnLoad
        fromRow=startCounter;
        toRow=fromRow+m_nOfRow-1;
        if(toRow>totRows-1)
            toRow=totRows-1;
        m_tables[0]->getColumn(colList,m_nOfCol, fromRow, toRow, m_fArray);
        
        if(m_ngp)
        {
            int ind;
            int nCell=m_sampleDimensions[0]*m_sampleDimensions[1]*m_sampleDimensions[2];
            // NGP on points
            for (int ptId=0; ptId < toRow-fromRow+1; ptId++)
            {
                float px[3];
                px[0]=m_fArray[0][ptId];
                px[1]=m_fArray[1][ptId];
                px[2]=m_fArray[2][ptId];
                
                float pos1 = (float) (px[0] - m_origin[0]) / m_spacing[0];
                float pos2 = (float) (px[1] - m_origin[1]) / m_spacing[1];
                float pos3 = (float) (px[2] - m_origin[2]) / m_spacing[2];
                
                // find nearest grid points
                int i1 = floor(pos1);
                if(fabs((pos1-i1))>0.5)
                    i1++;
                int i2 = floor(pos2);
                if(fabs((pos2-i2))>0.5)
                    i2++;
                int i3 = floor(pos3);
                if(fabs((pos3-i3))>0.5)
                    i3++;
                if(periodic)
                {
                    
                    if(i1<0)
                        i1=m_sampleDimensions[0]+i1;
                    if(i2<0)
                        i2=m_sampleDimensions[1]+i2;
                    if(i3<0)
                        i3=m_sampleDimensions[2]+i3;
                    if(i1>m_sampleDimensions[0]-1)
                        i1=i1-m_sampleDimensions[0];
                    if(i2>m_sampleDimensions[1]-1)
                        i2=i2-m_sampleDimensions[1];
                    if(i3>m_sampleDimensions[2]-1)
                        i3=i3-m_sampleDimensions[2];
                }
                if (i1 <m_sampleDimensions[0] && i1 >= 0 && i2 < m_sampleDimensions[1] && i2 >= 0 && i3 < m_sampleDimensions[2] && i3 >= 0)
                {
               		ind = i3*jkFactor + i2*jFactor + i1;
                    if(ind>nCell-1) continue;
                    
                    // calculate density
                    
                    if(ind<gridIndex[0] || ind>gridIndex[1]) //ind1 NOT in cache!
                    {
                        tableGrid.putColumn(gridList, nOfField, gridIndex[0],gridIndex[1],m_grid); //QUI controlla:estremi INCLUSI
                        gridIndex[0]=ind;
                        gridIndex[1]=gridIndex[0]+m_numNewPts-1;
                        if(gridIndex[1]>=gridPts)gridIndex[1]=gridPts-1;
                        tableGrid.getColumn(gridList, tableGrid.getNumberOfColumns(),gridIndex[0],gridIndex[1],m_grid);
                    }
                    for(int j=0;j<fieldList.size();j++)//Note: fieldList.size() MUST be equal 1 if m_avg
                    {
                        if(m_useConstant)
                            norm=m_constValue;
                        else
                            norm=m_fArray[3+j][ptId];
                        if(m_avg)
                        {
                            m_grid[0][ind-gridIndex[0]]+=1.0;
                            m_grid[1][ind-gridIndex[0]]+=norm;
                            m_grid[2][ind-gridIndex[0]]=m_grid[1][ind-gridIndex[0]]/m_grid[0][ind-gridIndex[0]];
                        }
                        else
                            m_grid[j][ind-gridIndex[0]]+=norm/cellVolume;
                    }
                } //if(pos...)
                
                // end main loop
                
            } //for(... ptId..)
        } // close if ngp
        
        
        if(m_cic)
        {
            float wc=0.;
            int nCell=m_sampleDimensions[0]*m_sampleDimensions[1]*m_sampleDimensions[2];
            // CIC on points
            for (int ptId=0; ptId < toRow-fromRow+1; ptId++)
            {
                //		std::clog<<ptId<<std::endl;
                wc=0.;
                float px[3];
                px[0]=m_fArray[0][ptId];
                px[1]=m_fArray[1][ptId];
                px[2]=m_fArray[2][ptId];
                
                float pos1 = (float) (px[0] - m_origin[0]) / m_spacing[0];
                float pos2 = (float) (px[1] - m_origin[1]) / m_spacing[1];
                float pos3 = (float) (px[2] - m_origin[2]) / m_spacing[2];
                
                int i1 = floor(pos1);
                /*       		if(fabs((pos1-i1))>0.5)
                 i1++;*/
                int i2 = floor(pos2);
                /*       		if(fabs((pos2-i2))>0.5)
                 i2++;*/
                int i3 = floor(pos3);
                /*       		if(fabs((pos3-i3))>0.5)
                 i3++;*/
                int i11=i1+1;
                int i21=i2+1;
                int i31=i3+1;
               	float dd1=pos1-(float)i1;
               	float dd2=pos2-(float)i2;
               	float dd3=pos3-(float)i3;
               	float de1=1.0-dd1;
               	float de2=1.0-dd2;
               	float de3=1.0-dd3;
                if(periodic)
                {
                    if(i1<0)
                        i1=m_sampleDimensions[0]+i1;
                    if(i2<0)
                        i2=m_sampleDimensions[1]+i2;
                    if(i3<0)
                        i3=m_sampleDimensions[2]+i3;
                    if(i1>m_sampleDimensions[0]-1)
                        i1=i1-m_sampleDimensions[0];
                    if(i2>m_sampleDimensions[1]-1)
                        i2=i2-m_sampleDimensions[1];
                    if(i3>m_sampleDimensions[2]-1)
                        i3=i3-m_sampleDimensions[2];
                    if(i11<0)
                        i11=m_sampleDimensions[0]+i11;
                    if(i21<0)
                        i21=m_sampleDimensions[1]+i21;
                    if(i31<0)
                        i31=m_sampleDimensions[2]+i31;
                    if(i11>m_sampleDimensions[0]-1)
                        i11=i11-m_sampleDimensions[0];
                    if(i21>m_sampleDimensions[1]-1)
                        i21=i21-m_sampleDimensions[1];
                    if(i31>m_sampleDimensions[2]-1)
                        i31=i31-m_sampleDimensions[2];
                }
                
                // calculate weights
                float d[8];
                d[0]= de1*de2*de3;
                d[1]= dd1*de2*de3;
                d[2]= de1*dd2*de3;
                d[3]= dd1*dd2*de3;
                d[4]= de1*de2*dd3;
                d[5]= dd1*de2*dd3;
                d[6]= de1*dd2*dd3;
                d[7]= dd1*dd2*dd3;
                
                // linearize coordinates
                unsigned long long int ind[8];
                ind[0] = i3*jkFactor + i2*jFactor + i1;
                ind[1] = i3*jkFactor + i2*jFactor + i11;
                ind[2] = i3*jkFactor + i21*jFactor + i1;
                ind[3] = i3*jkFactor + i21*jFactor + i11;
                ind[4] = i31*jkFactor + i2*jFactor + i1;
                ind[5] = i31*jkFactor + i2*jFactor + i11;
                ind[6] = i31*jkFactor + i21*jFactor + i1;
                ind[7] = i31*jkFactor + i21*jFactor + i11;
                
                
                // calculate density
                
                for(int n=0;n<8;n++)
                {
                    if(ind[n]<0  || ind[n]>=nCell)
                        continue;
                    if(ind[n]<gridIndex[0] || ind[n]>gridIndex[1]) //ind1 NOT in cache!
                    {
                        tableGrid.putColumn(gridList, nOfField, gridIndex[0],gridIndex[1],m_grid); //QUI controlla:estremi INCLUSI
                        gridIndex[0]=ind[n];
                        gridIndex[1]=gridIndex[0]+m_numNewPts-1;
                        if(gridIndex[1]>=gridPts)gridIndex[1]=gridPts-1;
                        tableGrid.getColumn(gridList, tableGrid.getNumberOfColumns(),gridIndex[0],gridIndex[1],m_grid);
                    }
                    for(int j=0;j<fieldList.size();j++)
                    {
                        if(m_useConstant)
                            norm=m_constValue;
                        else
                            norm=m_fArray[3+j][ptId];
                        m_grid[j][ind[n]-gridIndex[0]]+=d[n]*norm/cellVolume;
                        wc+=d[n]*norm;
                        //		outpippo<<"ptId="<<ptId<<" GRID j="<<j<<" i="<<ind[n]-gridIndex[0] <<" curr val="<<d[n]*norm<<" acc="<<m_grid[j][ind[n]-gridIndex[0]]<<std::endl; //AA
                    }
                }
                
                // end main loop
                if(wc>1.1*norm*fieldList.size())
                {
                    std::cerr<<"Error 2 on cic schema. Operation Aborted"<<std::endl;
                    if(colList!=NULL) delete [] colList;
                    if(gridList!=NULL)delete [] gridList;
                    return false;
                }
                
            } //for(... ptId..)
        } // close if cic
        if(m_tsc) //TSC
        {
            float wc=0;
            int nCell=m_sampleDimensions[0]*m_sampleDimensions[1]*m_sampleDimensions[2];
            for (int ptId=0; ptId < toRow-fromRow+1; ptId++)
            {
      	        float px[3];
                int ind_x,ind_y,ind_z;
                int xList[4],yList[4],zList[4];
                /*		m_spacing[0]=70./128;
                 m_spacing[1]=70./128;
                 m_spacing[2]=70./128;*/
                double node_x, node_y, node_z, dist_x, dist_y, dist_z, w_x, w_y, w_z, w;
                
                px[0]=m_fArray[0][ptId]-m_origin[0];
                if(periodic && px[0]<0.) px[0]+=m_sampleDimensions[0]*m_spacing[0];
                if(periodic && px[0]>=(m_sampleDimensions[0]*m_spacing[0])) px[0]-=m_sampleDimensions[0]*m_spacing[0];
                
                px[1]=m_fArray[1][ptId]-m_origin[1];
                if(periodic && px[1]<0.) px[1]+=m_sampleDimensions[1]*m_spacing[1];
                if(periodic && px[1]>=(m_sampleDimensions[1]*m_spacing[1])) px[1]-=m_sampleDimensions[1]*m_spacing[1];
                
                px[2]=m_fArray[2][ptId]-m_origin[2];
                if(periodic && px[2]<0.) px[2]+=m_sampleDimensions[2]*m_spacing[2];
                if(periodic && px[2]>=(m_sampleDimensions[2]*m_spacing[2])) px[2]-=m_sampleDimensions[2]*m_spacing[2];
                ; 		
                ind_x = (int)floor(px[0]/m_spacing[0]);
                ind_y = (int)floor(px[1]/m_spacing[1]);
                ind_z = (int)floor(px[2]/m_spacing[2]);
                wc=0.;
                for(int j=0;j<4;j++)
                {
                    xList[j]=ind_x+j-1;
                    yList[j]=ind_y+j-1;    
                    zList[j]=ind_z+j-1;
                }
                
                for(int ix = 0; ix <4;ix++)
                {
                    //		    std::clog<<"ix= "<<ix<<" wc= "<<wc<<std::endl;
                    for(int iy = 0; iy <4;iy++)
                    {
                        //			std::clog<<"iy= "<<iy<<" wc= "<<wc<<std::endl;
                        for(int iz = 0; iz <4;iz++)
                        {
                            dist_x = px[0]/m_spacing[0] - xList[ix];
                            dist_y = px[1]/m_spacing[1] - yList[iy];
                            dist_z = px[2]/m_spacing[2] - zList[iz];
                            if(dist_x<0) dist_x=-1*dist_x;
                            if(dist_y<0) dist_y=-1*dist_y;
                            if(dist_z<0) dist_z=-1*dist_z;
                            
                            if(dist_x <= .5)   w_x = .75 - dist_x*dist_x;
                            else if (dist_x <= 1.5) w_x = .5*(1.5-dist_x)*(1.5-dist_x);
                            else w_x = 0;
                            
                            if(dist_y <= .5)   w_y = .75 - dist_y*dist_y;
                            else if (dist_y <= 1.5) w_y = .5*(1.5-dist_y)*(1.5-dist_y);
                            else w_y = 0;
                            
                            if(dist_z <= .5)   w_z = .75 - dist_z*dist_z;
                            else if (dist_z <= 1.5) w_z = .5*(1.5-dist_z)*(1.5-dist_z);
                            else w_z = 0;                
                            
                            if( (w_x < 0. || w_y  < 0.) ||  w_z < 0.)
                                std::cerr<<"WARNING: particle "<<ptId<<" xyzList "<<xList[ix]<<" "<< yList[iy]<<" "<<zList[iz]<<" AND w_x,y,z "<<w_x<<" "<<w_y<<" "<<w_z<<std::endl;
                            
                            w = w_x * w_y * w_z;
                            //				std::clog<<"iz= "<<iz<<" wc= "<<wc<<" w="<<w<<std::endl;
                            if (w != 0)
                            {
                                int facX= xList[ix];
                                int facY= yList[iy];
                                int facZ= zList[iz];
                                if(w>1.01) 
                                {
                                    std::cerr<<"Error 1 on tsc schema. operation Aborted."<<std::endl;
                                    if(colList!=NULL) delete [] colList;
                                    if(gridList!=NULL)delete [] gridList;
                                    return false;
                                }
                                if(periodic && facX<0)
                                    facX+=m_sampleDimensions[0];
                                if(periodic && facX>=m_sampleDimensions[0]) facX-=m_sampleDimensions[0];
                                if(periodic && facY<0) facY+=m_sampleDimensions[1];
                                if(periodic && facY>=m_sampleDimensions[1]) facY-=m_sampleDimensions[1];
                                if(periodic && facZ<0) facZ+=m_sampleDimensions[2];
                                if(periodic && facZ>=m_sampleDimensions[2]) facZ-=m_sampleDimensions[2];
                                
                                if(facX<0  || facX>= m_sampleDimensions[0])
                                    continue;
                                if(facY<0  || facY>= m_sampleDimensions[1])
                                    continue;
                                if(facZ<0  || facZ>= m_sampleDimensions[2])
                                    continue;
                                
                                int ind_test = facX + m_sampleDimensions[0]*facY + m_sampleDimensions[0]*m_sampleDimensions[1]*facZ;
                                if(ind_test<0 || ind_test>(nCell-1))
                                {
                                    std::cerr<<"Error on particle "<<ptId<<" xyzList "<<xList[ix]<<" "<< yList[iy]<<" "<<zList[iz]<<" AND w w_x,y,z "<<w<<" "<<w_x<<" "<<w_y<<" "<<w_z<<std::endl;
                                }
                                if(ind_test<gridIndex[0] || ind_test>gridIndex[1]) //ind1 NOT in cache!
                                {
                                    std::cerr<<"Particle "<<ptId<<" xyzList "<<xList[ix]<<" "<< yList[iy]<<" "<<zList[iz]<<" ind_test= "<<ind_test<<" "<<gridIndex[0]<<" "<<gridIndex[1]<<std::endl;
                                    tableGrid.putColumn(gridList, nOfField, gridIndex[0],gridIndex[1],m_grid); //QUI controlla:estremi INCLUSI
                                    gridIndex[0]=ind_test;
                                    gridIndex[1]=gridIndex[0]+m_numNewPts-1;
                                    if(gridIndex[1]>=gridPts)gridIndex[1]=gridPts-1;
                                    tableGrid.getColumn(gridList,nOfField,gridIndex[0],gridIndex[1],m_grid); 	
                                }
                                
                                
                                for(int j=0;j<fieldList.size();j++)
                                {
                                    if(m_useConstant) 
                                        norm=m_constValue;
                                    else
                                        norm=m_fArray[3+j][ptId];
                                    m_grid[j][ind_test-gridIndex[0]]+=w*norm/cellVolume;
                                    wc+=w*norm;
                                    
                                }
                            }
                            
                        }// closes for iz
                    }// closes for iy
                }// closes for ix
                if(wc>1.1*norm)
                {						
                    std::cerr<<"Error 2 on tsc schema. Operation Aborted"<<std::endl;
                    if(colList!=NULL) delete [] colList;
                    if(gridList!=NULL)delete [] gridList;
                    return false;
                }
                
            }// closes for ptId
        } // close if TSC
        //
        startCounter=toRow+1;
        totEle=totEle-(toRow-fromRow+1);
        if(totEle<0) totEle=0;
    }  // while(totEle!=0)
    // final write of tableGrid
    if(m_avg) nOfField=3;
    tableGrid.putColumn(gridList,nOfField,gridIndex[0],gridIndex[1],m_grid);  
    
    m_executeDone=true;
    /*      	std::ofstream outtFile("/home/ube/test/1.txt",std::ios::out);
     outtFile<<"index_x  index_y  index_z  density"<<std::endl;
     int nCell=m_sampleDimensions[0]*m_sampleDimensions[1]*m_sampleDimensions[2];
     
     for(int ind_x = 0; ind_x <m_numNewPts; ++ind_x)
     {
     double den_w = m_grid[0][ind_x];
     outtFile<<ind_x<<" "<<m_grid[0][ind_x]<<std::endl;
     }// closes x cycle 
     outtFile.close();*/
    if(colList!=NULL) delete [] colList;
    if(gridList!=NULL)delete [] gridList;
    return true;
}


