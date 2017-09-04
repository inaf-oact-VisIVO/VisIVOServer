//
// C++ Interface: vsgrid2pointdistr
//
// Description: 
//
//
// Author:  <>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef VSGRID2POINTDISTR_H
#define VSGRID2POINTDISTR_H
#include "vstableop.h"

/**
	@author 
*/
class VSGrid2PointDistr: public VSTableOp
{
	float** m_fArray;
	float **m_grid;
	unsigned int m_nOfCol;
	unsigned long long int m_nOfRow;
	unsigned long long int m_nOfGridRow;

	
	int m_nOfEle;
 	unsigned long long int m_nOfGridEle; // numer of points grid input value and adjusted after allocation
	static const unsigned int MIN_NUMBER_OF_ROW;
	static const unsigned int MAX_NUMBER_TO_REDUCE_ROW;
 	bool allocateArray();
  bool m_tsc;
  bool m_cic;
  bool m_ngp;
   	float m_origin[3];
  	float m_spacing[3];
    bool setOrigin();
    bool setSpacing();
     const unsigned int * m_gridCellNumber;
	unsigned int m_sampleDimensions[3];
public:
    VSGrid2PointDistr();

    ~VSGrid2PointDistr();
    bool setOrigin(float x0, float y0,float z0);
    bool setSpacing(float xs, float ys,float zs);
    void printHelp();
    bool execute();

};

#endif
