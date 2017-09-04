//
// C++ Interface: vsvbt2ahf
//
// Description: 
//
//
// Author:  <>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef VSVBT2AHF_H
#define VSVBT2AHF_H
#include "vstableop.h"

/**
	@author 
*/
class VSvbt2ahf: public VSTableOp
{
	static const unsigned int MAX_NUMBER_TO_REDUCE_ROW;
	static const unsigned int MIN_NUMBER_OF_ROW;
	static const unsigned int NUM_AMIGA_OUT;
 	bool allocateArray();
	float **m_fArray;
	float *m_Amiga;
	int m_nOfCol;
	unsigned long long  int m_nOfRow;
	int m_nOfEle;
	float m_vMean;
	float m_velFactor;
public:
    VSvbt2ahf();
    ~VSvbt2ahf();
    void printHelp();
    bool execute();

};

#endif
