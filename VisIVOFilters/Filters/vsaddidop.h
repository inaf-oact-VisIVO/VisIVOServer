//
// C++ Interface: vsaddidop
//
// Description: 
//
//
// Author:  <>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef VSADDIDOP_H
#define VSADDIDOP_H
#include "vstableop.h"

/**
	@author 
*/
class VSAddIdOp: public VSTableOp
{
static const unsigned int MAX_NUMBER_TO_REDUCE_ROW;
static const unsigned int MIN_NUMBER_OF_ROW;

float **m_fArray;
  bool allocatefArray();
  unsigned  int m_nOfCol;
  unsigned  long long int m_nOfRow;
  int m_nOfEle;

public:
    VSAddIdOp();
    ~VSAddIdOp();
    void printHelp();
    bool execute();

};

#endif
