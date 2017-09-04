//
// C++ Interface: vsahfhalolistop
//
// Description: 
//
//
// Author:  <>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef VSAHFHALOLISTOP_H
#define VSAHFHALOLISTOP_H
#include "vstableop.h"

/**
	@author 
*/
class VSAhfHaloListOp: public VSTableOp
{
static const unsigned int MAX_NUMBER_TO_REDUCE_ROW;
static const unsigned int MIN_NUMBER_OF_ROW;

float **m_fArray;
  bool allocatefArray();
  bool createNewList();
  unsigned  int m_nOfCol;
  unsigned  long long int m_nOfRealList;
  unsigned  long long int m_nOfRow;
  int m_nOfEle;
std::ifstream m_fileInput;

std::ofstream m_fileOutput;

public:
    VSAhfHaloListOp();

    ~VSAhfHaloListOp();
    void printHelp();
    bool execute();

};

#endif
