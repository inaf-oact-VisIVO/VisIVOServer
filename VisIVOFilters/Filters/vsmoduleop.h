#ifndef __VSMODULEOP_H
#define __VSMODULEOP_H

#include "vstableop.h"


class VSModuleOp : public VSTableOp
{
static const unsigned int MIN_NUMBER_OF_ROW;
static const unsigned int MAX_NUMBER_TO_REDUCE_ROW;
 

 float **m_fArray;
  bool allocateArray();	
  unsigned int m_nOfCol;
  unsigned long long int m_nOfRows;
  int m_nOfEle;

  public:
    VSModuleOp();
    ~VSModuleOp();
    void printHelp();
    bool execute();
};

#endif

