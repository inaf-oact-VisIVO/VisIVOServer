#ifndef __VSPOLAROP_H 
#define __VSPOLAROP_H


#include "vstableop.h"


class VSPolarOp : public VSTableOp
{
 static const unsigned int MIN_NUMBER_OF_ROW;
 static const unsigned int MAX_NUMBER_TO_REDUCE_ROW; 

 float **m_fArray; 
 unsigned int m_nOfCol; 
 unsigned int m_nOfColIn; 
 unsigned long long int m_nOfRow; 
 int m_nOfElements;

 bool allocateArray(); 
 
 public:
 VSPolarOp(); 
 ~VSPolarOp(); 
 bool execute(); 
 void printHelp(); 
};

#endif
