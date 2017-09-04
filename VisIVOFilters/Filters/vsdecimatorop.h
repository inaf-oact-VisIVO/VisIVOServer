#ifndef __VSDECIMATORTABLEOP_H
#define __VSDECIMATORTABLEOP_H

#include "vstableop.h"


class VSDecimatorTableOp : public VSTableOp
{
 static const unsigned int MAX_NUMBER_TO_REDUCE_ROW;
 static const unsigned int MIN_NUMBER_OF_ROW;
 int m_step; 

 float **m_fArray;
 float **m_f1Array;
  bool m_firstRanSet;
  unsigned int m_nOfCol;
  unsigned long long int m_nOfRow;
  unsigned long long int m_nOfDecimatedRow;
  unsigned long long int  m_ranFromRow, m_ranToRow;
  bool allocatefArray();
  bool subset(unsigned long long int fromRow,unsigned long long int toRow);
 void computeDecimatedPoints();
 unsigned long long int  m_valueToExtract;
  unsigned long long int m_numberOfExtractedValue;//! number of already writed values
  unsigned long long int m_listElements;//! total number of write values
  public:
    VSDecimatorTableOp();
    ~VSDecimatorTableOp();
    void printHelp();
    bool execute();
};

#endif

