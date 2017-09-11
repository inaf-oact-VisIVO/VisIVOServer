#ifndef __VSRANDOMIZERTABLEOP_H
#define __VSRANDOMIZERTABLEOP_H

#include "vstableop.h"

#ifdef VSMPI
#include "mpi.h"
#endif

class VSRandomizerTableOp : public VSTableOp
{
 static const unsigned int MAX_NUMBER_TO_REDUCE_ROW;
 static const unsigned int MIN_NUMBER_OF_ROW;
 int m_lastRandomValue;
 float **m_fArray;
 float *m_tmpArray;
  bool m_firstRanSet;
  unsigned int m_nOfCol;
  unsigned long long int m_nOfRow;
  unsigned long long int  m_ranFromRow, m_ranToRow;
  bool allocatefArray();
  bool subset(unsigned long long int fromRow,unsigned long long int toRow);
 unsigned long long int  m_valueToExtract;
  unsigned long long int m_numberOfExtractedValue;//! number of already writed values
  unsigned long long int m_listElement;//! total number of write values


public:
    VSRandomizerTableOp();
    ~VSRandomizerTableOp();
    void printHelp();
    bool execute();
};

#endif

