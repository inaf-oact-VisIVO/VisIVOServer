//
// C++ Interface: vsextractlistrowsop
//
// Description: 
//
//
// Author:  <>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef VSEXTRACTLISTROWSOP_H
#define VSEXTRACTLISTROWSOP_H

#include "vstableop.h"
#include <fstream>
/**
	@author 
*/
class VSExtractListRowsOp: public VSTableOp
{
static const unsigned int MAX_NUMBER_TO_REDUCE_ROW;
static const unsigned int MIN_NUMBER_OF_ROW;
static const long int MAX_SEEK;
  float **m_fArray;
  float **m_wfArray;
  unsigned   int *m_list;
  unsigned  long long int *m_lllist;
  unsigned  int m_nOfCol;
  unsigned  long long int m_nOfRow;
  unsigned  long long int m_totalListEle;
  int m_nOfEle;
  int m_startPoint;
  int m_nOfEleList;
  int m_residualLoadList;
  bool allocatefArray();
  bool m_oneList;
  bool m_numberLists;
  bool m_listElemnts;
  int m_elementsInLists;
  bool m_listTotal;

bool readListsElements();
bool m_multlistFormatBinary;
int m_numberOfLists;
std::ifstream m_fileInput;
std::vector<unsigned long long int> m_numberOfElements;
public:
    VSExtractListRowsOp();
    ~VSExtractListRowsOp();
    void printHelp();
    bool execute();

};

#endif
