#ifndef __VSSIGMACONTOURS_H
#define __VSSIGMACONTOURS_H

#include "vstableop.h"

class VSSigmaContoursOp : public VSTableOp{

	float** m_fArray;
	unsigned int m_nOfCol;
	unsigned int m_nOfAllocatedCol;
	unsigned long long int* m_goodEventList;
	float m_nOfSigma;    
	
	unsigned long long int m_nOfRow;
	int m_nOfEle;
	static const unsigned int MIN_NUMBER_OF_ROW;
	static const unsigned int MAX_NUMBER_TO_REDUCE_ROW;
	bool m_exclude;  
 	bool allocateArray();
 	
public:

    VSSigmaContoursOp();
    ~VSSigmaContoursOp();
    bool execute();
    void printHelp();
    
protected:
    
    void SetNumberOfSigma(float value) {m_nOfSigma= value;};
 	  float GetNumberOfSigma() {return m_nOfSigma;};

    		
};

#endif
