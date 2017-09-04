//
// C++ Interface: vschangecolnameop
//
// Description: 
//
//
// Author:  <>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef VSCHANGECOLNAMEOP_H
#define VSCHANGECOLNAMEOP_H
#include "vstableop.h"

/**
	@author 
*/
class VSChangeColNameop: public VSTableOp
{
public:
    VSChangeColNameop();
    ~VSChangeColNameop();
    void printHelp();
    bool execute();

};

#endif
