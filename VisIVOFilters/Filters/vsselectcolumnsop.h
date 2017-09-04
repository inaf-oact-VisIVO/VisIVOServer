#ifndef __VSSELECTCOLUMNS_H
#define __VSSELECTCOLUMNS_H

#include "vstableop.h"


class VSSelectColumnsOp : public VSTableOp
{
  
  public:
    VSSelectColumnsOp();
    ~VSSelectColumnsOp();
    void printHelp();
    bool execute();
};

#endif

