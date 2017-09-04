//
// C++ Implementation: vschangecolnameop
//
// Description: 
//
//
// Author:  <>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <string>
#include <sstream>
#include <set>
#include "vstable.h"
#include "vschangecolnameop.h"

//---------------------------------------------------------------------
VSChangeColNameop::VSChangeColNameop()
{
}
//---------------------------------------------------------------------

//---------------------------------------------------------------------
VSChangeColNameop::~VSChangeColNameop()
{
}
//---------------------------------------------------------------------
//---------------------------------------------------------------------
void VSChangeColNameop::printHelp()
//---------------------------------------------------------------------
{
std::cout<<"Chage the column name in a VBT"<<std::endl<<std::endl;
std::cout<<"Usage: VisIVOFilters --op changecolname --field column_names --newnames new_names [--help] [--file] inputFile.bin"<<std::endl<<std::endl;

std::cout<<"Example: VisIVOFilters --op changecolname --field X_r Y_r --newnames X_n Y_n --file inputFile.bin"<<std::endl;

std::cout<<"Note: "<<std::endl;
std::cout<<"--field Valid columns names."<<std::endl;
std::cout<<"--newnames New column names."<<std::endl;
std::cout<<"--file  Input table filename."<<std::endl;
std::cout<<"--help produce this output "<<std::endl;

return;
 }
//---------------------------------------------------------------------

//---------------------------------------------------------------------
bool VSChangeColNameop::execute()
//---------------------------------------------------------------------
{
if(!isParameterPresent("field") || getParameterAsString("field").empty() || getParameterAsString("field")=="unknown" )
{
	std::cerr<<"Invalid parameter field. Operation aborted"<<std::endl;
	return false;
}
if(!isParameterPresent("newnames") || getParameterAsString("newnames").empty() || getParameterAsString("newnames")=="unknown" )
{
	std::cerr<<"Invalid parameter newnames. Operation aborted"<<std::endl;
	return false;
}
std::string colName;
std::vector<unsigned int> colNumberSet;
std::vector<std::string> newNameSet;
std::stringstream ssListparameters;
ssListparameters.str(getParameterAsString("field"));
while (!ssListparameters.eof())
{
	ssListparameters>>colName;
	if(m_tables[0] -> getColId(colName)>=0)
		colNumberSet.push_back(m_tables[0] -> getColId(colName));
}
std::stringstream ssNewName;
ssNewName.str(getParameterAsString("newnames"));
int totalCol=0;
while (!ssNewName.eof())
{
	std::string tmp;
	ssNewName>>tmp;
	newNameSet.push_back(tmp);
}
totalCol=colNumberSet.size();
if(totalCol>newNameSet.size()) 
	totalCol=newNameSet.size();
for(int i=0;i<totalCol;i++)
	m_tables[0]->setColName(colNumberSet[i],newNameSet[i]);
m_tables[0] -> writeHeader();
m_realOutFilename.push_back(m_tables[0]->getLocator());

return true;	
}
//---------------------------------------------------------------------

