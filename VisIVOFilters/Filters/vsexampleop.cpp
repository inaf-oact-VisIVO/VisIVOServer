//
// C++ Implementation: vsexampleop
//
// Description: 
//
//
// Author:  <>, (C) 2008
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
#include <fstream>
#include "vsexampleop.h"
#include "vstable.h"

//---------------------------------------------------------------------
bool VSExampleOp::execute()
//---------------------------------------------------------------------
{
//!This example method read data of the table and write on stdout and on file defined with --out in the command line the data of the header file, the first 20 values of the table and the following three values : the first, the second  and the last values of the of the first and the third columns in the table. The  example.bin table contains 10 columns and 1000 rows. The row 0 (first row) contain values from 0 to 999. the row 1 from 1000 to 1999, the row 2 from 2000 to 2999, and so on, the row 9 from 9000 to 9999

	/*** Output Filename  **/
	std::string fileNameOutput=getParameterAsString("out"); //! read the value as a string in the VSTableOP map that contain the key "out". The map is filled by the setParameters method called in the main.cpp program.

        std::ofstream fileOutput(fileNameOutput.c_str(), std::ios::out); //! open file for output

	unsigned long long int totRows=m_tables[0]->getNumberOfRows(); //! return the rows of the table

	unsigned int totColumns=m_tables[0]->getNumberOfColumns();

	std::cout <<"Rows in the table "<<totRows<<std::endl; 
	std::cout <<"Columns in the table "<<totColumns<<std::endl; 

	for(unsigned int i=0; i< totColumns;i++)
		std::cout<<m_tables[0]->getColName(i)<<std::endl; //! write columns names

//! Read and  Write the first 20  elements (if exist)  of all the columns in the table

	int numberOfElements=20; // twenty values to be read.
	int numberOfColumns=totColumns; // all the columns will be read.

	if(totRows<numberOfElements)numberOfElements=totRows; //! if the table has less than 20 values we will read all the values

	unsigned int *colList;
	colList=new unsigned int [numberOfColumns]; //! this array contain the numbers of the columns that I want to obtain from the table.
	
	for(int i=0;i<numberOfColumns;i++)
		colList[i]=i;

	float **fArray;  //! data read from file are saved here
	fArray=new  float*[numberOfColumns];
	for(int i=0;i<numberOfColumns;i++)
		fArray[i]= new float[numberOfElements];


	unsigned long long int fromRow=0; //! first row number to be read
	unsigned long long int toRow=numberOfElements-1; //! last row number to be read


	m_tables[0]->getColumn(colList, numberOfColumns,  fromRow, toRow, fArray); //! read data from table and put the value in fArray

//! Write data on file

	int width=15;
	int precision=8;

  	fileOutput.width(width);
 	fileOutput.precision(precision);
 	fileOutput.setf(std::ios::left);
	for(unsigned int i=0;i<numberOfColumns;i++)
	{ 
		fileOutput<<m_tables[0]->getColName(i).c_str(); //write on file column names
  		fileOutput.width(width);  
  		fileOutput.precision(precision);
	}
	fileOutput<<std::endl;


	for(unsigned int j=0;j<numberOfElements;j++)
	{
		for(unsigned int i=0;i<numberOfColumns;i++)
		{
				fileOutput<<fArray[i][j]; //write on file the table values
  				fileOutput.width(width); 
  				fileOutput.precision(precision);
		}
		fileOutput<<std::endl;
	}
	fileOutput<<std::endl<<std::endl<<std::endl;

	delete [] colList;
	for(int i=0; i<	numberOfColumns;i++)
		delete [] fArray[i];
	delete [] fArray;

//! Read and  Write 3  elements (if exist)  of the columns 2 and 3 in the table

	numberOfElements=3; // 3 values to be read.

	if(totRows<numberOfElements)numberOfElements=totRows; //! if the table has less than 3 values we will read all the values

	numberOfColumns=2; //! two columns to be read
	if(totColumns<2)numberOfColumns=totColumns;//! if the table has less than 2 columns we will read all the columns

	colList= new unsigned int [numberOfColumns]; //! this array contain the numbers of the columns I want to obtain from the table.
	
	colList[0]=0; //! first column to be read (numbered from 0)

	if(totColumns>=3)
		colList[1]=2;  //! third column to be read
	if(totColumns==2)
		colList[1]=1;  //! second column to be read


	fArray=new  float*[numberOfColumns];
	for(int i=0;i<numberOfColumns;i++)
		fArray[i]= new float[numberOfElements];

	unsigned long long int *list; //! contain the list of elements I want to read from table
	list= new unsigned long long int[numberOfElements];
	
	list[0]=0; //! first element to be read
	if(numberOfElements>=2)
		list[1]=1;  //! second element to be read
	if(numberOfElements>2)
		list[2]=totRows-1;  //! last element to be read
	

	m_tables[0]->getColumnList(colList,numberOfColumns, list, numberOfElements, fArray); //! read elements listed in list from columns listed in colList  and put the value in fArray

//! write data on file
	
	for(unsigned int i=0;i<numberOfColumns;i++)
	{ 
		fileOutput<<m_tables[0]->getColName(colList[i]).c_str(); //write on file column names
  		fileOutput.width(width);  
  		fileOutput.precision(precision);
	}
	fileOutput<<std::endl;

	for(unsigned int j=0;j<numberOfElements;j++)
	{
		for(unsigned int i=0;i<numberOfColumns;i++)
		{
				fileOutput<<fArray[i][j]; //write on file the table values
  				fileOutput.width(width); 
  				fileOutput.precision(precision);
		}
		fileOutput<<std::endl;
	}

	fileOutput.close();
   return true;
}
