/***************************************************************************
 *   Copyright (C) 2008 by Gabriella Caniglia, Roberto Munzone *
 *  gabriella.caniglia@oact.inaf.it *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifndef GADGETSOURCE_H
#define GADGETSOURCE_H

#include "abstractsource.h"

#include <vector>
#include <string>


struct headerType2
{
  char     first_fortran_legacy[4]; //!snapformat=2 : you have to "jump" 4+2*4 =12 bytes and you'll find the the Block size
  int      boh[1];
  int	   nBlock[1];	
  int      size[1];
  unsigned int      npart[6];
  double   mass[6];
  double   time[1];
  double   redshift[1];
  int      flag_sfr[1];
  int      flag_feedback[1];
  int      npartTotal[6];
  int      foling[1];
  int      num_files[1];
  double   BoxSize[1];
  double   Omega0[1];
  double   OmegaLambda[1];
  double   HubbleParam[1];
  int      FlagAge[1];
  int      FlagMetals[1];
  int      NallWH[6];
  int      flag_entr_ics[1]; 
  char     fill[256- 6*sizeof(int)- 6*sizeof(double)- 2*sizeof(double)- 2*sizeof(int)- 6*sizeof(int)- 2*sizeof(int)- 
      4*sizeof(double)- 9*sizeof(int)]; /*! fills to 256 Bytes */
  int      final_boh[1];
  int	   final_nBlock[1];
  
  char     tagFirstBlock[4];
  int      first_boh[1];
  int      first_nBlock[1];
  int      sizeFirstBlock[1]; 
   
};


struct headerType1
{
	//char     first_fortran_legacy[4]; //snapformat=2 : you have to "jump" 4+2*4 =12 bytes and you'll find the the Block size
	//int      boh[2];	
	int      size[1];
	unsigned int      npart[6];
	double   mass[6];
	double   time[1];
	double   redshift[1];
	int      flag_sfr[1];
	int      flag_feedback[1];
	int      npartTotal[6];
	int      foling[1];
	int      num_files[1];
	double   BoxSize[1];
	double   Omega0[1];
	double   OmegaLambda[1];
	double   HubbleParam[1];
	int      FlagAge[1];
	int      FlagMetals[1];
	int      NallWH[6];
	int      flag_entr_ics[1]; 
	char     fill[256- 6*sizeof(int)- 6*sizeof(double)- 2*sizeof(double)- 2*sizeof(int)- 6*sizeof(int)- 2*sizeof(int)- 
			4*sizeof(double)- 9*sizeof(int)]; /* fills to 256 Bytes */
	int      final_boh[1];
	int      sizeFirstBlock[1];
   
};

class GadgetSource : public AbstractSource
   
{
  public: //! Read the headerType2 file and set the basic table parameters
    int readHeader();
    int readData();
        
  private:
    std::vector <std::string> m_fieldsNames;
    int m_nRows; 
    char m_dataType, m_Endian;
    
     int m_snapformat;
     char tmpType[4]; int numBlock; int m_sizeBlock[1];
    
    std::vector<std::string> checkType;
    std::string tagType;
    	
    void swapHeaderType2();
    void swapHeaderType1();
    
    struct headerType2 m_pHeaderType2;
    struct headerType1 m_pHeaderType1;
  
};
  

#endif
