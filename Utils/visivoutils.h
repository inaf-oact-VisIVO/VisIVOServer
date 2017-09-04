/***************************************************************************
 *   Copyright (C) 2008 by Gabriella Caniglia, Marco Comparato *
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

#ifndef VISIVOUTILS_H
#define VISIVOUTILS_H



#include <string>
#include <vector>


double masToRad(double mas);
void masToRad(double *mas, int n);
void masToRad(float *mas, int n);
double degToRad(double deg);
void degToRad(double *deg, int n);
void degToRad(float *deg, int n);
double hmsToRad(char *hms);
void hmsToRad(char **hmsIn, double *radOut, int n);
void hmsToRad(char **hmsIn, float *radOut, int n);
double dmsToRad(char *dms);
void dmsToRad(char **dmsIn, double *radOut, int n);
void dmsToRad(char **dmsIn, float *radOut, int n);

int iCompare(std::string str1, std::string str2);

std::string getTempFileName(std::string suffix);
std::string getDir(std::string path);
std::string getExt(std::string path);
std::string getName(std::string path);

std::string trimRight(const std::string & source, const std::string & t = " ");
std::string trimLeft(const std::string & source, const std::string & t = " ");
std::string trim(const std::string & source, const std::string & t = " ");

void findAndReplace(std::string &str, char find, char replace);

 void makeHeader( unsigned long long int rows,std::string path,const std::vector<std::string>,double size[],double comp[],std::string file);
  void makeHeader( int rows,std::string path,const std::vector<std::string>,double size[],double comp[],std::string file);
 double tryToSetDimension(int nRows);
 
//int VisIVORemoteLoadToFile(std::string path, std::string tempFileName);
 double doubleSwap(char *value);
 long double longdoubleSwap(char *value);
 float floatSwap(char *value);
 int intSwap(char *value);
 long int longintSwap(char *value);
 long long int longlongintSwap(char *value);
 bool remoteDownloadFiles(std::string remoteFile,std::string login,std::string downloadedFile);
void HSVtoRGB( double *r, double *g, double *b, double h, double s, double v );
  void sortarray(int *vect, int elements);
 bool remoteLfn(std::string local,std::string se, std::string vo,std::string lfnFile,bool vbt=false);

#endif

