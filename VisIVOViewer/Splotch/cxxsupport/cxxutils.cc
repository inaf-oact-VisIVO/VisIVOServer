/*
 *  This file is part of libcxxsupport.
 *
 *  libcxxsupport is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxsupport is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxsupport; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  This file contains the implementation of various convenience functions
 *  used by the Planck LevelS package.
 *
 *  Copyright (C) 2002 - 2010 Max-Planck-Society
 *  Authors: Martin Reinecke, Reinhard Hell
 */

// if we are using g++, check for version 3.0 or higher
#ifdef __GNUC__
#if (__GNUC__<3)
#error your C++ compiler is too old. g++ version 3.0 or higher is required.
#endif
#endif

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include "cxxutils.h"
#include "datatypes.h"
#include "openmp_support.h"
#include "mpi_support.h"
#include "sse_utils.h"

using namespace std;

string SPLtrim (const string &orig)
  {
  string::size_type p1=orig.find_first_not_of(" \t");
  if (p1==string::npos) return "";
  string::size_type p2=orig.find_last_not_of(" \t");
  return orig.substr(p1,p2-p1+1);
  }

template<typename T> string dataToString (const T &x)
  {
  ostringstream strstrm;
  strstrm << x;
  return SPLtrim(strstrm.str());
  }

template<> string dataToString (const bool &x)
  { return x ? "T" : "F"; }
template<> string dataToString (const string &x)
  { return SPLtrim(x); }
template<> string dataToString (const float &x)
  {
  ostringstream strstrm;
  strstrm << setprecision(8) << x;
  return SPLtrim(strstrm.str());
  }
template<> string dataToString (const double &x)
  {
  ostringstream strstrm;
  strstrm << setprecision(16) << x;
  return SPLtrim(strstrm.str());
  }

template string dataToString (const signed char &x);
template string dataToString (const unsigned char &x);
template string dataToString (const short &x);
template string dataToString (const unsigned short &x);
template string dataToString (const int &x);
template string dataToString (const unsigned int &x);
template string dataToString (const long &x);
template string dataToString (const unsigned long &x);
template string dataToString (const long long &x);
template string dataToString (const unsigned long long &x);

string intToString(int64 x, tsize width)
  {
  ostringstream strstrm;
  strstrm << setw(width) << setfill('0') << x;
  return SPLtrim(strstrm.str());
  }

namespace {

void end_stringToData (const string &x, const char *tn, istringstream &strstrm)
  {
  string error = string("conversion error in stringToData<")+tn+">(\""+x+"\")";
  planck_assert (strstrm,error);
  string rest;
  strstrm >> rest;
//  rest=SPLtrim(rest);
  planck_assert (rest.length()==0,error);
  }

} // unnamed namespace

template<typename T> void stringToData (const string &x, T &value)
  {
  istringstream strstrm(x);
  strstrm >> value;
  end_stringToData (x,type2typename<T>(),strstrm);
  }

template<> void stringToData (const string &x, string &value)
  { value = SPLtrim(x); }

template<> void stringToData (const string &x, bool &value)
  {
  const char *x2 = x.c_str();
  const char *fval[] = {"F","f","n","N","false",".false.","FALSE",".FALSE." };
  const char *tval[] = {"T","t","y","Y","true",".true.","TRUE",".TRUE." };
  for (tsize i=0; i< sizeof(fval)/sizeof(fval[0]); ++i)
    if (strcmp(x2,fval[i])==0) { value=false; return; }
  for (tsize i=0; i< sizeof(tval)/sizeof(tval[0]); ++i)
    if (strcmp(x2,tval[i])==0) { value=true; return; }
  planck_fail("conversion error in stringToData<bool>(\""+x+"\")");
  }

template void stringToData (const string &x, signed char &value);
template void stringToData (const string &x, unsigned char &value);
template void stringToData (const string &x, short &value);
template void stringToData (const string &x, unsigned short &value);
template void stringToData (const string &x, int &value);
template void stringToData (const string &x, unsigned int &value);
template void stringToData (const string &x, long &value);
template void stringToData (const string &x, unsigned long &value);
template void stringToData (const string &x, long long &value);
template void stringToData (const string &x, unsigned long long &value);
template void stringToData (const string &x, float &value);
template void stringToData (const string &x, double &value);

bool equal_nocase (const string &a, const string &b)
  {
  if (a.size()!=b.size()) return false;
  for (tsize m=0; m<a.size(); ++m)
    if (tolower(a[m])!=tolower(b[m])) return false;
  return true;
  }

string tolower(const string &input)
  {
  string result=input;
  for (tsize m=0; m<result.size(); ++m)
    result[m]=char(tolower(result[m]));
  return result;
  }

namespace {

void openmp_status()
  {
  int threads = openmp_max_threads();
  if (threads>1)
    cout << "OpenMP active: max. " << threads << " threads. " << endl;
  }

void MPI_status()
  {
  int tasks = mpiMgr.num_ranks();
  if (tasks>1)
    cout << "MPI active with " << tasks << " tasks. " << endl;
  }

void SSE_status()
  {
#if(defined(PLANCK_HAVE_SSE)||defined(PLANCK_HAVE_SSE2))
  cout << "Processor features detected: ";
#if(defined(PLANCK_HAVE_SSE)&&defined(PLANCK_HAVE_SSE2))
  cout << "SSE, SSE2" << endl;
#elif(defined(PLANCK_HAVE_SSE))
  cout << "SSE" << endl;
#else
  cout << "SSE2" << endl;
#endif
#endif
  }

} //unnamed namespace

void announce (const string &name)
  {
  cout << endl << "+-";
  for (tsize m=0; m<name.length(); ++m) cout << "-";
  cout << "-+" << endl;
  cout << "| " << name << " |" << endl;
  cout << "+-";
  for (tsize m=0; m<name.length(); ++m) cout << "-";
  cout << "-+" << endl << endl;
  SSE_status();
  openmp_status();
  MPI_status();
  cout << endl;
  }

void module_startup (const string &name, int argc, const char **,
  int argc_expected, const string &argv_expected, bool verbose)
  {
  if (verbose) announce (name);
  if (argc==argc_expected) return;
  cerr << "Usage: " << name << " " << argv_expected << endl;
  planck_fail_quietly ("Incorrect usage");
  }

void parse_file (const string &filename, map<string,string> &dict)
  {
  int lineno=0;
  dict.clear();
  ifstream inp(filename.c_str());
  planck_assert (inp,"Could not open parameter file '"+filename+"'.");
  while (inp)
    {
    string line;
    getline(inp, line);
    ++lineno;
    // remove potential carriage returns at the end of the line
    line=line.substr(0,line.find("\r"));
    line=line.substr(0,line.find("#"));
    line=SPLtrim(line);
    if (line.size()>0)
      {
      string::size_type eqpos=line.find("=");
      if (eqpos!=string::npos)
        {
        string key=SPLtrim(line.substr(0,eqpos)),
               value=SPLtrim(line.substr(eqpos+1,string::npos));
        if (key=="")
          cerr << "Warning: empty key in '" << filename << "', line "
               << lineno << endl;
        else
          {
          if (dict.find(key)!=dict.end())
            cerr << "Warning: key '" << key << "' multiply defined in '"
                 << filename << "', line " << lineno << endl;
          dict[key]=value;
          }
        }
      else
        cerr << "Warning: unrecognized format in '" << filename << "', line "
             << lineno << ":\n" << line << endl;
      }
    }
  }

void calcShareGeneral (int64 glo, int64 ghi, int64 nshares, int64 myshare,
  int64 &lo, int64 &hi)
  {
  int64 nwork = ghi-glo;
  int64 nbase = nwork/nshares;
  int64 additional = nwork%nshares;
  lo = glo+myshare*nbase + ((myshare<additional) ? myshare : additional);
  hi = lo+nbase+(myshare<additional);
  }

template<typename T> void split (istream &stream, vector<T> &list)
  {
  while (stream)
    {
    string word;
    stream >> word;
    planck_assert (stream||stream.eof(),
      string("error while splitting stream into ") + type2typename<T>()
      + "components");
    if (stream) list.push_back(stringToData<T>(word));
    }
  }

template<typename T> void split (const string &inp, vector<T> &list)
  {
  istringstream stream(inp);
  split (stream,list);
  }

template void split (const string &inp, vector<string> &list);
template void split (const string &inp, vector<float> &list);
template void split (const string &inp, vector<double> &list);
template void split (const string &inp, vector<int> &list);
template void split (const string &inp, vector<long> &list);

void tokenize (const string &inp, char delim, vector<string> &list)
  {
  istringstream stream(inp);
  string token;
  list.clear();
  while (getline(stream,token,delim))
    list.push_back(token);
  }
