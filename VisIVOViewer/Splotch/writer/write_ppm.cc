#include <algorithm>
#include <sstream>

#include "cxxsupport/bstream.h"
#include "writer/writer.h"

using namespace std;

void write_ppm_ascii(paramfile &params, const arr2<COLOUR> &pic,
  const string &frame_name)
  {
  cout << " writing ASCII PPM file '" << frame_name << "'" << endl;

  ofstream file(frame_name.c_str());

  file << "P3" << endl << pic.size1() << endl << pic.size2() << endl << 255 << endl; 
  for (tsize y=0; y<pic.size2(); ++y)
    {
    for (tsize x=0; x<pic.size1(); ++x)
      file << min(255,int(256*pic[x][y].r)) << " "
           << min(255,int(256*pic[x][y].g)) << " "
           << min(255,int(256*pic[x][y].b)) << "   ";
    file << endl;
    }
  }

void write_ppm_bin(paramfile &params, const arr2<COLOUR> &pic,
  const string &frame_name)
  {
  cout << " writing binary PPM file '" << frame_name << "'" << endl;

  bofstream file(frame_name.c_str(),file_is_natural);

  ostringstream header;
  header << "P6" << endl << pic.size1() << endl << pic.size2() << endl << 255
         << endl;
  string hdrdata = header.str();
  file.put(hdrdata.c_str(),hdrdata.size());

  for (tsize y=0; y<pic.size2(); ++y)
    for (tsize x=0; x<pic.size1(); ++x)
      {
      uint8 pix[3];
      pix[0] = uint8(min(255,int(256*pic[x][y].r)));
      pix[1] = uint8(min(255,int(256*pic[x][y].g)));
      pix[2] = uint8(min(255,int(256*pic[x][y].b)));
      file.put(&pix[0],3);
      }
  }
