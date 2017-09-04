#include <algorithm>

#include "cxxsupport/bstream.h"
#include "writer/writer.h"

using namespace std;

namespace {

class pixel
  {
  public:
    uint8 r,g,b;

    pixel () : r(0),g(0),b(0) {}
    pixel (const COLOUR &c)
      : r(uint8(min(255,int(256*c.r)))),
        g(uint8(min(255,int(256*c.g)))),
        b(uint8(min(255,int(256*c.b)))) {}

    bool operator== (const pixel &other) const
      { return (r==other.r)&&(g==other.g)&&(b==other.b); }
    bool operator!= (const pixel &other) const
      { return (r!=other.r)||(g!=other.g)||(b!=other.b); }
  };

void write_equal_range (const arr<pixel> &px, tsize begin, tsize end,
  bofstream &file)
  {
  chunkMaker cm (end-begin,128);
  uint64 cbeg, csz;
  while (cm.getNext(cbeg,csz))
    {
    file << uint8(csz-1+128);
    file << px[begin].b << px[begin].g << px[begin].r;
    }
  }
void write_unequal_range (const arr<pixel> &px, tsize begin, tsize end,
  bofstream &file)
  {
  chunkMaker cm (end-begin,128);
  uint64 cbeg, csz;
  while (cm.getNext(cbeg,csz))
    {
    file << uint8(csz-1);
    for (tsize cnt=begin+cbeg; cnt< begin+cbeg+csz; ++cnt)
      file << px[cnt].b << px[cnt].g << px[cnt].r;
    }
  }

} // unnamed namespace

void write_tga(paramfile &params, const arr2<COLOUR> &pic,
  const string &frame_name)
  {
  cout << " writing tga file '" << frame_name << "'" << endl;
  tsize xres=pic.size1(), yres=pic.size2();
  const uint8 header[18] = { 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    xres%256, xres/256, yres%256, yres/256, 24, 32 };

  bofstream file(frame_name.c_str(),file_is_natural);

  file.put(&header[0],18);
  for (tsize y=0; y<yres; ++y)
    for (tsize x=0; x<xres; ++x)
      {
      uint8 pix[3];
      pix[0] = uint8(min(255,int(256*pic[x][y].b)));
      pix[1] = uint8(min(255,int(256*pic[x][y].g)));
      pix[2] = uint8(min(255,int(256*pic[x][y].r)));
      file.put(&pix[0],3);
      }
  }

void write_tga_rle(paramfile &params, const arr2<COLOUR> &pic,
  const string &frame_name)
  {
  cout << " writing tga file '" << frame_name << "'" << endl;
  tsize xres=pic.size1(), yres=pic.size2();
  const uint8 header[18] = { 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    xres%256, xres/256, yres%256, yres/256, 24, 32 };

  bofstream file(frame_name.c_str(),file_is_natural);

  file.put(&header[0],18);
  for (tsize y=0; y<yres; ++y)
    {
    arr<pixel> px(xres);
    for (tsize x=0; x<xres; ++x) px[x] = pic[x][y];
    tsize xstart=0;
    while (xstart<xres)
      {
      if (xstart==xres-1)
        {
        write_unequal_range (px,xstart,xstart+1,file);
        xstart=xres;
        }
      else
        {
        if (px[xstart+1]==px[xstart]) // range of equal pixels
          {
          tsize xend=xstart+2;
          while ((xend<xres) && (px[xend]==px[xstart])) ++xend;
          write_equal_range (px,xstart,xend,file);
          xstart=xend;
          }
        else
          {
          tsize xend=xstart+2;
          while ((xend<xres) && (px[xend]!=px[xend-1])) ++xend;
          write_unequal_range (px,xstart,xend,file);
          xstart=xend;
          }
        }
      }
    }
  }
