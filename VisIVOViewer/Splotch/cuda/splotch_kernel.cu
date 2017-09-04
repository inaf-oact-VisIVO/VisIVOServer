/*
Try accelerating splotch with CUDA. July 2009.
Copyright things go here.
*/

//#include "splotch_kernel.h"
#include "splotch_cuda.h"

//MACROs
#define Pi 3.14159265358979323846264338327950288
#define get_xy_from_sn(sn, xmin, ymin, ymax, x, y)\
        {int x1 =sn/(ymax-ymin); int y1 =sn-x1*(ymax-ymin);\
         x  =x1 +xmin; y  =y1 +ymin;}
#define get_sn_from_xy(x,y,maxy,miny, sn)\
    {sn =x*(maxy-miny) +y;}

#define get_minmax(minv, maxv, val) \
         minv=min(minv,val); \
         maxv=max(maxv,val);
#define MAXSIZE 1000

/////////constant memory declaration /////////////////////
__constant__ cu_color_map_entry dmap[MAXSIZE];
__constant__ int ptype_points[10];
__constant__ cu_param dparams;

/////////help functions///////////////////////////////////
__device__ float    my_asinh(float val)
  {
  return log(val+sqrt(1.+val*val));
  }

__device__ void my_normalize(float minv, float maxv, float &val)
  {
  if (minv!=maxv) val =  (val-minv)/(maxv-minv);
  }

__device__ void clamp (float minv, float maxv, float &val)
  {
  val = min(maxv, max(minv, val));
  }

//fetch a color from color table on device
__device__ cu_color get_color(int ptype, float val, int mapSize, int ptypes)
  {
  __shared__ int map_size;
  __shared__ int map_ptypes;

  map_size = mapSize;
  map_ptypes = ptypes;
  //first find the right entry for this ptype
  int     start, end;
  start =ptype_points[ptype];
  if ( ptype == map_ptypes-1)//the last type
    end =map_size-1;
  else
    end =ptype_points[ptype+1]-1;

  //search the section of this type to find the val
  int i=start;
  while ((val>dmap[i+1].val) && (i<end)) ++i;

  const float fract = (val-dmap[i].val)/(dmap[i+1].val-dmap[i].val);
  cu_color clr1=dmap[i].color, clr2=dmap[i+1].color;
  cu_color        clr;
  clr.r =clr1.r + fract*(clr2.r-clr1.r);
  clr.g =clr1.g + fract*(clr2.g-clr1.g);
  clr.b =clr1.b + fract*(clr2.b-clr1.b);

  return clr;
  }

__global__ void k_post_process(cu_color *pic, int n)
  {
  //first get the index m of this thread
  int m=blockIdx.x *blockDim.x + threadIdx.x;
  if (m >=n)
    m =n;

  //each pic[m] should do the same calc, so sequence does not matter!
  pic[m].r =1.0 - exp( pic[m].r);
  pic[m].g =1.0 - exp( pic[m].g);
  pic[m].b =1.0 - exp( pic[m].b);
  }

__global__ void k_combine
  (int minx, int miny, int maxx, int maxy, int xres, int yres,
  cu_particle_splotch *p, int pStart, int pEnd, cu_fragment_AeqE *fbuf, cu_color *pic)
  {
  int m =blockIdx.x *blockDim.x + threadIdx.x;
  int n =(maxx-minx)*(maxy-miny);
  if (m >=n)
    m =n;

  //get global coordinate point(x,y) of this thread
  int point_x, point_y;
  get_xy_from_sn(m, minx, miny, maxy, point_x, point_y);

  //go through all particles, for each particle p if point(x,y) is in its region
  //p(minx,miny, maxx,maxy) do the following.
  //find the sequencial number sn1 in p(minx,miny, maxx,maxy), the fragment we are looking
  //for in fragment buffer is fragBuf[ sn1+p.posInFBuf ]
  //grab the fragment f(deltaR,deltaG,deltaB)
  //find the sequencial number sn2 of point(x,y) in the output pic.
  //pic[sn2] += f
  int sn1, sn2, local_x, local_y, fpos;
  for (int i=pStart; i<=pEnd; i++)
    {
    if ( point_x >=p[i].minx && point_x<p[i].maxx &&
         point_y >=p[i].miny && point_y<p[i].maxy)
      {
      local_x =point_x -p[i].minx;
      local_y =point_y -p[i].miny;
      get_sn_from_xy(local_x, local_y, p[i].maxy, p[i].miny,sn1);
      fpos =sn1 +p[i].posInFragBuf;

      get_sn_from_xy(point_x, point_y, yres,0, sn2);
      pic[sn2].r +=fbuf[fpos].aR;
      pic[sn2].g +=fbuf[fpos].aG;
      pic[sn2].b +=fbuf[fpos].aB;
      }
    }
  }

//device render function k_render1
__global__ void k_render1
  (cu_particle_splotch *p, int nP,
  void *buf, bool a_eq_e, float grayabsorb, int mapSize, int types)
  {
  //first get the index m of this thread
  int m;

  m =blockIdx.x *blockDim.x + threadIdx.x;
  if (m >=nP)//m goes from 0 to nP-1
    return;

  // coloring
  int ptype = p[m].type;
  float col1=p[m].e.r,col2=p[m].e.g,col3=p[m].e.b;
  clamp (0.0000001,0.9999999,col1);
  if (dparams.col_vector[ptype])
    {
    clamp (0.0000001,0.9999999,col2);
    clamp (0.0000001,0.9999999,col3);
    }
  float intensity=p[m].I;
  clamp (0.0000001,0.9999999,intensity);
  intensity *= dparams.brightness[ptype];

  cu_color e;
  if (dparams.col_vector[ptype])   // color from file
    {
    e.r=col1*intensity;
    e.g=col2*intensity;
    e.b=col3*intensity;
    }
  else   // get color, associated from physical quantity contained in e.r, from lookup table
    {
  //first find the right entry for this ptype
      if (ptype<types)
      {
        e = get_color(ptype, col1, mapSize, types);
        e.r *= intensity;
        e.g *= intensity;
        e.b *= intensity;
      }
      else
      { e.r =e.g =e.b =0.0; }
    }

  //make fbuf the right type
  cu_fragment_AeqE        *fbuf;
  cu_fragment_AneqE       *fbuf1;
  if (a_eq_e)
    fbuf =(cu_fragment_AeqE*) buf;
  else
    fbuf1 =(cu_fragment_AneqE*)buf;

  //now do the rendering
  const float powtmp = pow(Pi,1./3.);
  const float sigma0 = powtmp/sqrt(2*Pi);

  const float r = p[m].r;
  const float radsq = 2.25*r*r;
  const float stp = -0.5/(r*r*sigma0*sigma0);

  cu_color q; //e=p[m].e;
  if (!a_eq_e)
   {
     q.r = e.r/(e.r+grayabsorb);
     q.g = e.g/(e.g+grayabsorb);
     q.b = e.b/(e.b+grayabsorb);
   }
  const float intens = -0.5/(2*sqrt(Pi)*powtmp);
  e.r*=intens; e.g*=intens; e.b*=intens;

  const float posx=p[m].x, posy=p[m].y;
  unsigned int fpos =p[m].posInFragBuf;

  if (a_eq_e)
  {
    for (int x=p[m].minx; x<p[m].maxx; ++x)
    {
     float dxsq=(x-posx)*(x-posx);
     for (int y=p[m].miny; y<p[m].maxy; ++y)
      {
        float dsq = (y-posy)*(y-posy) + dxsq;
        if (dsq<radsq)
        {
          float att = __expf(stp*dsq);
          fbuf[fpos].aR = att*e.r;
          fbuf[fpos].aG = att*e.g;
          fbuf[fpos].aB = att*e.b;
        }
        else
        {
          fbuf[fpos].aR =0.0;
          fbuf[fpos].aG =0.0;
          fbuf[fpos].aB =0.0;
        }
      //for each (x,y)
      fpos++;
      }//y
    }//x
  }
  else
  {
    for (int x=p[m].minx; x<p[m].maxx; ++x)
    {
     float dxsq=(x-posx)*(x-posx);
     for (int y=p[m].miny; y<p[m].maxy; ++y)
      {
        float dsq = (y-posy)*(y-posy) + dxsq;
        if (dsq<radsq)
        {
          float att = __expf(stp*dsq);
          float   expm1;
          expm1 =__expf(att*e.r)-1.0;
          fbuf1[fpos].aR = expm1;
          fbuf1[fpos].qR = q.r;
          expm1 =__expf(att*e.g)-1.0;
          fbuf1[fpos].aG = expm1;
          fbuf1[fpos].qG = q.g;
          expm1 =__expf(att*e.b)-1.0;
          fbuf1[fpos].aB = expm1;
          fbuf1[fpos].qB = q.b;
        }
        else
        {
          fbuf1[fpos].aR =0.0;
          fbuf1[fpos].aG =0.0;
          fbuf1[fpos].aB =0.0;
          fbuf1[fpos].qR =1.0;
          fbuf1[fpos].qG =1.0;
          fbuf1[fpos].qB =1.0;
        }
      //for each (x,y)
      fpos++;
      }//y
    }//x
  }
 }


//Transform by kernel
__global__ void k_transform
  (cu_particle_sim *p, cu_particle_splotch *p2, int n)
  {

  //first get the index m of this thread
  int m=blockIdx.x *blockDim.x + threadIdx.x;
  if (m >n) m =n;

  //now do x,y,z
  float x,y,z;
  x =p[m].x*dparams.p[0] + p[m].y*dparams.p[1] + p[m].z*dparams.p[2] + dparams.p[3];
  y =p[m].x*dparams.p[4] + p[m].y*dparams.p[5] + p[m].z*dparams.p[6] + dparams.p[7];
  z =p[m].x*dparams.p[8] + p[m].y*dparams.p[9] + p[m].z*dparams.p[10]+ dparams.p[11];

  //do r
  float xfac = dparams.xfac;
  const float   res2 = 0.5*dparams.xres;
  const float   ycorr = .5f*(dparams.yres-dparams.xres);
  if (!dparams.projection)
    {
    x = res2 * (x+dparams.fovfct*dparams.dist)*xfac;
    y = res2 * (y+dparams.fovfct*dparams.dist)*xfac + ycorr;
    }
  else
    {
    xfac=1./(dparams.fovfct*z);
    x = res2 * (x+dparams.fovfct*z)*xfac;
    y = res2 * (y+dparams.fovfct*z)*xfac + ycorr;
    }

  float r = p[m].r;
  p[m].I /= r;
  r *= res2*xfac;

  const float rfac= sqrt(r*r + 0.25*dparams.minrad_pix*dparams.minrad_pix)/r;
  r *= rfac;
  p2[m].I = p[m].I/rfac;

  p2[m].isValid = false;

  // compute region occupied by the partile
  const float rfacr=dparams.rfac*r;
  int minx=int(x-rfacr+1);
  if (minx>=dparams.xres) return;
  minx=max(minx,0);

  int maxx=int(x+rfacr+1);
  if (maxx<=0) return;
  maxx=min(maxx,dparams.xres);
  if (minx>=maxx) return;

  int miny=int(y-rfacr+1);
  if (miny>=dparams.yres) return;
  miny=max(miny,0);

  int maxy=int(y+rfacr+1);
  if (maxy<=0) return;
  maxy=min(maxy,dparams.yres);
  if (miny>=maxy) return;

  p2[m].minx =minx;  p2[m].miny =miny;
  p2[m].maxx =maxx;  p2[m].maxy =maxy;

  p2[m].isValid = true;
  p2[m].x = x;
  p2[m].y = y;
  p2[m].r = r;
  p2[m].e.r = (float) p[m].e.r;
  p2[m].e.g = (float) p[m].e.g;
  p2[m].e.b = (float) p[m].e.b;
  p2[m].type = p[m].type;

  }


//colorize by kernel
__global__ void k_colorize
  (int n, cu_particle_splotch *p2, int mapSize, int types)
  {

  //first get the index m of this thread
  int m=blockIdx.x *blockDim.x + threadIdx.x;
  if (m >n) m =n;
 
  int ptype = p2[m].type;
  float col1=p2[m].e.r,col2=p2[m].e.g,col3=p2[m].e.b;
  clamp (0.0000001,0.9999999,col1);
  if (dparams.col_vector[ptype])
    {
    clamp (0.0000001,0.9999999,col2);
    clamp (0.0000001,0.9999999,col3);
    }
  float intensity=p2[m].I;
  clamp (0.0000001,0.9999999,intensity);
  intensity *= dparams.brightness[ptype];

  cu_color e;
  if (dparams.col_vector[ptype])   // color from file
    {
    e.r=col1*intensity;
    e.g=col2*intensity;
    e.b=col3*intensity;
    }
  else   // get color, associated from physical quantity contained in e.r, from lookup table
    {
  //first find the right entry for this ptype
      if (ptype<types)
      {
        e = get_color(ptype, col1, mapSize, types);
        e.r *= intensity;
        e.g *= intensity;
        e.b *= intensity;
      }
      else
      { e.r =e.g =e.b =0.0; }
    }
  p2[m].e=e;

  }



