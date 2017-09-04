#ifndef H_SPLOTCH_DU_DATA
#define	H_SPLOTCH_DU_DATA

const int   nParticle =512*512; //const causes unroll!

struct  G_VARS{
    float   rfac, bfak,i00,sigma0;
    float   brightness, grayabsorb;	//only for now, type all same
    int     res, ycut0, ycut1;
	float	zmaxval,zminval;
    //for xexp function
};

struct PARTICLE
{
  float x,y,z,r,ro,I,T;
  int type;
  //jin: for CUDA parallel
  //p_vars which is the result of calculation
//  int posFBuf; //start pos in fragment buffer. No use in Plan A
  int	minx, maxx, miny, maxy;
  float prefac1, prefac2;
//  COLOUR	q, a;
  float	q[3], a[3];
};

/*No use in Plan A
struct PIXEL
{
	float	r,g,b;
};
*/

struct	FRAGMENT
{
	float	k_red, c_red;
	float	k_green, c_green;
	float	k_blue, c_blue;
};

//typedef FRAGMENT* FRAGMENT_PTR;
#endif