# include "Galaxy.h"
#include"cxxsupport/arr.h"
#include"kernel/colourmap.h"
#include"cxxsupport/vec3.h"

float box_muller(float m, float s);
float box_uniform(float m, float s);
float pi=3.141592654;

// this is used to create spheroids of points: OPTION 1
long GaussRFunc (paramfile &params, string ComponentName, long number_of_points, 
		 float * coordx, float * coordy, float * coordz) 
{
  float sigma[3];
	
  srand(time(NULL));

  sigma[0] = params.find<float>("Sigmax"+ComponentName,0);
  sigma[1] = params.find<float>("Sigmay"+ComponentName,0);
  sigma[2] = params.find<float>("Sigmaz"+ComponentName,0);

  printf("    spheroid component with sigma_[x,y,z] = %f,%f,%f\n", sigma[0],sigma[1],sigma[2]);

// x/y/z coord

  for (long i=0; i<number_of_points; i++)
    {
      coordx[i] = box_muller(0.0, sigma[0]);
      coordy[i] = box_muller(0.0, sigma[1]);
      coordz[i] = box_muller(0.0, sigma[2]);
    }

  return number_of_points;

}

//this is used for the stars distribution in a face-on spiral galaxy: OPTION 3
/* 
MAIN PARAMETERS:
sigma_z = thickness of the disk (dispersion of the seed points)
Sigma_g = dispersion of point clouds around a "seed point"
n_per_pixel = number of points around each seed
*/
long GaussRDiscFunc (paramfile &params, string ComponentName, long number_of_points, long ntot, 
		 float * coordx, float * coordy, float * coordz, long nx, long ny) 
{
  srand(time(NULL));

  float sigma_z = params.find<float>("Sigmaz"+ComponentName,0.01);
  float sigma_g = params.find<float>("Sigmag"+ComponentName,0.001);
  long n_per_pixel = params.find<float>("NperPixel"+ComponentName,1);

  printf("      disk with sigma_z = %f, %d particles per pixel\n", sigma_z, n_per_pixel);

  long count = 0;

  float* coordx_save;
  float* coordy_save;
  coordx_save = new float [number_of_points];
  coordy_save = new float [number_of_points];

  for (long i=0; i<number_of_points; i++)
    {
      coordx_save[i] = coordx[i];
      coordy_save[i] = coordy[i];
    }

  float pixsizex = 4.0/(float)nx;
  float pixsizey = 4.0/(float)ny;
  float coordz_aux;
  for (long i=0; i<number_of_points; i+=1)
    {
      float x0=coordx_save[i];
      float y0=coordy_save[i];
      float distance = sqrt(x0*x0+y0*y0);
      float sigma_z_eff = sigma_z*(1.0-distance)*(1.0-distance);
      coordz_aux = box_muller(0.0, sigma_z_eff);
      
      ///count++;
      if(count >= ntot)
	{
	  printf("Generating more particles than allowed (%d>=%d)\n",count,ntot);
	  exit(3);
	}

      ///for(long k=1; k<n_per_pixel; k++)
      for(long k=0; k<3*n_per_pixel; k++)
	{
	  coordx[count] = box_uniform(x0, pixsizex);
	  coordy[count] = box_uniform(y0, pixsizey);
	  //coordx[count] = box_muller(x0, pixsizex);
	  //coordy[count] = box_muller(y0, pixsizey);
	  //coordx[count] = box_muller(x0, sigma_g);
	  //coordy[count] = box_muller(y0, sigma_g);
	  coordz[count] = 0.0;
          if(k<n_per_pixel)
	     coordz[count] = box_muller(coordz_aux, sigma_g);
/*
          if(coordz[count]-coordz_aux > sigma_g)
            {
              k--;
              continue;
            }
*/
	  count++;
	  if(count >= ntot)
	    {
	      printf("Generating more particles than allowed (%d>=%d)\n",count,ntot);
	      exit(3);
	    }
	}
    }

  delete [] coordx_save;
  delete [] coordy_save;

  return count;

}


//This is used for the gas for any galaxy: OPTION 4
/* 
MAIN PARAMETERS:
sigma = parameter to increase/decreas gas thickness in third dimension (with respect to the velocity dispersion)
sigma_fixed = dispersion of point clouds around the corresponding seed point
npergroup = number of points per seed
rx, ry = rescaled resolutions (in pixels)
*/

long RDiscFunc (paramfile &params, string ComponentName, long number_of_points, long ntot, 
		 float * coordx, float * coordy, float * coordz,
                 float * III, long nx, long ny) 
{
  srand(time(NULL));

  float sigma = params.find<float>("Sigmaz"+ComponentName,0);
  float sigma_fixed = params.find<float>("Sigmazfixed"+ComponentName,0.1);
  long npergroup = params.find<long>("NperGroup"+ComponentName,0);
  long rx = params.find<long>("Scaledxres"+ComponentName,1);
  long ry = params.find<long>("Scaledyres"+ComponentName,1);
  float * resolution;

  resolution = new float [rx*ry];
  long countin = 0;
  long countout = 0;

  for(long i=0; i<rx*ry; i++)
    resolution[i] = 0;

  for (long i=0; i<number_of_points; i++)    // find unique points related to low res image description
    {
      int ix = (int) round((0.5*(coordx[i]+1.0))*(rx));
      int iy = (int) round((0.5*(coordy[i]+1.0))*(ry));
      int irx = (int) round((0.5*(coordx[i]+1.0))*nx);
      int iry = (int) round((0.5*(coordy[i]+1.0))*ny);

      long index = ix+iy*rx;
      long rindex = irx+iry*nx;
      
      if(resolution[index] == 0.0)
	{
	  if(III[rindex] > 0.0)
	    {
	      coordz[i] = III[rindex];
	      resolution[index] = 1.0;
	      countin++;
	    }
	  else
	    {
	      coordz[i] = 0.0;
	      countout++;
	    }
	} 
      else 
	{
	  coordz[i] = 0.0;
	  countout++;
	}
    }

  delete [] resolution;

  long iaux = 0;

  float xmin=1e10;
  float xmax=-1e10;
  float ymin=1e10;
  float ymax=-1e10;
  for (long i=0; i<number_of_points; i++)    // copy pixel positions to the leading of the array
    {
      if(coordz[i] == 0.0)
	continue;
      coordx[iaux] = coordx[i];
      coordy[iaux] = coordy[i];
      coordz[iaux] = coordz[i];
      xmax = max(coordx[iaux],xmax);
      ymax = max(coordy[iaux],ymax);
      xmin = min(coordx[iaux],xmin);
      ymin = min(coordy[iaux],ymin);
      iaux++;
    }

  float ref_size = ((xmax-xmin)+(ymax-ymin))/2.0;
  cout << "    CHARACHTERISTIC SIZE = " << ref_size << "\n";

  long abscounter = countin;

  long ii = 0;
  float pixsizex = 4./rx;
  float pixsizey = 4./ry;
  while (abscounter < ntot-npergroup && ii < countin)
    {
      //original by Klaus
      //float gsigmaz = sigma * (sigma_fixed + coordy[ii]);
      //Claudio:
      float weight = coordz[ii];
      float thick = sqrt(weight)*ref_size;
      //float gsigmaz = sqrt(sigma) * 0.25* thick;
      float gsigmaz = sigma * thick;

      coordz[ii] = box_muller(0.0, gsigmaz);
      for(long jj=1; jj<(long)(weight*npergroup); jj++)
	{
	  //coordx[abscounter] = box_uniform(coordx[ii], pixsizex);
	  //coordy[abscounter] = box_uniform(coordy[ii], pixsizey);
	  coordx[abscounter] = box_muller(coordx[ii], pixsizex*1.0);
	  coordy[abscounter] = box_muller(coordy[ii], pixsizey*1.0);
	  float zaux;

	  do {zaux = box_muller(coordz[ii], sigma_fixed*gsigmaz);}
	  //////////do {zaux = box_muller(coordz[ii], sigma_fixed);}
	  //while (zaux*zaux > 0.25*ref_size*ref_size);
	  while (zaux*zaux > 0.25*ref_size*ref_size);
          coordz[abscounter] = zaux;
	  abscounter++;
	}
      ii++;
    }

  return abscounter;

}

//this is used for stars distribution in irregular galaxies: OPTION 5
/* 
MAIN PARAMETERS:
sigma_z - third dimension thickness
npergroup = number of points around each seed (group points)
ndiffuse_factor = number of diffused point per (group point)
sigma_fixed = dispersion of point clouds around the corresponding seed point: NOTE - IN PIXELS!!!!!!!
sigma_g = sigma_r*sigma_fixed is the dispersion of diffused points
*/

long GaussRGlobFunc (paramfile &params, string ComponentName, long number_of_points, long ntot, 
		 float * coordx, float * coordy, float * coordz, float * III, long nx, long ny) 
{
  srand(time(NULL));
  float * xx;
  float * yy;
  float * zz;

  float sigma_z = params.find<float>("Sigmaz"+ComponentName,0.8);
  long npergroup = params.find<long>("NperGroup"+ComponentName,0);
  float ndiffuse_factor = params.find<float>("NdiffuseFactor"+ComponentName,0);
  float sigma_fixed = params.find<float>("Sigmazfixed"+ComponentName,0.1);
  float sigma_g = params.find<float>("Sigmag"+ComponentName,1);


  float norm = 0.5*(float)nx;
  FILE * pFile;

// find local maxima in the pixels distribution and the charachteristic size of the galaxy

        float xmin=1e10;
        float xmax=-1e10;
        float ymin=1e10;
        float ymax=-1e10;
        float xm;
        float ym;
        float ref_size;
        float * max_mask = new float[nx*ny];
        int maxcount=0;
        for(long i; i<nx*ny; i++)max_mask[i]=0.0;
        for (int iy=1; iy<ny-1; iy++)
        for (int ix=1; ix<nx-1; ix++)
           {
             long index = ix+iy*nx;
             float IIImax=-0.5;
             for(int iiy=0; iiy<3; iiy++)
             for(int iix=0; iix<3; iix++)
               {
                   long index_aux = (ix-iix+1)+(iy-iiy+1)*nx;
                   IIImax = max(IIImax,III[index_aux]);
               }
             if(IIImax == III[index] && IIImax > 0.65)
               {
                   max_mask[index] = 1.0; 
                   //cout << ix<< " " << iy << " " << max_mask[index] << endl;
                   maxcount++;
               }

             if(III[index] > 0.1)
               {
                   xmax = max(float(ix),xmax);
                   ymax = max(float(iy),xmax);
                   xmin = min(float(ix),xmin);
                   ymin = min(float(iy),xmin);
               }
               xm=xmax-xmin;
               ym=ymax-ymin; 
               ref_size = min(fabs(xm),fabs(ym));    
           }

	cout << "Number of MAXIMA = " << maxcount << endl;
// clean the maxima distribution
       
        int stencil = 3;
        for (int iy=stencil; iy<ny-stencil; iy++)
        for (int ix=stencil; ix<nx-stencil; ix++)
           {

             long index = ix+iy*nx;
             if (max_mask[index] == 1.0)
               {
                 for(int iiy=0; iiy<2*stencil+1; iiy++)
                 for(int iix=0; iix<2*stencil+1; iix++)
                   {
                      long index_aux = (ix-iix+stencil)+(iy-iiy+stencil)*nx;
                      if (max_mask[index_aux] == 1.0 && index != index_aux)
                        {
                           if(III[index] >= III[index_aux])max_mask[index_aux]=0;
                           if(III[index] < III[index_aux])max_mask[index]=0;
                        }
                   }
               }
           }
// create random (x,y,z) point distributions around each point:
// create clusters around III local maxima
// create diffuse component around all bright pixels

        //long nfinal = npergroup*maxcount;
        long nfinal = npergroup*nx*ny;
        xx = new float [nfinal]; 
        yy = new float [nfinal]; 
        zz = new float [nfinal]; 
        float r_dist;
        long pcounter = 0;
        for (int iy=0; iy<ny; iy++)
        for (int ix=0; ix<nx; ix++)
           {
             long index = ix+iy*nx;
             float xcenter = float(ix-nx/2);
             float ycenter = float(iy-ny/2);
             float radius_min = sqrt(xcenter*xcenter+ycenter*ycenter);
        
             //float beta=1.0;
             //float sigma_L = -(1.0/(beta*beta*ref_size))*radius_min*radius_min+ref_size;
             //if(sigma_L < 0.0)sigma_L=0.0;
             //////float radius = box_uniform(0,4*sigma_L);
             //////float radius = sigma_L;

             long IIIth=0;
             float r_0;
             float x_ref;
             float x_x;
             float sigma_L;
             if (max_mask[index] == 1.0) 
               {
                 IIIth = long((npergroup/5)*III[index]*III[index]*III[index]+10); 
                 r_0 = ref_size*sigma_z;
                 //r_0 = ref_size*0.8;
                 x_ref = 0.0;
                 r_dist=sigma_fixed;
                 x_x = radius_min / r_0;
                 if(x_x < x_ref) x_x = x_ref; 
                 sigma_L = x_x * exp(-x_x) * r_0;

                 float zcenter = box_muller(0,sigma_L);
                 long nrandom = IIIth;   
                 for (int ii=0;ii<nrandom;ii++)
                 {
                   xx[pcounter]=box_muller(xcenter,r_dist);
                   yy[pcounter]=box_muller(ycenter,r_dist);
                   zz[pcounter]=box_muller(zcenter,r_dist);
                   pcounter++;
                 }
                 nrandom = long(ndiffuse_factor*npergroup);
                 for (int ii=0;ii<nrandom;ii++)
                 {
                   xx[pcounter]=box_muller(xcenter,sigma_g*r_dist);
                   yy[pcounter]=box_muller(ycenter,sigma_g*r_dist);
                   zz[pcounter]=box_muller(zcenter,sigma_g*r_dist);
                   pcounter++;
                 }
                 if(pcounter >= ntot)
                 {
                  printf("Generating more particles than allowed (%ld>=%ld)\n",pcounter,ntot);
                  exit(3);
                 }

               }
           }

           for (long ii=0; ii<pcounter; ii++) 
           {  
                coordx[ii] = xx[ii]/norm;
                coordy[ii] = yy[ii]/norm;
                coordz[ii] = zz[ii]/norm;
           }


// write data
/*
         pFile = fopen("test.dat", "w");
	 for (long ii=0; ii<pcounter; ii++)
             fprintf(pFile,"%f %f %f\n",xx[ii],
                                        yy[ii],
                                        zz[ii]);
         fclose(pFile);
*/

         delete [] xx;
         delete [] yy;
         delete [] zz;
         return pcounter;
}

//This is used for the gas for any tilted galaxy: OPTION 6
/* 
MAIN PARAMETERS:
TirificModel      = file containg parameters of the TiRiFiC model
PixelToTirific    = Conversion factor between pixel size and TiRiFiC (radial) units
TirificPartReduce = Reduction factor between region covered by TiRiFiC model and image size 
NperPixel         = fraction of particles with height > 0 (1 =all, N =1 every N)
Sigmaz            = disk thickness if not read from the model file (in pixels)
*/

long RDiscFuncTirific (paramfile &params, string ComponentName, long number_of_pixels, long ntot, 
		 float * coordx, float * coordy, float * coordz,
                 float * III, long nx, long ny) 
{
  srand(time(NULL));

  float pixeltotirific = params.find<float>("PixelToTirific"+ComponentName,1.0);
  float npartfix = params.find<float>("TirificPartReduce"+ComponentName,0.75);
  float extenddisk = params.find<float>("TirificExtendDisk"+ComponentName,1.0);
  float rmax = params.find<float>("RmaxMask"+ComponentName,nx/2.);
  long n_per_pixel = params.find<float>("NperPixel"+ComponentName,1);
  float sigma = params.find<float>("Sigmaz"+ComponentName,0);
  COLOURMAP model;

  ifstream infile (params.find<string>("TirificModel"+ComponentName).c_str());
  planck_assert (infile,"could not open palette file  <" + params.find<string>("TirificModel"+ComponentName) + ">");
  string dummy;
  int nModel;
  infile >> nModel;
  infile >> dummy >> dummy >> dummy >> dummy >> dummy;
  cout << "      loading " << nModel << " entries of tirific model table " << endl;
  float rrr,vvv,zzz,iii,ppp;
  for (int i=0; i<nModel; i++)
    {
      infile >> rrr >> vvv >> zzz >> iii >> ppp;
      model.addVal(rrr,COLOUR(zzz,iii,ppp));
      //      cout << rrr << " , " << zzz << " , " << iii << " , " << ppp << endl;
    }

  long ntrial = ntot * nx * ny / number_of_pixels * npartfix;
  float particlesize = nx * ny / (float)ntrial;
  cout << "ntot, nx, ny, number_of_pixels, npartfix, ntrial, particlesize = " <<
          ntot << " " <<
          nx << " " <<
          ny << " " <<
          number_of_pixels << " " <<
          npartfix << " " <<
          ntrial << " " <<
          particlesize << " " << endl ;

  cout << "M_PI = " << M_PI << endl;

  long icount = 0;
  float r0 = 0.0;

  cout << "NTRIAL " << ntrial << endl;
  for (long i=0; i<ntrial; i++)    // loop over all possible particles
    {
      //if(!(i%1000))cout << "R0 = " << r0 << endl;
      float Rnorm = r0/rmax;
      float r1m = sqrt( particlesize / 2 / M_PI + r0 * r0);     // radius associated to the particle 
      float phi = box_uniform(1.0, 2.0) * M_PI;                 // random phase of particle
  
      // Interpolate the tirific model to this radial distance
      COLOUR ring = model.getVal(r1m * pixeltotirific);
      float thick = ring.r;
      float pa = (ring.b + 180) / 180 * M_PI;
      float inc = ring.g / 180 * M_PI;

            //cout << i << " " << r1m << " " << thick << " " << pa << " " << inc << " " << phi << endl;


      // Calculate the normalized normal vector of the ring from the tirific model
      vec3 nn(0,0,1);
      vec3 m1_1(1,0,0),m1_2(0,cos(inc),-sin(inc)),m1_3(0,sin(inc),cos(inc));
      vec3 m2_1(cos(pa),-sin(pa),0),m2_2(sin(pa),cos(pa),0),m2_3(0,0,1);
      vec3 n1(dotprod(m1_1,nn),dotprod(m1_2,nn),dotprod(m1_3,nn));
      vec3 n(dotprod(m2_1,n1),dotprod(m2_2,n1),dotprod(m2_3,n1));
      n.Normalize();

      // Find the radius vector within the x/y plane
      float x0 = sqrt(r1m*r1m/(1+(n.x/n.y)*(n.x/n.y)));
      float y0 = -(n.x/n.y)*x0;
      vec3 rr(x0,y0,0);

      // Rotate the radius vector arround the normal vector by phi
      float sp = sin(phi);
      float cp = cos(phi);
      vec3 m3_1(n.x*n.x*(1-cp)+cp    ,n.x*n.y*(1-cp)-n.z*sp,n.x*n.z*(1-cp)+n.y*sp);
      vec3 m3_2(n.y*n.x*(1-cp)+n.z*sp,n.y*n.y*(1-cp)+cp    ,n.y*n.z*(1-cp)-n.x*sp);
      vec3 m3_3(n.z*n.x*(1-cp)-n.y*sp,n.z*n.y*(1-cp)+n.x*sp,n.z*n.z*(1-cp)+cp);
      vec3 xx(dotprod(m3_1,rr),dotprod(m3_2,rr),dotprod(m3_3,rr));

      // Use the thickness of the disk to displace the particle randomly withi the disk
      float height = 0.0;
      float height_aux = thick / pixeltotirific * extenddisk;
      if(sigma != 0.0)
      {
        height_aux=sigma;
        float sigma_cut = 1.5*sigma;
        float sigma_z_eff = height_aux*(1.0-Rnorm)*(1.0-Rnorm);
        if(!(icount%n_per_pixel)) 
        {
           do {
             height = box_muller(0, sigma_z_eff);
           } while (height > 1.5*sigma_z_eff);
        }

      } else {
        height = box_muller(0, height_aux);
      } 
      vec3 xxx = xx + n * height;

      // Convert coordinates into ranges [-1,1] ...
      vec3 xxxx(2*xxx.x/nx,2*xxx.y/ny,2*xxx.z/(0.5*(nx+ny)));

      // Find index in image
      long ix = 0.5 * (xxxx.x + 1.0) * nx;
      long iy = 0.5 * (xxxx.y + 1.0) * ny;
      long index = ix + iy * nx;

      if (index < nx * ny)
	if(III[index] > 0.0)
	  //	if(r1m < 100)
	  {
	    coordx[icount] = xxxx.x;
	    coordy[icount] = xxxx.y;
	    coordz[icount] = xxxx.z;
	    icount++;
	    if(icount >= ntot-1)
	      {
		cout << "More particles produced than intended, have to stop !!" << endl;
		cout << "Reduce the value of <TirificPartReduce"+ComponentName+"> (currently " << npartfix << ")" << endl; 
		cout << i << "/" << icount << "/" << ntot << " : " << r0 << " r1m = " << r1m << " t = " << thick << " pa = " << pa << " inc = " << inc << endl;
		exit(0);
	      }
	  }

      // This was the old idl part
      //      nn=[0,0,1.]
      //      m1=[[1,0,0],[0,cos(ic),-sin(ic)],[0,sin(ic),cos(ic)]]
      //      m2=[[cos(pa),-sin(pa),0],[sin(pa),cos(pa),0],[0,0,1]]
      //      n1=m1##nn
      //      n=m2##n1
      //      x0 = SQRT(r^2/(1+(n[0]/n[1])^2))
      //      y0 = SQRT(r^2 - x0^2)
      //      r0 = [x0,y0,0]
      //      sp=sin(phi)
      //      cp=cos(phi)
      //      m3=[[n[0]*n[0]*(1-cp)+cp     ,n[0]*n[1]*(1-cp)-n[2]*sp,n[0]*n[2]*(1-cp)+n[1]*sp],$
      //          [n[1]*n[0]*(1-cp)+n[2]*sp,n[1]*n[1]*(1-cp)+cp     ,n[1]*n[2]*(1-cp)-n[0]*sp],$
      //          [n[2]*n[0]*(1-cp)-n[1]*sp,n[2]*n[1]*(1-cp)+n[0]*sp,n[2]*n[2]*(1-cp)+cp     ]]
      //      xx=m3##r0

      // Update covered radii
      r0 = sqrt( particlesize / M_PI + r0 * r0);
      if(r0 > rmax){
        cout << "FINAL RADIUS" << r0 << endl; 
	break;
      }
    }
  return icount;

}


