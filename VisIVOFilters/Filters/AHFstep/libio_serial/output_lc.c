/* The lightcone output routine.
 * We write the same header as the normal output, 
 * to be able to use the standard AMIGA tools
 * (note that the buffer is larger anyway, to
 * include z). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void dummy_output_lc()
{
}

#ifdef LIGHTCONE

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "io_serial.h"
#include "../libutility/utility.h"

#ifdef MULTIMASS
#define BUFSIZE 8
#else
#define BUFSIZE 7
#endif

/* write particles to a lightcone file, with a redshift extension */
void output_lc()
{
  char    file_ncpy[100];     /* copy of file name                      */
  char    file_no[10];        /* file number                            */
  char   *f_name;             /* file name pointer                      */
  int     i;		      /* loop index                             */
  int     NC,NCB;	      /* cube index range and number of cubes   */
  double  r0,r1,r0sq,r1sq,dr; /* previous and current lightcone radii   */
  double  z0,z1;	      /* redshifts: prev and curr               */
  cubeptr cur_cube;	      /* current cube                           */
  int     npart;	      /* no of particles written                */
  partptr cur_part;           /* current particle                       */
  bckptr  backup_part;        /* its coordinates at the previous moment */
  double  x,y,z;	      /* particle coordinates                   */
  double  xb,yb,zb;	      /* previous coordinates                   */
  double  dx,dy,dz;	      /* coordinate diffs and coords in cubes   */
  double  dxb,dyb,dzb;	      /* previous coords in cubes               */
  double  ex,ey,ez,eysq,ersq; /* rotated coords                         */
  double  prsq,prsqb,pr,prb;  /* distances from the observer            */
  double  a,a1;		      /* interpolation variables                */
  char    machine_sizeof_long;
  flouble buf[BUFSIZE];	      /* particle data for the lightcone        */

  if(io.lightcone==NULL) {		/* file not open */
    f_name=strcpy(file_ncpy,io.lightcone_prefix);
    sprintf(file_no,"z%.4f",fabs(global.z));
    f_name=strcat(f_name, file_no);
    if((io.lightcone=fopen(f_name,"wb"))==NULL) {
      fprintf(io.logfile,"output: could not open file %s\n", f_name);
      exit(1); 
    }
		
    /* write a simple "1" to file (for BYTESWAP testing) */
    machine_sizeof_long = sizeof(long);
    fwrite(&machine_sizeof_long, sizeof(int), 1, io.lightcone);
		
    /* write the IO header */
    io.header.K_current   = energy.K_current;
    io.header.U_current   = energy.U_current;
    io.header.Eintegral   = energy.integral;
    io.header.Econst      = energy.econst;
    io.header.no_timestep = global.no_timestep;
    io.header.a_current = global.a;
    if(fwrite(&(io.header),sizeof(io.header),1,io.lightcone)!=1) {
      fprintf(io.logfile,"\n\noutput: could not write header\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
    }
  }
  /* the real work */
  r0=io.rcone0; r0sq=r0*r0;
  z0=io.z0;
  r1=r_cone(global.a);
  r1sq=r1*r1; 
  z1=global.z;
  dr=r0-r1;
  NC=r1+2;		/* r1 is in internal units, 2 for safety */

  /* create the list of relevant cubes */
  /* allocate memory, if necessary (the first shell) */
  if(io.fst_cube==NULL) {
    NCB=8*NC*NC*NC;		/* cube, -NC..NC-1 */
    if((io.fst_cube=(cubeptr)calloc(NCB,sizeof(cube)))==NULL) {
      fprintf(stderr,"Cannot allocate memory for the cube list\n");
      exit(1);
    }
  }
  /* find the cubes, which cross the lightcone */
  switch(io.conetype) {
  case -1:
    create_slice(NC,r0sq,r1sq); 	/* sky patch with args in struct io */
    break;
  case 1:
    create_cone(-NC,-NC,-NC,NC,r0sq,r1sq);/* full-sky cone */
    break;
  case 2:
    create_cone(-NC,-NC,0,NC,r0sq,r1sq);	/* half-sky cone */
    break;
  case 8:
    create_cone(0,0,0,NC,r0sq,r1sq);		/* 1/8 of the sky */
    break;
  default:
    fprintf(io.logfile,"unsupported lightcone request %d in output_lc\n",
	    io.conetype);
    fclose(io.logfile);
    exit(1);
  }

  /* loop over particles -> for each loop over cubes and write, if neccessary */
  cur_part=global.fst_part;
  backup_part=io.fst_backup_part;
  for(i=npart=0;i<global.no_part;i++,cur_part++,backup_part++) {
    x=cur_part->pos[0]; y=cur_part->pos[1]; z=cur_part->pos[2];
    xb=backup_part->pos[0]; yb=backup_part->pos[1]; zb=backup_part->pos[2];
    /* periodicity */
    if((dx=x-xb)<-0.5) xb-=1.; else if(dx>0.5) xb+=1.;	
    if((dy=y-yb)<-0.5) yb-=1.; else if(dy>0.5) yb+=1.;	
    if((dz=z-zb)<-0.5) zb-=1.; else if(dz>0.5) zb+=1.;
    buf[1]=cur_part->mom[0]; /* the leapfrog v is constant from xb to x */
    buf[3]=cur_part->mom[1]; buf[5]=cur_part->mom[2];
#ifdef MULTIMASS
    buf[7]=cur_part->weight;
#endif
    for(cur_cube=io.fst_cube;cur_cube->used;cur_cube++) {
      dx=x+cur_cube->I; dy=y+cur_cube->J; dz=z+cur_cube->K;
      dxb=xb+cur_cube->I; dyb=yb+cur_cube->J; dzb=zb+cur_cube->K;
      prsq=dx*dx+dy*dy+dz*dz;
      prsqb=dxb*dxb+dyb*dyb+dzb*dzb;
      if(prsqb<r0sq && prsq>=r1sq) {	/* crosses the lightcone */
	/* linear interpolation in all lengths! */
	pr=sqrt(prsq); prb=sqrt(prsqb);
	a=(r0-prb)/(pr-prb+dr);				
	a1=1.-a;
	dx=a*dx+a1*dxb;			/* reusing dx,dxb */
	dy=a*dy+a1*dyb;
	dz=a*dz+a1*dzb;
	if(io.conetype==-1) {	/* check euler() for explanation */
	  ey=-(io.R[0][0]*dx+io.R[0][1]*dy+io.R[0][2]*dz);
	  ex=-(io.R[1][0]*dx+io.R[1][1]*dy+io.R[1][2]*dz);
	  ez=io.R[2][0]*dx+io.R[2][1]*dy+io.R[2][2]*dz;
	  eysq=ey*ey; ersq=ex*ex+eysq+ez*ez;
	  /* spherical geometry! */
	  if(ez<-ZERO || ex*ex>(ersq-eysq)*io.xcoef || eysq>ersq*io.ycoef) 
	    continue;		/* the point is not in the patch */
	  buf[0]=ex;
	  buf[2]=ey;
	  buf[4]=ez;
	}
	else {
	  buf[0]=dx;
	  buf[2]=dy;
	  buf[4]=dz;
	}
	buf[6]=a*z1+a1*z0;	/* redshift */
	fwrite(buf,BUFSIZE*sizeof(flouble),1,io.lightcone);
	npart++;
      }
    }
  }
  fflush(io.lightcone);
  if(npart) fprintf(io.logfile,"\nwrote %d particles to lightcone\n",npart);
  fflush(io.logfile);
}

/* this function is called from main() */
void store_backup()
{
  partptr cur_part;           /* current particle  */
  bckptr backup_part;         /* its coordinates at the previous moment */
  int i,j;

  backup_part=io.fst_backup_part;
  cur_part=global.fst_part;
  for(i=0;i<global.no_part;i++) {
    for(j=X;j<=Z;j++) 
      backup_part->pos[j] = cur_part->pos[j];
    backup_part++; cur_part++;
  }
}

#endif

