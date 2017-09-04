/* The routines to select the cubes
 * which contain parts of the lightcone. */
#include <stdio.h>
#include <stdlib.h>

void dummy_gen_cubes()
{
}

#ifdef LIGHTCONE 

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "utility.h"
#include "../libio_serial/io_serial.h"


/* creates the list of cubes, which
 * intersect a simple lightcone */
void create_cone(X1,Y1,Z1,XYZ2,rsq1,rsq2)
int X1,Y1,Z1;	        /* minimum cell indices */
int XYZ2;		/* the common max cell index */
double rsq1,rsq2;	/* the lightcone shell radii squared */
{
        cubeptr cur_cube;			/* current cube */
        int i,j,k;					/* indices */
        double x,y,z;				/* cube vertices */

	cur_cube=io.fst_cube;
	for(i=X1;i<XYZ2;i++)
	for(j=Y1;j<XYZ2;j++)
	for(k=Z1;k<XYZ2;k++) {
		x=i;y=j;z=k;		
		if(!radius_check(x,y,z,rsq1) && !radius_check(x,y,z,rsq2))
		  continue;	 /* the cube does not cover the lightcone shell */
		cur_cube->I=i; cur_cube->J=j; cur_cube->K=k;
		cur_cube->used=1;
		cur_cube++;
	}
	cur_cube->used=0;
}


/* creates the list of cubes, which intersect a slice-type lightcone */
#define EPS 0.0001	/* to avoid -0 situations */
/* these static structures are needed for several subroutines */
double vert[8][3];			/* cube vertices, rotated coords */
double *face[6][4];			/* pointers to vertices, facewise */
double *edge[12][2];		        /* pointers to vertices, edgewise */

void create_slice(N,rsq1,rsq2)
int N;				/* the max cell index */
double rsq1,rsq2;	        /* the lightcone shell radii squared */
{
        cubeptr cur_cube;		/* current cube */
        int i,j,k;			/* indices */
        double x,y,z;			/* vertices */

	cur_cube=io.fst_cube;
	for(i=-N;i<N;i++)
	for(j=-N;j<N;j++)
	for(k=-N;k<N;k++) {
		x=i;y=j;z=k;		
		/* check for the radius, first, as in create_cone */
		if(!radius_check(x,y,z,rsq1) && !radius_check(x,y,z,rsq2))
		  continue;	 /* the cube does not cover the lightcone shell */
		/* now check if any vertex is inside the sky patch */
		if(vertex_check(x,y,z)) goto found; 	
		/* check if the patch edges cross the faces of the cube */
		if(face_check()) goto found; 		
		/* check if cube edges intersect with patch faces */
		if(!edge_check()) continue;
	found: 		/* found a cube */
		cur_cube->I=i; cur_cube->J=j; cur_cube->K=k;
		cur_cube->used=1;
		cur_cube++;
	}
	cur_cube->used=0;
}

/* check if the cube starting at x,y,z intersects with
 * the lightcone shell between the previous and present radii. */
int radius_check(x,y,z,rsq)
double x,y,z;
double rsq;
{
	int s;
	double xsq,ysq,zsq,x1sq,y1sq,z1sq;

	xsq=x*x; ysq=y*y; zsq=z*z;
	x1sq=xsq+2*x+1; y1sq=ysq+2*y+1; z1sq=zsq+2*z+1;
	s=  (xsq+ysq+zsq<rsq)? 0:1;
	s+= (xsq+ysq+z1sq<rsq)? 0:1;
	s+= (xsq+y1sq+zsq<rsq)? 0:1;
	s+= (xsq+y1sq+z1sq<rsq)? 0:1;
	s+= (x1sq+ysq+zsq<rsq)? 0:1;
	s+= (x1sq+ysq+z1sq<rsq)? 0:1;
	s+= (x1sq+y1sq+zsq<rsq)? 0:1;
	s+= (x1sq+y1sq+z1sq<rsq)? 0:1;
	if(s==0 || s==8) return 0;	/* the cube does not cover the lightcone */
	else return 1;
}

/* check if the vertices of a cube starting at x,y,z
 * are inside the sky patch. Also fills the vertices array. */
int vertex_check(x,y,z)
double x,y,z;
{
	int i,j,k,iv;
	double x1,y1,z1,ysq,rsq;

	for(i=iv=0;i<=1;i++)
	for(j=0;j<=1;j++)
	for(k=0;k<=1;k++) {
		euler(x+i,y+j,z+k,&x1,&y1,&z1);
		ysq=y1*y1;
		rsq=x1*x1+ysq+z1*z1;
		if(z1>-EPS && x1*x1<=(rsq-ysq)*io.xcoef && ysq<=rsq*io.ycoef) 
		  return 1;			/* vertex in the sky patch */
		vert[iv][0]=x1; vert[iv][1]=y1; vert[iv][2]=z1;
		iv++;
	}
	return 0;
}

/* checks if the edges of the patch pyramid
 * (assuming flat faces, might allow extra cubes)
 * intersect the cube faces */
int face_check()
{
	int i,j,k,k1;
	double x1,y1,z1,x2,y2,z2,x3,y3,z3;
	double A,B,C;
	double xc,yc,zc;
	double dx1,dy1,dz1,dx2,dy2,dz2;
	double nk,u,rsq1,rsq2,tota,cosa;

	/* set up the faces, first */
	face[0][0]=vert[0]; face[0][1]=vert[1]; 
	face[0][2]=vert[3]; face[0][3]=vert[2];
	face[1][0]=vert[2]; face[1][1]=vert[3]; 
	face[1][2]=vert[7]; face[1][3]=vert[6]; 
	face[2][0]=vert[0]; face[2][1]=vert[2]; 
	face[2][2]=vert[6]; face[2][3]=vert[4]; 
	face[3][0]=vert[0]; face[3][1]=vert[1]; 
	face[3][2]=vert[5]; face[3][3]=vert[4];
	face[4][0]=vert[4]; face[4][1]=vert[5];
	face[4][2]=vert[7]; face[4][3]=vert[6];
	face[5][0]=vert[1]; face[5][1]=vert[3];
	face[5][2]=vert[7]; face[5][3]=vert[5];
	for(i=0;i<6;i++) {	/* loop over faces */
	  /* n=(A,B,C) is the normal to a face */
		x1=face[i][0][0]; y1=face[i][0][1]; z1=face[i][0][2];
		x2=face[i][1][0]; y2=face[i][1][1]; z2=face[i][1][2];
		x3=face[i][2][0]; y3=face[i][2][1]; z3=face[i][2][2];
		A=y1*(z2-z3)+y2*(z3-z1)+y3*(z1-z2);
		B=z1*(x2-x3)+z2*(x3-x1)+z3*(x1-x2);
		C=x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2);
		for(j=0;j<4;j++) { /* loop over sample pyramid edges */
		  /* n.(x-v0)=0, x=uk -> u=(n.v0)/(n.k) */
			if((nk=A*io.kpatch[j][0]+B*io.kpatch[j][1]+C*io.kpatch[j][2])<EPS)
			  continue;	/* parallel, will intersect another face */
			u=(A*x1+B*y1+C*z1)/nk;
			xc=u*io.kpatch[j][0]; yc=u*io.kpatch[j][1]; zc=u*io.kpatch[j][2];
			if(zc<-EPS) continue;		/* slices have z>0 */
			/* if sum angle(v_i-x,v_j-x)== 2pi, x is inside the face */
			for(k=0,tota=0.;k<4;k++) { /* loop over face vertex pairs */
				if((k1=k+1)==4) k1=0;
				dx1=face[i][k][0]-xc; dy1=face[i][k][1]-yc; 
				dz1=face[i][k][2]-zc; 
				dx2=face[i][k1][0]-xc; dy2=face[i][k1][1]-yc; 
				dz2=face[i][k1][2]-zc; 
				rsq1=dx1*dx1+dy1*dy1+dz1*dz1;
				rsq2=dx2*dx2+dy2*dy2+dz2*dz2;
				if(rsq1<EPS || rsq2<EPS) return 1;	/* found a vertex */
				cosa=(dx1*dx2+dy1*dy2+dz1*dz2)/sqrt(rsq1*rsq2);
				tota+=acos(cosa);
			}
			if(tota>=2.*PI-EPS) return 1;		/* xc inside the face */
		}
	}
	return 0;
}

/* checks if the cube edges intersect the patch pyramid planes
 * (assuming a flat-face pyramid, might allow extra cubes) */
int edge_check()
{
	int i,j;
	double x1,y1,z1,x2,y2,z2;
	double ne1,ne2,u,u1,xc,yc,zc,rsq;

	/* set up the edges, first */
	edge[0][0]=vert[0]; edge[0][1]=vert[2];
	edge[1][0]=vert[2]; edge[1][1]=vert[6];
	edge[2][0]=vert[6]; edge[2][1]=vert[4];
	edge[3][0]=vert[4]; edge[3][1]=vert[0];
	edge[4][0]=vert[1]; edge[4][1]=vert[3];
	edge[5][0]=vert[3]; edge[5][1]=vert[7];
	edge[6][0]=vert[7]; edge[6][1]=vert[5];
	edge[7][0]=vert[5]; edge[7][1]=vert[1];
	edge[8][0]=vert[0]; edge[8][1]=vert[1];
	edge[9][0]=vert[2]; edge[9][1]=vert[3];
	edge[10][0]=vert[4]; edge[10][1]=vert[5];
	edge[11][0]=vert[6]; edge[11][1]=vert[7];
	/* n.x=0, x=e1+u*(e2-e1) -> u=n.e1/(n.e1-n.e2) */
	for(i=0;i<12;i++) { 		/* loop over cube edges */
		x1=edge[i][0][0]; y1=edge[i][0][1]; z1=edge[i][0][2];
		x2=edge[i][1][0]; y2=edge[i][1][1]; z2=edge[i][1][2];
		for(j=0;j<4;j++) {		/* loop over pyramid planes */
			ne1=x1*io.npatch[j][0]+y1*io.npatch[j][1]+z1*io.npatch[j][2];
			ne2=x2*io.npatch[j][0]+y2*io.npatch[j][1]+z2*io.npatch[j][2];
			if(ne1*ne2>=0) continue; 	/* at the same side of the plane */
			u=ne1/(ne1-ne2); u1=1.-u;
			xc=u1*x1+u*x2; yc=u1*y1+u*y2; zc=u1*z1+u*z2;
			rsq=xc*xc+yc*yc+zc*zc;
			if(zc>-EPS && xc*xc<=rsq*io.xcoef && yc*yc<=rsq*io.ycoef)
			  return 1;		/* intersection point in the patch */
		}
	}
	return 0;
}
		
/* rotates the coordinates to the patch system:
 * the z-axis is directed towards the centre of the sky patch,
 * the y-axis lies along the shortest axis of the patch
 * and the x-axis --along the longest. */
void euler(x,y,z,px,py,pz)
double x,y,z;			/* old coords */
double *px,*py,*pz;		/* pointers to new coords */
{
  /* Euler rotation +
   * (rot + refl) in xy-plane: *px=-y1; *py=-x1; */
	*py=-(io.R[0][0]*x+io.R[0][1]*y+io.R[0][2]*z);
	*px=-(io.R[1][0]*x+io.R[1][1]*y+io.R[1][2]*z);
	*pz=io.R[2][0]*x+io.R[2][1]*y+io.R[2][2]*z;
}
#endif /*LIGHTCONE*/



