/*
Copyright (c) 2012, MB
All rights reserved.
*/

#include <iostream>
#include <math.h>
#include "vsLine.h"
#include "vsVector.h"

#define SMALL_NUM  0.00000000001 // anything that avoids division overflow

using namespace std;

/**********************************
Line::Line
costruttore
**********************************/

Line::Line(){}
Line::Line(Vector p1, Vector p2){
	
	P0=p1;
	P1=p2;
	u=P1.diffP(&P0); //calcolo di u=P1-P0	
	u_n=u.norm();
	
}

Line::~Line(){};

/**********************************
Line::poca
input:due rette
output: un punto di scattering
**********************************/

Vector Line::poca(Line *l){

	Vector Psc, Qtc, Usc, Vtc, w0,scattP;
	
	w0=P0.diffP(&l->P0); //w0=P0-Q0
	double a= u.dotProduct(&u); //a=u*u
	double b= u.dotProduct(&l->u); //b=u*v
	double c= l->u.dotProduct(&l->u); //c=v*v
	double d= u.dotProduct(&w0); //d=u*w0
	double e= l->u.dotProduct(&w0); //e=v*w0
 	double D = a*c - b*b;       
	double sc, tc;

        double moduloU= u.module();
	double moduloV=l->u.module();

	/*
	calcolo del cosTheta senza segno	
	double cosTheta=b/(moduloU*moduloV);
	//if(cosTheta>1) cosTheta=1; //controllo sul valore numerico di cosTheta 
	//if (cosTheta<-1) cosTheta=-1;
	*theta=acos(cosTheta);
	*/
	
	 //se le rette sono quasi parallele ---> non c'è scattering
	if (D<SMALL_NUM) 
	{      	
		//restituisce come pto di scattering P0
		sc=0;
		tc=e/c;
	
	}
    	else 
	{
		sc = (b*e - c*d)/D;
		tc = (a*e - b*d)/D;
    	}

	
	Usc=u.multP(sc); //usc=u*sc
	Vtc=l->u.multP(tc); //vtc=v*tc	
	
	Psc=P0.addP(&Usc); //Psc=P0+usc
	Qtc=l->P0.addP(&Vtc); //Qtc=Q0+vtc
	
	scattP=Psc.mediumP(&Qtc); //scattering Vector= punto medio Psc e Qtc
	return scattP;
}

/**********************************
Line::scatteringAngle
input:due rette
output:angolo di scattering in deg
***********************************/
double Line::scatteringAngle(Line *l){
    
    
    
    
    
    Vector v_r1, v_r2;
    double theta=0,dumb=0,dotUV=0;
    //double cosTheta=0,senTheta=0;
    
    //1.
    //ricavo i versori (normalizzati) delle due rette r1 e r2
    
    
    
    v_r1=u_n;
    v_r2=l->u_n;
    dotUV=v_r1.dotProduct(&v_r2);
    dumb=dotUV/(v_r1.module()*v_r2.module());
    if(dumb>=1) theta=0;
    else theta=acos(dumb);
    if (dotUV<0){
        if(theta>0) theta+=M_PI; else theta-=M_PI;
    }
    
    theta=theta*180/M_PI;
    return theta;
}









