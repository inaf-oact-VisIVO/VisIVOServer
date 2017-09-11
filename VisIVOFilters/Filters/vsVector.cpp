/*
Copyright (c) 2012, MB
All rights reserved. DA FARE: creare la classe Point come ereditata da Vector!!!
*/

#include "vsVector.h"
#include <iostream>
#include <math.h>
#ifndef PASSO_DISC
#define PASSO_DISC 0.5
#endif



using namespace std;


Vector::Vector()
	{
		x=0;
		y=0;
		z=0;
	}
Vector::Vector(float c1, float c2, float c3)
	{
	x=c1;
	y=c2;
	z=c3;		
	}

Vector::~Vector()
{}

/**********************************
Vector::dotProduct= prodotto scalare
input: 2 Vector
output: prodotto scalare
**********************************/
float Vector::dotProduct(Vector *v){

	float prod;
	prod=x*v->x+y*v->y+z*v->z;
	return prod;

}

/**********************************
Vector::module= modulo di un vettore
input: 1 vettore
output: modulo
**********************************/
float Vector::module(){
	float mod;
	mod=sqrt(x*x+y*y+z*z);
	return mod;
}

/**********************************
Vector::norm= calcolo del versore di un vettore
input: 1 vettore
output: versore
**********************************/

Vector Vector::norm(){
	Vector p(0, 0, 0);
	float d=sqrt(x*x+y*y+z*z);
	if(d!=0)
	{
		p.x=x/d;
		p.y=y/d;
		p.z=z/d;
	}
	return p;
}




Vector Vector::diffP(Vector *p1){ //differenza tra vettori

	Vector p;
	p.x=x-p1->x;
	p.y=y-p1->y;
	p.z=z-p1->z;
	return p;
}

Vector Vector::addP(Vector* p1){
	Vector p;
	p.x=x+p1->x;
	p.y=y+p1->y;
	p.z=z+p1->z;
	return p;
}

Vector Vector::multP(double c){
	Vector p;	
	p.x=x*c;
	p.y=y*c;
	p.z=z*c;
	return p;
}	

Vector Vector::divP(double c){
	Vector p(0,0,0);
	if(c!=0)
	{	
		p.x=x/c;
		p.y=y/c;
		p.z=z/c;
	} return p;
	
}

void Vector::Eulerodiscretization(Vector *p){ //Eulero in avanti
	
	if(fabs((x-(int)x))>PASSO_DISC)
	{
		if(x>0)	
		p->x=(int)x+PASSO_DISC;
		else
		p->x=(int)x-PASSO_DISC;
	}
	else
	p->x=(int)x;
	
	if(fabs((y-(int)y))>PASSO_DISC)
	{
		if(y>0)
		p->y=(int)y+PASSO_DISC;
		else
		p->y=(int)y-PASSO_DISC;
	}
	else
	p->y=(int)y;
	
	if(fabs((z-(int)z))>PASSO_DISC)
	{
		if(z>0)
		p->z=(int)z+PASSO_DISC;
		else 
		p->z=(int)z-PASSO_DISC;
	}	
	else
	p->z=(int)z;

}


void Vector::discretization(Vector *p, int dim_x, int dim_y, int dim_z){

	//discretizzazione con passo 0.5 centrata in un intorno del punto di approssimazione
	//p->x=x+300;
	//p->y=y+150;
	int xi, yi, zi;
	xi=(int)(x+dim_x/2);
	yi=(int)(y+dim_y/2);
	zi=(int)(z+dim_z/2);
	
	//discretizzazione asse x
	if(fabs(x)<(fabs(xi)+0.25)) p->x=xi-0.5;
	else	
	if(fabs(x)>(fabs(xi)+0.25)&&fabs(x)<(fabs(xi)+0.75)) p->x=xi;	
	else
	if(fabs(x)>(fabs(xi)+0.75)) p->x=xi+0.5;
	
	//discretizzazione asse y
	if(fabs(y)<(fabs(yi)+0.25)) p->y=yi-0.5;
	else	
	if(fabs(y)>(fabs(yi)+0.25)&&fabs(y)<(fabs(yi)+0.75)) p->y=yi;	
	else
	if(fabs(y)>(fabs(yi)+0.75)) p->y=yi+0.5;
		
	//discretizzazione asse z
	if(fabs(z)<(fabs(zi)+0.25)) p->z=zi-0.5;
	else	
	if(fabs(z)>(fabs(zi)+0.25)&&fabs(z)<(fabs(zi)+0.75)) p->z=zi;	
	else
	if(fabs(z)>(fabs(zi)+0.75)) p->z=zi+0.5;
}

int Vector::Point2Voxel(int dim_x, int dim_y, int dim_z, int dim_voxel)
{
		
	int i;	
	//aggiungo alla coordinata x, DIM_X/2 per avere solo punti con coordinate positive [0, 599];
	int a= int((x+dim_x/2)/dim_voxel); 
	int b= int((y+dim_y/2)/dim_voxel)*dim_x/dim_voxel;
	int c= int((z+dim_z/2)/dim_voxel)*dim_x*dim_y/(dim_voxel*dim_voxel);
	
	i=a+b+c;
	return i;	
}

Vector Vector::mediumP(Vector *p1){
	
	Vector p;
	p.x=(x+p1->x)/2;
	p.y=(y+p1->y)/2;
	p.z=(z+p1->z)/2;
	return p;
}



/*int main()

{
	Vector p1=Vector(1, 1, 1);
	std::cout<<"P1: "<<p1.x<<" "<<p1.y<<" " <<p1.z<<" !\n";	
	cout<<"ok!\n";
}

*/



