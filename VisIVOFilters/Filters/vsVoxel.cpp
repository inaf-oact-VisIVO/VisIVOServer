/*

Copyright (c) 2012, MB
All rights reserved.

Un voxel è caratterizzato da un indice nel vettore, ha dimensioni fisiche DIM_VOXEL*DIM_VOXEL*DIM_VOXEL
ed è caratterizzato da un valore double che rappresenta al somma dei ThetaQuadri dei punti che 
cadono nel volume
*/

#include <iostream>
#include "vsVector.h"
#include "vsVoxel.h"
#include <math.h>
#include <vector>



using namespace std;

Voxel::Voxel(){

	nPts=0;
	thetaM=0;
	sumTheta=0;
	thetaQuadroM=0;
	sumThetaQuadro=0;
	sigma=0;
	error=0;

}

Voxel::~Voxel(){}

//funzione che aggiunge il theta corrispondente ad un punto al voxel
void Voxel::addTheta(double theta)
{
	sumTheta=sumTheta+theta;
}

void Voxel::addThetaQuadro(double tQ)
{
	sumThetaQuadro=sumThetaQuadro+tQ;
}

void Voxel::memTheta(double c){

	thetaVector.push_back(c); //inserisce in coda
}
	
void Voxel::setThetaM(){
	if(nPts==0)
	thetaM=0;
	else
	thetaM=sumTheta/nPts;
	}

void Voxel::setThetaQuadroM(){ 
	if(nPts==0)
	thetaQuadroM=0;
	else
	thetaQuadroM=sumThetaQuadro/nPts;
	}

void Voxel::setSigma(){ //calcolo della deviazione standard: sigma=sqrt((sum(x_i-x_medio)^2)/N_eventi)
	
	if(nPts==0)
	sigma=0;
	else
	{
		while(!thetaVector.empty())
		{
		sigma+=pow((thetaVector.back()-thetaM), 2);
		thetaVector.pop_back();
		}
		sigma=sqrt(sigma/nPts);
	}
}
void Voxel::setError(){
	if(nPts==1 || nPts==0)
	error=1;
	else	
	error=1/sqrt(2*(nPts-1));
	}


/*
int main()
{
float x, y, z;
Vector p1;
cout<<"Inserisci le tre coordinate del punto:\n";
cin>>x;
cin>>y;
cin>>z;

p1=Vector(x, y, z);
int index;

index=p1.Vector2Voxel();

cout<<"Al punto corrisponde il voxel n. "<<index<< "\n";



}*/
