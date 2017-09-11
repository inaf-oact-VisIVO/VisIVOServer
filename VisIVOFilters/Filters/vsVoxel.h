/*

Copyright (c) 2012, MB
All rights reserved.

Un voxel è caratterizzato da un indice nel vettore, ha dimensioni fisiche DIM_VOXEL*DIM_VOXEL*DIM_VOXEL
ed è caratterizzato da un valore double (sumTheta) che rappresenta al somma dei ThetaQuadri dei punti che 
cadono nel volume
*/
#include <iostream>
#include <vector>
#include <math.h>

#ifndef VOXEL_H
#define VOXEL_H

using namespace std;

class Voxel{

	
	private:
	
	double thetaM, sumTheta, thetaQuadroM, sumThetaQuadro,sigma, error;
	vector<double> thetaVector;

	public:
	int nPts;
	
	
	Voxel();
	~Voxel();
	
	void addTheta(double theta);
	void addThetaQuadro(double tQ);
	void memTheta(double theta);	
	
	double getSumTheta(){return sumTheta;}
	double getSumThetaQuadro(){return sumThetaQuadro;}
	double getThetaM(){return thetaM;}
	double getThetaQuadroM(){return thetaQuadroM;}
	double getSigma(){return sigma;}
	double getError(){return error;}

	void setThetaM();
	//void setThetaQuadroM();
	void setThetaQuadroM();
	void setSigma();
	void setError();
	
	//int Vector2Voxel(Vector *p);

};

#endif
