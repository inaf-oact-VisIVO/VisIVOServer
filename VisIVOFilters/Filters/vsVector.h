/*
Copyright (c) 2012, MB
All rights reserved.
*/
#ifndef Vector_H
#define Vector_H

class Vector{

	public:
		float x, y, z; //coordinate del punto-vettore
		Vector(); //costruttore
		Vector(float, float, float); //costruttore
		~Vector();

		Vector mediumP(Vector*); //punto medio tra due punti
		int Point2Voxel(int dim_x, int dim_y, int dim_z, int dim_voxel); //trasformare un punto in un voxel
		float dotProduct(Vector *); //prodotto scalare Vettore-Vettore
		float module(); //modulo di un vettore
		Vector norm(); //calcolo del versore di un vettore		
		Vector diffP(Vector*); //differenza tra il punto/vettore e un altro punto/vettore
		Vector addP(Vector*); //somma di vettori
		Vector multP(double); //moltiplicazione di un double per un vettore
		Vector divP(double); //divisione di un vettore per un double
		void Eulerodiscretization(Vector *); //discretizzazione dei punti del piano secondo il metodo di Eulero in avanti
		void discretization(Vector *, int dim_x, int dim_y, int dim_z); //discretizzazione dei punti del piano
		
};

#endif
