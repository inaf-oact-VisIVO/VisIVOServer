/*
Copyright (c) 2012, MB
All rights reserved.
*/

#ifndef LINE_H
#define LINE_H

#include "vsVector.h"

class Line{

	public:
	Vector P0;
	Vector P1;
	Vector u; //vettore direzione
	Vector u_n; //versore	

	Line(Vector p1, Vector p2);
	Line();
	~Line();
	//float dist3D_Line_to_Line(Line l);
	Vector poca(Line *);
	double scatteringAngle(Line *l);//restituisce l'angolo di scattering
	


};

#endif
