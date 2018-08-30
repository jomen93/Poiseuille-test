#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LatticeBoltzmann.h"

const int  Nx = 64;
const int Ny = 64;
const int  Nz = 16;
const int  Nx1 = Nx+1;
const int Ny1 = Ny+1;
const int Nz1 = Nz+1;
const int  L = (Ny1+1);
const int  Q = 19;
const double rho0 = 1.0;
const double ux0 = 0.0;
const double uy0 = 0.0;
const double uz0 = 0.0;
const double g = 0.000001;
const double dt = 1;

double Fx,Fy,Fz;
double tau = 0.9;
int n = 0;
int M = Nx/2;
int N = Ny/2;
int O = Nz/2;

int main(){
	LB Fluido;
	Fluido.Inicializacion();
	Fluido.Imprimase();

}

