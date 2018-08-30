#ifndef LATTICEBOLTZMANN_H
#define LATTICEBOLTZMANN_H


class LB{
private:
	double w[Q];
	int cx[Q], cy[Q], cz[Q];
	double f[Nx1][Ny1][Nz1][Q], f_post[Nx1][Ny1][Nz1][Q];
	bool Esfrontera[Nx1][Ny1][Nz1];
	bool Inlet[Nx1][Ny1][Nz1];
	bool Outlet[Nx1][Ny1][Nz1];
public:
	LB(void);
	void Inicializacion(void);
	void Colision(void);
	void Streaming(void);
	void Macroscopicas(void);
	void BBOS(void);
	void Datos(void);


#endif

	double feq(double RHO, double U, double V, double W, int q);
	double Fi(double RHO, double U, double V, double W, int q);
	double Fx, Fy, Fz;
	double rho[Nx1][Ny1][Nz1];
	double ux[Nx1][Ny1][Nz1];
	double uy[Nx1][Ny1][Nz1];
};


