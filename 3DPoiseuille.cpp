#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
using namespace std;

#define Nx 128
#define Ny 128
#define Nz 1
#define Nx1 Nx+1
#define Ny1 Ny+1 
#define Nz1 Nz+1
#define Q 19
#define rho0 1.0
#define ux0 0.0
#define uy0 0.0
#define uz0 0.0
#define g 0.000001
#define dt 1

int cx[Q] = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0};
int cy[Q] = {0, 0, 0, 1,-1, 0, 0, 1,-1, 1, 0, 1,-1,-1, 1, 0, 0, 1,-1};
int cz[Q] = {0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 1, 0,-1, 1,-1, 1};
//			 0  1  2  3. 4. 5. 6. 7. 8. 9. 10 11 12 13 14 15 16 17 18

bool EsFrontera[Ny1][Nx1][Nz1];
double Fx, Fy, Fz;
double f[Ny1][Nx1][Nz1][Q];
double f_post[Ny1][Nx1][Nz1][Q];
double rho[Ny1][Nx1][Nz1], ux[Ny1][Nx1][Nz1]; 
double uy[Ny1][Nx1][Nz1], uz[Ny1][Nx1][Nz1] ;
double tau;
double w[Q] = {1.0/3 ,1.0/18,1.0/18,1.0/18,1.0/18,1.0/18,
			   1.0/18,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,
			   1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36};
  			
			//   0      1       2        3      4       5      
  			//   6      7       8        9     10      11      
			//   12     13     14       15     16      17  	18 


void Init_eq(void);
double feq(double RHO, double U, double V, double W, int q);
void Coll_BGK(void);
void Streaming(void);
void Den_Vel(void);
void Bounce_back(void);
double u0[Ny1][Nx1][Nz1],v0[Ny1][Nx1][Nz1],w0[Ny1][Nx1][Nz1];
void Data_Output(void);
double Si(double RHO, double U, double V, double W, int q);

//double zrho[Nz1];


int main()
{
	int k = 0;
	int M2, N2, O2;
	tau = 0.5;
	M2=Ny/2; N2=Nx/2; O2=Nz/2;
	Init_eq();
	while(k <= 100)
	{
		k++;
		Coll_BGK();
		Streaming();
		//Bounce_back();
		Den_Vel();
		printf("rho=%e ux_c=%e uy_c=%e uz_c=%e k=%d\n",rho[M2][N2][0],ux[M2][N2][0],uy[M2][N2][0],uz[M2][N2][0], k);
	}
	Data_Output();
}



void Init_eq()
{
	int j , i, k, q;
	for(j = 0; j <= Ny; j++) 
		for(i = 0; i <= Nx; i++) 
			for(k = 0; k <= Nz; k++)
	{
		rho[j][i][k] = rho0;
		ux[j][i][k] = ux0;
		uy[j][i][k] = uy0;
		uz[j][i][k] = uz0;
		EsFrontera[j][i][k] = false;
		EsFrontera[0][i][k] = true;
		EsFrontera[j][0][k] = true;
		EsFrontera[j][i][0] = true;
		EsFrontera[Ny][i][k] = true;
		EsFrontera[j][Nx][k] = true;
		EsFrontera[j][i][Nz] = true;
		for (q = 0; q < Q; q++)
			f[j][i][k][q] = feq(rho[j][i][k],ux[j][i][k],uy[j][i][k],uz[j][i][k],q);
	}


}


double feq(double RHO, double U, double V, double W,int q)
{
	double cu, U2;
	cu = cx[q]*U + cy[q]*V + cz[q]*W; // c k*u
	U2= U*U + V*V + W*W; // u*u; norma al cuadrado
	return w[q]*RHO*(1.0+3.0*cu+4.5*cu*cu-1.5*U2);
}


double Si(double RHO, double U, double V, double W, int q)
{
	double Fx = 0.0, Fy = -RHO*g, Fz = 0.0;
	double t1 = cx[q]*Fx + cy[q]*Fy + cz[q]*Fz;
	double t2 = cx[q]*cx[q]*U*Fx + cx[q]*cy[q]*(V*Fx + U*Fy) + 
				cx[q]*cz[q]*(W*Fz + U*Fz) + cy[q]*cy[q]*V*Fy +
				cy[q]*cz[q]*(V*Fz + W*Fy) + cz[q]*cz[q]*W*Fz;
	double t3 = U*Fx + V*Fy + W*Fz;
	return (1-((dt)/(2*tau)))*w[q]*(3.0*t1+ 9.0*t2- 3.0*t3);
}


void Coll_BGK()
{
	int i, j, k, q;
	double FEQ;
	for(j = 0; j <= Ny1; j++) 
		for(i = 0; i <= Nx1; i++) 
			for(k = 0; k <= Nz1; k++) 
				for (q = 0; q < Q; q++)
	{
		FEQ = feq(rho[j][i][k],ux[j][i][k],uy[j][i][k],uz[j][i][k],q);
		f_post[j][i][k][q] = (1 - (dt/tau))*f[j][i][k][q];// + (dt/tau)*FEQ;
		 		//+dt*Si(rho[j][i][k],ux[j][i][k],uy[j][i][k],uz[j][i][k],q);
	}
}


void Streaming()
{
	int j, i, k, jd, id, kd, q;
	for(j = 0; j <= Ny1; j++) 
		for(i = 0; i <= Nx1; i++) 
			for(k = 0; k <= Nz1; k++) 
				for(q = 0; q < Q; q++)
	{
		jd = j-cy[q]; id = i-cx[q]; kd = k-cz[q];
		if(jd >= 0 && jd <= Ny && id >= 0 && id <= Nx && kd >= 0 && kd <= Nz)
			f[j][i][k][q] = f_post[jd][id][kd][q]; 
	}
}









void Den_Vel()
{
	int j, i, k;
	double Ax = 0.0;
	double Ay = -g;
	double Az = 0.0;
	for(j = 0; j <= Ny1; j++) 
		for(i = 0; i <= Nx1; i++) 
			for(k = 0; k <= Nz1; k++)
	{
		rho[j][i][k] = f[j][i][k][0]+f[j][i][k][1]+f[j][i][k][2]+f[j][i][k][3]+
		f[j][i][k][4]+f[j][i][k][5]+f[j][i][k][6]+f[j][i][k][7]+f[j][i][k][8]+
		f[j][i][k][9]+f[j][i][k][10]+f[j][i][k][11]+f[j][i][k][12]+f[j][i][k][13]+
		f[j][i][k][14]+f[j][i][k][15]+f[j][i][k][16]+f[j][i][k][17]+f[j][i][k][18];
		
		ux[j][i][k] = (f[j][i][k][1]-f[j][i][k][2]+f[j][i][k][7]-f[j][i][k][8] +
			f[j][i][k][9]-f[j][i][k][10]+f[j][i][k][13]-f[j][i][k][14]
			+f[j][i][k][15]-f[j][i][k][16])/rho[j][i][k] + 0.5*dt*Ax;
		
		uy[j][i][k] = (f[j][i][k][3]-f[j][i][k][4]+f[j][i][k][7]-f[j][i][k][8]+
			f[j][i][k][9]+f[j][i][k][11]-f[j][i][k][12]-f[j][i][k][13]+f[j][i][k][14]
			+f[j][i][k][17]-f[j][i][k][18])/rho[j][i][k] + 0.5*dt*Ay;
		
		uz[j][i][k] = (f[j][i][k][5]-f[j][i][k][6]+f[j][i][k][9]-f[j][i][k][10]+
			f[j][i][k][11]-f[j][i][k][12]+f[j][i][k][13]-f[j][i][k][15]+f[j][i][k][16]
			-f[j][i][k][17]+f[j][i][k][18])/rho[j][i][k] + 0.5*dt*Az;
	}
}




void Data_Output() // Datos de salida
{
	int i,j,k;
	FILE *fp;
	fp=fopen("x.dat","w+");
	for(i=0;i<=Nx;i++) fprintf(fp,"%e \n", float(i)/Nx);
	fclose(fp);
	fp=fopen("y.dat","w+");
	for(j=0;j<=Ny;j++) fprintf(fp,"%e \n", float(j)/Ny);
	fclose(fp);
	fp=fopen("z.dat","w+");
	for(k=0;j<=Nz;k++) fprintf(fp,"%e \n", float(k)/Nz);
	fclose(fp);

	fp=fopen("vx.dat","w");
	for(j=0;j<=Ny;j++) {
	for (i=0; i<=Nx; i++) fprintf(fp,"%e ",ux[j][i][0]);
	fprintf(fp,"\n");
	}
	fclose(fp);
	
	fp=fopen("vy.dat","w");
	for(j=0;j<=Ny;j++){
	for (i=0; i<=Nx; i++) fprintf(fp,"%e ",uy[j][i][0]);
	fprintf(fp,"\n");
	}
	fclose(fp);

	fp=fopen("vz.dat","w");
	for(j=0;j<=Ny;j++){
	for (i=0; i<=Nx; i++) fprintf(fp,"%e ",uz[j][i][0]);
	fprintf(fp,"\n");
	}
	fclose(fp);
	
	fp=fopen("denrho.dat","w");
	for(j=0;j<=Ny;j++){
	for (i=0; i<=Nx; i++) fprintf(fp,"%e ",rho[j][i][0]);
	fprintf(fp,"\n");
	}
	fclose(fp);

}






















