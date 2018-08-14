#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define Nx 128 
#define Ny 128 
#define Nz 10 
#define Nx1 (Nx+1)
#define Ny1 (Ny+1)
#define Nz1 (Nz+1)
#define L (Ny+1)
#define Q 19 		
#define rho0 1.0  
#define ux0 0.0 
#define uy0 0.0   
#define uz0 0.0   
#define g 0.000001 
#define dt 1



int cx[Q] = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0};
int cy[Q] = {0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1};
int cz[Q] = {0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1};
//			 0  1  2  3. 4. 5. 6. 7. 8. 9. 10 11 12 13 14 15 16 17 18



bool EsFrontera[Nx1][Ny1][Nz1];
bool Inlet[Nx1][Ny1][Nz1];
bool Outlet[Nx1][Ny1][Nz1];
double Fx, Fy, Fz;
double f[Nx1][Ny1][Nz1][Q]; 
double f_post[Nx1][Ny1][Nz1][Q]; 
double rho[Nx1][Ny1][Nz1], ux[Nx1][Ny1][Nz1], uy[Nx1][Ny1][Nz1],uz[Nx1][Ny1][Nz1];
double tau; 
double w[Q] = {1.0/3 ,1.0/18,1.0/18,1.0/18,1.0/18,1.0/18,
			   1.0/18,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,
			   1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36};
			//   0      1       2        3      4       5      
  			//   6      7       8        9     10      11      
			//   12     13     14       15     16      17  	18 


void Init_Eq(void);
double feq(double RHO, double U, double V, double W,int q);
void Coll_BGK(void);
void Streaming(void); 
void Den_Vel(void); 
void BBOS(void);
double u0[Nx1][Ny1][Nz1],v0[Nx1][Ny1][Nz1],w0[Nx1][Ny1][Nz1];
void Data_Output(void);
double Fi(double RHO, double U, double V, double W,int q);



int main()
{
	int n,M2,N2,O2;

	M2=Nx/2; N2=Ny/2; O2 = Nz/2;
	n=0;
	tau=0.9;
	Init_Eq();
	while(n <=500)
	{
		n++;
		Coll_BGK(); 
		Streaming(); 
		//BBOS();
		Den_Vel(); 
		printf("rho=%e ux_c=%e uy_c=%e uz_c=%e k=%d\n",rho[M2][N2][O2],ux[M2][N2][O2],uy[M2][N2][O2],uz[M2][N2][O2], n); 	
	}
	Data_Output(); 
}



void Init_Eq()
{
	int j, i, k, q;
	for(i=0;i<=Nx;i++) for(j=0;j<=Ny;j++) for(k = 0; k<=Nz; k++)
	{
		rho[i][j][k] = rho0;
		ux[i][j][k] = ux0;
		uy[i][j][k] = uy0;
		uz[i][j][k] = uz0;
		EsFrontera[i][j][k] = false;
		EsFrontera[0][j][k] = true;
		Outlet[i][0][k] = true;
		EsFrontera[i][j][0] = true;
		EsFrontera[Nx][j][k] = true;
		Inlet[i][Ny][k] = true;
		EsFrontera[i][j][Nz] = true;
		for(q=0;q<Q;q++)
		f[i][j][k][q]=feq(rho[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],q);
	}
}
 
double feq(double RHO, double U, double V, double W,int q)
{
	double cu, U2;
	cu=cx[q]*U+cy[q]*V+cz[q]*W; 
	U2=U*U+V*V+W*W; 
	return w[q]*RHO*(1.0+3.0*cu+4.5*cu*cu-1.5*U2);
}
double Si(double RHO, double U, double V, double W,int q)
{
	double Fx = 0.0, Fy = -RHO*g, Fz = 0.0;
	double t1 = cx[q]*Fx + cy[q]*Fy + cz[q]*Fz;
	double t2 = cx[q]*cx[q]*U*Fx + cx[q]*cy[q]*(V*Fx + U*Fy) + 
				cx[q]*cz[q]*(W*Fx + U*Fz) + cy[q]*cy[q]*V*Fy +
				cy[q]*cz[q]*(V*Fz + W*Fy) + cz[q]*cz[q]*W*Fz;
	double t3 = U*Fx + V*Fy + W*Fz;
	return (1-((dt)/(2*tau)))*w[q]*(3.0*t1+ 9.0*t2- 3.0*t3);	



}
void Coll_BGK()
{
	int j, i, k, q;
	double FEQ;
	for (i=0;i<=Nx1;i++) for(j=0;j<=Ny1;j++) for(k=0;k<=Nz1;k++) for(q=0;q<Q;q++)
	{
		FEQ=feq(rho[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],q); 
		f_post[i][j][k][q] = (1 - (dt/tau))*f[i][j][k][q] + (dt/tau)*FEQ
		+dt*Si(rho[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],q);
	}
}



void Streaming()
{
	int j, i, k, jd, id, kd, q;
	for (i=0;i<Nx1;i++) for(j=0;j<Ny1;j++) for(k=0;k<Nz1;k++) for(q=0;q<Q;q++){
	jd=j-cy[q]; id=i-cx[q]; kd=k-cz[q]; 
	if(jd>=0 && jd<=Ny && id>=0 && id<=Nx && kd>=0 && kd<=Nz){
		f[i][j][k][q]=f_post[id][jd][kd][q]; // streaming
	}
	}
}

void BBOS(){
	int i,j,k,q;
	for (i=0;i<=Nx;i++) for(j=0;j<=Ny;j++) for(k=0;k<=Nz;k++){ 
		
		if(EsFrontera[i][j][k]==true){
			for (int q = 1; q < Q; q++)
				f[i][j][k][q]=f_post[i][j][k][int(q+pow(-1,q))];
		
		if(Inlet[i][j][k] == true)
			//plano arriba 
			f[i][j][k][3]=f_post[i-int(cx[3]*dt)][0-int(cy[3]*dt)][k-int(cz[3]*dt)][3];
			f[i][j][k][7]=f_post[i-int(cx[7]*dt)][0-int(cy[3]*dt)][k-int(cz[7]*dt)][7];
			f[i][j][k][11]=f_post[i-int(cx[11]*dt)][0-int(cy[3]*dt)][k-int(cz[11]*dt)][11];
			f[i][j][k][14]=f_post[i-int(cx[14]*dt)][0-int(cy[3]*dt)][k-int(cz[14]*dt)][14];
			f[i][j][k][17]=f_post[i-int(cx[17]*dt)][0-int(cy[3]*dt)][k-int(cz[17]*dt)][17];
		
		if (Outlet[i][j][k] == true)
			//plano abajo
			f[i][j][k][3]=f_post[i-int(cx[3]*dt)][Ny-int(cy[3]*dt)][k-int(cz[3]*dt)][3];
			f[i][j][k][7]=f_post[i-int(cx[7]*dt)][Ny-int(cy[7]*dt)][k-int(cz[7]*dt)][7];
			f[i][j][k][11]=f_post[i-int(cx[11]*dt)][Ny-int(cy[11]*dt)][k-int(cz[11]*dt)][11];
			f[i][j][k][14]=f_post[i-int(cx[14]*dt)][Ny-int(cy[14]*dt)][k-int(cz[14]*dt)][14];
			f[i][j][k][17]=f_post[i-int(cx[17]*dt)][Ny-int(cy[17]*dt)][k-int(cz[17]*dt)][17];
		
		}
	}
}

void Den_Vel()
{
	int j, i, k;
	double Ax = 0.0;
	double Ay = -g;
	double Az = 0.0;
	for(i = 0; i <= Nx; i++) 
		for(j = 0; j <= Ny; j++) 
			for(k = 0; k <= Nz; k++)
	{
		rho[i][j][k] = f[i][j][k][0]+f[i][j][k][1]+f[i][j][k][2]+f[i][j][k][3]+
		f[i][j][k][4]+f[i][j][k][5]+f[i][j][k][6]+f[i][j][k][7]+f[i][j][k][8]+
		f[i][j][k][9]+f[i][j][k][10]+f[i][j][k][11]+f[i][j][k][12]+f[i][j][k][13]+
		f[i][j][k][14]+f[i][j][k][15]+f[i][j][k][16]+f[i][j][k][17]+f[i][j][k][18];
		
		ux[i][j][k] = (f[i][j][k][1]+f[i][j][k][7]+f[i][j][k][9]+f[i][j][k][13]+f[i][j][k][15]
			-f[i][j][k][2]-f[i][j][k][8]-f[i][j][k][10]-f[i][j][k][14]-f[i][j][k][16])/rho[i][j][k];

		uy[i][j][k] = (f[i][j][k][3] + f[i][j][k][7] + f[i][j][k][11] + 
			f[i][j][k][14] + f[i][j][k][17] - f[i][j][k][4] - f[i][j][k][8]
			-f[i][j][k][12] - f[i][j][k][13] - f[i][j][k][18])/rho[i][j][k];

		uz[i][j][k] = (f[i][j][k][5] + f[i][j][k][9] + f[i][j][k][11]
			+f[i][j][k][16] + f[i][j][k][18] - f[i][j][k][6] - f[i][j][k][10]
				-f[i][j][k][12] - f[i][j][k][15] -f[i][j][k][17])/rho[i][j][k];
		 
	}
}







void Data_Output() 
{
	int i,j,k,z;
	z = 0;
	FILE *fp;
	fp=fopen("x.dat","w+");
	for(i=0;i<=Nx;i++) fprintf(fp,"%e \n", float(i)/L);
	fclose(fp);
	fp=fopen("y.dat","w+");
	for(j=0;j<=Ny;j++) fprintf(fp,"%e \n", float(j)/L);
	fclose(fp);
	fp=fopen("vx.dat","w");
	for(i=0;i<=Nx;i++) {
	for (j=0; j<=Ny; j++) fprintf(fp,"%e ",ux[i][j][z]);
	fprintf(fp,"\n");
	}
	fclose(fp);
	
	fp=fopen("vy.dat","w");
	for(i=0;i<=Nx;i++){
	for (j=0; j<=Ny; j++) fprintf(fp,"%e ",uy[i][j][z]);
	fprintf(fp,"\n");
	}
	fclose(fp);

	fp=fopen("rho.dat","w");
	for(i=0;i<=Nx;i++){
	for (j=0; j<=Ny; j++) fprintf(fp,"%e ",rho[j][i][z]);
	fprintf(fp,"\n");
	}
	fclose(fp);
}