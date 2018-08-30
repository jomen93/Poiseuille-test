#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "LatticeBoltzmann.h"


LB::LB(void){
	w[Q] = {1.0/3 ,1.0/18,1.0/18,1.0/18,1.0/18,1.0/18,
			   1.0/18,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,
			   1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36};
	cx[Q] = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0};
	cy[Q] = {0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1};
	cz[Q] = {0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1};
}


double LB::feq(double RHO, double U, double V, double W,int q){
	double cu, U2;
	cu=cx[q]*U+cy[q]*V+cz[q]*W; 
	U2=U*U+V*V+W*W; 
	return w[q]*RHO*(1.0+3.0*cu+4.5*cu*cu-1.5*U2);
}


double LB::Si(double RHO, double U, double V, double W,int q){
	double Fx = 0.0, Fy = -RHO*g, Fz = 0.0;
	double t1 = cx[q]*Fx + cy[q]*Fy + cz[q]*Fz;
	double t2 = cx[q]*cx[q]*U*Fx + cx[q]*cy[q]*(V*Fx + U*Fy) + 
				cx[q]*cz[q]*(W*Fx + U*Fz) + cy[q]*cy[q]*V*Fy +
				cy[q]*cz[q]*(V*Fz + W*Fy) + cz[q]*cz[q]*W*Fz;
	double t3 = U*Fx + V*Fy + W*Fz;
	return (1-((dt)/(2*tau)))*w[q]*(3.0*t1+ 9.0*t2- 3.0*t3);
}


void LB::Inicializacion(void){
	int i,j,k,q:
	for(i=0;i<=Nx;i++)for(j=0;j<=Ny;j++)for(k=0;k<=Nz;k++)
	{
		rho[i][j][k] = rho0;
		ux[i][j][k] = ux0;
		uy[i][j][k] = uy0;
		uz[i][j][k] = uz0;
		EsFrontera[i][j][k] = false;
		EsFrontera[0][j][k] = true;
		EsFrontera[Nx][j][k] = true;
		EsFrontera[i][j][0] = true;
		EsFrontera[i][j][Nz] = true;
		Inlet[i][Ny][k] = true;
		Outlet[i][0][k] = true;
		for(q=0;q<Q;q++)
			f[i][j][k][q]=feq(rho[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],q);
	}
}


void LB::Colision(void){
	int j, i, k, q;
	double FEQ;
	for (i=0;i<=Nx1;i++) for(j=0;j<=Ny1;j++) for(k=0;k<=Nz1;k++) for(q=0;q<Q;q++)
	{
		FEQ=feq(rho[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],q); 
		f_post[i][j][k][q] = (1 - (dt/tau))*f[i][j][k][q] + (dt/tau)*FEQ
		+dt*Si(rho[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],q);
	}
}


void LB::Streaming(void){
	int j, i, k, jd, id, kd, q;
	for (i=0;i<Nx1;i++) for(j=0;j<Ny1;j++) for(k=0;k<Nz1;k++) for(q=0;q<Q;q++){
	jd=j-cy[q]; id=i-cx[q]; kd=k-cz[q]; 
	if(jd>=0 && jd<=Ny && id>=0 && id<=Nx && kd>=0 && kd<=Nz){
		f[i][j][k][q]=f_post[id][jd][kd][q]; // streaming
		}
	}
}


void LB::BBOS(void){
	int i,j,k,q;
	for (i=0;i<=Nx;i++) for(j=0;j<=Ny;j++) for(k=0;k<=Nz;k++){ 
		
		if(Inlet[i][j][k] == true){
			//plano arriba 
			f[i][j][k][4]=f_post[i-cx[4]*dt][0][k-cz[4]*dt][4];
			f[i][j][k][8]=f_post[i-cx[8]*dt][0][k-cz[8]*dt][8];
			f[i][j][k][12]=f_post[i-cx[12]*dt][0][k-cz[12]*dt][12];
			f[i][j][k][13]=f_post[i-cx[13]*dt][0][k-cz[13]*dt][13];
			f[i][j][k][18]=f_post[i-cx[18]*dt][0][k-cz[18]*dt][18];
		}

		if (Outlet[i][j][k] == true){
			//plano abajo
			f[i][j][k][3]=f_post[i-cx[3]*dt][Ny][k-cz[3]*dt][3];
			f[i][j][k][7]=f_post[i-cx[7]*dt][Ny][k-cz[7]*dt][7];
			f[i][j][k][11]=f_post[i-cx[11]*dt][Ny][k-cz[11]*dt][11];
			f[i][j][k][14]=f_post[i-cx[14]*dt][Ny][k-cz[14]*dt][14];
			f[i][j][k][17]=f_post[i-cx[17]*dt][Ny][k-cz[17]*dt][17];
			}
		
		if(EsFrontera[i][j][k]==true){
			//for (int q = 1; q < Q; q++)
				//f[i][j][k][q]=f_post[i][j][k][int(q+pow(-1,q))];
				
				f[i][j][k][1]=f_post[i][j][k][2];
				f[i][j][k][2]=f_post[i][j][k][1];
				f[i][j][k][3]=f_post[i][j][k][4];
				f[i][j][k][4]=f_post[i][j][k][3];
				f[i][j][k][5]=f_post[i][j][k][6];
				f[i][j][k][6]=f_post[i][j][k][5];
				f[i][j][k][7]=f_post[i][j][k][8];
				f[i][j][k][8]=f_post[i][j][k][7];
				f[i][j][k][9]=f_post[i][j][k][10];
				f[i][j][k][10]=f_post[i][j][k][9];
				f[i][j][k][11]=f_post[i][j][k][12];
				f[i][j][k][12]=f_post[i][j][k][11];
				f[i][j][k][13]=f_post[i][j][k][14];
				f[i][j][k][14]=f_post[i][j][k][13];
				f[i][j][k][15]=f_post[i][j][k][16];
				f[i][j][k][16]=f_post[i][j][k][15];
				f[i][j][k][17]=f_post[i][j][k][18];
				f[i][j][k][18]=f_post[i][j][k][17];
			}
		}	
}


void LB::Datos(void){
	int i,j,k,z;
	z = 1;
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