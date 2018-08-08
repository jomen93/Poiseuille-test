#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define Nx 256 // Numero de cuadro en la direccion x
#define Ny 256 // Numero de cuadro en la direccion y
#define Nz 10 // Numero de cuadro en la direccion y
#define Nx1 (Nx+1)
#define Ny1 (Ny+1)
#define Nz1 (Nz+1)
#define L (Ny+1)
#define Q 19 		// Numero de velocidades discretas
#define rho0 1.0  // Densidad Inicial
#define ux0 0.0 // Velocidad inicial en la componente x
#define uy0 0.0   // Velocidad inicial en la componente y 
#define uz0 0.0   // Velocidad inicial en la componente y 
#define g 0.000001 // verificar para conservar compresibilidad !!!
#define dt 1

//int cx[Q]={0, 1, 0, -1, 0, 1, -1, -1, 1};
//int cy[Q]={0, 0, 1, 0, -1, 1, 1, -1, -1};

int cx[Q] = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0};
int cy[Q] = {0, 0, 0, 1,-1, 0, 0, 1,-1, 1, 0, 1,-1,-1, 1, 0, 0, 1,-1};
int cz[Q] = {0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 1, 0,-1, 1,-1, 1};
//			 0  1  2  3. 4. 5. 6. 7. 8. 9. 10 11 12 13 14 15 16 17 18



bool EsFrontera[Ny1][Nx1][Nz1];
bool EsFronteraPeriodica[Ny1][Nx1][Nz1];
double Fx, Fy, Fz;
double f[Ny1][Nx1][Nz1][Q]; //Arreglo de las funciones de distribucion
double f_post[Ny1][Nx1][Nz1][Q]; // Arreglo de las funciones de distribucion luego de la colision
double rho[Ny1][Nx1][Nz1], ux[Ny1][Nx1][Nz1], uy[Ny1][Nx1][Nz1],uz[Ny1][Nx1][Nz1];
// Arreglo de la densidad, velocidad en x e y 
double tau; // Tiempo de relajacion en el modelo BGK
double w[Q] = {1.0/3 ,1.0/18,1.0/18,1.0/18,1.0/18,1.0/18,
			   1.0/18,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,
			   1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36};
			//   0      1       2        3      4       5      
  			//   6      7       8        9     10      11      
			//   12     13     14       15     16      17  	18 

//double w[Q]={4.0/9 ,1.0/9 ,1.0/9 ,1.0/9 ,1.0/9 ,1.0/36 ,1.0/36,
//1.0/36,1.0/36}; // Pesos 
//int rc[Q]={0,3,4,1,2,7,8,5,6}; // index of reversed velocity
void Init_Eq(void); //Funcion de Initializacion
double feq(double RHO, double U, double V, double W,int q);
// Funcionde distribucion de equilibrio 
void Coll_BGK(void); // BGK colision
void Streaming(void); // Streaming (Transmision)
void Den_Vel(void); // Variables macroscopicas
void BBOS(void);
double u0[Ny1][Nx1][Nz1],v0[Ny1][Nx1][Nz1],w0[Ny1][Nx1][Nz1]; // definicion de matrices de condiciones iniciales
void Data_Output(void); // Funcion que escribe los datos
double Fi(double RHO, double U, double V, double W,int q); // funcion para agregar el forzamiento 
//=========================================================
//=========================================================
int main()
{
	int n,M2,N2,O2;

	M2=Ny/2; N2=Nx/2; O2 = Nz/2;
	n=0;
	tau=50;// Tiempo de relajacion para BGK
	Init_Eq();
	while(n <=100)
	{
		n++;
		Coll_BGK(); //BGK colision
		Streaming(); // Streaming (Transmision)
		BBOS();
		Den_Vel(); // Variables macroscopicas de fluido 
		printf("rho=%e ux_c=%e uy_c=%e uz_c=%e k=%d\n",rho[M2][N2][O2],ux[M2][N2][O2],uy[M2][N2][O2],uz[M2][N2][O2], n); 	
	}
	Data_Output(); //Escribir los pasos cuando acabe la iteracion
}


// Funcion que inicializa las matrices
void Init_Eq()
{
	int j, i, k, q;
	for(j=0;j<=Ny;j++) for(i=0;i<=Nx;i++) for(k = 0; k<=Nz; k++)
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
		for(q=0;q<Q;q++)
		f[j][i][k][q]=feq(rho[j][i][k],ux[j][i][k],uy[j][i][k],uz[j][i][k],q);
	}

}


// Calculo de la distribucion de equilibrio 
double feq(double RHO, double U, double V, double W,int q)
{
	double cu, U2;
	cu=cx[q]*U+cy[q]*V+cz[q]*W; // c k*u
	U2=U*U+V*V+W*W; // u*u; norma al cuadrado
	return w[q]*RHO*(1.0+3.0*cu+4.5*cu*cu-1.5*U2);
}

double Si(double RHO, double U, double V, double W,int q)
{
	double Fx = 0.0, Fy = -RHO*g, Fz = 0.0;
	double t1 = cx[q]*Fx + cy[q]*Fy + cz[q]*Fz;
	double t2 = cx[q]*cx[q]*U*Fx + cx[q]*cy[q]*(V*Fx + U*Fy) + 
				cx[q]*cz[q]*(W*Fz + U*Fz) + cy[q]*cy[q]*V*Fy +
				cy[q]*cz[q]*(V*Fz + W*Fy) + cz[q]*cz[q]*W*Fz;
	double t3 = U*Fx + V*Fy + W*Fz;
	return (1-((dt)/(2*tau)))*w[q]*(3.0*t1+ 9.0*t2- 3.0*t3);	

	//double Fx = 0, Fy = -RHO*g;
	//double t1 = cx[q]*Fx + cy[q]*Fy;
	//double t2 = cx[q]*cx[q]*U*Fx +cx[q]*cy[q]*(V*Fx+U*Fy)+cy[q]*cy[q]*V*Fy;
	//double t3 = U*Fx + V*Fy;
	//return (1-((dt)/(2*tau)))*w[q]*(3.0*t1+ 9.0*t2- 3.0*t3);
}


// Funcion que hace la colision BGK
void Coll_BGK()
{
	int j, i, k, q;
	double FEQ;
	for (j=0;j<=Ny1;j++) for(i=0;i<=Nx1;i++) for(k=0;k<=Nz1;k++) for(q=0;q<Q;q++)
	{
		FEQ=feq(rho[j][i][k],ux[j][i][k],uy[j][i][k],uz[j][i][k],q); // EDF
		f_post[j][i][k][q] = (1 - (dt/tau))*f[j][i][k][q] + (dt/tau)*FEQ +
		 dt*Si(rho[j][i][k],ux[j][i][k],uy[j][i][k],uz[j][i][k],q) ;// Post-collision funciones de distribucion
	}
}


// Streaming (Transmision de la informacion)

void Streaming()
{
	int j, i, k, jd, id, kd, q;
	for (j=0;j<Ny1;j++) for(i=0;i<Nx1;i++) for(k=0;k<Nz1;k++) for(q=0;q<Q;q++){
	jd=j-cy[q]; id=i-cx[q]; kd=k-cz[q]; 
	if(jd>=0 && jd<=Ny && id>=0 && id<=Nx){
		f[j][i][k][q]=f_post[jd][id][kd][q]; // streaming
	}
	}
}


//Bounce Back
void BBOS(){
	int i,j,k;
	for (j=0;j<=Ny1;j++) for(i=0;i<=Nx1;i++) for(k=0;k<=Nz1;k++){ 
		if(EsFrontera[j][i][k]==true){
			//plano derecho
			f[j][i][k][1]=f_post[j][i][k][2];
			f[j][i][k][7]=f_post[j][i][k][8];
			f[j][i][k][9]=f_post[j][i][k][10];
			f[j][i][k][13]=f_post[j][i][k][14];
			f[j][i][k][15]=f_post[j][i][k][16];
			//plano izquierdo
			f[j][i][k][2]=f_post[j][i][k][1];
			f[j][i][k][8]=f_post[j][i][k][7];
			f[j][i][k][10]=f_post[j][i][k][9];
			f[j][i][k][14]=f_post[j][i][k][13];
			f[j][i][k][16]=f_post[j][i][k][15];
			//plano atras
			f[j][i][k][6]=f_post[j][i][k][5];
			f[j][i][k][10]=f_post[j][i][k][9];
			f[j][i][k][12]=f_post[j][i][k][11];
			f[j][i][k][15]=f_post[j][i][k][16];
			f[j][i][k][17]=f_post[j][i][k][18];
			//plano adelante
			f[j][i][k][5]=f_post[j][i][k][6];
			f[j][i][k][9]=f_post[j][i][k][10];
			f[j][i][k][11]=f_post[j][i][k][12];
			f[j][i][k][16]=f_post[j][i][k][15];
			f[j][i][k][18]=f_post[j][i][k][17];

		if(EsFrontera[j][i][k] == true)
			//plano arriba

			f[0][i][k][3]=f_post[Ny-int(cy[3]*dt)][i-int(cx[3]*dt)][k-int(cz[3]*dt)][3];
			f[0][i][k][7]=f_post[Ny-int(cy[7]*dt)][i-int(cx[7]*dt)][k-int(cz[7]*dt)][7];
			f[0][i][k][11]=f_post[Ny-int(cy[11]*dt)][i-int(cx[11]*dt)][k-int(cz[11]*dt)][11];
			f[0][i][k][14]=f_post[Ny-int(cy[14]*dt)][i-int(cx[14]*dt)][k-int(cz[14]*dt)][14];
			f[0][i][k][17]=f_post[Ny-int(cy[17]*dt)][i-int(cx[17]*dt)][k-int(cz[17]*dt)][17];

			//f[j][i][k][3]=f_post[j][i][k][4];
			//f[j][i][k][7]=f_post[j][i][k][8];
			//f[j][i][k][11]=f_post[j][i][k][12];
			//f[j][i][k][14]=f_post[j][i][k][13];
			//f[j][i][k][17]=f_post[j][i][k][18];
			//plano abajo

			f[Ny][i][k][4]=f_post[0-int(cy[4]*dt)][i-int(cx[4]*dt)][k-int(cz[4]*dt)][4];
			f[Ny][i][k][8]=f_post[0-int(cy[8]*dt)][i-int(cx[8]*dt)][k-int(cz[8]*dt)][8];
			f[Ny][i][k][12]=f_post[0-int(cy[12]*dt)][i-int(cx[12]*dt)][k-int(cz[12]*dt)][12];
			f[Ny][i][k][13]=f_post[0-int(cy[13]*dt)][i-int(cx[13]*dt)][k-int(cz[13]*dt)][13];
			f[Ny][i][k][18]=f_post[0-int(cy[18]*dt)][i-int(cx[18]*dt)][k-int(cz[18]*dt)][18];
			
			//f[j][i][k][4]=f_post[j][i][k][3];
			//f[j][i][k][8]=f_post[j][i][k][7];
			//f[j][i][k][12]=f_post[j][i][k][11];
			//f[j][i][k][13]=f_post[j][i][k][14];
			//f[j][i][k][18]=f_post[j][i][k][17];
		}
	}
}


//Calculo de las variables macroscopicas
void Den_Vel()
{
	int j, i, k;
	double Ax = 0.0;
	double Ay = -g;
	double Az = 0.0;
	for(j = 0; j < Ny1; j++) 
		for(i = 0; i < Nx1; i++) 
			for(k = 0; k < Nz1; k++)
	{
		rho[j][i][k] = f[j][i][k][0]+f[j][i][k][1]+f[j][i][k][2]+f[j][i][k][3]+
		f[j][i][k][4]+f[j][i][k][5]+f[j][i][k][6]+f[j][i][k][7]+f[j][i][k][8]+
		f[j][i][k][9]+f[j][i][k][10]+f[j][i][k][11]+f[j][i][k][12]+f[j][i][k][13]+
		f[j][i][k][14]+f[j][i][k][15]+f[j][i][k][16]+f[j][i][k][17]+f[j][i][k][18];
		
		ux[j][i][k] = (f[j][i][k][1]-f[j][i][k][2]+f[j][i][k][7]-f[j][i][k][8]-
		f[j][i][k][9]-f[j][i][k][10]+f[j][i][k][13]-f[j][i][k][14]+f[j][i][k][15]-
		f[j][i][k][16])/rho[j][i][k] + 0.5*dt*Ax;

		uy[j][i][k] = (f[j][i][k][3]-f[j][i][k][4]+f[j][i][k][7]-f[j][i][k][8]+
			+f[j][i][k][11]-f[j][i][k][12]-f[j][i][k][13]+f[j][i][k][14]
			+f[j][i][k][17]-f[j][i][k][18])/rho[j][i][k] + 0.5*dt*Ay;
		
		uz[j][i][k] = (f[j][i][k][5]-f[j][i][k][6]+f[j][i][k][9]-f[j][i][k][10]+
			f[j][i][k][11]-f[j][i][k][12]-f[j][i][k][15]+f[j][i][k][16]
			-f[j][i][k][17]+f[j][i][k][18])/rho[j][i][k] + 0.5*dt*Az;
	}
}



void Data_Output() // Datos de salida
{
	int i,j;
	FILE *fp;
	fp=fopen("x.dat","w+");
	for(i=0;i<=Nx;i++) fprintf(fp,"%e \n", float(i)/L);
	fclose(fp);
	fp=fopen("y.dat","w+");
	for(j=0;j<=Ny;j++) fprintf(fp,"%e \n", float(j)/L);
	fclose(fp);
	fp=fopen("vx.dat","w");
	for(j=0;j<=Ny;j++) {
	for (i=0; i<=Nx; i++) fprintf(fp,"%e ",ux[j][i][5]);
	fprintf(fp,"\n");
	}
	fclose(fp);
	
	fp=fopen("vy.dat","w");
	for(j=0;j<=Ny;j++){
	for (i=0; i<=Nx; i++) fprintf(fp,"%e ",uy[j][i][5]);
	fprintf(fp,"\n");
	}
	fclose(fp);
	
	fp=fopen("rho.dat","w");
	for(j=0;j<=Ny;j++){
	for (i=0; i<=Nx; i++) fprintf(fp,"%e ",rho[j][i][5]);
	fprintf(fp,"\n");
	}
	fclose(fp);
}