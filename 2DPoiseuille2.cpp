#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define Nx 256 // Numero de cuadro en la direccion x
#define Ny 256 // Numero de cuadro en la direccion y
#define Nx1 (Nx+1)
#define Ny1 (Ny+1)
#define L (Ny+1) // Ancho de la cavidad
#define Q 9 		// Numero de velocidades discretas
#define rho0 1.0  // Densidad Inicial
#define ux0 0.0 // Velocidad inicial en la componente x
#define uy0 0.0   // Velocidad inicial en la componente y 
#define g 0.000001 // verificar para conservar compresibilidad !!!
#define dt 1

int cx[Q]={0, 1, 0, -1, 0, 1, -1, -1, 1};
int cy[Q]={0, 0, 1, 0, -1, 1, 1, -1, -1};

bool Esfrontera[Ny1][Nx1];
bool EsfronteraPeriodica[Ny1][Nx1];
double Fx, Fy;
double f_post[Ny1][Nx1][Q]; // Arreglo de las funciones de distribucion luego de la colision
double f[Ny1][Nx1][Q];
double rho[Ny1][Nx1], ux[Ny1][Nx1], uy[Ny1][Nx1];
// Arreglo de la densidad, velocidad en x e y 
double tau; // Tiempo de relajacion en el modelo BGK
double w[Q]={4.0/9 ,1.0/9 ,1.0/9 ,1.0/9 ,1.0/9 ,1.0/36 ,1.0/36,
1.0/36,1.0/36}; // Pesos 
//int rc[Q]={0,3,4,1,2,7,8,5,6}; // index of reversed velocity
void Init_Eq(void); //Funcion de Initializacion
double feq(double RHO, double U, double V, int k);
// Funcionde distribucion de equilibrio 
void Coll_BGK(void); // BGK colision
void Streaming(void); // Streaming (Transmision)
void Den_Vel(void); // Variables macroscopicas
void BBOS(void);
double Err(void); // Funcion del error para parar el ciclo 
double u0[Ny1][Nx1],v0[Ny1][Nx1]; // definicion de matrices de condiciones iniciales
void Data_Output(void); // Funcion que escribe los datos
double Fi(double RHO, double U, double V, int k); // funcion para agregar el forzamiento 

//=========================================================
//=========================================================
int main()
{
	int k,M2,N2;
	double err;
	M2=Ny/2; N2=Nx/2;
	k=0;
	tau=0.9;// Tiempo de relajacion para BGK
	Init_Eq();
	while(k <=100)
	{
		k++;
		Coll_BGK(); //BGK colision
		Streaming(); // Streaming (Transmision)
		BBOS();
		Den_Vel(); // Variables macroscopicas de fluido 
		printf("rho=%e ux_center=%e uy_center=%e k=%d\n",rho[M2][N2],ux[M2][N2],uy[M2][N2], k); 	
	}
	Data_Output(); //Escribir los pasos cuando acabe la iteracion
}


// Funcion que inicializa las matrices
void Init_Eq()
{
	int j, i, k;
	for(j=0;j<=Ny;j++) for(i=0;i<=Nx;i++)
	{
		rho[j][i]=rho0;
		ux[j][i]=ux0;
		uy[j][i]=uy0;
		Esfrontera[j][i] = false; 
		EsfronteraPeriodica[0][i] = true;
		EsfronteraPeriodica[Ny][i] = true;
		Esfrontera[j][0] = true;
		Esfrontera[j][Nx] = true;
		for(k=0;k<Q;k++)
		f[j][i][k]=feq(rho[j][i],ux[j][i],uy[j][i],k);
	}

}


// Calculo de la distribucion de equilibrio 
double feq(double RHO, double U, double V, int k)
{
	double cu, U2;
	cu=cx[k]*U+cy[k]*V; // c k*u
	U2=U*U+V*V; // u*u; norma al cuadrado
	return w[k]*RHO*(1.0+3.0*cu+4.5*cu*cu-1.5*U2);
}


double Si(double RHO, double U, double V, int k)
{
	double Fx = 0, Fy = -RHO*g;
	double t1 = cx[k]*Fx + cy[k]*Fy;
	double t2 = cx[k]*cx[k]*U*Fx +cx[k]*cy[k]*(V*Fx+U*Fy)+cy[k]*cy[k]*V*Fy;
	double t3 = U*Fx + V*Fy;
	return (1-((dt)/(2*tau)))*w[k]*(3.0*t1+ 9.0*t2- 3.0*t3);
}


// Funcion que hace la colision BGK
void Coll_BGK()
{
	int j, i, k;
	double FEQ;
	for (j=0;j<=Ny1;j++) for(i=0;i<=Nx1;i++) for(k=0;k<Q;k++)
	{
		FEQ=feq(rho[j][i],ux[j][i],uy[j][i],k); // EDF
		f_post[j][i][k] = (1 - (dt/tau))*f[j][i][k] + (dt/tau)*FEQ +
		 dt*Si(rho[j][i],ux[j][i],uy[j][i],k) ;// Post-collision funciones de distribucion
	}
}


// Streaming (Transmision de la informacion)

void Streaming()
{
	int j, i, jd, id, k;
	for (j=0;j<Ny1;j++) for(i=0;i<Nx1;i++) for(k=0;k<Q;k++){
	jd=j-cy[k]; id=i-cx[k]; 
	if(jd>=0 && jd<=Ny && id>=0 && id<=Nx){
		f[j][i][k]=f_post[jd][id][k]; // streaming
	}
	}
}

//Bounce Back
void BBOS(){
	int i,j;
	for (j=0;j<=Ny;j++) for(i=0;i<=Nx;i++){ 
		
		if(Esfrontera[j][i] == true){
			f[j][i][1]=f_post[j][i][3];
			f[j][i][5]=f_post[j][i][7];
			f[j][i][8]=f_post[j][i][6];

			f[j][i][7]=f_post[j][i][5];
			f[j][i][3]=f_post[j][i][1];
			f[j][i][6]=f_post[j][i][8];
		}
		if(EsfronteraPeriodica[j][i] ==true){
			f[Ny][i][4]=f_post[0-int(cy[4]*dt)][i-int(cx[4]*dt)][4];
			f[Ny][i][7]=f_post[0-int(cy[7]*dt)][i-int(cx[7]*dt)][7];
			f[Ny][i][8]=f_post[0-int(cy[8]*dt)][i-int(cx[8]*dt)][8];
			f[0][i][2]=f_post[Ny-int(cy[2]*dt)][i-int(cx[2]*dt)][2];
			f[0][i][5]=f_post[Ny-int(cy[5]*dt)][i-int(cx[5]*dt)][5];
			f[0][i][6]=f_post[Ny-int(cy[6]*dt)][i-int(cx[6]*dt)][6];
		}

	}
}



//Calculo de las variables macroscopicas
void Den_Vel()
{
	int j, i;
	double Fx = 0.0;
	double Fy = -g;
	for(j=0;j<=Ny;j++) for(i=0;i<=Nx;i++)
	{

	rho[j][i]=f[j][i][0]+f[j][i][1]+f[j][i][2]+f[j][i][3]
	+f[j][i][4]+f[j][i][5]+f[j][i][6]+f[j][i][7]+
	f[j][i][8];

	ux[j][i]=(f[j][i][1]+f[j][i][5]+f[j][i][8]-f[j][i][3]-
	f[j][i][6]-f[j][i][7])/rho[j][i] + 0.5*dt*Fx;
	
	uy[j][i]=(f[j][i][5]+f[j][i][6]+f[j][i][2]-f[j][i][7]-
	f[j][i][8]-f[j][i][4])/rho[j][i] + 0.5*dt*Fy;
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
	for (i=0; i<=Nx; i++) fprintf(fp,"%e ",ux[j][i]);
	fprintf(fp,"\n");
	}
	fclose(fp);
	
	fp=fopen("vy.dat","w");
	for(j=0;j<=Ny;j++){
	for (i=0; i<=Nx; i++) fprintf(fp,"%e ",uy[j][i]);
	fprintf(fp,"\n");
	}
	fclose(fp);
	
	fp=fopen("rho.dat","w");
	for(j=0;j<=Ny;j++){
	for (i=0; i<=Nx; i++) fprintf(fp,"%e ",rho[j][i]);
	fprintf(fp,"\n");
	}
	fclose(fp);
}