/*________________________________________________________________________________________________________
 
programa C Simulacion desarrollo: difusion y red genetica
___________________________________________________________________________________________________________
*/

#include <cstdlib>
#include <pngwriter.h>
#include <iostream.h>
#include <string>
//using namespace std;
#include <stdio.h>

#include <time.h>

/* definicion de las estructuras de datos a emplear:
	- Genes, estructura primitiva de las células, serán de tres tipos
	- Posicion, será una terna de números que localice a la célula
	- Prop Física, actualmente no se implementan
	- Espacio difusivo es un array de estructuras espDifusivo
	- Espacio fisico, todavía no implementado.
	- Células, estructura que tiene anidados genes, posición y prop físicas.
*/

#define 	LARGO		400	
#define 	ALTO		1
#define 	ANCHO		400
#define 	FACTOR		1
#define 	MAX_CEL		1000
#define 	MIN_CEL		200
#define 	MAXGEN_REG	5
#define 	MAXGEN_SEN	2
#define 	MAXGEN_EST	2

/* Ahora debemos definir unos limites razonables para inicializar aleatoriamente los valores de los parámetros
del modelo. Proceso arduo. Se deben ajustar los rangos de las afinidades, umbrales, tasas de degradación y difusion
de todas las sustancias.
*/
#define AFF		8
#define UMB		12
#define DEGR		9
#define DIF		90
#define FIS_K		10
#define FIS_MU		8
#define SENAL		2
#define T_SIMULACION	2000
#define INTERVALO	0.1

int configura(void);

int inicializa(void);

int simula(void);

int difusion_engine(void);

int dibujar(void);





int celulas_TOT=0;


//int espDifusivo2D[LARGO][ANCHO][MAXGEN_SEN];
int espDifusivo2D_old[LARGO][ANCHO][MAXGEN_SEN];
int espDifusivo2D[LARGO][ANCHO][MAXGEN_SEN]; //Esta matriz contendrá los valores donde sí existe célula


float t_simulacion;	

int main(void)
{
	
	configura();
	inicializa();
	simula();
	//visualiza();
	return(0);
}

int configura(void)
{
	srand( (unsigned int)time( NULL ) );
	return(0);
}

int inicializa(void)
{
	int i,j,s;

	for (i=0; i<LARGO; i++)
	{	for (j=0; j<ANCHO; j++) 
		{
			//for (s=0;s<MAXGEN_SEN; s++)
				espDifusivo2D[i][j][0]= i*(rand() % (DIF));
				espDifusivo2D[i][j][1]= (ANCHO - i)*(rand() % (DIF));
		}//espDifusivo2D[i][j] =  i*DIF+5;

	} 
//fclose(archivo);
return(0);
}


//==============================================================================================//
//				ALGORITMO DE SIMULACION						//
//==============================================================================================//

int simula(void) //Esta función lleva a cabo un bucle llamando al simulador de física las veces que le digamos
{	//Si te viene mejor que la iteración sea en la otra función pues chachi. Lo quitas de aquí y listo.
	
	for(t_simulacion=0; t_simulacion<T_SIMULACION; t_simulacion = t_simulacion + INTERVALO)
	{
		//red_genetica_engine();
		difusion_engine();
		//fisica_engine();
	}
	return(0);
}




int difusion_engine(void)
{
	static int cont = 0;
	int i,j,s;
	cont++;

	//espDifusivo2D[15][15][0] = espDifusivo2D[15][15][0] + 1000;
	
	//espDifusivo2D[15][25][0] = espDifusivo2D[15][25][0] + 1000;
	
	espDifusivo2D[200][350][0] = espDifusivo2D[200][350][0] + 35000*sin(t_simulacion);

	espDifusivo2D[200][255][0] = espDifusivo2D[200][255][0] + 35000*sin(t_simulacion+1.62);

//	espDifusivo2D[200][200][1] = espDifusivo2D[200][200][0] + 25000;

//	espDifusivo2D[200+2][350+2][0] = espDifusivo2D[200][350][0] + 25000;

//	espDifusivo2D[200+2][255+2][0] = espDifusivo2D[200][255][0] + 25000;

//	espDifusivo2D[200-1][200-1][1] = espDifusivo2D[200][200][0] + 35000;
	

	espDifusivo2D[200][150][0] = espDifusivo2D[200][150][1] + 35000*sin(t_simulacion);
	
	espDifusivo2D[200][55][1] = espDifusivo2D[200][55][1] + 35000*sin(t_simulacion+1.62);
	
//	espDifusivo2D[15][35][1] = espDifusivo2D[15][35][1] + 25000;

//	espDifusivo2D[15][55][0] = espDifusivo2D[15][55][0] + 25000;

//	espDifusivo2D[15][75][1] = espDifusivo2D[15][75][1] + 25000;

	



	for (i=0; i<LARGO; i++)
	{
		
		for(j=0; j<ANCHO; j++)
		{
			for (s=0;s<MAXGEN_SEN; s++)
				espDifusivo2D_old[i][j][s]=espDifusivo2D[i][j][s];
		}
	}
	
	if ((cont % 8) == 1)
		dibujar();

	for (i=0; i<LARGO; i++)
		{
			for(j=0; j<ANCHO; j++)
			{ 
				for(s=0;s<MAXGEN_SEN; s++)
				{
				if ((i == 0) | (j == 0) | (i == (LARGO - 1)) | (j == (ANCHO - 4)))
					espDifusivo2D[i][j][s] = espDifusivo2D[i][j][s] - 4;
				else if ( (espDifusivo2D[i][j][s] =  espDifusivo2D_old[i][j][s] + (int) (0.125*(espDifusivo2D_old[i-1][j][s] + espDifusivo2D_old[i+1][j][s] + espDifusivo2D_old[i][j-1][s] + espDifusivo2D_old[i][j+1][s] - 4*espDifusivo2D_old[i][j][s]))) < 0)
						espDifusivo2D[i][j][s] = 0;
					
				}
			
			}
		}

		
	return(0);
}

int dibujar()
{
	static float contador = 0;
	int dec;
	int sing;
	int len, q, i, j,s;
	int rangocolor = (int) (65535 * 1.5 / (LARGO*DIF));
	//FILE *archivo1, *archivo2;
	//archivo1=fopen("datos_difusion_sen0", "a");
	//archivo2=fopen("datos_difusion_sen1", "a");
	char nombre[30]="difusion";
	char numero[8]="0";
	char *numero_tmp;
	char *numero_fin;
	contador++;
	numero_tmp = fcvt(contador,0,&dec,&sing);
	strcat(numero, numero_tmp);
	numero_tmp = numero;
	numero_fin = strchr(numero,'\0');
	len = (int) (numero_fin - numero_tmp) / sizeof(char);
	if (len>7)
		return (0);
	for (q=0; q<(7-len); q++)
		strcat(nombre,"0");
		
	strcat(nombre, numero);
	strcat(nombre, ".png");
//	printf("\t%s\n", nombre);
	pngwriter one((ANCHO*FACTOR),(LARGO*FACTOR),0,nombre);
	one.setcompressionlevel(9);
/*	
      s=0;
	fprintf(archivo1, "\n==================================================================\n");
	fprintf(archivo1, "\titeracion nº:  %f\t Señal:  %d\n", t_simulacion,s);
	fprintf(archivo1, "==================================================================\n");
	

      s=1;
	fprintf(archivo2, "\n==================================================================\n");
	fprintf(archivo2, "\titeracion nº:  %f\t Señal:  %d\n", t_simulacion,s);
	fprintf(archivo2, "==================================================================\n");
	
*/	
	
	for (i=0; i<LARGO; i++)
	{
		
		for(j=0; j<ANCHO; j++)
		{
			//for(s=0; s<MAXGEN_SEN;s++)
			//{
			one.filledsquare((i * FACTOR),(j * FACTOR),(i * FACTOR + FACTOR),(j * FACTOR + FACTOR),(int) (rangocolor*espDifusivo2D_old[i][j][0]),0,(int) (rangocolor*espDifusivo2D_old[i][j][1]));
			//fprintf(archivo1, "%3d ", espDifusivo2D[i][j][0]);
			//fprintf(archivo2, "%3d ", espDifusivo2D[i][j][1]);
			//if (j == (ANCHO - 1)) { fprintf(archivo1,"\n"); fprintf(archivo2, "\n"); }
			//}
		}
	}
	//fclose(archivo1);
	//fclose(archivo2);
	one.close();
	return (0);
}
















