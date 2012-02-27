///========================================================================================================
///________________________________________________________________________________________________________
/// 
///	Simutron2006
///	programa escrito en C 
///	Simulacion desarrollo: difusion y red genetica
///	Autores: Alvaro, Alberto y Jorge.
///________________________________________________________________________________________________________
///
///========================================================================================================


#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <pngwriter.h>

#include "conf.h"
#include "data_struct.h"
#include "funcionescpp.h"




/* definicion de las estructuras de datos a emplear:
	- Genes, estructura primitiva de las células, serán de tres tipos
	- Posicion, será una terna de números que localice a la célula
	- Prop Física, actualmente no se implementan
	- Espacio difusivo es un array de estructuras espDifusivo
	- Espacio fisico, todavía no implementado.
	- Células, estructura que tiene anidados genes, posición y prop físicas.
*/




main()
{
	
	configura();
	inicializa();
	simula();
	visualiza();
	return(0);
}


//==============================================================================================//
//	FUNCION QUE CONFIGURA EL MODO EN QUE TRABAJARÁ EL SIMULADOR. POR AHORA, NO PERMITE 
//      OPCIONES, AUNQUE SE PUEDEN IMPLEMENTAR MUCHOS MÉTODOS DE PARAMETRIZACION Y PERSONALIZACION//
//==============================================================================================//


int configura(void)
{
	srand( (unsigned int)time( NULL ) );
	return(0);
}

//==============================================================================================//
//	FUNCION QUE INICIALIZA LAS VARIABLES DEL PROGRAMA Y GUARDA ESOS DATOS EN conf_inicial	//
//==============================================================================================//


int inicializa(void)
{
	int i,j,k=0,r,s,cont;
	for (i=0; i<LARGO; i++)
	{	for (j=0; j<ANCHO; j++, k++) 
		{ //	Llamada a crear_celula recursiva hasta completar la poblacion inicial
			if ((r = (rand() % 10)) > 5)
			{ 
				espDifusivo2D[i][j][MAXGEN_SEN] = 0;
				k--;
				continue;
			}

			if ( (celulas[k] = crear_celula(i,j,0)) == NULL ) // Si el proceso de division falla...
				{ k--; espDifusivo2D[i][j][MAXGEN_SEN] = 0; }
			else	
				espDifusivo2D[i][j][MAXGEN_SEN] = k;
	
	 	for (s=0; s<MAXGEN_SEN; s++)  
			espDifusivo2D[i][j][s] = 0.5*(i+1)*(rand() % DIF);  //Inicializacion del espacio difusivo
		}
	}

	print_init_conf();
	return(0);
}




//==============================================================================================//
//	FUNCION MAESTRA COORDINADORA DEL ALGORITMO DE SIMULACION Y DEL TEMPO DE LA		//
//					 SIMULACION						//
//==============================================================================================//

int simula(void)
{
	
	for(t_simulacion=0; t_simulacion<T_SIMULACION; t_simulacion = t_simulacion + INTERVALO)
	{
		
		difusion_engine();
		fisica_engine();
		red_genetica_engine();
		vidaymuerte();
		if (celulas_TOT<MIN_CEL)
			break;
	}
	return(0);
}





//==============================================================================================//
//	FUNCIONES PARA VISUALIZAR LOS DATOS DE DIVERSAS FORMAS					//
//==============================================================================================//





int visualiza(void) 
{
return(0);

}


//==============================================================================================//
//	FUNCIONES QUE CONTROLAN LOS DIVERSOS ALGORITMOS DE SIMULACION				//
//												//
//==============================================================================================//



//==============================================================================================//
//	FUNCIONES QUE CONTROLAN EL MODELO FISICO DE MUELLES PARA LA DINAMICA DEL TEJIDO 	//
//==============================================================================================//


int fisica_engine(void)
{
	return (0);
}



//==============================================================================================//
//	FUNCION GENERADORA DEL GRADIENTE DIFUSIVO DE LOS DISTINTOS MORFOGENES			//
//		GUARDA SUS RESULTADOS EN RESPECTIVOS ARCHIVOS datos_difusionX			//
//==============================================================================================//
		
int difusion_engine(void)
{
	int i,j,s;

	espDifusivo2D[10][12][0] += 10; //Fuente de morfogen en [10,12]
	

	for(s=0; s<MAXGEN_SEN; s++)
	{
		print_espDifusivo(s);
		for (i=0; i<LARGO; i++)
		{
			for(j=0; j<ANCHO; j++)
			{ 
				if (espDifusivo2D_old[i][j][s] < 1)
					espDifusivo2D[i][j][s] = 0;
				else
				{
					if (i == 0 | j == 0 | i == (LARGO - 1) | j == (ANCHO - 1))
						espDifusivo2D[i][j][s] -= 9;
					else if ( (espDifusivo2D[i][j][s] += (int) 0.15*(espDifusivo2D_old[i-1][j][s] + espDifusivo2D_old[i+1][j][s] + espDifusivo2D_old[i][j-1][s] + espDifusivo2D_old[i][j+1][s] - 4*espDifusivo2D_old[i][j][s])) < 0)
espDifusivo2D[i][j][s] = 0;
					
				}
			}
		}
	}
	return(0);
}


//==============================================================================================//
//	FUNCIONES QUE CONTROLAN EL MODELO SIMPLE DE RED GENETICA REGULATORIA			//
//==============================================================================================//

int red_genetica_engine(void)
{
	int i;

	dibujar_red_genetica();

	for (i=0; i<celulas_TOT; i++)
	{
		if (genes_estructurales(i) == APOPTOSIS)
			continue;
		genes_reguladores(i);
		genes_senales(i);
		
	}
	return(0);
}

int vidaymuerte(void)
{
	int coordenada[3],a,i;
	for (i=0; i<celulas_TOT; i++)
	{
		if ((*celulas[i]).vida_muerte==APOPTOSIS)
		{
			muerte(i);
			i--;
		}
		
		else if ((*celulas[i]).vida_muerte==DIVISION)
		{
			(*celulas[i]).vida_muerte=0;
			if ( search(coordenada, i) )	continue;
			a=retornar_indice(coordenada[0], coordenada[1], coordenada[2]);
			if((celulas[a] = crear_celula(coordenada[0], coordenada[1], coordenada[2])) == NULL) continue;
			espDifusivo2D[*coordenada][*(coordenada+1)][MAXGEN_SEN] = a;
		}
	}
	return(0);
}

//==============================================================================================//
//	FUNCIONES SIMULADORAS DE LA DINÁMICA GENÉTICA EN LAS CÉLULAS				//
//==============================================================================================//

int genes_estructurales(int i)
{
	int ret=0,k,l,act_tot, X[MAXGEN_REG];

	for (l=MAXGEN_SEN; l<(MAXGEN_EST+MAXGEN_SEN);l++)
	{                                                   
		act_tot=0;
		for (k=0;k<MAXGEN_REG;k++)
		{
			X[k] =  (*celulas[i]).Reguladores[k].afinidad[l][0] * (*celulas[i]).Reguladores[k].actividad - (*celulas[i]).Reguladores[k].afinidad[l][1];
			//LINEA LARGA DIVIDIDA
			if (X[k] > 0)
				X[k] = (*celulas[i]).Reguladores[k].afinidad[l][2];
			else 
				X[k] = 0;
			act_tot = act_tot + X[k];
		}
		
		(*celulas[i]).Estructurales[l-MAXGEN_SEN].produccion = 3*(1 + tanh(2*act_tot - 5));

		if ( ((*celulas[i]).Estructurales[l-MAXGEN_SEN].actividad = (*celulas[i]).Estructurales[l-MAXGEN_SEN].actividad + (*celulas[i]).Estructurales[l-MAXGEN_SEN].produccion - (*celulas[i]).Estructurales[l-MAXGEN_SEN].degradacion) < 0)
			(*celulas[i]).Estructurales[l-MAXGEN_SEN].actividad = 0;
		
		if ((*celulas[i]).Estructurales[l-MAXGEN_SEN].actividad < UMB_ACT)
			continue;
		//                                                fprintf(archivo, "%3d\t", act_tot);
// Ahora toca decidir en funcion del tipo de gen estructural, las acciones a realizar	
///*
		switch ((*celulas[i]).Estructurales[l-MAXGEN_SEN].tipo) 
		{	
			case	GENEST_DIVISION		:
				if(celulas_TOT > MAX_CEL)	break;
				(*celulas[i]).vida_muerte=DIVISION;
				break;

			case	GENEST_APOPTOSIS	:
				//if ((*celulas[i]).Estructurales[l-MAXGEN_SEN].actividad > UMB_ACT + 0.5)
				{
					(*celulas[i]).vida_muerte=APOPTOSIS;
					ret=APOPTOSIS;
				}
				break;
			case	GENEST_FIS_MU		:

			case	GENEST_FIS_K		:

			case	GENEST_FIS_GAMMA	:
				break;
		}
	} 
	return (ret);
}





//////////////////////////////////////////////////////////////////////////////////////////////

int genes_reguladores(int i)
{
	int k,l,act_tot, X[MAXGEN_SEN];
	for (l=0; l<MAXGEN_REG;l++)
	{
		act_tot=0;
		for (k=0;k<MAXGEN_SEN;k++)
		{
			X[k] = (*celulas[i]).Senales[k].afinidad[l][2] * (*celulas[i]).Senales[k].afinidad[l][0] * espDifusivo2D[(*celulas[i]).posicion.x][(*celulas[i]).posicion.y][k]/500 - (*celulas[i]).Senales[k].afinidad[l][1];
			//LINEA LARGA DIVIDIDA
			
			if (X[k] > 0)
				X[k] = (*celulas[i]).Senales[k].afinidad[l][2];
			else 
				X[k] = 0;
			act_tot = act_tot + X[k];
		}

		(*celulas[i]).Reguladores[l].produccion = (4 + (1 + 5*tanh(2*act_tot -5)));

		if( ((*celulas[i]).Reguladores[l].actividad = (*celulas[i]).Reguladores[l].actividad + (*celulas[i]).Reguladores[l].produccion -(*celulas[i]).Reguladores[l].degradacion) < 0 )
			(*celulas[i]).Reguladores[l].actividad = 0;
		//LINEA LARGA DIVIDIDA
	}	
	return (0);
}

//////////////////////////////////////////////////////////////////////////////////////////////

int genes_senales(int i)
{

	int k,l,act_tot, X[MAXGEN_REG];
	for (l=0; l<MAXGEN_SEN;l++)
	{
		act_tot=0;
		for (k=0;k<MAXGEN_REG;k++)
		{
			X[k] = (*celulas[i]).Reguladores[k].afinidad[l][2] * (*celulas[i]).Reguladores[k].afinidad[l][0] * (*celulas[i]).Reguladores[k].actividad - (*celulas[i]).Reguladores[k].afinidad[l][1];
			//LINEA LARGA DIVIDIDA
			if (X[k] > 0)
				X[k] = (*celulas[i]).Reguladores[k].afinidad[l][2];
			else 
				X[k] = 0;
			act_tot = act_tot + X[k];
		}
		if( ( (*celulas[i]).Senales[l].produccion = 500*(1 + tanh(2*act_tot -5)) ) < 0)
			(*celulas[i]).Senales[l].produccion = 0;
		else
			espDifusivo2D[(*celulas[i]).posicion.x][(*celulas[i]).posicion.y][l] = espDifusivo2D[(*celulas[i]).posicion.x][(*celulas[i]).posicion.y][l] + (*celulas[i]).Senales[l].produccion;
		//LINEA LARGA DIVIDIDA
	}

	return (0);
}


//==============================================================================================//
//	FUNCIONES PARA DESARROLLAR FUNCIONES CELULARES						//
//==============================================================================================//




//==============================================================================================//
//	FUNCION QUE ELIMINA LA I-ÉSIMA CÉLULA DEL ARREGLO *Celulas[Celulas_TOT]			//
//==============================================================================================//

int muerte(int i)
{
	int l;
	espDifusivo2D[(*celulas[i]).posicion.x][(*celulas[i]).posicion.y][MAXGEN_SEN]=0;
	for(l=i;l<celulas_TOT-1;l++)
	{
		celulas[l]=(celulas[l+1]);
		espDifusivo2D[(*celulas[l]).posicion.x][(*celulas[l]).posicion.y][MAXGEN_SEN]=l;
	}
	//free((void *) celulas[celulas_TOT - 1]);
	celulas_TOT--;	
	return (0);
}






//==============================================================================================//
//	FUNCION QUE BUSCA LA CÉLULA MÁS CERCANA AL PUNTO DONDE QUEREMOS CREAR UNA CELULA NUEVA	//
//==============================================================================================//

int retornar_indice(int x, int y, int z)
{
	int i, j, r,bandera=0;
	
		
		for(j = y; j >= 0; j--)
		{	
			if(r=espDifusivo2D[x][j][MAXGEN_SEN])
			{
				bandera=1;
				r += 1;
				break;
			}
			
		}
		
		
	for(i = x-1; i >= 0; i--)
	{
		if (bandera) break;
		for(j = ANCHO-1; j >= 0; j--)
		{	
			if(r=espDifusivo2D[i][j][MAXGEN_SEN])
			{
				bandera=1;
				r += 1;
				break;
			}
			
		}
		
	}
	for(i=celulas_TOT;i>r;i--)
	{
		celulas[i]=(celulas[i-1]);
		espDifusivo2D[(*celulas[i]).posicion.x][(*celulas[i]).posicion.y][MAXGEN_SEN]=i;
	}
	//printf("el indice buscado es : %d", r);
	return (r);

}



//==============================================================================================//
//	FUNCION CREADORA DE CÉLULAS A PETICION, REQUIERE UNAS COORDENADAS (x, y, z)		//
//==============================================================================================//


struct celula *crear_celula(int x,int y,int z)
{
	
	int i,j,r;
	struct celula *cel=NULL;
	cel = (struct celula *) malloc(sizeof(struct celula));
	if (inicializado == 0)
	{
		inicializado = 1;
		for(i = 0; i < MAXGEN_SEN; i++)
		{
			for(j=0; j < MAXGEN_REG; j++)
			{	
				aff_senales[i][j][0]=rand() % AFF;
				aff_senales[i][j][1]=rand() % UMB;
				aff_senales[i][j][2] = ((rand() % EFF) < 2) ?  -1 :  1;
			}
		}
	
		for(i = 0; i<MAXGEN_REG; i++)
		{
			for(j=0; j < (MAXGEN_SEN + MAXGEN_EST); j++)
			{	
				aff_reguladores[i][j][0]=rand() % AFF;
				aff_reguladores[i][j][1]=rand() % UMB;
				aff_reguladores[i][j][2] = ((rand() % EFF) < 2) ?  -1 : 1;
			}
		}

		
	}



/////////////////////////////////////////////////////////////////////////////

	for(i = 0; i < MAXGEN_SEN; i++)
	{
		for(j=0; j < MAXGEN_REG; j++)
		{	
			(*cel).Senales[i].afinidad[j][0] = aff_senales[i][j][0];
			(*cel).Senales[i].afinidad[j][1] = aff_senales[i][j][1];
			(*cel).Senales[i].afinidad[j][2] = aff_senales[i][j][2];
		}

		(*cel).Senales[i].difusion = rand() % DIF;
		(*cel).Senales[i].coef_difusion = rand() % DIF;
		(*cel).Senales[i].produccion = rand() % 5;	

	}

	for(i = 0; i<MAXGEN_REG; i++)
	{
		for(j=0; j < (MAXGEN_SEN + MAXGEN_EST); j++)
		{	
			(*cel).Reguladores[i].afinidad[j][0] = aff_reguladores[i][j][0];
			(*cel).Reguladores[i].afinidad[j][1] = aff_reguladores[i][j][1];
			(*cel).Reguladores[i].afinidad[j][2] = aff_reguladores[i][j][2];
		} 
		
		(*cel).Reguladores[i].actividad = rand() % 4;
		(*cel).Reguladores[i].produccion = rand() % 4;
		(*cel).Reguladores[i].degradacion = rand() % DEGR;
	}

	for(i = 0; i<MAXGEN_EST; i++)
	{
		(*cel).Estructurales[i].degradacion = rand() % DEGR;
		(*cel).Estructurales[i].actividad = rand() % 3;
		(*cel).Estructurales[i].produccion = rand() % 3;
	}
	(*cel).Estructurales[0].tipo = GENEST_APOPTOSIS	;
	(*cel).Estructurales[1].tipo = GENEST_DIVISION	;
	

	(*cel).posicion.x = x;
	(*cel).posicion.y = y;
	(*cel).posicion.z = 0;

	(*cel).fisica.k = rand() % AFF;
	(*cel).fisica.mu = rand() % AFF;
	(*cel).fisica.vecinos = 4;

	(*cel).vida_muerte=0;
	celulas_TOT++;
	
	return(cel);
}





//==============================================================================================//
//   FUNCION QUE BUSCA ESPACIO LIBRE DONDE CREAR UNA CELULA, DEVUELVE  COORDENADAS (x, y, z)	//
//==============================================================================================//


int search(int *coor, int i)
{
int j,k,l,r, ccc=10;
j = (*celulas[i]).posicion.x;
k = (*celulas[i]).posicion.y;
l = (*celulas[i]).posicion.z;
while (1)
{	r = rand() % 8;
	
	if ( (j == 0) | (k == 0) | (j == (ANCHO -1)) | (k == (LARGO -1)) ) return 1;
	if(!(ccc--)) return 1;

	switch (r)
	{
			case 0	:
			if (espDifusivo2D[j+1][k][MAXGEN_SEN] == 0)
				{ *coor = j+1; *(coor + 1) = k; *(coor + 2) = 0; return 0;}
			break;
			case 1	:
			if (espDifusivo2D[j-1][k][MAXGEN_SEN] == 0)
				{ *coor = j-1; *(coor + 1)=k; *(coor + 2)=0;  return 0;}
			break;
			case 2	:
			if (espDifusivo2D[j+1][k+1][MAXGEN_SEN] == 0)
				{ *coor  =j+1; *(coor + 1)=k+1; *(coor +2)=0; return 0;}
			break;
			case 3	:
			if (espDifusivo2D[j-1][k-1][MAXGEN_SEN] == 0)
				{ *coor =j-1; *(coor + 1)=k-1; *(coor +2)=0; return 0;}
			break;
			case 4	:
			if (espDifusivo2D[j][k-1][MAXGEN_SEN] == 0)
				{ *coor =j; *(coor + 1)=k-1; *(coor +2)=0; return 0;}
			break;
			case 5	:
			if (espDifusivo2D[j][k+1][MAXGEN_SEN] == 0)
				{ *coor =j; *(coor + 1)=k+1; *(coor +2)=0; return 0;}
			break;
		 	case 6	:
			if (espDifusivo2D[j+1][k-1][MAXGEN_SEN] == 0)
				{ *coor =j+1; *(coor + 1)=k-1; *(coor +2)=0; return 0;}
			break;
			case 7	:
			if (espDifusivo2D[j-1][k+1][MAXGEN_SEN] == 0)
				{ *coor =j-1; *(coor + 1)=k+1; *(coor +2)=0; return 0;}
			break;
			
			
			
	}
}
//return 1;
}




//==============================================================================================//
//   FUNCIONES PARA TAREAS AUXILIARES DE PRESENTACION DE DATOS	Y TAL				//
//==============================================================================================//




	//-------------------------------------------------------
	//	IMPRESION DE LA CONFIGURACION INICIAL EN conf_inicial 


int print_init_conf(void)
{
/*

	int k,s,cont;
	FILE *file=NULL;
	if ( !(file = fopen("datos_inicial", "w")) )
		printf("No puedo abrir datos\n");
	
	for (k=0; k<celulas_TOT; k++) 	
	{
		fprintf(file, "\n\ncélula nº:  %d\t Dirección puntero:  %d\n", k, &celulas[k]);
		fprintf(file, "ESTADO DE SUS GENES\n");
		fprintf(file, "\tSEÑAL 1\n");
		fprintf(file, "\tDifusion \tAfinidades \n");
		fprintf(file, "\t%f", (*celulas[k]).Senales[0].difusion);
		for (cont=0; cont<4; cont++) 
			fprintf(file, "\t%f", (*celulas[k]).Senales[0].afinidad[cont][0]);
		fprintf(file, "\n\tSEÑAL 2\n");
		fprintf(file, "\tDifusion \tAfinidades \n");
		fprintf(file, "\t%f", (*celulas[k]).Senales[1].difusion);
		for (cont=0; cont<4; cont++) 
			fprintf(file, "\t%f", (*celulas[k]).Senales[1].afinidad[cont][0]);
		
		fprintf(file, "\n\tREGULADOR 1\n");
		fprintf(file, "\tDegradacion \tAfinidades \n");
		fprintf(file, "\t%f", (*celulas[k]).Reguladores[0].degradacion);
		for (cont=0; cont<4; cont++) 
			fprintf(file, "\t%f", (*celulas[k]).Reguladores[0].afinidad[cont][0]);
		fprintf(file, "\n\tREGULADOR 2\n");
		fprintf(file, "\tDegradacion \tAfinidades \n");
		fprintf(file, "\t%f\n", (*celulas[k]).Reguladores[1].degradacion);
		for (cont=0; cont<4; cont++) 
			fprintf(file, "\t%f", (*celulas[k]).Reguladores[1].afinidad[cont][0]);
		fprintf(file, "\n\tREGULADOR 3\n");
		fprintf(file, "\tDegradacion \tAfinidades \n");
		fprintf(file, "\t%f\n", (*celulas[k]).Reguladores[2].degradacion);
		for (cont=0; cont<4; cont++) 
			fprintf(file, "\t%f", (*celulas[k]).Reguladores[2].afinidad[cont][0]);
		fprintf(file, "\nPOSICION\n");
fprintf(file, "\tx=%d \ty=%d \tz=%d\n",(*celulas[k]).posicion.x, (*celulas[k]).posicion.y, (*celulas[k]).posicion.z);
		for (s=0; s<MAXGEN_SEN;s++)
		{
		fprintf(file, "\n\nVALOR DE LAS CONCENTRACIONES DE MORFOGEN\n");
fprintf(file, "\tSEÑAL=%d \tCONCENTRACION=%d", s, espDifusivo2D[(*celulas[k]).posicion.x][(*celulas[k]).posicion.y][s]);
		}
		fprintf(file, "\n\n\n");
		
	}
	fclose(file);*/
	return(0);
}



int print_red_genetica(void)
{
/*	int i,j;
	FILE *archivo1=NULL;
	FILE *archivo2=NULL;
	if ( !(archivo1 = fopen("datos_celulas", "a")) )
		printf("No puedo abrir datos\n");
	if ( !(archivo2 = fopen("datos_produccion", "a")) )
		printf("No puedo abrir datos\n");
	
	fprintf(archivo1, "\n==================================================================\n");
	fprintf(archivo1, "\titeracion nº:  %f\t Celulas Totales:  %d\n", t_simulacion, celulas_TOT);
	fprintf(archivo1, "==================================================================\n");
		
	
	for (i=0;i<ANCHO;i++)
	{
		for (j=0;j<LARGO;j++)
		{
			fprintf(archivo1, "%4d", espDifusivo2D[i][j][MAXGEN_SEN]);
		}
		fprintf(archivo1, "\n");
	}
		fprintf(archivo2, "\n==================================================================\n");
		fprintf(archivo2, "\titeracion nº:  %f\t Celula nº:  %d\n", t_simulacion, i);
		fprintf(archivo2, "==================================================================\n");
		fprintf(archivo2, "Estructurales\t\t%4f\t\t Act %4f\n", (*celulas[7]).Estructurales[0].produccion, (*celulas[5]).Estructurales[0].actividad);	
		fprintf(archivo2, "Estructurales\t\t%4f\t\t Act %4f\n", (*celulas[7]).Estructurales[1].produccion,(*celulas[5]).Estructurales[1].actividad);	
		fprintf(archivo2, "Reguladores\t\t%4f\t\t Act %4f\n", (*celulas[7]).Reguladores[0].produccion,(*celulas[5]).Reguladores[0].actividad);	
		fprintf(archivo2, "Reguladores\t\t%4f\t\t Act %4f\n", (*celulas[7]).Reguladores[1].produccion,(*celulas[5]).Reguladores[1].actividad);	
		fprintf(archivo2, "Reguladores\t\t%4f\t\t Act %4f\n", (*celulas[7]).Reguladores[2].produccion,(*celulas[5]).Reguladores[2].actividad);
		fprintf(archivo2, "Reguladores\t\t%4f\t\t Act %4f\n", (*celulas[7]).Reguladores[3].produccion,(*celulas[5]).Reguladores[3].actividad);	
		fprintf(archivo2, "Señales\t\t%4f\t\t Act %4f\n", (*celulas[7]).Senales[0].produccion);
		fprintf(archivo2, "Señales\t\t%4f\t\t Act %4f\n", (*celulas[7]).Senales[1].produccion);
		
	//(*celulas[5]).Reguladores[0].actividad
	//(*celulas[5]).Reguladores[1].actividad
		
	fclose(archivo1);	
	fclose(archivo2);
	return (0);

}

int print_genestructurales()
{
int l;
FILE *archivo=NULL;
for (l=MAXGEN_SEN; l<(MAXGEN_EST+MAXGEN_SEN);l++)
{
	if ( !(archivo = fopen("datos_estructurales", "a")) )
		printf("No puedo abrir datos\n");
	fprintf(archivo, "\n%5f\t\t%3f\n", (*celulas[7]).Estructurales[0].produccion, (*celulas[7]).Estructurales[l-MAXGEN_SEN].actividad);
	fprintf(archivo, "\n%5f\t\t%3f\n", (*celulas[7]).Estructurales[1].produccion, (*celulas[7]).Estructurales[l-MAXGEN_SEN].actividad);
	fprintf(archivo, "\n%5f\t\t%3f\n", (*celulas[7]).Estructurales[0].produccion, (*celulas[7]).Estructurales[l-MAXGEN_SEN].actividad);
	fprintf(archivo, "\n%5f\t\t%3f\n", (*celulas[7]).Estructurales[1].produccion, (*celulas[7]).Estructurales[l-MAXGEN_SEN].actividad);
}
fclose(archivo);*/
return(0);
}

int print_espDifusivo(int s)
{
	int i,j;
	FILE *archivo1, *archivo2, *archivo3;
	archivo1 = fopen("datos_difusion1", "a");
	archivo2 = fopen("datos_difusion2", "a");
	archivo3 = fopen("datos_difusion3", "a");
	
	
	switch(s)
	{
	case 0:
	{
		fprintf(archivo1, "\n==================================================================\n");
		fprintf(archivo1, "\titeracion nº:  %f\t Señal:  %d\n", t_simulacion, s);
		fprintf(archivo1, "==================================================================\n");
		for (i=0; i<LARGO; i++)
		{
			for(j=0; j<ANCHO; j++)
			{	
				fprintf(archivo1, "%3d ", espDifusivo2D_old[i][j][s] = espDifusivo2D[i][j][s]);
				if (j == (ANCHO - 1))
					fprintf(archivo1,"\n");
			}
			
		}
	break;
	}
	
	case 1:
	{
		fprintf(archivo2, "\n==================================================================\n");
		fprintf(archivo2, "\titeracion nº:  %f\t Señal:  %d\n", t_simulacion, s);
		fprintf(archivo2, "==================================================================\n");
		for (i=0; i<LARGO; i++)
		{
			for(j=0; j<ANCHO; j++)
			{	
				fprintf(archivo2, "%3d ", espDifusivo2D_old[i][j][s] = espDifusivo2D[i][j][s]);
				if (j == (ANCHO - 1))
					fprintf(archivo2,"\n");
			}
			
		}
	break;
	}	

	case 2:
	{
		fprintf(archivo3, "\n==================================================================\n");
		fprintf(archivo3, "\titeracion nº:  %f\t Señal:  %d\n", t_simulacion, s);
		fprintf(archivo3, "==================================================================\n");
		for (i=0; i<LARGO; i++)
		{
			for(j=0; j<ANCHO; j++)
			{	
				fprintf(archivo3, "%3d ", espDifusivo2D_old[i][j][s] = espDifusivo2D[i][j][s]);
				if (j == (ANCHO - 1))
					fprintf(archivo3,"\n");
			}
			
		}
	break;
	}
	}
	fclose(archivo1);
	fclose(archivo2);
	fclose(archivo3);
	return(0);
}
/////////////////////////////////////////////////////////////////////////////


int dibujar_red_genetica(void)
{

	static float contador = 0;
	int dec;
	int sing;
	int len, q, i, j,s;
	int rangocolor = (int) (65535 * 1.5 / (LARGO*DIF));
	char nombre[30]="celulas";
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

	pngwriter one((ANCHO*FACTOR),(LARGO*FACTOR),0,nombre);
	one.setcompressionlevel(9);

	
	for (i=0; i<LARGO; i++)
	{
		
		for(j=0; j<ANCHO; j++)
		{
			if (espDifusivo2D[i][j][MAXGEN_SEN])
				one.filledsquare((i * FACTOR),(j * FACTOR),(i * FACTOR + FACTOR),(j * FACTOR + FACTOR),1.0,1.0,1.0);
		}
	}
	
	one.close();
	return (0);


}
