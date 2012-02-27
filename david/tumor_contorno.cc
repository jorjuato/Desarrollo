/*________________________________________________________________________________________________________
 
programa C Simulacion desarrollo: difusion y red genetica
___________________________________________________________________________________________________________
*/

//#include <cstdlib>	
#include <pngwriter.h>
#include <iostream.h>
//#include <string>
//#include <math.h>
#include <stdio.h>

//#include <time.h>


#define 	LARGO		201	
#define 	ALTO		1
#define 	ANCHO		200
#define 	FACTOR		1


int condiciones_iniciales[LARGO][ANCHO];
int leer_imagen(void);
int print_archivo();	

int main(void)
{
	
	leer_imagen();
	print_archivo();
	return(0);
}



int leer_imagen()
{
   int i,j;

   pngwriter tumor(1,1,0,"copia_tumor.png"); 
  
   /* readfromfile()
    * Ahora especificamos la ubicacion del archivo que deseamos colocar
    * dentro de la instancia de PNGwriter llamada tumor.
    * */
   std::cout << "Opening tumor.png...";   
   tumor.readfromfile("tumor.png"); //It really is that easy.
   std::cout << " done." << std::endl;
   

   int tumorwidth = tumor.getwidth();
   int tumorheight = tumor.getheight();
  
   
   std::cout << "The image that has just been read from disk (tumor.png) is " << tumor.getheight();
   std::cout << " pixels high and " << tumor.getwidth()<<" pixels wide."<<std::endl;
   std::cout << "Bit depth is " << tumor.getbitdepth()<<std::endl;
   std::cout << "Image gamma is: " << tumor.getgamma() << std::endl; 
   for (i=1;i<ANCHO;i++)
	{
	for (j=1;j<LARGO;j++)
		condiciones_iniciales[i][j]=tumor.dread(i,j);
	}
   tumor.close();
return (0);
}

int print_archivo()
{

	int i,j,num;
	FILE *archivo1=NULL;
	if ( !(archivo1 = fopen("matriz_inicial", "w")) )
		printf("No puedo abrir datos\n");
	
	
	fprintf(archivo1, "\n==================================================================\n");
	fprintf(archivo1, "Matriz de la imagen \n");
	fprintf(archivo1, "==================================================================\n");
	for (i=1;i<ANCHO;i++)
	{
		for (j=1;j<LARGO;j++)
		{	num=(condiciones_iniciales[i][j]); //<ceil
			fprintf(archivo1, "\t%f",num);
		}
		fprintf(archivo1, "\n");
	}		
fclose(archivo1);

return(0);
}
