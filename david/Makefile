############ MAKEFILE FOR PNGWRITER ######################################
#
#
#       Website: Main:             http://pngwriter.sourceforge.net/
#                Sourceforge.net:  http://sourceforge.net/projects/pngwriter/
#                Freshmeat.net:    http://freshmeat.net/projects/pngwriter/
#  
#       Author:                    Paul Blackburn
#
#       Email:                     individual61@users.sourceforge.net
#
#       Version:                   0.5.2   (10 / I / 2005)
#
#       Description:               Library that allows plotting a 48 bit
#                                  PNG image pixel by pixel, which can 
#                                  then be opened with a graphics program.
#  
#       License:                   GNU General Public License
#                                  © 2002, 2003, 2004, 2005 Paul Blackburn
#
###############################################################################

######################## IMPORTANT - IMPORTANTE ###############################
# ENGLISH
# This makefile is meant to help beginning programmers learn about simple 
# Makefiles. The compilation uses the object file, and not the library,
# because I belive it to be a more 'hands-on' approach.
#
# CASTELLANO/ESPANIOL
# Este makefile se hizo teniendo en mente a los programadores principiantes,
# para ayudarles a aprender acerca de Makefiles simples. La compilacion usa
# el archivo objeto, no la libreria, porque creo que es una forma mas
# aplicada de hacer las cosas.
# La explicacion en castellano se encuentra al final de este Makefile.
###############################################################################


include ../make.include
CXX= g++
EXAMPLES= leer_arch

all: $(EXAMPLES)


leer_arch: tumor_contorno.cc
	@ echo "#"
	@ echo "#"
	@ echo "#  COMPILACION DEL SIMULADOR DE DIFUSION"
	@ echo "#  GRAFICOS DE SALIDA EN FORMATO PNG"
	@ echo "#  "
	@ echo "#  LIBRERIA GRAFICA EMPLEADA"
	@ echo "#  pngWriter C++, linux/windows"
	@ echo "#"
	@ echo "#  http://pngwriter.sourceforge.net/"
	@ echo "#"
	
	$(CXX) $(CXXFLAGS) $(INC) tumor_contorno.cc -g -o tumor_contorno $(LIBS)

clean	:    
	rm -f $(EXAMPLES) *~ arcoiris.png copiaburro.png one.png two.png
	rm -f lyapunov.cc~ lyapunov.espaniol.cc~ Makefile~ pngtest.espaniol.cc~
	rm -f .DS_Store out.png triangles.png
	



################################### CASTELLANO / ESPANIOL ###############################################
# INTRODUCCION
# Este makefile te puede ayudar si no sabes como usar uno. 
#
# Lo importante que
# tienes que saber es que lo basico para usar PNGwriter es que se pueda
# generar pngwriter.o. Necesita de libpng, una libreria de funciones para
# facilitar la vida a la gente que trabaja con lectores/creadores de imagenes
# PNG. Debes tener esta libreria ya instalada. Si no es asi, anda a
# www.libpng.org y bajala.
#
# Nota: Este documento se ha escrito sin acentos para que sea mas portatil.
#
# OPCIONES
# Que Opciones De Compilacion Puedo Elegir? *IMPORTANTE*
#
# PNGwriter puede ser compilado en muchas configuraciones distintas de sistemas,
# pero tu tendras que crear un Makefile personalizado.
#
#
# ASPECTOS BASICOS DE LOS MAKEFILE
# Que significan las opciones que se usan aqui?
#
# -O3 (raya oh tres). Optimizacion del codigo. 1 es leve, 3 es maximo.
#
# -o "Lo que produzcas, crealo con este nombre".
#
# -c Compila el codigo, pero no lo linkees. Piensa en esto como un carpintero
# que hace una silla pero todavia no pega las distintas partes. El respaldo es
# un respaldo, pero como silla, es inutil, necesita del resto de las partes.
# Eso hace el linker. Un .o es como el respaldo de la silla, o una pata. La gracia es que
# PNGwriter es como una de las patas, y uno puede usar la pata de PNGwriter ya
# tallada (ya compilada) en otra silla (otro proyecto).
#
# -I (i mayuscula): Si tu, como yo, tienes algunos "include" (.h) instalados en un lugar no
# estandard, especifica donde con -I. 
#
# -L: Lo mismo, pero para librerias.
#
# -l (L minuscula): Para usar una libreria. Para usar libpng, usa -lpng. Para Zlib,
# usa -lz.
#######################################################################################

