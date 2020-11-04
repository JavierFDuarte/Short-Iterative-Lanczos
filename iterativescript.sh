#!/bin/bash

#NOTA: Hay que estar en un enviroment con conda y numpy instalados... ACTIVARLOS A MANOS



#Loop-Info
	#Se pueden definir variables en un bashscript, como por ejemplo:
	#A tener en cuenta las variables cuando son llamadas, se hace con el signo $

ruta=/home/jfd/Escritorio/Tesis/Activ/T5-Lanczos;
	#Tambien se pueden hacer loops. Este por ejemplo crea archivos
	#Nota en el seq se usa la tilde al revés, la tecla que está arriba del shift y a la derecha del acento circunflejo

cp $ruta/Respaldo/initialstate.csv $ruta
#Compilation
gfortran Lanczos-TempEvol.f90 -l lapack -l blas -o Lanczos-TempEvol.x;

#Iterations
rm -r $ruta/Iteraciones
statename=initialstate;
mkdir $ruta/Iteraciones;
for k in `seq 1 3`
do
	python Lanczos-Algth.py;
	./Lanczos-TempEvol.x;
	python Lanczos-InitStateFormer.py;
	#Me guardo el resultado final para hacer plots después, el cual es el estado inicial de la siguiente iteracion ---> Me lo guardo en CSV así nomás, después graficaré en mathemática, Fortran+GNUPlot = HOSHIBLE
	cp $ruta/$statename.csv $ruta/$statename$k.csv;
	mv $ruta/$statename$k.csv $ruta/Iteraciones/;
done


#Por alguna razón cuando pongo asó el seq sí me quiere tomar el límite superior... Bueno, en vez de poner el n de los códigos, poner n-1 acá
n=3

for kk in `seq 0 $n`
do
	rm temporalfile_v$kk.csv
done
rm coeficientesBeta.txt
rm coeficientesAlpha.txt
rm CoefEvTemp.txt
rm Lanczos-TempEvol.x
rm initialstate.csv


