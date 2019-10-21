close all
clear all
clc

%% ingresando datos
elasticidad = input('Ingrese el modulo de eslasticidad: ');
area = input('Ingrese el area: ');
longitud = input('Ingrese la longitud: ');
tramos = input('Ingrese la cantidad de tramos: ');
altura = input('Ingrese la altura del puente: ');

%% grados de libertad globales
if(tramos ==0)
    n = 4
elseif(tramos==1)
    
elseif(tramos>1)
    
end
for i=1:n
   gdl(i,:) = [2*i-1 2*i]; 
end







