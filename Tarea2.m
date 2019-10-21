close all
clear all
clc

%% ingresando datos
elasticidad = input('Ingrese el modulo de eslasticidad: ');
area = input('Ingrese el area: ');
longitud = input('Ingrese la longitud: ');
altura = input('Ingrese la altura del puente: ');
tramos = input('Ingrese la cantidad de tramos: ');

%% grados de libertad globales
if(tramos ==0)
    n = 4;
elseif(tramos>0)
    n = 4+2*tramos;
end

for i=1:n
   gdl(i,:) = [2*i-1 2*i]; 
end
%% coordenadas del sistema

xy= zeros(n, 2);

m= n/2+2;

for i=2:m
    xy(i,1)= longitud*(i-1) 
end

for j= m:n
    xy(j,1)= longitud*()
end

for j= m:n
    xy(j,2)= altura
end





