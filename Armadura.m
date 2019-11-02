
% Como ejemplo, consideremos la armadura met�lica mostrada en la figura anterior
%sometida a la acci�n % de dos cargas verticales y una horizontal. 
%Todas las barras tienen un m�dulo de elasticidad E = 2 �10^8kN/m2 
%y un �rea seccional A = 0.005m2
%%
close all
clear all
clc
%% Datos

E=2.08e7; A=0.0013;

l_1=1.25; 
l_2=2.5; 
l_3=sqrt(1.25^2+2.5^2);
l_4=2.5; 
l_5=1.25; 

%% Encontramos los grados de libertad globales segun su nodo
n=4;
for i=1:n
gdl(i,:)=[2*i-1 2*i];
end
%% Coordenadas del Sistema
xy=[0 0
    0 1.25
    2.5 0
    2.5 1.25];
%% Angulos de las barras de la armadura
Ang=[90
    0
    -26.57
    0
    90];
 %% Calculamos la matriz de rigidez por elemento
 k_1=Barra_Axial2(E,A,l_1);
 k_2=Barra_Axial2(E,A,l_2);
 k_3=Barra_Axial2(E,A,l_3);
 k_4=Barra_Axial2(E,A,l_4);
 k_5=Barra_Axial2(E,A,l_5);
 
%% Matriz de Transformacion por Elemento
%Elemoento 1
Ang_1=Ang(1);
eta=cosd(Ang_1);
mu=sind(Ang_1);

T_1= [eta mu 0 0;-mu eta 0 0;
    0 0 eta mu; 0 0 -mu eta];

K_1=T_1'*k_1*T_1;

%Elemento 2
Ang_2=Ang(2);
eta=cosd(Ang_2);
mu=sind(Ang_2);

T_2= [eta mu 0 0;-mu eta 0 0;
    0 0 eta mu; 0 0 -mu eta];

K_2=T_2'*k_2*T_2;

%Elemento 3
Ang_3=Ang(3);
eta=cosd(Ang_3);
mu=sind(Ang_3);

T_3= [eta mu 0 0;-mu eta 0 0;
    0 0 eta mu; 0 0 -mu eta];

K_3=T_3'*k_3*T_3;
%Elemento 4
Ang_4=Ang(4);
eta=cosd(Ang_4);
mu=sind(Ang_4);

T_4= [eta mu 0 0;-mu eta 0 0;
    0 0 eta mu; 0 0 -mu eta];

K_4=T_4'*k_4*T_4;
%Elemento 5
Ang_5=Ang(5);
eta=cosd(Ang_5);
mu=sind(Ang_5);

T_5= [eta mu 0 0;-mu eta 0 0;
    0 0 eta mu; 0 0 -mu eta];

K_5=T_5'*k_5*T_5;

%% Ensamblaje de Matriz de Rigidez
n=4;
K=zeros(n*2,n*2);
%Elemento 1
g_1=[1 2 3 4];
DeltaK_1=zeros(8,8);
DeltaK_1(g_1,g_1)=K_1;
K=K+DeltaK_1;
%Elemento 2
g_2=[1 2 5 6];
DeltaK_2=zeros(8,8);
DeltaK_2(g_2,g_2)=K_2;
K=K+DeltaK_2;
%Elemento 3
g_3=[3 4 5 6];
DeltaK_3=zeros(8,8);
DeltaK_3(g_3,g_3)=K_3;
K=K+DeltaK_3;
%Elemento 4
g_4=[3 4 7 8];
DeltaK_4=zeros(8,8);
DeltaK_4(g_4,g_4)=K_4;
K=K+DeltaK_4;
%Elemento 5
g_5=[5 6 7 8];
DeltaK_5=zeros(8,8);
DeltaK_5(g_5,g_5)=K_5;
K=K+DeltaK_5;

%% Grados de Libertad Libres y Restringidos
a=[1 2 5 7 8]; b=[3 4 6];

K_aa=K(a,a); K_ab=K(a,b); K_ba=K(b,a); K_bb=K(b,b);

P=[0 0 0 0 0 -568 0 0]';
%%1.93e-14 fuerza para desplazar 2.5cm hacia abajo
P_b=P(b);

D_b=K_bb\P_b;
P_a=K_ab*D_b;

D=zeros(8,1);
D(b)=D_b
D(6,1)
%% Esfuerzos Internos
beta=90;
eta=cosd(beta); mu=sind(beta);
D_1=D(g_1);
sigma_1=E*[-eta -mu eta mu]*D_1/l_1

beta=0;
eta=cosd(beta); mu=sind(beta);
D_2=D(g_2);
sigma_2=E*[-eta -mu eta mu]*D_2/l_2

beta=-26.57;
eta=cosd(beta); mu=sind(beta);
D_3=D(g_3);
sigma_3=E*[-eta -mu eta mu]*D_3/l_3

beta=0;
eta=cosd(beta); mu=sind(beta);
D_4=D(g_4);
sigma_4=E*[-eta -mu eta mu]*D_4/l_4

beta=90;
eta=cosd(beta); mu=sind(beta);
D_5=D(g_5);
sigma_5=E*[-eta -mu eta mu]*D_5/l_5

%% Fuerzas 
N_1=A*sigma_1
N_2=A*sigma_2
N_3=A*sigma_3
N_4=A*sigma_4
N_5=A*sigma_5

%% Plotear Deformacion
XY=zeros(4,2);
XY(1,:)=[0 0];
XY(2,:)=[0 1.25];
XY(3,:)=[2.5 0];
XY(4,:)=[2.5 1.25];

XYdef=zeros(size(XY));
fac=5;
c=0;
for i=1:4
c=c+1;
XYdef(i,1)=XY(i,1)+fac*D(c);
c=c+1;
XYdef(i,2)=XY(i,2)+fac*D(c);
end

IJ=zeros(5,2);
IJ(1,:)=[1 2];
IJ(2,:)=[1 3];
IJ(3,:)=[2 3];
IJ(4,:)=[2 4];
IJ(5,:)=[3 4];

%% figura
for e=1:5
Q=[XY(IJ(e,1),1) XY(IJ(e,1),2);...
XY(IJ(e,2),1) XY(IJ(e,2),2)];
Qdef=[XYdef(IJ(e,1),1) XYdef(IJ(e,1),2);...
XYdef(IJ(e,2),1) XYdef(IJ(e,2),2)];
plot(Q(:,1),Q(:,2),'--b',Qdef(:,1),Qdef(:,2),'-r')
hold on
end
xlabel('x')
ylabel('y')
axis equal

%% funcion

function [k]=Barra_Axial2(E,A,Le)

k=E*A/Le*[1 0 -1 0
          0 0 0 0
          -1 0 1 0
          0 0 0 0];
end