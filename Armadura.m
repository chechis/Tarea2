
% Como ejemplo, consideremos la armadura met�lica mostrada en la figura anterior
%sometida a la acci�n % de dos cargas verticales y una horizontal. 
%Todas las barras tienen un m�dulo de elasticidad E = 2 �10^8kN/m2 
%y un �rea seccional A = 0.005m2
%%
close all
clear all
clc
%% Datos

E=2e8; A=0.005;

l_1=4; 
l_2=sqrt(4^2+3^2); 
l_3=3; 
l_4=4; 
l_5=sqrt(4^2+3^2);
l_6=4; 
l_7=3; 
l_8=4; 
l_9=sqrt(4^2+3^2);
%% Encontramos los grados de libertad globales segun su nodo
n=6;
for i=1:n
gdl(i,:)=[2*i-1 2*i];
end
%% Coordenadas del Sistema
xy=[0 0
    4 0
    4 3
    8 0
    8 3
    12 0];
%% Angulos de las barras de la armadura

Ang=[0
    36.87
    90
    0
    -36.87
    0
    90
    0
    -36.87];

 % Ang=rad2deg(asin(xy(

 %% Calculamos la matriz de rigidez por elemento
 k_1=Barra_Axial2(E,A,l_1);
 k_2=Barra_Axial2(E,A,l_2);
 k_3=Barra_Axial2(E,A,l_3);
 k_4=Barra_Axial2(E,A,l_4);
 k_5=Barra_Axial2(E,A,l_5);
 k_6=Barra_Axial2(E,A,l_6);
 k_7=Barra_Axial2(E,A,l_7);
 k_8=Barra_Axial2(E,A,l_8);
 k_9=Barra_Axial2(E,A,l_9);

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
beta=90;
eta=cosd(beta); mu=sind(beta);
%Elemento 3
T_3= [ eta mu 0 0; -mu eta 0 0;
0 0 eta mu; 0 0 -mu eta];
K_3=T_3'*k_3*T_3;
beta=0;
eta=cosd(beta); mu=sind(beta);
%Elemento 4
T_4= [ eta mu 0 0; -mu eta 0 0;
0 0 eta mu; 0 0 -mu eta];
K_4=T_4'*k_4*T_4;
beta=-36.87;
eta=cosd(beta); mu=sind(beta);
%Elemento 5
T_5= [ eta mu 0 0; -mu eta 0 0;
0 0 eta mu; 0 0 -mu eta];
K_5=T_5'*k_5*T_5;
beta=0;
eta=cosd(beta); mu=sind(beta);

%Elemento 6
T_6= [ eta mu 0 0; -mu eta 0 0;
0 0 eta mu; 0 0 -mu eta];
K_6=T_6'*k_6*T_6;
beta=90;
eta=cosd(beta); mu=sind(beta);

%Elemento 7
T_7= [ eta mu 0 0; -mu eta 0 0;
0 0 eta mu; 0 0 -mu eta];
K_7=T_7'*k_7*T_7;
beta=0;
eta=cosd(beta); mu=sind(beta);
%Elemento 8
T_8= [ eta mu 0 0; -mu eta 0 0;
0 0 eta mu; 0 0 -mu eta];
K_8=T_8'*k_8*T_8;
beta=-36.87;
eta=cosd(beta); mu=sind(beta);
%Elemento 9
T_9= [ eta mu 0 0; -mu eta 0 0;
0 0 eta mu; 0 0 -mu eta];
K_9=T_9'*k_9*T_9;
%% Ensamblaje de Matriz de Rigidez
n=6;
K=zeros(n*2,n*2);
%Elemento 1
g_1=[1 2 3 4];
DeltaK_1=zeros(12,12);
DeltaK_1(g_1,g_1)=K_1;
K=K+DeltaK_1;
%Elemento 2
g_2=[1 2 5 6];
DeltaK_2=zeros(12,12);
DeltaK_2(g_2,g_2)=K_2;
K=K+DeltaK_2;
%Elemento 3
g_3=[3 4 5 6];
DeltaK_3=zeros(12,12);
DeltaK_3(g_3,g_3)=K_3;
K=K+DeltaK_3;
%Elemento 4
g_4=[3 4 7 8];
DeltaK_4=zeros(12,12);
DeltaK_4(g_4,g_4)=K_4;
K=K+DeltaK_4;
%Elemento 5
g_5=[5 6 7 8];
DeltaK_5=zeros(12,12);
DeltaK_5(g_5,g_5)=K_5;
K=K+DeltaK_5;
%Elemento 6
g_6=[5 6 9 10];
DeltaK_6=zeros(12,12);
DeltaK_6(g_6,g_6)=K_6;
K=K+DeltaK_6;
%Elemento 7
g_7=[7 8 9 10];
DeltaK_7=zeros(12,12);
DeltaK_7(g_7,g_7)=K_7;
K=K+DeltaK_7;
%Elemento 8
g_8=[7 8 11 12];
DeltaK_8=zeros(12,12);
DeltaK_8(g_8,g_8)=K_8;
K=K+DeltaK_8;
%Elemento 9
g_9=[9 10 11 12];
DeltaK_9=zeros(12,12);
DeltaK_9(g_9,g_9)=K_9;
K=K+DeltaK_9;

%% Grados de Libertad Libres y Restringidos
a=[1 2 11 12]; b=[3 4 5 6 7 8 9 10];

K_aa=K(a,a); K_ab=K(a,b); K_ba=K(b,a); K_bb=K(b,b);

P=[0 0 0 0 4 -20 0 0 0 -20 0 0]';

P_b=P(b);

D_b=K_bb\P_b;
P_a=K_ab*D_b;

D=zeros(12,1);
D(b)=D_b;
%% Esfuerzos Internos
beta=0;

eta=cosd(beta); mu=sind(beta);
D_1=D(g_1);
sigma_1=E*[-eta -mu eta mu]*D_1/l_1;

beta=36.87;
eta=cosd(beta); mu=sind(beta);
D_2=D(g_2);
sigma_2=E*[-eta -mu eta mu]*D_2/l_2;
beta=90;
eta=cosd(beta); mu=sind(beta);
D_3=D(g_3);
sigma_3=E*[-eta -mu eta mu]*D_3/l_3;
beta=0;
eta=cosd(beta); mu=sind(beta);
D_4=D(g_4);
sigma_4=E*[-eta -mu eta mu]*D_4/l_4;
beta=-36.87;
eta=cosd(beta); mu=sind(beta);
D_5=D(g_5);
sigma_5=E*[-eta -mu eta mu]*D_5/l_5;
beta=0;
eta=cosd(beta); mu=sind(beta);
D_6=D(g_6);
sigma_6=E*[-eta -mu eta mu]*D_6/l_6;
beta=90;
eta=cosd(beta); mu=sind(beta);
D_7=D(g_7);
sigma_7=E*[-eta -mu eta mu]*D_7/l_7;
beta=0;
eta=cosd(beta); mu=sind(beta);
D_8=D(g_8);
sigma_8=E*[-eta -mu eta mu]*D_8/l_8;
beta=-36.87;
eta=cosd(beta); mu=sind(beta);
D_9=D(g_9);
sigma_9=E*[-eta -mu eta mu]*D_9/l_9;
%% Fuerzas 
N_1=A*sigma_1
N_2=A*sigma_2
N_3=A*sigma_3
N_4=A*sigma_4
N_5=A*sigma_5
N_6=A*sigma_6
N_7=A*sigma_7
N_8=A*sigma_8
N_9=A*sigma_9
%% Plotear Deformacion
XY=zeros(6,2);
XY(1,:)=[0 0];
XY(2,:)=[4 0];
XY(3,:)=[4 3];
XY(4,:)=[8 0];
XY(5,:)=[8 3];
XY(6,:)=[12 0];

XYdef=zeros(size(XY));
fac=500;
c=0;
for i=1:6
c=c+1;
XYdef(i,1)=XY(i,1)+fac*D(c);
c=c+1;
XYdef(i,2)=XY(i,2)+fac*D(c);
end

IJ=zeros(9,2);
IJ(1,:)=[1 2];
IJ(2,:)=[1 3];
IJ(3,:)=[2 3];
IJ(4,:)=[2 4];
IJ(5,:)=[3 4];
IJ(6,:)=[3 5];
IJ(7,:)=[4 5];
IJ(8,:)=[4 6];
IJ(9,:)=[5 6];
figure
for e=1:9
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