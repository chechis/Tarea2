
% Como ejemplo, consideremos la armadura metálica mostrada en la figura anterior
%sometida a la acción % de dos cargas verticales y una horizontal. 
%Todas las barras tienen un módulo de elasticidad E = 2 ×10^8kN/m2 
%y un área seccional A = 0.005m2
%%
close all
clear all
clc
%% Datos

E=2e10; A=0.04;

l_x=zeros(65,1);

l_1=3.5;
l_2=4;
l_3=sqrt(3.5^2+4^2);
l_4=3.5;

for i =1:16
    l_x(4*i-3,1)= l_1;
    l_x(4*i+1-3,1)= l_2;
    l_x(4*i+2-3,1)= l_3;
    l_x(4*i+3-3,1)= l_4;
end

l_x(65,1)=4;


%% Encontramos los grados de libertad globales segun su nodo
n=34;
gdl=zeros(n,2);
for i=1:n
gdl(i,:)=[2*i-1 2*i];
end
g_x = zeros(65,4);

g_x(1,:)= [gdl(1,:) gdl(3,:)];     g_x(2,:)= [gdl(1,:) gdl(2,:)];      g_x(3,:)= [gdl(2,:) gdl(3,:)];      g_x(4,:)= [gdl(2,:) gdl(4,:)];
g_x(5,:)= [gdl(3,:) gdl(5,:)];     g_x(6,:)= [gdl(3,:) gdl(4,:)];      g_x(7,:)= [gdl(4,:) gdl(5,:)];      g_x(8,:)= [gdl(4,:) gdl(6,:)];
g_x(9,:)= [gdl(5,:) gdl(7,:)];     g_x(10,:)= [gdl(5,:) gdl(6,:)];     g_x(11,:)= [gdl(6,:) gdl(7,:)];     g_x(12,:)= [gdl(6,:) gdl(8,:)];
g_x(13,:)= [gdl(7,:) gdl(9,:)];    g_x(14,:)= [gdl(7,:) gdl(8,:)];     g_x(15,:)= [gdl(8,:) gdl(9,:)];     g_x(16,:)= [gdl(8,:) gdl(10,:)];
g_x(17,:)= [gdl(9,:) gdl(11,:)];   g_x(18,:)= [gdl(9,:) gdl(10,:)];    g_x(19,:)= [gdl(10,:) gdl(11,:)];   g_x(20,:)= [gdl(10,:) gdl(12,:)];
g_x(21,:)= [gdl(11,:) gdl(13,:)];  g_x(22,:)= [gdl(11,:) gdl(12,:)];   g_x(23,:)= [gdl(12,:) gdl(13,:)];   g_x(24,:)= [gdl(12,:) gdl(14,:)];
g_x(25,:)= [gdl(13,:) gdl(15,:)];  g_x(26,:)= [gdl(13,:) gdl(14,:)];   g_x(27,:)= [gdl(14,:) gdl(15,:)];   g_x(28,:)= [gdl(14,:) gdl(16,:)];
g_x(29,:)= [gdl(15,:) gdl(17,:)];  g_x(30,:)= [gdl(15,:) gdl(16,:)];   g_x(31,:)= [gdl(16,:) gdl(17,:)];   g_x(32,:)= [gdl(16,:) gdl(18,:)];

g_x(33,:)= [gdl(17,:) gdl(19,:)];  g_x(34,:)= [gdl(17,:) gdl(18,:)];   g_x(35,:)= [gdl(17,:) gdl(20,:)];   g_x(36,:)= [gdl(18,:) gdl(20,:)];
g_x(37,:)= [gdl(19,:) gdl(21,:)];  g_x(38,:)= [gdl(19,:) gdl(20,:)];   g_x(39,:)= [gdl(19,:) gdl(22,:)];   g_x(40,:)= [gdl(20,:) gdl(22,:)];
g_x(41,:)= [gdl(21,:) gdl(23,:)];  g_x(42,:)= [gdl(21,:) gdl(22,:)];   g_x(43,:)= [gdl(21,:) gdl(24,:)];   g_x(44,:)= [gdl(22,:) gdl(24,:)];
g_x(45,:)= [gdl(23,:) gdl(25,:)];  g_x(46,:)= [gdl(23,:) gdl(24,:)];   g_x(47,:)= [gdl(23,:) gdl(26,:)];   g_x(48,:)= [gdl(24,:) gdl(26,:)];
g_x(49,:)= [gdl(25,:) gdl(27,:)];  g_x(50,:)= [gdl(25,:) gdl(26,:)];   g_x(51,:)= [gdl(25,:) gdl(28,:)];   g_x(52,:)= [gdl(26,:) gdl(28,:)];
g_x(53,:)= [gdl(27,:) gdl(29,:)];  g_x(54,:)= [gdl(27,:) gdl(28,:)];   g_x(55,:)= [gdl(27,:) gdl(30,:)];   g_x(56,:)= [gdl(28,:) gdl(30,:)];
g_x(57,:)= [gdl(29,:) gdl(31,:)];  g_x(58,:)= [gdl(29,:) gdl(30,:)];   g_x(59,:)= [gdl(29,:) gdl(32,:)];   g_x(60,:)= [gdl(30,:) gdl(32,:)];
g_x(61,:)= [gdl(31,:) gdl(34,:)];  g_x(62,:)= [gdl(31,:) gdl(32,:)];   g_x(63,:)= [gdl(31,:) gdl(34,:)];   g_x(64,:)= [gdl(32,:) gdl(34,:)];
                                   g_x(65,:)= [gdl(33,:) gdl(34,:)];
%% Coordenadas del Sistema
xy= zeros(34,2);

for i = 1:17
   xy(2*i-1,1) = 3.5*i-3.5;
   xy(2*i,1) = 3.5*i-3.5;

   xy(2*i,2) = 4;
   xy(2*i-1,2) = 0;
end

%% Angulos de las barras de la armadura

Ang= zeros(65,1);

for i= 1:8
   Ang(4*i-3,1)= 0;
   Ang(4*i-2,1)= 90;
   Ang(4*i-1,1)= -48.8;
   Ang(4*i,1)= 0;
end
for i= 9:16
   Ang(4*i-3,1)= 0;
   Ang(4*i-2,1)= 90;
   Ang(4*i-1,1)= 48.8;
   Ang(4*i,1)= 0;
end

Ang(65,1)= 90;

 %% Calculamos la matriz de rigidez por elemento
 
 K= zeros (n*2, n*2);
 for i= 1:65
    ang_x=Ang(i);

    t_x=t_xx(ang_x);

    k_x=barraAxial(E,A,l_x(i,1));
    K_x=t_x'*k_x*t_x;
    
    deltaK_x= zeros(n*2, n*2);
    
    deltaK_x(g_x(i,:), g_x(i,:))=K_x;
    K= K + deltaK_x;
 end
 
%% Grados de Libertad Libres y Restringidos
a=[2 65 66]; b=zeros(1, 65);
for i= 3:64
    b(1,1)= 1;    b(1,64)= 67;    b(1,65)= 68;
    if (i == 2)
        continue
    elseif(i == 65)
        continue
    elseif(i == 66)
        continue
    else
    b(1, i-1)=i;
    end
end

K_aa=K(a,a); K_ab=K(a,b); K_ba=K(b,a); K_bb=K(b,b);

condCarga = input ('Condicion de carga 1 o 2: ');
P = zeros(n*2,1);
    if(condCarga==1)
            for i = 1:17
                P(4*i-2,1)= -1500; 
            end
    elseif(condCarga==2)
            for i = 1:34
                P(2*i-1,1)= -250; 
            end
    else
            for i = 1:34
                P(2*i-1,1)= -250; 
            end
    end

P_b=P(b);

D_b=K_bb\P_b;
P_a=K_ab*D_b;

D=zeros(n*2,1);
D(b)=D_b;
%% Esfuerzos Internos
sigmas = zeros(65, 1);
for i= 1:65
    D_x = D(g_x(i,:));
    sigmas(i,1) = sigma_x(Ang(i,1),E, D_x, l_x(i,1));
end

%% Fuerzas 
N_x = zeros(65, 1);
for i = 1:65
N_x(i,1) = A*sigmas(i,1);
end

%% Plotear Deformacion
XY = xy;
xydef = zeros(size(XY));
fac = input('Ingrese el factor, para observar la deformacion: ');
c=0;

for i =1:n
    c=c+1;
    xydef(i,1)= XY(i,1)+fac*D(c);
    c=c+1;
    xydef(i,2)= XY(i,2)+fac*D(c);   
end

IJ= zeros(65,2);

for i= 1:(65-1)/2
    IJ(2*i-1,:)= [i i+1];
    IJ(2*i,:)=[i i+2];
end
IJ(65,:)=[n-1 n];

%% plotear

for e=1:65
Q=[XY(IJ(e,1),1) XY(IJ(e,1),2);...
XY(IJ(e,2),1) XY(IJ(e,2),2)];
Qdef=[xydef(IJ(e,1),1) xydef(IJ(e,1),2);...
xydef(IJ(e,2),1) xydef(IJ(e,2),2)];
plot(Q(:,1),Q(:,2),'--b',Qdef(:,1),Qdef(:,2),'-r')
hold on
end
xlabel('x')
ylabel('y')
axis equal

%% funciones

function [k]=barraAxial(E,A,Le)

k=E*A/Le*[1 0 -1 0
          0 0 0 0
          -1 0 1 0
          0 0 0 0];
end
function [t]=t_xx(ang_x)

    eta=cosd(ang_x);
    mu = sind(ang_x);
    t=[eta mu 0 0; -mu eta 0 0;
        0 0 eta mu; 0 0 -mu eta];
end
 function [sigma_xx]=sigma_x(beta, E, D_x, L)

    eta=cosd(beta);
    mu = sind(beta);
    sigma_xx= E*[-eta -mu eta mu]*D_x/L;
 end