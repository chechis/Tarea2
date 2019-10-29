close all
clear all
clc

%% ingresando datos
% elasticidad = input('Ingrese el modulo de eslasticidad: ');
% area = input('Ingrese el area: ');
% longitud = input('Ingrese la longitud: ');
% altura = input('Ingrese la altura del puente: ');

elasticidad = 2e8;
area = 0.005;
longitud = 4;
altura = 3;

tramos = input('Ingrese la cantidad de tramos: ');

%% grados de libertad globales y cargas
if(tramos ==0)
    n = 4;
elseif(tramos>0)
    n = 4+2*tramos;
end

p=[0 0 0 0 4 -20 0 0 0 -20 0 0]';

% cantidadCargas = input('Ingrese la cantidad de cargas: ');
% p= zeros(1, n*2);
% for i = 1:cantidadCargas
%    
%     carga = input ('Ingrese la carga: ');
%     cargaPos=input ('Ingrese la posicion de la carga en el gdl: ');
%     p(1, cargaPos)= carga;
%     
% end

gdl= zeros(n,2);

for i=1:n
   gdl(i,:) = [2*i-1 2*i]; 
end
%% coordenadas del sistema

xy= zeros(n, 2);

m= tramos+1;

for i=2:(n-m)
    xy(i,1)= longitud*(i-1);
end

for j= (n-m+1):n
    xy(j,1)= (j-(n-m))*longitud;
end

for j= (n-m+1):n
    xy(j,2)= altura;
end

%% angulos de las barras y longitd barra inclinada
angulo1 = 0;
angulo2 = atand(altura/longitud);
angulo3 = 90;
inclinada= sqrt(altura^2+longitud^2);

%% Angulos de las barras de la armadura

if(tramos ==0)
    barras = 5;
elseif(tramos>0)
    barras = 5+4*tramos;
end

ang = zeros(barras, 1);

for i = 1:3
    if(i==1)
        ang(2*i-1,1)=angulo1;
        ang(2*i,1)=angulo2;
        
    elseif(i==2)
        
        if(tramos==0)
            ang(i+1,1)=angulo3;
        elseif(tramos>0)
            for t = 1:tramos
                ang(4*t-1,1)=angulo3;
                ang(4*t,1)=angulo1;
                ang(4*t+1,1)=-1*angulo2;
                ang(4*t+2,1)=angulo1;
            end
        end
        
    elseif(i==3)
        ang(barras-2,1)=angulo3;
        ang(barras-1,1)=angulo1;
        ang(barras,1)=-1*angulo2;
    end
        
end

%% Longitdud de cada barra

if(tramos ==0)
    barras = 5;
elseif(tramos>0)
    barras = 5+4*tramos;
end

g_x = zeros(barras, 4);
longitudes = zeros(barras, 1);

for i = 1:3
    if(i==1)
        longitudes(2*i-1,1)=longitud;
        g_x(2*i-1,:)= [gdl(1,:) gdl(2,:)];
        
        longitudes(2*i,1)=inclinada;
        g_x(2*i,:)= [gdl(1,:) gdl(3,:)];
        
        longitudes(2*i+1,1)=altura;
        g_x(2*i+1,:)= [gdl(2,:) gdl(3,:)];
        
    elseif(i==2)
        
        if(tramos==0)
            continue
        elseif(tramos>0)

            for t = 1:tramos

            longitudes(4*t,1)=longitud;
            g_x(4*t,:)= [gdl(t*2,:) gdl(t*2+2,:)];
            
            longitudes(4*t+1,1)=inclinada;
            g_x(4*t+1,:)= [gdl(t*2+1,:) gdl(t*2+2,:)];
            
            longitudes(4*t+2,1)=longitud;
            g_x(4*t+2,:)= [gdl(t*2+1,:) gdl(t*2+3,:)];
            
            longitudes(4*t+3,1)=altura;
            g_x(4*t+3,:)= [gdl(t*2+2,:) gdl(t*2+3,:)];

            end
        end
        
    elseif(i==3)

        longitudes(barras-1,1)=longitud;
        g_x(barras-1,:)= [gdl(n-2,:) gdl(n,:)];
        
        longitudes(barras,1)=inclinada;
        g_x(barras,:)= [gdl(n-1,:) gdl(n,:)];

    end
        
end

 %% Matriz de rigidez por elemento
 
K=zeros(n*2, n*2);

for i = 1:barras
    
    ang_x=ang(i);

    t_x=t_xx(ang_x);

    k_x=barraAxial(elasticidad,area,longitudes(i,1));
    K_x=t_x'*k_x*t_x;
    
    deltaK_x= zeros(n*2, n*2);
    
    deltaK_x(g_x(i,:), g_x(i,:))=K_x;
    K= K + deltaK_x;
end

%% Grados de libertar libres y restringidos
a=[1 2 gdl(n,:)]; 
b= zeros(1, n-4);
for i = 1:2*n-4
    b(1,i)= i+2;
end

k_aa = K(a,a); k_ab=K(a,b); k_ba=K(b,a); k_bb=K(b,b); 

p_b= p(b);

D_b=k_bb\p_b;
P_a=k_ab*D_b;

D=zeros(n*2,1);
D(b)=D_b; 

%% Esfuerzos internos
sigmas = zeros(barras, 1);
for i= 1:barras

    D_x = D(g_x(i,:));
    sigmas(i,1) = sigma_x(ang(i,1),elasticidad, D_x, longitudes(i,1));
end

%% Fuerzas 
N_x = zeros(barras, 1);
for i = 1:barras
N_x(i,1) = area*sigmas(i,1);
end

%% Plotear deformacion

xydef = zeros(size(xy));
fac = 500;
c=0;

for i =1:6
    c=c+1;
    xydef(i,1)= xy(i,1)+fac*D(c);
    c=c+1;
    xydef()= xy()+fac*D(c);   
end
%% Funciones

 function [sigma_xx]=sigma_x(beta, E, D_x, L)

    eta=cosd(beta);
    mu = sind(beta);
    sigma_xx= E*[-eta -mu eta mu]*D_x/L;
 end

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
 