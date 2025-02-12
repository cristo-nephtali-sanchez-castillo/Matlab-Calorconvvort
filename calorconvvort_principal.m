%% Ecuación de calor convectivo (con función de corriente).
% Ecuación de calor convectivo:
%
% /H\2           H /         \  H     H
% |-| *T  + T  - -*|Ѱ T -Ѱ T |= -*T - -*F
% \L/   xx   yy  L \ y x  x y/  L  t  L
%
% ------------------------------------------------------------------------
% Considerese una lámina de tamaño H*L.
% Condiciones Neumann en x = 0 (x'=0) y x = L (x'=1):
%
%   ∂T                 ∂T
% - -- (0,y,t) = 0 ;   -- (L,y,t) = 0
%   ∂x                 ∂x
% 
% Se supone una lamina aislada lateralmente (flujo de calor cero)
% Condiciones Dirichlet en y = 0 (y'=0) y y = H (y'=1):
%
% T(x,0,t) = (b-a)/(b-a)=1 ; T(x,H,t) = (a-a)/(b-a)=0
% ------------------------------------------------------------------------
% ========================================================================
% Ecuación de vorticidad (función de corriente):
% H     L
% -*Ѱ  +-*Ѱ  =-Ra*T  
% L  xx H  yy      x  
% 
% donde:
% 
%    (ρc)_f*ρ_r*g*β_th*ΔT*H*k
% Ra=------------------------
%               η*λ
% ------------------------------------------------------------------------
% Condiciones Dirichlet en todas las fronteras del dominio Ω:
%
% Ѱ = 0 para todo ∂Ω (siendo ∂Ω la frontera del dominio Ω)
% 
% ------------------------------------------------------------------------

%% Limpiar ventana de comandos

clc, clear, close all

%% Datos del modelo
disp('*******************************************************************')
disp('*MODELADOR DE TRANSPORTE DE MASA Y CALOR CON CONVECCIÓN MEDIANTE  *')
disp('*ECUACIÓN DE VORTICIDAD (FUNCIÓN CORRIENTE)                       *')
disp('*******************************************************************')

disp(' ')
disp('Cargando datos del modelo...')
disp(' ')

M = 80; % número de nodos en la direccion X
N = 60; % número de nodos en la dirección Y
L = 8000; % longitud del dominio en [m]
H = 4000; % altura del dominio en [m]
tfin = 5.3611e+12;%10e5*130000; % tiempo final en [s]
NT = 300; % % número de pasos de tiempo
Ta = 150; % Temperatura en la base del dominio en [°C]
Tb = 15; % Temperatura en la cima del dominio en [°C]
condu = [1.65; 1.865; 2.79; 2.085; 2.68]; % Conductividad térmica del medio poroso en [W/(m.°C)]; cada elemento del vector corresponde al valor para cada litologia.
fdens = 1000; % densidad del fluido en [kg/m^3]
fdensr = 968; % densidad del fluido a 25°C en [kg/m^3]
fcalespec = 4184; % calor especifico del fluido a presión constante en [J/(kg.°C)]
rdens = [2341; 2270.5; 2500; 2650; 2761]; % densidad del medio poroso en [kg/m^3]; cada elemento del vector corresponde al valor para cada litologia.
rcalespec = [840; 840; 1035; 960; 950]; % calor especifico del medio poroso a presión constante en [J/(kg.°C)]; cada elemento del vector corresponde al valor para cada litologia.
poros = [0.12; 0.1; 0.075; 0.055; 0.05]; % porosidad del medio poroso [fracción]; cada elemento del vector corresponde al valor para cada litologia.
perm = [1e-17; 5.05e-16; 5.05e-14; 5e-14; 5.08e-17]; % permeabilidad intrinseca del medio poroso en [m^2]
grav = 9.81; % valor de la gravedad en [m/s^2]
fthexp = 8.275e-4; % coeficiente de expansión volumétrica del fluido en [1/°C]
fvisc = 8.71e-4; % viscosidad del fluido en [Pa*s]
f = 0; % término fuente/sumidero [W/m^3]
% Ra = 80; % Número de Rayleigh
alfa = H/L; % Relación H/L
gam = 0.001; % Coeficiente de variación de la viscosidad con la temperatura
q = -0.3; % Flujo de calor en la frontera [W/m^2]
p = 300; % Valor constante correspondiente a la condicion de frontera en la función corriente, flujo de agua en la frontera [m^3/s] 
start = 1; % Bandera para indicar inicio de simulación (1 indica primera etapa de la simulación, otro valor numérico indica etapas posteriores de la simulación)

disp('Hecho.')
disp(' ')
%% Generar parámetros adimensionales y auxiliares

disp('Generando parámetros adimensionales y auxiliares...')
disp(' ')

[Ladim,hx,Hadim,hy,MT,kappa,rfdensCS,k,csi,dtemp,fu,tau,aux,aux2,Ra,conduprom] = VarAux(L,M,H,N,poros,rdens,rcalespec,fdens,...
    fdensr,fcalespec,condu,perm,grav,fthexp,fvisc,tfin,NT,alfa,Ta,Tb,f,q);

lito=load('Heterogeneidad.mat');

disp('Hecho.')
disp(' ')

%% Implementar condición inicial

disp('Generando condición inicial...')
disp(' ')

if start == 1
    p = 0;
    [x,y,T0,Psi0] = CondInic(Ladim,M,Hadim,N,Ta,Tb,start,conduprom,fdens,fcalespec);
    contour_levels = 25;
else
    [x,y,T0,Psi0] = CondInic(Ladim,M,Hadim,N,Ta,Tb,start,conduprom,fdens,fcalespec);
    contour_levels = 60; % Este valor se puede ajustar a fin de hacer la gráfica "presentable"
end

disp('Hecho.')
disp(' ')

%% Matriz A de Ax=b

disp('Generando matriz A...')
disp(' ')

% Vorticidad

disp('Vorticidad')
disp(' ')

A1 = MatAVort(M,N,MT,kappa,alfa);

disp('OK')

disp('Calor')
disp(' ')

A2 = MatACal(M,N,MT,Psi0,csi,kappa,alfa,aux);

disp('OK')
disp(' ')

%% Vector de coeficientes b de Ax=b

disp('Generando vector b...')
disp(' ')

% Vorticidad


disp('Vorticidad')
disp(' ')

b1 = VecbVort(M,N,MT,tau,T0,alfa,p,hx,lito.Flag);

disp('OK')

disp('Calor')
disp(' ')

b2 = VecbCal(M,N,MT,csi,T0,fu,Psi0,Ta,Tb,kappa,aux,aux2,hy);

disp('OK')
disp(' ')

%% Vectores auxiliares para la solución

TNew = zeros(1,MT);
TI = T0';
TOld = T0;
PsiNew = zeros(1,MT);
PsiI = Psi0';
PsiOld = Psi0;
tiempo = 0;

%%

% Expresar resultados de variables adimensionales a variables dimensionales

[PsiDim,TDim,xdm,ydm,tiempo_dim,ConduRes] = GeneraVariableDim(PsiOld,TOld,dtemp,Tb,...
    x,y,L,H,tiempo,conduprom,rfdensCS,MT,M,N,gam,fdens,fcalespec);

disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
disp(['t = ' num2str(tiempo_dim) ': Graficando estado inicial'])

SolT = reshape(TDim,[M,N]);
MaxTem = max(max(TDim));
SolS = reshape(PsiDim,[M,N]);

subplot(1,2,1)
contourf(xdm,ydm,SolT',20);
view(2);
colorbar('location','eastoutside','fontsize',12);
clim([Tb MaxTem])
xlabel('X [m]','fontSize',12);
ylabel('Y [m]','fontSize',12);
title('Comportamiento Térmico en [°C]','fontsize',12);

subplot(1,2,2)
contourf(xdm,ydm,SolS',25);
view(2);
colorbar('location','eastoutside','fontsize',12);
xlabel('X [m]','fontSize',12);
ylabel('Y [m]','fontSize',12);
title('Función Corriente [m^2/s]','fontsize',12);

fh = figure(1);
set(fh, 'color', 'white');

pause;

%% Resolver sistema lineal y graficar los resultados, para cada paso de tiempo.

disp('*******************************************************************')
disp('Solución de las ecuaciones de calor y vorticidad')

tol = 1e-3;
maxiter = 2000;

for t = 1:NT

    disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ')
    tiempo = tiempo + k;
    disp(['Paso de tiempo ' num2str(t)])
    disp(['Tiempo adimensional = ' num2str(tiempo) ':'])

    [PsiI,TI] = Solver(A1,b1,PsiI,A2,b2,TI,M,N,tol,maxiter,Ta,Tb,aux2,hy,p,hx);

    disp('-------------------------------')

    A2 = MatACal(M,N,MT,PsiI,csi,kappa,alfa,aux);
    b1 = VecbVort(M,N,MT,tau,TI,alfa,p,hx,lito.Flag);
    b2 = VecbCal(M,N,MT,csi,TI,fu,PsiI,Ta,Tb,kappa,aux,aux2,hy);

    [PsiDim,TDim,xdm,ydm,tiempo_dim,ConduRes] = GeneraVariableDim(PsiI,TI,dtemp,Tb,...
    x,y,L,H,tiempo,conduprom,rfdensCS,MT,M,N,gam,fdens,fcalespec);

    disp(['t = ' num2str(tiempo_dim) ' [s] = ' num2str(tiempo_dim/(60*60*24*365)) ' [años]'])
    disp(' ')

    SolT = reshape(TDim,[M,N]);
    MaxTem = max(max(TDim));
    SolS = reshape(PsiDim,[M,N]);

    PsiI = PsiI';
    TI = TI';
    
    subplot(1,2,1)
    contourf(xdm,ydm,SolT',20);
    view(2);
    colorbar('location','eastoutside','fontsize',12);
    clim([Tb MaxTem])
    xlabel('X [m]','fontSize',12);
    ylabel('Y [m]','fontSize',12);
    title('Comportamiento Térmico en [°C]','fontsize',12);
    
    subplot(1,2,2)
    contourf(xdm,ydm,SolS',contour_levels);
    view(2);
    colorbar('location','eastoutside','fontsize',12);
    xlabel('X [m]','fontSize',12);
    ylabel('Y [m]','fontSize',12);
    title('Función Corriente [m^2/s]','fontsize',12);
    
    fh = figure(1);
    set(fh, 'color', 'white');

    disp('------------------------------------------------------------')

    if t == NT/2
        pause;
    end
end

save("SaveIncon.mat","SolS","SolT",'-mat');

disp('SIMULACIÓN COMPLETA!!!!!')
disp('****************************************************************')