%% Función CondInic
%
% Función CondInic: Genera el estado inicial para el modelo, ya sea desde
% cero o a partir de un estado generado de una simulación previa 
% ------------------------------------------------------------------------
% Datos de entrada:
%   Ladim    : Longitud adimensional
%   M        : Número de nodos en la dirección x
%   Hadim    : Altura adimensional
%   N        : Número de nodos en la dirección y
%   TempBot  : Temperatura en la base del dominio
%   TempTop  : Temperatura en la cima del dominio
%   start    : Bandera de inicio
%   Conduprom: Comductividad térmica promedio (en el caso de un modelo
%              heterogéneo)
%   fDens    : Densidad del fluido
%   fCalSpc  : Calor específico del fluido
% ------------------------------------------------------------------------
% Datos de salida:
%
%   x        : vector de valores de "x"
%   y        : vector de valores de "y"
%   Temp0    : Estado de temperatura inicial (como vector)
%   Psi0     : Estado de función corriente inicial (como vector)
% ------------------------------------------------------------------------

function [x,y,Temp0,Psi0] = CondInic(Ladim,M,Hadim,N,TempBot,TempTop,start,Conduprom,fDens,fCalSpc)

x = linspace(0,Ladim,M);
y = linspace(0,Hadim,N);

% ====== TEMPERATURA =====================================================

if start == 1
    CIT = sin(pi*x/Ladim).'*sin(pi*y/Ladim);
    CIT = CIT';
    CIT(1,:) = (TempBot-TempTop)/(TempBot-TempTop);
    CIT(N,:) = (TempTop-TempTop)/(TempBot-TempTop);
    for j = 2:N-1
        CIT(j,1) = CIT(j,2);
        CIT(j,M) = CIT(j,M-1);
    end
    Temp0 = reshape(CIT',[1,M*N]);
else
    EstInicial = load("SaveIncon.mat");
    Tdim0 = (EstInicial.SolT)';
    CIT(1,:) = (Tdim0(1,:)-TempTop)/(TempBot-TempTop);
    CIT(N,:) = (Tdim0(N,:)-TempTop)/(TempBot-TempTop);
    for j = 2:N-1
        CIT(j,:) = (Tdim0(j,:)-TempTop)/(TempBot-TempTop);
    end
    Temp0 = reshape(CIT',[1,M*N]);
end

% ====== VORTICIDAD ======================================================

if start == 1
    CIST = sin(pi*x/Ladim).'*sin(pi*y/Ladim);
    CIST = CIST';
    CIST(1,:) = 0;
    CIST(N,:) = 0;
    for j = 2:N-1
        CIST(j,1) = 0;
        CIST(j,M) = 0;
    end
    Psi0 = reshape(CIST',[1,M*N]);
else
    alpha_f = Conduprom/(fDens*fCalSpc);
    Psidim0 = (EstInicial.SolS)';
    CIST(1,:) = Psidim0(1,:)/alpha_f;
    CIST(N,:) = Psidim0(N,:)/alpha_f;
    for j = 2:N-1
        CIST(j,:) = Psidim0(j,:)/alpha_f;
    end
    Psi0 = reshape(CIST',[1,M*N]);
end

end