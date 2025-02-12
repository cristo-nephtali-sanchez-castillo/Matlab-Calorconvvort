%% Función GeneraVariableDim

% Función GeneraVariableDim: Genera los valores en las variables
% dimensionales a partir de valores en las variables adimensionales.
%
% ------------------------------------------------------------------------
% Datos de entrada:
%   Psi      : Estado de función corriente adimensional
%   T        : Estado de temperatura adimensional
%   DTemp    : Diferencia de temperaturas en la cima y base del dominio
%   Ttop     : Temperatura en la cima del dominio
%   x        : Vector de valores en "x" (adimensional)
%   y        : Vector de valores en "y" (adimensional)
%   L        : Longitud del dominio
%   H        : Altura del dominio
%   tiempo   : Tiempo adimensional
%   Conduprom: Valor promedio de la condcutividad térmica del dominio
%   rfdensCS : Cantidad de calor especifico de roca + fluido
%   MT       : Número total de nodos en el modelo
%   M        : Número de nodos en la dirección "x"
%   N        : Número de nodos en la dirección "y"
%   gam      : Coeficiente de variación de la viscosidad con respecto a la
%              temperatura
%   fDens    : Densidad del fluido
%   fCalSpc  : Calor específico del fluido
% ------------------------------------------------------------------------
% Datos de salida:
%   PsiDim    : Estado de funcion corriente dimensional
%   TempDim   : Estado de temperatura dimensional
%   xdm       : Vector de valores en la dimensión "x" dimensional
%   ydm       : Vector de valores en la dimensión "y" dimensional
%   tiempo_dim: tiempo dimensional
%   ConduRes  : Estado de conductividad térmica resultante
% ------------------------------------------------------------------------
% Nota: Para la variable PsiDim, la expresión correspondiente a dicha
% variable viene del proceso de adimensionalización de la ecuación de
% vorticidad en su forma dimensional, en donde dicha expresión se usó para
% transformar en el término derecho de la ecuación de vorticidad en forma
% adimensional, expresada por
%
%                           -Ra*T_x
%
% donde Ra es el número de Rayleigh:
%
%    (ρc)_f*ρ_r*g*β_th*ΔT*H*k
% Ra=------------------------
%               η*λ
%
% y T_x es la derivada parcial de T con respecto a x.

function [PsiDim,TempDim,xdm,ydm,tiempo_dim,ConduRes] = GeneraVariableDim(Psi,T,DTemp,Ttop,x,y,L,H,tiempo,Conduprom,...
    rfdensCS,MT,M,N,gam,fDens,fCalSpc)

PsiDim = zeros(1,MT);
TempDim = zeros(1,MT);
ConduRes = zeros(1,MT);

for i = 1:MT
    PsiDim(i) = Psi(i)*(Conduprom/(fDens*fCalSpc));
    TempDim(i) = T(i)*DTemp + Ttop;
    ConduRes(i) = Conduprom/(1+gam*(T(i)-Ttop));
end

xdm = zeros(1,M);
ydm = zeros(1,N);

for i = 1:M
    xdm(i) = x(i)*L;
end

for i = 1:N
    ydm(i) = y(i)*H;
end

tiempo_dim = (L*H*rfdensCS*tiempo)/Conduprom;

end