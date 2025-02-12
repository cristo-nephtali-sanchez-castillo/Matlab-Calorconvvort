%% Función MatAVort
% Función MatAVort: Genera la matriz A del sistema de ecuaciones lineales
%                               A * x = b
% correspondiente a la ecuación de vorticidad
% ------------------------------------------------------------------------
% Datos de entrada:
% Mx    : Número de nodos en la dirección "x"
% My    : Número de nodos en la dirección "y"
% TotNod: Número total de nodos en el modelo
% kappa : Constante Kappa
% alfa  : Constante Alfa
% ------------------------------------------------------------------------
% Datos de salida:
% A     : Matriz de coeficientes A
% ------------------------------------------------------------------------

function A = MatAVort(Mx,My,TotNod,kappa,alfa)

% Construir las diagonales de la matriz A
C = ones(1,TotNod);
E = zeros(1,TotNod-1);
W = E;
N = zeros(1,TotNod-Mx);
S = N;

for j = 2:My-1
    C((j-1)*Mx+1) = 1;
    if j == 3 || j == 4 || j == 5 || j == 6 || j == 7
        C((j-1)*Mx+2) = -(alfa^2+2*kappa);
        C((j-1)*Mx+3:j*Mx-1) = -(2*alfa^2+2*kappa);
        E((j-1)*Mx+1) = -1;
    elseif j == My-2 || j == My-3 || j == My-4 || j == My-5 || j == My-6
        C((j-1)*Mx+2:j*Mx-2) = -(2*alfa^2+2*kappa);
        C(j*Mx-1) = -(alfa^2+2*kappa);
        W(j*Mx-1) = -1;
    else
        C((j-1)*Mx+2:j*Mx-1) = -(2*alfa^2+2*kappa);
    end
    C(j*Mx) = 1;
    E((j-1)*Mx+2:j*Mx-2) = alfa^2;
    W((j-1)*Mx+2:j*Mx-2) = alfa^2;
end

for j = 2:My-2
    N((j-1)*Mx+2:j*Mx-1) = kappa;
end

for j = 2:My-2
    S((j-1)*Mx+2:j*Mx-1) = kappa;
end

% Construir la matriz A
A = zeros(TotNod,TotNod);
A = A + diag(C,0) + diag(W,-1) + diag(E,1) + diag(N,Mx) + diag(S,-Mx);

end