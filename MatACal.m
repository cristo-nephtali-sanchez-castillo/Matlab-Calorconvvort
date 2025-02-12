%% Función MatACal
% Función MatACal: Genera la matriz A del sistema de ecuaciones lineales
%                           A * x = b
% para la ecuación de calor.
% ------------------------------------------------------------------------
% Datos de entrada:
% Mx    : Número de nodos en la dirección "x"
% My    : Número de nodos en la dirección "y"
% TotNod: Número total de nodos en el modelo
% Psi   : Estado de función corriente adimensional
% Csi   : Constante Csi
% kappa : Constante Kappa
% alfa  : Constante Alfa
% aux   : Constante auxiliar
% ------------------------------------------------------------------------
% Datos de salida:
% A     : Matriz de coeficientes A
% ------------------------------------------------------------------------

function A = MatACal(Mx,My,TotNod,Psi,Csi,kappa,alfa,aux)

% Construir las diagonales de la matriz A
C = ones(1,TotNod);
E = zeros(1,TotNod-1);
W = E;
N = zeros(1,TotNod-Mx);
S = N;

Epsx = zeros(1,TotNod);
Epsy = zeros(1,TotNod);

for I = Mx+1:(My-1)*Mx
    Epsx(I) = Psi(I+1) - Psi(I-1);
    Epsy(I) = Psi(I+Mx) - Psi(I-Mx);
end

for j = 2:My-1
    C((j-1)*Mx+1) = 1;
    if j == 2
        C((j-1)*Mx+2) = -(alfa^2 + kappa + Csi - aux*Epsy((j-1)*Mx+2)/4 + aux*Epsx((j-1)*Mx+2)/4);
        C((j-1)*Mx+3:j*Mx-2) = -(2*alfa^2 + kappa + Csi + aux*Epsx((j-1)*Mx+2)/4);
        C(j*Mx-1) = -(alfa^2 + kappa + Csi + aux*Epsy(j*Mx-1)/4 + aux*Epsx((j-1)*Mx+2)/4);
    else
        C((j-1)*Mx+2) = -(alfa^2 + 2*kappa + Csi - aux*Epsy((j-1)*Mx+2)/4);
        C((j-1)*Mx+3:j*Mx-2) = -(2*alfa^2 + 2*kappa + Csi);
        C(j*Mx-1) = -(alfa^2 +2*kappa + Csi + aux*Epsy(j*Mx-1)/4);
    end
    C(j*Mx) = 1;
    E((j-1)*Mx+1) = -1;
    E((j-1)*Mx+2:j*Mx-2) = alfa^2 - aux*Epsy((j-1)*Mx+2:j*Mx-2)/4;
    W((j-1)*Mx+2:j*Mx-2) = alfa^2 + aux*Epsy((j-1)*Mx+3:j*Mx-1)/4;
    W(j*Mx-1) = -1;
end

N(1:Mx)=-1;

for j = 2:My-2
    N((j-1)*Mx+2:j*Mx-1) = kappa + aux*Epsx((j-1)*Mx+2:j*Mx-1)/4;
end

for j = 2:My-2
    S((j-1)*Mx+2:j*Mx-1) = kappa - aux*Epsx(j*Mx+2:(j+1)*Mx-1)/4;
end

% Construir la matriz A
A = zeros(TotNod,TotNod);
A = A + diag(C,0) + diag(W,-1) + diag(E,1) + diag(N,Mx) + diag(S,-Mx);

end