%% Función VecbCal
% Función VecbCal: Genera el vector de constantes (lado derecho) de
%                                   A * x = b
% para la ecuación de calor
% ------------------------------------------------------------------------
% Datos de entrada:
%   M: Número de nodos en la dirección "x"
%   N: Número de nosod en la dirección "y"
%   TotNod: Número total de nodos en el dominio
%   Csi: Coeficiente Csi
%   Temp: Estado de temperaturas
%   F: Término fuente sumidero adimensional
%   Psi: Estado de función corriente
%   TempBot: Temperatura en la base del dominio
%   TempTop: Temperatura en la cima del dominio
%   kappa: Coeficiente Kappa
%   aux y aux2: Coeficientes auxiliares
%   hy: Delta y (diferencia espacial en "y")
% ------------------------------------------------------------------------
% Datos de salida:
%   b: Vector b del sistema de ecuaciones lineales (calor)
% ------------------------------------------------------------------------

function b = VecbCal(M,N,TotNod,Csi,Temp,F,Psi,TempBot,TempTop,kappa,aux,aux2,hy)

b = zeros(TotNod,1);
Tbotad = (TempBot-TempTop)/(TempBot-TempTop);
Ttopad = (TempTop-TempTop)/(TempBot-TempTop);

% Aquí se implementa la CF Neumann inferior
b(1:M) = -aux2*hy;

b(M+1) = 0; % CF Neuman izquierda
for i = M+2:2*M-1
    Epsx = Psi(i+1) - Psi(i-1);
    b(i) =-(Csi*Temp(i) + F - (kappa-aux*Epsx/4)*aux2*hy);
end
b(2*M) = 0; % CF Neuman derecha

for j = 3:N-2
    b((j-1)*M+1) = 0; % CF Neuman izquierda
    for i = (j-1)*M+2:j*M-1
        b(i) = -(Csi*Temp(i) + F);
    end
    b(j*M) = 0; % CF Neuman derecha
end

b((N-2)*M+1) = 0; % CF Neuman izquierda
for i = (N-2)*M+2:(N-1)*M-1
    Epsx = Psi(i+1) - Psi(i-1);
    b(i) = -(Csi*Temp(i) + F + (kappa+aux*Epsx/4)*Ttopad);
end
b((N-1)*M) = 0; % CF Neuman derecha

% Aquí se implementa la CF Dirichlet superior
b((N-1)*M+1:TotNod) = Ttopad;

end