%% Función VecbVort
% Función VecbVort: Genera el vector de constantes b (lado derecho) de
%                               A * x = b
% para la ecuación de vorticidad.
% ------------------------------------------------------------------------
% Datos de entrada:
%   M: Número de nodos en la dirección "x"
%   N: Número de nodos en la dirección "y"
%   TotNod: Número total de nodos en el dominio
%   Tau: Coeficiente Tau
%   Temp: Estado de tempeeraturas
%   alfa: Coeficiente Alfa
%   p: Valor constante correspondiente a la condicion de frontera en la 
%   función corriente
%   hx: Delta x (diferencia espacial en "x")
%   Flag: Bandera correspondiente a la litología deseada.
% ------------------------------------------------------------------------
% Datos de salida:
%   b: Vector de constantes b (lado derecho) para la ecuación de
%   vorticidad.
% ------------------------------------------------------------------------

function b = VecbVort(M,N,TotNod,Tau,Temp,alfa,p,hx,Flag)

b = zeros(TotNod,1);
b(1:M) = 0; % CF Dirichlet inferior
for j = 2:N-1
    if j == 3 || j == 4 || j == 5 || j == 6 || j == 7
        b((j-1)*M+1) = -p*hx; % CF Neumann izquierda
        for i = (j-1)*M+2:j*M-1
            if i == (j-1)*M+2
                b(i) = -Tau(Flag(i))*(Temp(i+1) - Temp(i-1)) + alfa^2*p*hx;
            else
                b(i) = -Tau(Flag(i))*(Temp(i+1) - Temp(i-1));
            end
        end
        b(j*M) = 0; % CF Dirichlet derecha
    elseif j == N-2 || j == N-3 || j == N-4 || j == N-5 || j == N-6
        b((j-1)*M+1) = 0; % CF Dirichlet izquierda
        for i = (j-1)*M+2:j*M-1
            if i == j*M-1
                b(i) = -Tau(Flag(i))*(Temp(i+1) - Temp(i-1)) - alfa^2*p*hx;
            else
                b(i) = -Tau(Flag(i))*(Temp(i+1) - Temp(i-1));
            end
        end
        b(j*M) = p*hx; % CF Neumann derecha
    else
        b((j-1)*M+1) = 0; % CF Dirichlet izquierda
        for i = (j-1)*M+2:j*M-1
            b(i) = -Tau(Flag(i))*(Temp(i+1) - Temp(i-1));
        end
        b(j*M) = 0; % CF Dirichlet derecha
    end
end
b((N-1)*M+1:TotNod) = 0; % CF Dirichlet superior

end