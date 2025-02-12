%% Función Solver
% Función Solver: Se encarga de resolver de manera numérica los sistemas de
% ecuaciones lineales correspondientes tanto para la ecuación de calor como
% para la ecuación de vorticidad, mediante el método de Gauss-Seidel
% -------------------------------------------------------------------------
% Datos de entrada:
% A1: Matriz de coeficientes A para el sistema de ecuaciones (vorticidad)
% B1: Vector de constantes (lado derecho de A*x=b), para la ecuación de
% vorticidad
% P1: Aproximación inicial (vorticidad)
% A2: Matriz de coeficientes A para el sistema de ecuaciones (calor)
% B1: Vector de constantes (lado derecho de A*x=b), para la ecuación de
% calor
% P2: Aproximación inicial (calor)
% Mx: Número de nodos en la dirección "x"
% My: Número de nodos en la dirección "y"
% delta: Tolerancia establecida
% max1: Número máximo de iteraciones establecida
% TempBot: Temperatura en la base del dominio
% TempTop: Temperatura en la cima del dominio
% aux2: Valor del flujo de calor en la fuente regional del dominio
% hy: Delta y (diferencia espacial en "y")
% p: Valor constante correspondiente a la condicion de frontera en la 
% función corriente
% hx: Delta x (diferencia espacial en "x")
% -------------------------------------------------------------------------
% Datos de salida:
% X1: Vector solución para el estado, a tiempo actual, de la función
% corriente
% X2: Vector solución para el estado, a tiempo actual, de la temperatura.
% -------------------------------------------------------------------------

function [X1,X2] = Solver(A1,B1,P1,A2,B2,P2,Mx,My,delta,max1,TempBot,TempTop,aux2,hy,p,hx)

% Inicializar error en 100%
err1 = 1;
err2 = 1;

% Se resolverá por el método de Gauss-Seidel.

N=length(B1);

for iter = 1:max1
    
    % Vorticidad
    for j=1:Mx
        if j==1
            X1(1)=(B1(1)-A1(1,[2,Mx+1])*P1([2,Mx+1]))/A1(1,1);
        else 
        %             X contiene la aproximacion k-esima
        %             y P la (k-1)-esima
            X1(j)=(B1(j)-A1(j,j-1)*X1(j-1)'-A1(j,[j+1,j+Mx])*P1([j+1,j+Mx]))/A1(j,j);
        end
    end

    for j=Mx+1:(My-1)*Mx
        X1(j) = (B1(j)-A1(j,[j-Mx,j-1])*X1([j-Mx,j-1])'-A1(j,[j+1,j+Mx])*P1([j+1,j+Mx]))/A1(j,j);
    end

    for j=(My-1)*Mx:N
        if j==N
            X1(N)=(B1(N)-A1(N,[N-Mx,N-1])*X1([N-Mx,N-1])')/A1(N,N);
        else
            X1(j)=(B1(j)-A1(j,[j-Mx,j-1])*X1([j-Mx,j-1])'-A1(j,j+1)*P1(j+1))/A1(j,j);
        end
    end

    err1=abs(norm(X1'-P1));
    relerr1=err1/(norm(X1)+eps);
    P1=X1';

    % Calor
    for j=1:Mx
        if j==1
            X2(1)=(B2(1)-A2(1,[2,Mx+1])*P2([2,Mx+1]))/A2(1,1);
        else 
        %             X contiene la aproximacion k-esima
        %             y P la (k-1)-esima
            X2(j)=(B2(j)-A2(j,j-1)*X2(j-1)'-A2(j,[j+1,j+Mx])*P2([j+1,j+Mx]))/A2(j,j);
        end
    end

    for j=Mx+1:(My-1)*Mx
        X2(j) = (B2(j)-A2(j,[j-Mx,j-1])*X2([j-Mx,j-1])'-A2(j,[j+1,j+Mx])*P2([j+1,j+Mx]))/A2(j,j);
    end

    for j=(My-1)*Mx:N
        if j==N
            X2(N)=(B2(N)-A2(N,[N-Mx,N-1])*X2([N-Mx,N-1])')/A2(N,N);
        else
            X2(j)=(B2(j)-A2(j,[j-Mx,j-1])*X2([j-Mx,j-1])'-A2(j,j+1)*P2(j+1))/A2(j,j);
        end
    end
    err2=abs(norm(X2'-P2));
    relerr2=err2/(norm(X2)+eps);
    P2=X2';

    if(err1<delta && err2<delta)||(relerr1<delta && relerr2<delta)
        disp(['Iteraciones alcanzadas: ' num2str(iter)])
        disp(['Error (Vorticidad): ' num2str(err1)])
        disp(['Error (Calor): ' num2str(err2)])
        break
    end
    
end

PsiI=X1';
TempI=X2'; 

%Asegurarse de que se cumplan la CF Dirichlet superior y la CF Neumann inferior;
tbadim = (TempBot-TempTop)/(TempBot-TempTop);
ttadim = (TempTop-TempTop)/(TempBot-TempTop);
TempI(1:Mx) = TempI(Mx+1:2*Mx) - aux2*hy;
TempI((My-1)*Mx+1:N) = ttadim;
%Asegurarse de que se cumplan las CF Neumann laterales
for j = 2:My-1
    TempI((j-1)*Mx+1) = TempI((j-1)*Mx+2);
    TempI(j*Mx) = TempI(j*Mx-1);
    if j == 3 || j == 4 || j == 5 || j == 6 || j == 7
        PsiI((j-1)*Mx+1) = PsiI((j-1)*Mx+2)-p*hx;
    elseif j == My-6 || j == My-5 || j == My-4 || j == My-3 || j == My-2
        PsiI(j*Mx) = PsiI(j*Mx-1)+p*hx;
    end
end

end