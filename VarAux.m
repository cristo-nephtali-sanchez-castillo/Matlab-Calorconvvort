%% Función VarAux
% Función VarAux: Genera las variables auxiliares y valores adimensionales
% necesarias para la ejecución de la simulación del modelo
% ------------------------------------------------------------------------
% Datos de entrada:
% L: Longitud del dominio
% M: Número de nodos en la direeción "x"
% H: Altura del dominio
% N: Número de nodos en la dirección "y"
% poros: Porosidad
% rDens: Densidad de la matriz de roca
% rCalSpc: Calor específico de la matriz de roca
% fDens: Densidad del fluido
% fDensr: Densidad del fluido a condiciones de referencia
% fCalSpc: Calor específico del fluido
% Condu: Conductividad térmica del medio
% Perm: Permeabilidad
% Grav: Valor de la aceleración por gravedad
% fthexp: Expansión térmica del fluido
% fvisc: Viscosidad del fluido
% tfin: Tiempo final de simulación
% PasTemp: Número de pasos de tiempo
% alfa: Coeficiente Alfa
% TempBot: Temperatura en la base del dominio
% TempTop: Temperatura en la cima del dominio
% f: Valor del término fuente sumidero
% q: Valor del flujo de calor en la frontera
% ------------------------------------------------------------------------
% Datos de salida:
% Ladim: Longitud adimensional
% DeltaEspx: Delta x (diferencia espacial en "x")
% Hadim: Altura adimensional
% DeltaEspy: Delta x (diferencia espacial en "y")
% TotNod: Número total de nodos en el modelo
% kappa: Coeficiente Kappa
% rfdensCS: densidad de la roca + fluido
% DTiempo: diferencia temporal
% Csi: Constante Csi
% DeltaTemper: Diferencia de temperaturas en la cima y base del dominio
% Fuen: Término fuente/sumidero adimensional
% Tau: Coeficiente Tau
% aux y aux2: Coeficientes auxiliares
% Ra: Número de Rayleigh
% Conduprom: Valor promedio de la conductividad térmica en el medio
% poroso.
% ------------------------------------------------------------------------

function [Ladim,DeltaEspx,Hadim,DeltaEspy,TotNod,kappa,rfdensCS,DTiemp,Csi,DeltaTemper,Fuen,Tau,aux,aux2,Ra,Conduprom] = ...
    VarAux(L,M,H,N,poros,rDens,rCalSpc,fDens,fDensr,fCalSpc,Condu,Perm,Grav,fthexp,fvisc,tfin,PasTemp,alfa,TempBot,TempTop,f,q)

Ladim = L/L;
DeltaEspx = Ladim/M;
Hadim = H/H;
DeltaEspy = Hadim/N;
TotNod = M*N;
kappa = (DeltaEspx)^2/((DeltaEspy)^2);
rfdensCS1 = (1-poros(1))*(rDens(1)*rCalSpc(1)) + poros(1)*(fDens*fCalSpc);
rfdensCS2 = (1-poros(2))*(rDens(2)*rCalSpc(2)) + poros(2)*(fDens*fCalSpc);
rfdensCS3 = (1-poros(3))*(rDens(3)*rCalSpc(3)) + poros(3)*(fDens*fCalSpc);
rfdensCS4 = (1-poros(4))*(rDens(4)*rCalSpc(4)) + poros(4)*(fDens*fCalSpc);
rfdensCS5 = (1-poros(5))*(rDens(5)*rCalSpc(5)) + poros(5)*(fDens*fCalSpc);
rfdensCS = (rfdensCS1 + rfdensCS2 + rfdensCS3 + rfdensCS4 + rfdensCS5)/5;
Conduprom = sum(Condu)/5;
tfadim = (Conduprom*tfin)/(L*H*rfdensCS);
DTiemp = tfadim/PasTemp;
DifTerm = Condu./(fDens*fCalSpc);
Csi = (alfa*(DeltaEspx)^2)/DTiemp;
DeltaTemper = TempBot - TempTop;
FuenAux = (f*H*L)/(Conduprom*DeltaTemper);
Fuen = alfa*FuenAux*(DeltaEspx)^2;
Ra = (fDens*fCalSpc*fDensr*Grav*fthexp*DeltaTemper*(H/5)*Perm)./(fvisc*Condu);
Tau = (alfa*Ra*DeltaEspx)/2;
aux = (alfa*DeltaEspx)/DeltaEspy;
%aux2= (H*q)/(Conduprom*DeltaTemper);
aux2 = q;

end