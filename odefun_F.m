function F=odefun_F(~,Y,k,D_AB, V, S)
%varibales depentidentes 
C_A=Y(1); u=Y(2);
%Sistema de ecuaciones diferenciales
dC_Adz=u;
dudz=(k*V*C_A)/(-S*D_AB);
%salida
F=[dC_Adz; dudz];
end

