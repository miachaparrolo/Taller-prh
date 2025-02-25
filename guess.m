function g= guess (z,sigma,C_A_front,D_AB,k)
C_A_0=C_A_front (1); C_A_sig=C_A_front (2); 
g=[((C_A_sig-C_A_0)*z/sigma)+C_A_0; (k/D_AB)*(((C_A_sig-C_A_0)*z^2/(2*sigma))+(C_A_0*z))];
end