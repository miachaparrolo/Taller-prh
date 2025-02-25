function res = bcfun(ya,yb,C_A_front)
    C_A_0=C_A_front (1); C_A_sig=C_A_front (2); 
    res=[ya(1)-C_A_0; yb(1)-C_A_sig];
end