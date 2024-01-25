function C_I_e = CalcSagnacMatrix(range)
    Define_Constants;
    rep = omega_ie * range/c;
    C_I_e = [1,rep,0;-rep,1,0;0,0,1];
end