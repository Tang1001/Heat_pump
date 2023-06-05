function COP=COP1(T_in)
load("cop_P1.mat")
a1 = cop_P1(1);
a2 = cop_P1(2);
a3 = cop_P1(3);

COP = a1 + a2 * T_in + a3 * T_in.^2;
end