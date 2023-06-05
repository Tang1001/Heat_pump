function COP=COP2(T_in,T_amb)
load("cop_P2.mat")
a1 = cop_P2(1);
a2 = cop_P2(2);
a3 = cop_P2(3);
a4 = cop_P2(4);
a5 = cop_P2(5);

COP = a1 + a2 * T_in + a3 * T_in.^2 + a4 * T_amb + a5 * T_in .* T_amb;

end