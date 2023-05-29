function T_results=system_simulation_Euler_2and4(time_data,sample_time,T_in_interp,m_p_dot_interp,T_buttom_water_tank2_interp,T_upper_water_tank1_interp,T_upper_water_tank2_interp,m_s_dot_interp,params)

%params: 1R  diff_T_c T_s m_c_dot   4m

R12_tank=params(1);
m_c_dot= params(2); % flow rate of hot water circulation    [m^3/h]=1000*[L/h]=1000*[kg/h]
diff_T_c = params(3); % temperature difference  [K]
T_s = params(4); % cold water temperature      [°C]
m3 = params(5);
m4 = params(6);
m5 = params(7);
m6 = params(8);

% Parameters and initial conditions
cp = 4186; % specific heat  [J/kg·K]        cp = 4217 - 0.3 * T
k= 0.6; % the thermal conductivity of water  [W/m·K]
r=0.65/2; % radius of tank  [m]
A= pi*r^2; % cross-sectional area between layers    [m^2]
H=1.64; % height of tank    [m]



% mass layer   L~kg
% height of layer    [m]
m1 = 500/2;
m2 = 500/2;
h1 = H/2;
h2 = H/2;

% tank2
h3=H*(m3/500);
h4=H*(m4/500);
h5=H*(m5/500);
h6=H*(m6/500);


% thermal resistances
%R_i = h_i / (k_i * A)
R1 = h1/(k*A);
R2 = h2/(k*A);

R3=h3/(k*A);
R4=h4/(k*A);
R5=h5/(k*A);
R6=h6/(k*A);

%H/(k*A)
R_tank1_tot = R1+R2;
R_tank2_tot = R3+R4+R5+R6;


% thermal efficiency between layers
%k*A/h_i      [W/K]
R12 = 1/(2 * (R1 * R2) / (R1 + R2));

R34 = 1/(2 * (R3 * R4) / (R3 + R4));
R45 = 1/(2 * (R4 * R5) / (R4 + R5));
R56 = 1/(2 * (R5 * R6) / (R5 + R6));



%% Initialization
T1_initial=T_upper_water_tank1_interp(1);
T3_initial=T_upper_water_tank2_interp(1);
T6_initial=T_buttom_water_tank2_interp(1);

d_T1 = (T1_initial - T3_initial) * (R1 / R_tank1_tot);
T2_initial = T1_initial - d_T1;     % layer 2 initial temperature
d_T3 = (T3_initial - T6_initial) * (R3 / R_tank2_tot);
T4_initial= T3_initial - d_T3; % layer 4 initial temperature
d_T4 = (T3_initial - T6_initial) * (R4 / R_tank2_tot);
T5_initial= T4_initial - d_T4; % layer 4 initial temperature
T_initial = [T1_initial; T2_initial; T3_initial; T4_initial; T5_initial; T6_initial];

% T_results = system_of_equations_Euler_Nlayers_singlePipe(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, num_layer_tank1, num_layer_tank2, m_layer_tank1,m_layer_tank2, R_thermal_tank1, R_thermal_tank2, R12, cp, m_c_dot, diff_T_c, T_s);
T_results = system_of_equations_Euler_6layers_5R_singlePipe(time_data, T_initial, m_p_dot_interp, T_in_interp, m_s_dot_interp, m1, m2, m3, m4, m5, m6, cp, R12, R12_tank, R34, R45, R56, m_c_dot, diff_T_c, T_s);




end
