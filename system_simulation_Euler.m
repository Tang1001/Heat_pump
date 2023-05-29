%%
% The number of layers:
% num_layer_tank1
% num_layer_tank2

% Temperature interpolated data of inlet hot water
% T_in_data_interp

% Flow rate:
% FR_water_hp

% Temperatures interpolated data of tanks:
% T_buttom_water_tank2_interp
% T_upper_water_tank1_interp
% T_upper_water_tank2_interp

% The consumption data of water:
% m_s_dot_data_interp


function T_results=system_simulation_Euler(time_data,sample_time,num_layer_tank1,num_layer_tank2,T_in_interp,m_p_dot_interp,T_buttom_water_tank2_interp,T_upper_water_tank1_interp,T_upper_water_tank2_interp,m_s_dot_interp,params)
% Parameters and initial conditions
cp = 4186; % specific heat  [J/kg·K]        cp = 4217 - 0.3 * T
k= 0.6; % the thermal conductivity of water  [W/m·K]
r=0.65/2; % radius of tank  [m]
A= pi*r^2; % cross-sectional area between layers    [m^2]
H=1.64; % height of tank    [m]


R12=params(1);
m_c_dot= params(2); % flow rate of hot water circulation    [m^3/h]=1000*[L/h]=1000*[kg/h]
diff_T_c = params(3); % temperature difference  [K]
T_s = params(4); % cold water temperature      [°C]


% mass layer   L~kg
% height of layer    [m]
m_layer_tank1 = 500/num_layer_tank1;
h_layer_tank1 = H/num_layer_tank1;

m_layer_tank2 = 500/num_layer_tank2;
h_layer_tank2 = H/num_layer_tank2;

% thermal resistances
%R_i = h_i / (k_i * A)
R_layer_tank1 = h_layer_tank1/(k*A);
R_layer_tank2 = h_layer_tank2/(k*A);
R_tank1_tot = R_layer_tank1*num_layer_tank1;
R_tank2_tot = R_layer_tank2*num_layer_tank2;


% initial thermal efficiency between layers
%k*A/h_i      [W/K]
R_thermal_tank1 = 1/R_layer_tank1;
R_thermal_tank2 = 1/R_layer_tank2;




%% Initialization
T_initial_layer_tank1 = zeros(1,num_layer_tank1);
T_initial_layer_tank2 = zeros(1,num_layer_tank2);

% measurement values
T_initial_layer_tank1(1)=T_upper_water_tank1_interp(1);
T_initial_layer_tank2(1)=T_upper_water_tank2_interp(1);
T_initial_layer_tank2(end)=T_buttom_water_tank2_interp(1);

d_T_tank1 = (T_initial_layer_tank1(1) - T_initial_layer_tank2(1)) * (R_layer_tank1 / R_tank1_tot);
d_T_tank2 = (T_initial_layer_tank2(1) - T_initial_layer_tank2(end)) * (R_layer_tank2 / R_tank2_tot);

%Tank1
for i=2:num_layer_tank1
    T_initial_layer_tank1(i) = T_initial_layer_tank1(1) - d_T_tank1; 
end

%Tank2
for i=2:num_layer_tank2
    T_initial_layer_tank2(i) = T_initial_layer_tank2(1) - d_T_tank2;
end


T=[];
for i=1:num_layer_tank1
    T=[T;T_initial_layer_tank1(i)];
end
for i=1:num_layer_tank2
    T=[T;T_initial_layer_tank2(i)];
end


T_results = system_of_equations_Euler_Nlayers_singlePipe(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, num_layer_tank1, num_layer_tank2, m_layer_tank1,m_layer_tank2, R_thermal_tank1, R_thermal_tank2, R12, cp, m_c_dot, diff_T_c, T_s);




end
