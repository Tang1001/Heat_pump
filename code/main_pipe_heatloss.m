%% the Euler method
clear all
close all

%m cp dT/dt = m_dot * cp * (T1 - T2) = Q_dot 
%[kg] * [J/kg·K] [K/s] = [kg/s] * [J/kg·K] * [K] = [J/s] = [W]
%[kg] * [J/kg·K] [K/h] = [kg/h] * [J/kg·K] * [K] = [J/h] = 1/(60*60) * [J/s] = 1/(60*60) * [W] = 1/(60*60) * [J/s]

%R (T1-T2)
%[W/K] * [K] = [W] = [J/s] = (60*60) * [J/h] 

%% Read data
%1min
Temp_out_water_hp = readcell('../Data/Mon 1 May 2023 outlet water T in heat pump.csv');
Temp_out_water_hp=cell2table(Temp_out_water_hp(2:end,:));
[date_out_hp,T_out_hp]=read_csv(Temp_out_water_hp);

Temp_in_water_hp = readcell('../Data/Mon 1 May 2023 inlet water T in heat pump.csv');
Temp_in_water_hp=cell2table(Temp_in_water_hp(2:end,:));
[date_in_hp,T_in_hp]=read_csv(Temp_in_water_hp);


%% Choose data between 
%% 1 May
%5:19 ~ 7:13
%9:08 ~ 9:59
%11:53 ~ 13:16
%17:02 ~ 23:59

%17:02 ~ 23:59
%17+2/60~23+59/60
% % Time span for the simulation
% t_start = 17+2/60;
% t_end = 23+59/60;     % end time [h]
% sample_time = 2;    % [min]
% time_span = [t_start:sample_time/60:t_end];
% time_data = time_span;
% 
% %1min
% T_out_hp_data=T_out_hp(1+17*60+2:1+23*60+59);
% T_in_hp_data=T_in_hp(1+17*60+2:1+23*60+59);
% 
% 
% time_data_1min=t_start:1/60:t_end;
% 
% 
% % Interpolation
% T_out_hp_interp = interp1(time_data_1min,T_out_hp_data,time_data)';
% T_in_hp_interp = interp1(time_data_1min,T_in_hp_data,time_data)';




%11:53 ~ 13:16
%11+53/60~13+16/60
% % Time span for the simulation
% t_start = 11+53/60;
% t_end = 13+16/60;     % end time [h]
% sample_time = 2;    % [min]
% time_span = [t_start:sample_time/60:t_end];
% time_data = time_span;
% 
% %1min
% T_out_hp_data=T_out_hp(1+11*60+53:1+13*60+16);
% T_in_hp_data=T_in_hp(1+11*60+53:1+13*60+16);
% 
% 
% time_data_1min=t_start:1/60:t_end;
% 
% 
% % Interpolation
% T_out_hp_interp = interp1(time_data_1min,T_out_hp_data,time_data)';
% T_in_hp_interp = interp1(time_data_1min,T_in_hp_data,time_data)';






%9:08 ~ 9:59
%9+8/60~9+59/60
% % Time span for the simulation
% t_start = 9+8/60;
% t_end = 9+59/60;     % end time [h]
% sample_time = 2;    % [min]
% time_span = [t_start:sample_time/60:t_end];
% time_data = time_span;
% 
% %1min
% T_out_hp_data=T_out_hp(1+9*60+8:1+9*60+59);
% T_in_hp_data=T_in_hp(1+9*60+8:1+9*60+59);
% 
% 
% time_data_1min=t_start:1/60:t_end;
% 
% 
% % Interpolation
% T_out_hp_interp = interp1(time_data_1min,T_out_hp_data,time_data)';
% T_in_hp_interp = interp1(time_data_1min,T_in_hp_data,time_data)';




%5:19 ~ 7:13
%5+19/60~7+13/60
% Time span for the simulation
t_start = 5+19/60;
t_end = 7+13/60;     % end time [h]
sample_time = 2;    % [min]
time_span = [t_start:sample_time/60:t_end];
time_data = time_span;

%1min
T_out_hp_data=T_out_hp(1+5*60+19:1+7*60+13);
T_in_hp_data=T_in_hp(1+5*60+19:1+7*60+13);


time_data_1min=t_start:1/60:t_end;


% Interpolation
T_out_hp_interp = interp1(time_data_1min,T_out_hp_data,time_data)';
T_in_hp_interp = interp1(time_data_1min,T_in_hp_data,time_data)';



%%
% Properties
height = 1; % Height of the buffer in meters
radius = 0.1; % Radius of the buffer in meters
cp = 4186; % Specific heat capacity of water in J/kg/K
rho = 1000; % Density of water in kg/m^3
h = 10; % Heat transfer coefficient in W/m^2/K  10-100

% Derived properties
V = pi * radius^2 * height; % Volume of the buffer in cubic meters
A = 2 * pi * radius * height; % Surface area of the buffer in square meters
m_pipe = rho * V; % Mass of the water in kg
R_loss_pipe = h*A;  %[W/K]


%% Euler method
% T_results = system_of_equations_Euler_6layers_5R_singlePipe(time_data, T_initial, m_p_dot_interp, T_in_interp, m_s_dot_interp, m1, m2, m3, m4, m5, m6, cp, R12, R_tank12, R34, R45, R56, m_c_dot, diff_T_c, T_s);

T_ambient = 9*ones(length(time_data),1);

T_results_in_hp = system_of_equations_pipe_heatloss(time_data, T_out_hp_interp(1), T_ambient, m_pipe, R_loss_pipe);
plot_comparison_T_pipe_heatloss(time_data,T_results_in_hp,time_data_1min,T_out_hp_data)



T_results_in_hp = system_of_equations_pipe_heatloss(time_data, T_in_hp_interp(1), T_ambient, m_pipe, R_loss_pipe);
plot_comparison_T_pipe_heatloss(time_data,T_results_in_hp,time_data_1min,T_in_hp_data)





%% Identify optimal parameters
%% R_loss_pipe m_pipe
% pipe of outlet water
P0 = [R_loss_pipe,m_pipe];

% Define the objective function
objective_function = @(params) objFun_Euler_pipe_heatloss(time_data, T_out_hp_interp(1), T_ambient,  params, T_out_hp_interp)
    
% Call fminunc
% Aeq = [0, ones(1,num_layer_tank2), 0,0,0];
% beq = [500];
lb = [0,0];

ub = [1,10];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
options.MaxFunctionEvaluations = 3000;
P_out_hp_optimal = fmincon(objective_function, P0, [], [], [], [], lb, ub, [], options)
% save("R_loss.mat", "R_loss_optimal")

R_out_hp_loss=P_out_hp_optimal(1);
m_out_hp_pipe=P_out_hp_optimal(2);

T_results_out_hp = system_of_equations_pipe_heatloss(time_data, T_out_hp_interp(1), T_ambient, m_out_hp_pipe, R_out_hp_loss);
plot_comparison_T_pipe_heatloss(time_data,T_results_out_hp,time_data_1min,T_out_hp_data)

%%
% pipe of inlet water

% Define the objective function
objective_function = @(params) objFun_Euler_pipe_heatloss(time_data, T_in_hp_interp(1), T_ambient,  params, T_in_hp_interp)
    
% Call fminunc
% Aeq = [0, ones(1,num_layer_tank2), 0,0,0];
% beq = [500];
lb = [0,0];
ub = [1,10];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
options.MaxFunctionEvaluations = 3000;
P_in_hp_optimal = fmincon(objective_function, P0, [], [], [], [], lb, ub, [], options)
% save("R_loss.mat", "R_loss_optimal")

R_in_hp_loss=P_in_hp_optimal(1);
m_in_hp_pipe=P_in_hp_optimal(2);

T_results_in_hp = system_of_equations_pipe_heatloss(time_data, T_in_hp_interp(1), T_ambient, m_in_hp_pipe, R_in_hp_loss);
plot_comparison_T_pipe_heatloss(time_data,T_results_in_hp,time_data_1min,T_in_hp_data)


 
