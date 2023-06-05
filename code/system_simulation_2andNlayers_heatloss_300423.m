%% the Euler method
clear all
close all

%m cp dT/dt = m_dot * cp * (T1 - T2) = Q_dot 
%[kg] * [J/kg·K] [K/s] = [kg/s] * [J/kg·K] * [K] = [J/s] = [W]
%[kg] * [J/kg·K] [K/h] = [kg/h] * [J/kg·K] * [K] = [J/h] = 1/(60*60) * [J/s] = 1/(60*60) * [W] = 1/(60*60) * [J/s]

%R (T1-T2)
%[W/K] * [K] = [W] = [J/s] = (60*60) * [J/h] 

% Parameters and initial conditions
% cp = 4186; % specific heat  [J/kg·K]        cp = 4217 - 0.3 * T
% k= 0.6; % the thermal conductivity of water  [W/m·K]
% r=0.65/2; % radius of tank  [m]
% A= pi*r^2; % cross-sectional area between layers    [m^2]
% H=1.64; % height of tank    [m]

% set the number of layers
% num_layer_tank1=2;
num_layer_tank1=2;
% num_layer_tank2=6;
% num_layer_tank2=5;
num_layer_tank2=4;

% initial masses of layers   L~kg
m_layer_tank1=[250;250];
% m_layer_tank1=[200;200;100];
m_layer_tank2=[100;150;150;100];
% m_layer_tank2=[100;100;100;100;100];
% m_layer_tank2=[100;100;100;100;50;50];





R_tank12=0;
R_loss=0;
m_c_dot= 1*1000; % flow rate of hot water circulation    [m^3/h]=1000*[L/h]=1000*[kg/h]
diff_T_c = 2; % temperature difference  [K]
T_s = 12; % cold water temperature      [°C]


%% Read data
%1min
Temp_out_water_he = readcell('../Data/Sun 30 April 2023 outlet water T in heat exchanger.csv');
Temp_out_water_he=cell2table(Temp_out_water_he(2:end,:));
[date_out_he,T_out_he]=read_csv(Temp_out_water_he);

FlowRate_water_hp = readcell('../Data/Sun 30 April 2023 water flow rate in heat pump.csv');
FlowRate_water_hp=cell2table(FlowRate_water_hp(2:end,:));
[date_FR_hp,FR_hp]=read_csv(FlowRate_water_hp);

Temp_bottom_water_tank2 = readcell('../Data/Sun 30 April 2023 water T in bottom 2 tank.csv');
Temp_bottom_water_tank2=cell2table(Temp_bottom_water_tank2(2:end,:));
[date_bottom_tank2,T_bottom_tank2]=read_csv(Temp_bottom_water_tank2);

Temp_upper_water_tank1 = readcell('../Data/Sun 30 April 2023 water T in upper 1 tank.csv');
Temp_upper_water_tank1=cell2table(Temp_upper_water_tank1(2:end,:));
[date_upper_tank1,T_upper_tank1]=read_csv(Temp_upper_water_tank1);


%8min
Temp_upper_water_tank2 = readcell('../Data/Sun 30 April 2023 water T in upper 2 tank.csv');
Temp_upper_water_tank2=cell2table(Temp_upper_water_tank2(2:end,:));
[date_upper_tank2,T_upper_tank2]=read_csv(Temp_upper_water_tank2);

%1h
Consp_water = readcell('../Data/Sun 30 April 2023 water V consumption.csv');
Consp_water=cell2table(Consp_water(2:end,:));
[date_V_water,V_water]=read_csv(Consp_water);  %[m^3]




%% Choose data between 
%% 30 Apr
%3:55 ~ 4:33
%6:24 ~ 6:59
%11:52 ~ 12:32


%3:55 ~ 4:33
%3+55/60~4+33/60
% % Time span for the simulation
% t_start = 3+55/60;
% t_end = 4+33/60;     % end time [h]
% sample_time = 2;    % [min]
% time_span = [t_start:sample_time/60:t_end];
% time_data = time_span;
% 
% %1min
% m_p_dot_data=FR_hp(1+60*3+55:1+60*4+33)*1000;  %[L/h] 60*7+1+5~60*(7+11.5)+1
% T_in_data=T_out_he(1+60*3+55:1+60*4+33);
% T_upper_tank1_data = T_upper_tank1(1+60*3+55:1+60*4+33);
% T_bottom_tank2_data = T_bottom_tank2(1+60*3+55:1+60*4+33);
% T_upper_tank2_data = T_upper_tank2(1+60*3+55:1+60*4+33); %17:04~18:32
% 
% %1h
% m_s_dot_data=V_water(1+3:1+5)*1000;    %[L/h] 0:00~1:00
% 
% time_data_1min=t_start:1/60:t_end;
% time_data_1h=floor(time_data(1)):ceil(time_data(end));
% % time_data_8min=17+4/60:8/60:18+32/60;



%6:24 ~ 6:59
%6+24/60~6+59/60
% Time span for the simulation
t_start = 6+24/60;
t_end = 6+59/60;     % end time [h]
sample_time = 2;    % [min]
time_span = [t_start:sample_time/60:t_end];
time_data = time_span;

%1min
m_p_dot_data=FR_hp(1+60*6+24:1+60*6+59)*1000;  %[L/h] 60*7+1+5~60*(7+11.5)+1
T_in_data=T_out_he(1+60*6+24:1+60*6+59);
T_upper_tank1_data = T_upper_tank1(1+60*6+24:1+60*6+59);
T_bottom_tank2_data = T_bottom_tank2(1+60*6+24:1+60*6+59);
T_upper_tank2_data = T_upper_tank2(1+60*6+24:1+60*6+59); %17:04~18:32

%1h
m_s_dot_data=V_water(1+6:1+7)*1000;    %[L/h] 6:00~7:00

time_data_1min=t_start:1/60:t_end;
time_data_1h=floor(time_data(1)):ceil(time_data(end));
% time_data_8min=17+4/60:8/60:18+32/60;


%11:52 ~ 12:32
%11+52/60~12+32/60
% Time span for the simulation
% t_start = 11+52/60;
% t_end = 12+32/60;     % end time [h]
% sample_time = 2;    % [min]
% time_span = [t_start:sample_time/60:t_end];
% time_data = time_span;
% 
% %1min
% m_p_dot_data=FR_hp(1+60*11+52:1+60*12+32)*1000;  %[L/h] 60*7+1+5~60*(7+11.5)+1
% T_in_data=T_out_he(1+60*11+52:1+60*12+32);
% T_upper_tank1_data = T_upper_tank1(1+60*11+52:1+60*12+32);
% T_bottom_tank2_data = T_bottom_tank2(1+60*11+52:1+60*12+32);
% T_upper_tank2_data = T_upper_tank2(1+60*11+52:1+60*12+32); %17:04~18:32
% 
% %1h
% m_s_dot_data=V_water(1+11:1+13)*1000;    %[L/h] 11:00~13:00
% 
% time_data_1min=t_start:1/60:t_end;
% time_data_1h=floor(time_data(1)):ceil(time_data(end));
% % time_data_8min=17+4/60:8/60:18+32/60;

%% Interpolation
[m_p_dot_interp,T_in_interp,m_s_dot_interp,T_upper_tank1_interp,T_upper_tank2_interp,...
    T_bottom_tank2_interp]=interpolation_data_heatloss(time_data,sample_time,time_data_1min,time_data_1h,time_data_1min,m_p_dot_data,T_in_data,m_s_dot_data,T_upper_tank1_data,T_upper_tank2_data,T_bottom_tank2_data);


% Initial Temperature Setting

T_initial=initial_T_setting(num_layer_tank1,num_layer_tank2,m_layer_tank1,m_layer_tank2,T_upper_tank1_interp(1),T_upper_tank2_interp(1),T_bottom_tank2_interp(1));


% Define ground true data
T_ground_truth = [T_upper_tank1_interp',T_upper_tank2_interp',T_bottom_tank2_interp']; % Your historical data here






%% Euler method
% T_results = system_of_equations_Euler_6layers_5R_singlePipe(time_data, T_initial, m_p_dot_interp, T_in_interp, m_s_dot_interp, m1, m2, m3, m4, m5, m6, cp, R12, R_tank12, R34, R45, R56, m_c_dot, diff_T_c, T_s);

T_outdoor = 18.5*ones(length(time_data),1);


T_results = system_of_equations_Euler_Nlayers_heatloss(time_data, T_initial, T_outdoor, m_s_dot_interp, num_layer_tank1,num_layer_tank2, m_layer_tank1, m_layer_tank2, R_tank12,0, m_c_dot, diff_T_c, T_s);
% T_results = system_of_equations_Euler_Nlayers_heatloss(time_data, T_initial, T_outdoor, m_s_dot_interp, num_layer_tank1,num_layer_tank2, m_layer_tank1, m_layer_tank2, R_tank12,R_loss, m_c_dot, diff_T_c, T_s);
plot_comparison_T_heatloss(time_data,T_results,time_data_1min,time_data_1min,num_layer_tank1,num_layer_tank2,T_upper_tank1_data,T_upper_tank2_data,T_bottom_tank2_data)

%%
load("optimal_params.mat")
R_tank12_optimal = P_optimal(1);
m_layer_tank2_optimal=zeros(num_layer_tank2,1);
for i=1:num_layer_tank2
    m_layer_tank2_optimal(i) = P_optimal(i+1);
end

diff_T_c_optimal = P_optimal((1+num_layer_tank2)+1);
T_s_optimal = P_optimal((1+num_layer_tank2)+2);
m_c_dot_optimal = P_optimal((1+num_layer_tank2)+3);

T_results = system_of_equations_Euler_Nlayers_heatloss(time_data, T_initial, T_outdoor, m_s_dot_interp, num_layer_tank1,num_layer_tank2, m_layer_tank1, m_layer_tank2_optimal, R_tank12_optimal,R_loss, m_c_dot_optimal, diff_T_c_optimal, T_s_optimal);
plot_comparison_T_heatloss(time_data,T_results,time_data_1min,time_data_1min,num_layer_tank1,num_layer_tank2,T_upper_tank1_data,T_upper_tank2_data,T_bottom_tank2_data)



%% Identify optimal parameters
%% R_loss
R0 = [R_loss];
% R0 = [R_loss];
P0 = [R0];


% Define the objective function
objective_function = @(params) objFun_Euler_2andNlayers_heatloss(time_data, T_initial, T_outdoor, m_s_dot_interp,num_layer_tank2, m_layer_tank1,m_layer_tank2_optimal, R_tank12_optimal, m_c_dot_optimal, diff_T_c_optimal, T_s_optimal, params, T_ground_truth);
% Call fminunc
% Aeq = [0, ones(1,num_layer_tank2), 0,0,0];
% beq = [500];
lb = [1];
ub = [100000];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
options.MaxFunctionEvaluations = 3000;
R_loss_optimal = fmincon(objective_function, P0, [], [], [], [], lb, ub, [], options)
% save("R_loss.mat", "R_loss_optimal")


T_results = system_of_equations_Euler_Nlayers_heatloss(time_data, T_initial, T_outdoor, m_s_dot_interp, num_layer_tank1,num_layer_tank2, m_layer_tank1, m_layer_tank2_optimal, R_tank12_optimal,R_loss_optimal, m_c_dot_optimal, diff_T_c_optimal, T_s_optimal);
plot_comparison_T_heatloss(time_data,T_results,time_data_1min,time_data_1min,num_layer_tank1,num_layer_tank2,T_upper_tank1_data,T_upper_tank2_data,T_bottom_tank2_data)



%% 1R   m_layer_tank2  diff_T_c T_s m_c_dot
% R0 = [R_tank12];
% % M0 = [m3,m4,m5,m6];
% P0 = [R0,m_layer_tank2',diff_T_c,T_s,m_c_dot];
% 
% 
% % Define the objective function
% objective_function = @(params) objFun_Euler_2andNlayers(time_data, T_initial, m_p_dot_interp, T_in_interp, m_s_dot_interp,num_layer_tank2, m_layer_tank1, params, T_ground_truth);
% % Call fminunc
% Aeq = [0, ones(1,num_layer_tank2), 0,0,0];
% beq = [500];
% lb = [0, zeros(1,num_layer_tank2), 1,9,800];
% ub = [0.1, ones(1,num_layer_tank2)*500, 5,15,1200];
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% options.MaxFunctionEvaluations = 3000;
% P_optimal = fmincon(objective_function, P0, [], [], Aeq, beq, lb, ub, [], options);
% save("optimal_params.mat", "P_optimal")
% 
% R_tank12_optimal = P_optimal(1);
% m_layer_tank2_optimal=zeros(num_layer_tank2,1);
% for i=1:num_layer_tank2
%     m_layer_tank2_optimal(i) = P_optimal(i+1);
% end
% 
% diff_T_c_optimal = P_optimal((1+num_layer_tank2)+1);
% T_s_optimal = P_optimal((1+num_layer_tank2)+2);
% m_c_dot_optimal = P_optimal((1+num_layer_tank2)+3);
% 
% 
% T_results = system_of_equations_Euler_Nlayers(time_data, T_initial, m_p_dot_interp, T_in_interp, m_s_dot_interp, num_layer_tank1,num_layer_tank2, m_layer_tank1, m_layer_tank2_optimal, R_tank12_optimal, m_c_dot_optimal, diff_T_c_optimal, T_s_optimal);
% 
% 
% % Plot Figure
% plot_comparison_T(time_data,T_results,num_layer_tank1,num_layer_tank2,T_upper_tank1_data,T_upper_tank2_data,T_bottom_tank2_data)
% 
