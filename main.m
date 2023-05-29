clear all
close all

%% Read data
%1min
Temp_out_water_he = readcell('Data/Mon 1 May 2023 outlet water T in heat exchanger.csv');
Temp_out_water_he=cell2table(Temp_out_water_he(2:end,:));
[date_out_water_he,T_out_water_he]=read_csv(Temp_out_water_he);

FlowRate_water_hp = readcell('Data/Mon 1 May 2023 water flow rate in heat pump.csv');
FlowRate_water_hp=cell2table(FlowRate_water_hp(2:end,:));
[date_FR_water_hp,FR_water_hp]=read_csv(FlowRate_water_hp);

Temp_buttom_water_tank2 = readcell('Data/Mon 1 May 2023 water T in bottom 2 tank.csv');
Temp_buttom_water_tank2=cell2table(Temp_buttom_water_tank2(2:end,:));
[date_buttom_water_tank2,T_buttom_water_tank2]=read_csv(Temp_buttom_water_tank2);

Temp_upper_water_tank1 = readcell('Data/Mon 1 May 2023 water T in upper 1 tank.csv');
Temp_upper_water_tank1=cell2table(Temp_upper_water_tank1(2:end,:));
[date_upper_water_tank1,T_upper_water_tank1]=read_csv(Temp_upper_water_tank1);


% Temp_outdoor = readcell('Thu 13 April 2023 outdoor T.csv');
% Temp_outdoor=cell2table(Temp_outdoor(2:end,:));
% [date_outdoor,T_outdoor]=read_csv(Temp_outdoor);

%8min
Temp_upper_water_tank2 = readcell('Data/Mon 1 May 2023 water T in upper 2 tank.csv');
Temp_upper_water_tank2=cell2table(Temp_upper_water_tank2(2:end,:));
[date_upper_water_tank2,T_upper_water_tank2]=read_csv(Temp_upper_water_tank2);

%1h
Consp_water = readcell('Data/Mon 1 May 2023 water V consumption.csv');
Consp_water=cell2table(Consp_water(2:end,:));
[date_V_water,V_water]=read_csv(Consp_water);  %[m^3]


%% Choose data between 7:00~19:00
%1min
m_p_dot_data=FR_water_hp(60*7+1:60*19+1)*1000;  %[L/h] 60*7+1~60*19
T_in_data=T_out_water_he(60*7+1:60*19+1);
% T_outdoor_data=T_outdoor(60*7+1:60*19+1);
T_upper_water_tank1_data = T_upper_water_tank1(60*7+1:60*19+1);
T_buttom_water_tank2_data = T_buttom_water_tank2(60*7+1:60*19+1);

%8min
T_upper_water_tank2_data = T_upper_water_tank2(round(60/8*7):round(60/8*19)+1); %06:56~19:04



%1h
m_s_dot_data=V_water(1*8:1*20)*1000;    %[L/h] 07:00~19:00


%% Time span for the simulation
t_start = 0;
t_end = 12;     % end time [h]
sample_time = 2;    % [min]
time_span = [t_start:sample_time/60:t_end];
time_data = time_span;

time_data_1min=0:1/60:12;
time_data_1h=0:12;
time_data_8min=0-4/60:8/60:12+4/60;



%% Interpolation
[m_p_dot_interp,T_in_interp,m_s_dot_interp,...
    T_upper_water_tank1_interp,T_upper_water_tank2_interp,...
    T_buttom_water_tank2_interp]=...
interpolation_data(time_data, sample_time,m_p_dot_data,T_in_data,m_s_dot_data, ...
T_upper_water_tank1_data,T_upper_water_tank2_data,T_buttom_water_tank2_data);

% Define ground true data
data = [T_upper_water_tank1_interp',T_upper_water_tank2_interp',T_buttom_water_tank2_interp']; % Your historical data here



num_layer_tank1=2;
num_layer_tank2=4;


%Initial parameter
R12=0;  % Thermal transfer parameter between tank1 and tank2
m_c_dot= 1*1000; % flow rate of hot water circulation    [m^3/h]=1000*[L/h]=1000*[kg/h]
diff_T_c = 2; % temperature difference  [K]
T_s = 12; % cold water temperature      [°C]
params=[R12,m_c_dot,diff_T_c,T_s];


%%
T_results=system_simulation_Euler(time_data,sample_time,num_layer_tank1,num_layer_tank2,T_in_interp,m_p_dot_interp,T_buttom_water_tank2_interp,T_upper_water_tank1_interp,T_upper_water_tank2_interp,m_s_dot_interp,params);



%% Identify optimal parameters
% 1R  diff_T_c T_s m_c_dot
% Initial guess for R1 and R2
params0 = [R12,m_c_dot,diff_T_c,T_s];

params_optimal=Identify_params(time_data,sample_time,num_layer_tank1,num_layer_tank2,T_in_interp,m_p_dot_interp,T_buttom_water_tank2_interp,T_upper_water_tank1_interp,T_upper_water_tank2_interp,m_s_dot_interp,params0,data)

T_results=system_simulation_Euler(time_data,sample_time,num_layer_tank1,num_layer_tank2,T_in_interp,m_p_dot_interp,T_buttom_water_tank2_interp,T_upper_water_tank1_interp,T_upper_water_tank2_interp,m_s_dot_interp,params_optimal);


%% Identify optimal parameters of 2+4 models
% 1R  m_c_dot diff_T_c T_s  4m
num_layer_tank2=4;
m_tank2=ones(1,4)*500/num_layer_tank2;
params0 = [R12,m_c_dot,diff_T_c,T_s,m_tank2];

params_optimal=Identify_params_2and4(time_data,sample_time,T_in_interp,m_p_dot_interp,T_buttom_water_tank2_interp,T_upper_water_tank1_interp,T_upper_water_tank2_interp,m_s_dot_interp,params0,data)


m_c_dot_optimal = params_optimal(2);
diff_T_c_optimal = params_optimal(3);
T_s_optimal = params_optimal(4);
m3_optimal = params_optimal(5);
m4_optimal = params_optimal(6);
m5_optimal = params_optimal(7);
m6_optimal = params_optimal(8);



% % height of layer    [m]
% h3_optimal=H*(m3_optimal/500);
% h4_optimal=H*(m4_optimal/500);
% h5_optimal=H*(m5_optimal/500);
% h6_optimal=H*(m6_optimal/500);
% 
% % thermal resistances    R_i = h_i / (k_i * A)
% R3_optimal=h3_optimal/(k*A);
% R4_optimal=h4_optimal/(k*A);
% R5_optimal=h5_optimal/(k*A);
% R6_optimal=h6_optimal/(k*A);
% 
% % initial thermal efficiency between layers
% R23_optimal = params_optimal(1);
% R34_optimal = 1/(2 * (R3_optimal * R4_optimal) / (R3_optimal + R4_optimal));
% R45_optimal = 1/(2 * (R4_optimal * R5_optimal) / (R4_optimal + R5_optimal));
% R56_optimal = 1/(2 * (R5_optimal * R6_optimal) / (R5_optimal + R6_optimal));

 T_results=system_simulation_Euler_2and4(time_data,sample_time,T_in_interp,m_p_dot_interp,T_buttom_water_tank2_interp,T_upper_water_tank1_interp,T_upper_water_tank2_interp,m_s_dot_interp,params_optimal);

% %%
% m1=500/2;
% m2=500/2;
% R12_tank1=h_layer_tank1/(k*A);
% T_results = system_of_equations_Euler_6layers_5R_singlePipe(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, 500/2, 500/2, m3_optimal, m4_optimal, m5_optimal, m6_optimal, cp, R12_tank1, R23_optimal, R34_optimal, R45_optimal, R56_optimal, m_c_dot_optimal, diff_T_c_optimal, T_s_optimal);










