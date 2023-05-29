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
Temp_out_water_he = readcell('../Data/Mon 1 May 2023 outlet water T in heat exchanger.csv');
Temp_out_water_he=cell2table(Temp_out_water_he(2:end,:));
[date_out_he,T_out_he]=read_csv(Temp_out_water_he);

FlowRate_water_hp = readcell('../Data/Mon 1 May 2023 water flow rate in heat pump.csv');
FlowRate_water_hp=cell2table(FlowRate_water_hp(2:end,:));
[date_FR_hp,FR_hp]=read_csv(FlowRate_water_hp);

Temp_bottom_water_tank2 = readcell('../Data/Mon 1 May 2023 water T in bottom 2 tank.csv');
Temp_bottom_water_tank2=cell2table(Temp_bottom_water_tank2(2:end,:));
[date_bottom_tank2,T_bottom_tank2]=read_csv(Temp_bottom_water_tank2);

Temp_upper_water_tank1 = readcell('../Data/Mon 1 May 2023 water T in upper 1 tank.csv');
Temp_upper_water_tank1=cell2table(Temp_upper_water_tank1(2:end,:));
[date_upper_tank1,T_upper_tank1]=read_csv(Temp_upper_water_tank1);


%8min
Temp_upper_water_tank2 = readcell('../Data/Mon 1 May 2023 water T in upper 2 tank.csv');
Temp_upper_water_tank2=cell2table(Temp_upper_water_tank2(2:end,:));
[date_upper_tank2,T_upper_tank2]=read_csv(Temp_upper_water_tank2);

%1h
Consp_water = readcell('../Data/Mon 1 May 2023 water V consumption.csv');
Consp_water=cell2table(Consp_water(2:end,:));
[date_V_water,V_water]=read_csv(Consp_water);  %[m^3]


%%
% set the number of layers
num_layer_tank1=2;
num_layer_tank2=4;

% initial masses of layers   L~kg
m_layer_tank1=[250;250];
m_layer_tank2=[100;150;150;100];

R_tank12=0;
m_c_dot= 1*1000; % flow rate of hot water circulation    [m^3/h]=1000*[L/h]=1000*[kg/h]
diff_T_c = 2; % temperature difference  [K]
T_s = 12; % cold water temperature      [°C]


%% Choose data between 7:00 ~ 19:00
%7:00 ~ 19:00
%0~12
% Time span for the simulation
t_start = 0;
t_end = 12;     % end time [h]
sample_time = 2;    % [min]
time_span = [t_start:sample_time/60:t_end];
time_data = time_span;

%1min
m_p_dot_data=FR_hp(60*7+1:60*19+1)*1000;  %[L/h] 60*7+1~60*19
T_in_data=T_out_he(60*7+1:60*19+1);
T_upper_tank1_data = T_upper_tank1(60*7+1:60*19+1);
T_bottom_tank2_data = T_bottom_tank2(60*7+1:60*19+1);

%8min
T_upper_tank2_data = T_upper_tank2(round(60/8*7):round(60/8*19)+1); %06:56~19:04

%1h
m_s_dot_data=V_water(1*8:1*20)*1000;    %[L/h] 07:00~19:00



T_results=simulation_Euler_Nlayer(time_data,sample_time,num_layer_tank1,num_layer_tank2,m_layer_tank1, m_layer_tank2, R_tank12, m_c_dot, diff_T_c, T_s, T_in_data,m_p_dot_data,T_bottom_tank2_data,T_upper_tank1_data,T_upper_tank2_data,m_s_dot_data);


plot_comparison_T(time_data,T_results,num_layer_tank1,num_layer_tank2,T_upper_tank1_data,T_upper_tank2_data,T_bottom_tank2_data)


