%% 2023-04-13  [W]
clear all
close all

%m cp dT/dt = m_dot * cp * (T1 - T2) = Q_dot 
%[kg] * [J/kg·K] [K/s] = [kg/s] * [J/kg·K] * [K] = [J/s] = [W]
%[kg] * [J/kg·K] [K/h] = [kg/h] * [J/kg·K] * [K] = [J/h] = 1/(60*60) * [J/s] = 1/(60*60) * [W] = 1/(60*60) * [J/s]

%R (T1-T2)
%[W/K] * [K] = [W] = [J/s] = (60*60) * [J/h] 

% Parameters and initial conditions
m1 = 250; % mass layer 1    250L~250kg
m2 = 250; % mass layer 2    250L~250kg
m3 = 250; % mass layer 3    250L~250kg
m4 = 250; % mass layer 4    250L~250kg
cp = 4186; % specific heat  [J/kg·K]        cp = 4217 - 0.3 * T
k= 0.6; % the thermal conductivity of water  [W/m·K]
r=0.65/2; % radius of tank  [m]
A= pi*r^2; % cross-sectional area between layers    [m^2]
H=1.64; % height of tank    [m]
h= H/2; % height of layer ??      [m]
R1 = k*A/h;
R2 = k*A/h; %              [W/K]
m_c_dot= 1*1000; % flow rate of hot water circulation    [m^3/h]=1000*[L/h]=1000*[kg/h]
diff_T_c = 5; % temperature difference  [K]
T_s = 12; % cold water temperature      [°C]


%% 2023-03-27


%Read data
%1min
Temp_out_water_he_130423 = readcell('Thu 13 April 2023 outlet water T in heat exchanger.csv');
Temp_out_water_he_130423=cell2table(Temp_out_water_he_130423(2:end,:));
[date_out_water_he_130423,T_out_water_he_130423]=read_csv(Temp_out_water_he_130423);

FlowRate_water_hp_130423 = readcell('Thu 13 April 2023 water flow rate in heat pump.csv');
FlowRate_water_hp_130423=cell2table(FlowRate_water_hp_130423(2:end,:));
[date_FR_water_hp_13,FR_water_hp_130423]=read_csv(FlowRate_water_hp_130423);

Temp_buttom_water_tank2_130423 = readcell('Thu 13 April 2023 water T in bottom 2 tank.csv');
Temp_buttom_water_tank2_130423=cell2table(Temp_buttom_water_tank2_130423(2:end,:));
[date_buttom_water_tank2_130423,T_buttom_water_tank2_130423]=read_csv(Temp_buttom_water_tank2_130423);

Temp_upper_water_tank1_130423 = readcell('Thu 13 April 2023 water T in upper 1 tank.csv');
Temp_upper_water_tank1_130423=cell2table(Temp_upper_water_tank1_130423(2:end,:));
[date_upper_water_tank1_130423,T_upper_water_tank1_130423]=read_csv(Temp_upper_water_tank1_130423);


Temp_outdoor_130423 = readcell('Thu 13 April 2023 outdoor T.csv');
Temp_outdoor_130423=cell2table(Temp_outdoor_130423(2:end,:));
[date_outdoor_130423,T_outdoor_130423]=read_csv(Temp_outdoor_130423);


%8min
Temp_upper_water_tank2_130423 = readcell('Thu 13 April 2023 water T in upper 2 tank.csv');
Temp_upper_water_tank2_130423=cell2table(Temp_upper_water_tank2_130423(2:end,:));
[date_upper_water_tank2_130423,T_upper_water_tank2_130423]=read_csv(Temp_upper_water_tank2_130423);

%1h
Consp_water_130423 = readcell('Thu 13 April 2023 water V consumption.csv');
Consp_water_130423=cell2table(Consp_water_130423(2:end,:));
[date_V_water_130423,V_water_130423]=read_csv(Consp_water_130423);  %[m^3]



%7:00 %1min:60*7+1=421  8min:60/8*7+1=53.5
T1_initial=T_upper_water_tank1_130423(421); % layer 1 initial temperature
T3_initial=(T_upper_water_tank2_130423(53)+T_upper_water_tank2_130423(54))/2; % layer 3 initial temperature
T4_initial=T_buttom_water_tank2_130423(421); % layer 3 initial temperature
T2_initial=(T1_initial+T3_initial)/2; % layer 2 initial temperature
initial_conditions = [T1_initial; T2_initial; T3_initial; T4_initial]; % initial temperatures

%7:00 ~ 19:00
%0~12
% Time span for the simulation
t_start = 0;
t_end = 12; % end time
time_span = [t_start:1/60:t_end];

%%
% % Known input data
% % 2023-03-27 7:00 ~ 19:00
%1min
m_p_dot_data=FR_water_hp_130423(60*7+1:60*19+1)*1000;  %[L/h] 60*7+1~60*19
T_in_data=T_out_water_he_130423(60*7+1:60*19+1);
time_data_1min=0:1/60:12;
T_outdoor_data=T_outdoor_130423(60*7+1:60*19+1);

%1h
m_s_dot_data=V_water_130423(1*9:1*20)*1000;    %[L/h] 08:00~19:00
m_s_dot_data=[m_s_dot_data(1);m_s_dot_data];
time_data_1h=0:12;

time_data_8min=0:8/60:12;

% Solve the ODE
[t, T] = ode45(@(t, T) system_of_equations_4layers(t, T, m_p_dot_data, T_in_data, m_s_dot_data, time_data_1min, time_data_1h, m1, m2, m3, m4, cp, R1, R2, m_c_dot, diff_T_c, T_s), time_span, initial_conditions);

%% ID R1R2
% 
T_upper_water_tank1_130423_interp = interp1(time_data_1min,T_upper_water_tank1_130423(60*7+1:60*19+1),time_span);
T_upper_water_tank2_130423_interp = interp1(time_data_8min,T_upper_water_tank2_130423(round(60/8*7)+1:round(60/8*19)+1),time_span);
T_buttom_water_tank2_130423_interp = interp1(time_data_1min,T_buttom_water_tank2_130423(60*7+1:60*19+1),time_span);




data = [T_upper_water_tank1_130423_interp',T_upper_water_tank2_130423_interp',T_buttom_water_tank2_130423_interp']; % Your historical data here
% 

%R1R2
% Initial guess for the parameters
params0 = [R1; R2];
% params0 = [10; 10];

% Optimize the parameters
optimalParams = fminsearch(@(params) objFun(m_p_dot_data, T_in_data, m_s_dot_data, time_data_1min, time_data_1h, m1, m2, m3, m4, cp, m_c_dot, diff_T_c, T_s, params, time_span, initial_conditions, data), params0)

% Resolve the ODE
[t, T] = ode45(@(t, T) system_of_equations_4layers_IDR1R2(t, T, m_p_dot_data, T_in_data, m_s_dot_data, time_data_1min, time_data_1h, m1, m2, m3, m4, cp, m_c_dot, diff_T_c, T_s, optimalParams), time_span, initial_conditions);


% %R1R2 and heat loss
% % Initial guess for the parameters
% params0 = [R1; R2; zeros(4,1)];
% 
% % Optimize the parameters
% [optimalParams,fval,exitflag,output] = fminsearch(@(params) objFun_heatloss_outdoor(m_p_dot_data, T_in_data, T_outdoor_data, m_s_dot_data, time_data_1min, time_data_1h, m1, m2, m3, m4, cp, m_c_dot, diff_T_c, T_s, params, time_span, initial_conditions, data), params0)
% 
% % Resolve the ODE
% [t, T] = ode45(@(t, T) system_of_equations_4layers_IDR1R2Const_outdoor(t, T, m_p_dot_data, T_in_data, T_outdoor_data, m_s_dot_data, time_data_1min, time_data_1h, m1, m2, m3, m4, cp, m_c_dot, diff_T_c, T_s, optimalParams), time_span, initial_conditions);
% 

%% Genetic
% %R1R2
% fitnessfunction = @(params) objFun(m_p_dot_data, T_in_data, m_s_dot_data, time_data_1min, time_data_1h, m1, m2, m3, m4, cp, m_c_dot, diff_T_c, T_s, params, time_span, initial_conditions, data);
% 
% % Define the options
% options = optimoptions('ga','Display','iter','PlotFcn',@gaplotbestf);
% 
% % Run the genetic algorithm
% optimalParams = ga(fitnessfunction,2,[],[],[],[],[0 0],[Inf Inf],[],options);
% %%
% [t, T] = ode45(@(t, T) system_of_equations_4layers_IDR1R2(t, T, m_p_dot_data, T_in_data, m_s_dot_data, time_data_1min, time_data_1h, m1, m2, m3, m4, cp, m_c_dot, diff_T_c, T_s, optimalParams), time_span, initial_conditions);

%%
%R1R2R3
fitnessfunction = @(params) objFun_R1R2R3(m_p_dot_data, T_in_data, m_s_dot_data, time_data_1min, time_data_1h, m1, m2, m3, m4, cp, m_c_dot, diff_T_c, T_s, params, time_span, initial_conditions, data);

% Define the options
options = optimoptions('ga','Display','iter','PlotFcn',@gaplotbestf);

% Run the genetic algorithm
optimalParams = ga(fitnessfunction,3,[],[],[],[],[0 0],[Inf Inf],[],options);

%%
[t, T] = ode45(@(t, T) system_of_equations_4layers_IDR1R2R3(t, T, m_p_dot_data, T_in_data, m_s_dot_data, time_data_1min, time_data_1h, m1, m2, m3, m4, cp, m_c_dot, diff_T_c, T_s, optimalParams), time_span, initial_conditions);


%% Plot the results
figure;
plot(t, T(:, 1), 'r', 'LineWidth', 1.5);
hold on;
plot(t, T(:, 3), 'g', 'LineWidth', 1.5);
plot(t, T(:, 4), 'b', 'LineWidth', 1.5);
xlabel("Time (h; 0=07:00, 12=19:00)")
ylabel("Temp (°C)")
legend('Layer 1', 'Layer 3', 'Layer 4');
grid on;

%%
figure;
plot(date_upper_water_tank1_130423(60*7+1:60*19+1), T_upper_water_tank1_130423(60*7+1:60*19+1), 'r', 'LineWidth', 1.5);
hold on;
plot(date_upper_water_tank2_130423((15*3+8)+1:(9*15+8)+1), T_upper_water_tank2_130423((15*3+8)+1:(9*15+8)+1), 'g', 'LineWidth', 1.5);
plot(date_buttom_water_tank2_130423(60*7+1:60*19+1), T_buttom_water_tank2_130423(60*7+1:60*19+1), 'b', 'LineWidth', 1.5);
xlabel("Time (h; 0=07:00, 12=19:00)")
ylabel("Temp (°C)")
legend('upper water of tank 1', 'upper water of tank 2', 'buttom water of tank 2');
grid on;

%%
figure;
subplot(3,1,1)
plot(t, T(:, 1),time_data_1min, T_upper_water_tank1_130423(60*7+1:60*19+1));
legend("layer 1","upper water of tank 1")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
grid on
subplot(3,1,2)
plot(t, T(:, 3),time_data_8min, T_upper_water_tank2_130423((15*3+8)+1:(9*15+8)+1));
legend("layer 3","upper water of tank 2")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
grid on
subplot(3,1,3)
plot(t, T(:, 4),time_data_1min, T_buttom_water_tank2_130423(60*7+1:60*19+1));
legend("layer 4","buttom water of tank 2")
xlabel("Time (0=07:00, 12=19:00)")
ylabel("Temp (°C)")
grid on
set(gcf,'position',[400,300,1000,600])
sgtitle('Comparison of Tempertures (2023-04-13) ') 

% figure
% bar(time_data_1h,m_s_dot_data)
% grid on
% set(gcf,'position',[400,300,1000,600])
%% demo
% % figure
% % plot(time_data_1min, m_p_dot_data);
% % figure
% % plot(time_data_1min, T_in_data);
% figure
% plot(time_data_1h, m_s_dot_data);
% m_s_dot1 = interp1(time_data_1h, m_s_dot_data, 0:0.1:12,'next');
% figure
% plot(0:0.1:12, m_s_dot1);
% m_s_dot2 = interp1(time_data_1h, m_s_dot_data, 0:0.1:12,'previous');
% figure
% plot(0:0.1:12, m_s_dot2);
