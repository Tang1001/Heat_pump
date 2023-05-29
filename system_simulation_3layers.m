%% 2023-03-27 week 13
clear all
close all

%m cp dT/dt = m_dot * cp * (T1 - T2) = Q_dot 
%[kg] * [J/kg·K] [K/s] = [kg/s] * [J/kg·K] * [K] = [J/s] = [W]
%[kg] * [J/kg·K] [K/h] /(60*60) = [kg/h]/(60*60) * [J/kg·K] * [K] = [J/s]  = [W]
%[kg] * [J/kg·K] [K/h] = [kg/h] * [J/kg·K] * [K] = [J/h] = 1/(60*60) * [W] = 1/(60*60) * [J/s]

%R (T1-T2)
%[W/K]*[K]=[W]

% Parameters and initial conditions
m1 = 250; % mass layer 1    250L~250kg
m2 = 500; % mass layer 2    500L~500kg
m3 = 250; % mass layer 3    250L~250kg
cp = 4186; % specific heat  [J/kg·K]
k= 0.6; % the thermal conductivity of water  [W/m·K]
r=0.65/2; % radius of tank  [m]
A= pi*r^2; % cross-sectional area between layers    [m^2]
H=1.64; % height of tank    [m]
h= H/2; % height of layer ??      [m]
R1 = k*A/h;
R2 = k*A/h; %              [W/K]
m_c_dot= 1*1000; % flow rate of hot water circulation    [m^3/h]=1000*[L/h]=1000*[kg/h]
diff_T_c = 5; % temperature difference  [K]
T_s = 12; % cold water temperature      [^]


%% 2023-03-27
%7:00
T1_initial=62.5; % layer 1 initial temperature
T2_initial=55.5; % layer 2 initial temperature
T3_initial=53.1; % layer 3 initial temperature
initial_conditions = [T1_initial; T2_initial; T3_initial]; % initial temperatures

%7:00 ~ 19:00
%0~12
% Time span for the simulation
t_start = 0;
t_end = 12; % end time
time_span = [t_start, t_end];


%Read data
%1min
Temp_out_water_he_13 = readcell('week 13 outlet water T in heat exchanger.csv');
Temp_out_water_he_13=cell2table(Temp_out_water_he_13(2:end,:));
[date_out_water_he_13,T_out_water_he_13]=read_csv(Temp_out_water_he_13);

FlowRate_water_hp_13 = readcell('week 13 water flow rate in heat pump.csv');
FlowRate_water_hp_13=cell2table(FlowRate_water_hp_13(2:end,:));
[date_FR_water_hp_13,FR_water_hp_13]=read_csv(FlowRate_water_hp_13);

Temp_buttom_water_tank2_13 = readcell('week 13 water T in buttom 2 tank.csv');
Temp_buttom_water_tank2_13=cell2table(Temp_buttom_water_tank2_13(2:end,:));
[date_buttom_water_tank2_13,T_buttom_water_tank2_13]=read_csv(Temp_buttom_water_tank2_13);

%8min
Temp_upper_water_tank1_13 = readcell('week 13 water T in upper 1 tank.csv');
Temp_upper_water_tank1_13=cell2table(Temp_upper_water_tank1_13(2:end,:));
Temp_upper_water_tank2_13 = readcell('week 13 water T in upper 2 tank.csv');
Temp_upper_water_tank2_13=cell2table(Temp_upper_water_tank2_13(2:end,:));
[date_upper_water_tank1_13,T_upper_water_tank1_13]=read_csv(Temp_upper_water_tank1_13);
[date_upper_water_tank2_13,T_upper_water_tank2_13]=read_csv(Temp_upper_water_tank2_13);

%1h
Consp_water_13 = readcell('week 13 water V consumption.csv');
Consp_water_13=cell2table(Consp_water_13(2:end,:));
[date_V_water_13,V_water_13]=read_csv(Consp_water_13);  %[m^3]



%%
% % Known input data
% % 2023-03-27 7:00 ~ 19:00
%1min
m_p_dot_data=FR_water_hp_13(60*7+1:60*19+1)*1000;  %60*7+1~60*19
T_in_data=T_out_water_he_13(60*7+1:60*19+1);
time_data_1min=0:1/60:12;

%1h
m_s_dot_data=V_water_13(1*9:1*20)*1000;    %[L/h] 08:00~19:00
m_s_dot_data=[m_s_dot_data(1);m_s_dot_data];
time_data_1h=0:12;

time_data_8min=0:8/60:12;

% Solve the ODE
[t, T] = ode45(@(t, T) system_of_equations(t, T, m_p_dot_data, T_in_data, m_s_dot_data, time_data_1min, time_data_1h, m1, m2, m3, cp, R1, R2, m_c_dot, diff_T_c, T_s), time_span, initial_conditions);


%% Task 1
%sample change
%id parameters R1 R2


%% Plot the results
figure;
plot(t, T(:, 1), 'r', 'LineWidth', 1.5);
hold on;
plot(t, T(:, 2), 'g', 'LineWidth', 1.5);
plot(t, T(:, 3), 'b', 'LineWidth', 1.5);
xlabel("Time (h; 0=07:00, 12=19:00)")
ylabel("Temp (°C)")
legend('Layer 1', 'Layer 2', 'Layer 3');
grid on;

%%
figure;
plot(date_upper_water_tank1_13((15*3+8)+1:(9*15+8)+1), T_upper_water_tank1_13((15*3+8)+1:(9*15+8)+1), 'r', 'LineWidth', 1.5);
hold on;
plot(date_upper_water_tank2_13((15*3+8)+1:(9*15+8)+1), T_upper_water_tank2_13((15*3+8)+1:(9*15+8)+1), 'g', 'LineWidth', 1.5);
plot(date_buttom_water_tank2_13(60*7+1:60*19+1), T_buttom_water_tank2_13(60*7+1:60*19+1), 'b', 'LineWidth', 1.5);
xlabel("Time (h; 0=07:00, 12=19:00)")
ylabel("Temp (°C)")
legend('upper water of tank 1', 'upper water of tank 2', 'buttom water of tank 3');
grid on;

%%
figure;
subplot(3,1,1)
plot(t, T(:, 1),time_data_8min, T_upper_water_tank1_13((15*3+8)+1:(9*15+8)+1));
legend("layer 1","upper water of tank 1")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
grid on
subplot(3,1,2)
plot(t, T(:, 2),time_data_8min, T_upper_water_tank2_13((15*3+8)+1:(9*15+8)+1));
legend("layer 2","upper water of tank 2")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
grid on
subplot(3,1,3)
plot(t, T(:, 3),time_data_1min, T_buttom_water_tank2_13(60*7+1:60*19+1));
legend("layer 3","buttom water of tank 2")
xlabel("Time (0=07:00, 12=19:00)")
ylabel("Temp (°C)")
grid on
set(gcf,'position',[400,300,1000,600])
sgtitle('Comparison of Tempertures (2023-03-27)') 


% %% demo
% % figure
% plot(time_data_1min, m_p_dot_data);
% figure
% plot(time_data_1min, T_in_data);
% figure
% plot(time_data_1h, m_s_dot_data);
% m_s_dot = interp1(time_data_1h, m_s_dot_data, 0:0.1:12,'next');
% figure
% plot(0:0.1:12, m_s_dot);
