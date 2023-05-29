%% 2023-04-13  % the Euler method
clear all
close all

%m cp dT/dt = m_dot * cp * (T1 - T2) = Q_dot 
%[kg] * [J/kg·K] [K/s] = [kg/s] * [J/kg·K] * [K] = [J/s] = [W]
%[kg] * [J/kg·K] [K/h] = [kg/h] * [J/kg·K] * [K] = [J/h] = 1/(60*60) * [J/s] = 1/(60*60) * [W] = 1/(60*60) * [J/s]

%R (T1-T2)
%[W/K] * [K] = [W] = [J/s] = (60*60) * [J/h] 

% Parameters and initial conditions
cp = 4186; % specific heat  [J/kg·K]        cp = 4217 - 0.3 * T
k= 0.6; % the thermal conductivity of water  [W/m·K]
r=0.65/2; % radius of tank  [m]
A= pi*r^2; % cross-sectional area between layers    [m^2]
H=1.64; % height of tank    [m]

% mass layer   L~kg
m1 = 250; 
m2 = 250;
m3 = 100;
m4 = 150;
m5 = 150;
m6 = 100;

% height of layer    [m]
h1=H*(m1/500); 
h2=H*(m2/500);
h3=H*(m3/500);
h4=H*(m4/500);
h5=H*(m5/500);
h6=H*(m6/500);

% thermal resistances    R_i = h_i / (k_i * A)
R1=h1/(k*A);
R2=h2/(k*A);
R3=h3/(k*A);
R4=h4/(k*A);
R5=h5/(k*A);
R6=h6/(k*A);
R1_tot = R1 + R2;
R2_tot = R3 + R4 + R5 + R6;


% initial thermal efficiency between layers
R12 = 1/R1;         % k*A/h_i      [W/K]
R23 = 0;
R34 = 1/R3;
R45 = 1/R4;
R56 = 1/R5;

m_c_dot= 1*1000; % flow rate of hot water circulation    [m^3/h]=1000*[L/h]=1000*[kg/h]
diff_T_c = 2; % temperature difference  [K]
T_s = 12; % cold water temperature      [°C]


%% Read data 2023-03-27
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
T6_initial=T_buttom_water_tank2_130423(421); % layer 5 initial temperature


d_T1 = (T1_initial - T3_initial) * (R1 / R1_tot);
T2_initial = T1_initial - d_T1;     % layer 2 initial temperature
d_T3 = (T3_initial - T6_initial) * (R3 / R2_tot);
T4_initial= T3_initial - d_T3; % layer 4 initial temperature
d_T4 = (T3_initial - T6_initial) * (R4 / R2_tot);
T5_initial= T4_initial - d_T4; % layer 4 initial temperature



%% Known input data
%2023-03-27 7:00 ~ 19:00
%1min
m_p_dot_data=FR_water_hp_130423(60*7+1:60*19+1)*1000;  %[L/h] 60*7+1~60*19
T_in_data=T_out_water_he_130423(60*7+1:60*19+1);
T_outdoor_data=T_outdoor_130423(60*7+1:60*19+1);

%1h
m_s_dot_data=V_water_130423(1*8:1*20)*1000;    %[L/h] 07:00~19:00
% m_s_dot_data=V_water_130423(1*9:1*20)*1000;    %[L/h] 08:00~19:00
% m_s_dot_data=[m_s_dot_data(1);m_s_dot_data];

time_data_1min=0:1/60:12;
time_data_1h=0:12;
time_data_8min=0:8/60:12;

%7:00 ~ 19:00
%0~12
% Time span for the simulation
t_start = 0;
t_end = 12;     % end time [h]
sample_time = 2;    % [min]
time_span = [t_start:sample_time/60:t_end];
time_data = time_span;

% Interpolation
m_p_dot_interp = interp1(time_data_1min, m_p_dot_data, time_span);
T_in_interp = interp1(time_data_1min, T_in_data, time_span);

% estimate min-level flow rate of supplying water
% m_s_dot_interp = interp1(time_data_1h, m_s_dot_data, time_span,'next');
m_s_dot_interp = interp1(time_data_1h, m_s_dot_data, time_span,'spline');
% m_s_dot_interp = interp1(time_data_1h, m_s_dot_data, time_span,'nearest');


%%
% Adjust the minute-level data to make the total sums equal
total_hourly = sum(m_s_dot_data(2:end));
total_minute = sum(m_s_dot_interp(2:end))*(sample_time/60); % Don't divide by 60
%%
difference = total_hourly - total_minute;
m_s_dot_interp_adjusted = m_s_dot_interp;
m_s_dot_interp_adjusted(2:end) = m_s_dot_interp(2:end) + (difference/(sample_time/60))/ (length(m_s_dot_interp)-1);
% Calculate the new total for the minute-level data
total_minute_adjusted = sum(m_s_dot_interp_adjusted(2:end))*(sample_time/60); % Don't divide by 60
disp(['Total water consumption (Hourly data): ', num2str(total_hourly), ' liters']);
disp(['Total water consumption (Minute data, adjusted): ', num2str(total_minute_adjusted), ' liters']);


%%
m_s_dot_interp=m_s_dot_interp_adjusted;

% Initialization
T1 = T1_initial; % initial temperature T1
T2 = T2_initial; % initial temperature T2
T3 = T3_initial; % initial temperature T3
T4 = T4_initial; % initial temperature T4
T5 = T5_initial; % initial temperature T5
T6 = T6_initial; % initial temperature T6
T = [T1; T2; T3; T4; T5; T6];



%% Euler method
T_results = system_of_equations_Euler_6layers_5R(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, m1, m2, m3, m4, m5, m6, cp, R12, R23, R34, R45, R56, m_c_dot, diff_T_c, T_s);




%% Identify optimal parameters
% Define data
T_upper_water_tank1_130423_interp = interp1(time_data_1min,T_upper_water_tank1_130423(60*7+1:60*19+1),time_data);
T_upper_water_tank2_130423_interp = interp1(time_data_8min,T_upper_water_tank2_130423(round(60/8*7)+1:round(60/8*19)+1),time_data);
T_buttom_water_tank2_130423_interp = interp1(time_data_1min,T_buttom_water_tank2_130423(60*7+1:60*19+1),time_data);
data = [T_upper_water_tank1_130423_interp',T_upper_water_tank2_130423_interp',T_buttom_water_tank2_130423_interp']; % Your historical data here


%% 5R   4m  diff_T_c
% Initial guess for R1 and R2
R0 = [R12, R23, R34, R45, R56];
M0 = [m3,m4,m5,m6];
P0 = [R0,M0,diff_T_c];

% Define the objective function
objective_function = @(params) objFun_Euler_6layers_5R_4m_diffT(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, m1, m2, cp, m_c_dot, T_s, params, data);

% Call fminunc
Aeq = [0,0,0,0,0,1,1,1,1,0];
beq = 500;
lb = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
ub = [1,0.1,1,1,1, 500, 500, 500, 500, 5];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
R_optimal = fmincon(objective_function, P0, [], [], Aeq, beq, lb, ub, [], options)

R12_optimal = R_optimal(1);
R23_optimal = R_optimal(2);
R34_optimal = R_optimal(3);
R45_optimal = R_optimal(4);
R56_optimal = R_optimal(5);
m3_optimal = R_optimal(6);
m4_optimal = R_optimal(7);
m5_optimal = R_optimal(8);
m6_optimal = R_optimal(9);
diff_T_c_optimal = R_optimal(10);

%%
T_results = system_of_equations_Euler_6layers_5R(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, m1, m2, m3_optimal, m4_optimal, m5_optimal, m6_optimal, cp, R12_optimal, R23_optimal, R34_optimal, R45_optimal, R56_optimal, m_c_dot, diff_T_c_optimal, T_s);


%% 5R   4m  diff_T_c T_s
% Initial guess for R1 and R2
R0 = [R12, R23, R34, R45, R56];
M0 = [m3,m4,m5,m6];
P0 = [R0,M0,diff_T_c,T_s];

% Define the objective function
objective_function = @(params) objFun_Euler_6layers_5R_4m_diffT_Ts(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, m1, m2, cp, m_c_dot, params, data);

% Call fminunc
Aeq = [0,0,0,0,0,1,1,1,1,0,0];
beq = 500;
lb = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
ub = [1,0.1,1,1,1, 500, 500, 500, 500, 5,20];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
R_optimal = fmincon(objective_function, P0, [], [], Aeq, beq, lb, ub, [], options)

R12_optimal = R_optimal(1);
R23_optimal = R_optimal(2);
R34_optimal = R_optimal(3);
R45_optimal = R_optimal(4);
R56_optimal = R_optimal(5);
m3_optimal = R_optimal(6);
m4_optimal = R_optimal(7);
m5_optimal = R_optimal(8);
m6_optimal = R_optimal(9);
diff_T_c_optimal = R_optimal(10);
T_s_optimal = R_optimal(11);

%%
T_results = system_of_equations_Euler_6layers_5R(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, m1, m2, m3_optimal, m4_optimal, m5_optimal, m6_optimal, cp, R12_optimal, R23_optimal, R34_optimal, R45_optimal, R56_optimal, m_c_dot, diff_T_c_optimal, T_s_optimal);



%% 5R   4m  diff_T_c T_s m_c_dot
% Initial guess for R1 and R2
R0 = [R12, R23, R34, R45, R56];
M0 = [m3,m4,m5,m6];
P0 = [R0,M0,diff_T_c,T_s,m_c_dot];

% Define the objective function
objective_function = @(params) objFun_Euler_6layers_5R_4m_diffT_Ts_flowrate(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, m1, m2, cp, params, data);

% Call fminunc
Aeq = [0,0,0,0,0, 1,1,1,1, 0,0,0];
beq = 500;
lb = [0,0,0,0,0, 0,0,0,0, 0,0,800];
ub = [1,0.1,1,1,1, 500,500,500,500, 5,20,1200];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
R_optimal = fmincon(objective_function, P0, [], [], Aeq, beq, lb, ub, [], options)

R12_optimal = R_optimal(1);
R23_optimal = R_optimal(2);
R34_optimal = R_optimal(3);
R45_optimal = R_optimal(4);
R56_optimal = R_optimal(5);
m3_optimal = R_optimal(6);
m4_optimal = R_optimal(7);
m5_optimal = R_optimal(8);
m6_optimal = R_optimal(9);
diff_T_c_optimal = R_optimal(10);
T_s_optimal = R_optimal(11);
m_c_dot_optimal = R_optimal(12);
%%
T_results = system_of_equations_Euler_6layers_5R(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, m1, m2, m3_optimal, m4_optimal, m5_optimal, m6_optimal, cp, R12_optimal, R23_optimal, R34_optimal, R45_optimal, R56_optimal, m_c_dot_optimal, diff_T_c_optimal, T_s_optimal);




%% 1R   4m  diff_T_c T_s m_c_dot
% Initial guess for R1 and R2
R0 = [R23];
M0 = [m3,m4,m5,m6];
P0 = [R0,M0,diff_T_c,T_s,m_c_dot];

% Define the objective function
objective_function = @(params) objFun_Euler_6layers_1R_4m_diffT_Ts_flowrate(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, m1, m2, cp, R12, H, A, k, params, data);

% Call fminunc
Aeq = [0, 1,1,1,1, 0,0,0];
beq = [500];
lb = [0, 0,0,0,0, 0,0,800];
ub = [0.1, 500,500,500,500, 5,30,1200];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
R_optimal = fmincon(objective_function, P0, [], [], Aeq, beq, lb, ub, [], options)

R23_optimal = R_optimal(1);
m3_optimal = R_optimal(2);
m4_optimal = R_optimal(3);
m5_optimal = R_optimal(4);
m6_optimal = R_optimal(5);
diff_T_c_optimal = R_optimal(6);
T_s_optimal = R_optimal(7);
m_c_dot_optimal = R_optimal(8);


% height of layer    [m]
h3_optimal=H*(m3_optimal/500);
h4_optimal=H*(m4_optimal/500);
h5_optimal=H*(m5_optimal/500);
h6_optimal=H*(m6_optimal/500);

% thermal resistances    R_i = h_i / (k_i * A)
R3_optimal=h3_optimal/(k*A);
R4_optimal=h4_optimal/(k*A);
R5_optimal=h5_optimal/(k*A);
R6_optimal=h6_optimal/(k*A);


% initial thermal efficiency between layers
R34_optimal = 1/(2 * (R3_optimal * R4_optimal) / (R3_optimal + R4_optimal));
R45_optimal = 1/(2 * (R4_optimal * R5_optimal) / (R4_optimal + R5_optimal));
R56_optimal = 1/(2 * (R5_optimal * R6_optimal) / (R5_optimal + R6_optimal));


%%
T_results = system_of_equations_Euler_6layers_5R(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, m1, m2, m3_optimal, m4_optimal, m5_optimal, m6_optimal, cp, R12, R23_optimal, R34_optimal, R45_optimal, R56_optimal, m_c_dot_optimal, diff_T_c_optimal, T_s_optimal);



% %% R1R2R3
% % Initial guess for R1, R2 and R3
% R0 = ones(3,1)*R1;
% 
% % Define the objective function
% objective_function = @(params) objFun_Euler_R1R2R3(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, m1, m2, m3, m4, cp, m_c_dot, diff_T_c, T_s, params, data);
% 
% % Call fminunc
% lb = [0, 0 ,0];
% ub = [50,50, 10];
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% R_optimal = fmincon(objective_function, R0, [], [], [], [], lb, [], [], options)
% 
% R1_optimal = R_optimal(1);
% R2_optimal = R_optimal(2);
% R3_optimal = R_optimal(3);
% 
% %%
% T_results = system_of_equations_Euler_4layers_R1R2R3(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, m1, m2, m3, m4, cp, R1_optimal, R2_optimal, R3_optimal, m_c_dot, diff_T_c, T_s);



%% Plot the results
% figure;
% plot(time_data, T_results(1, :), 'r', 'LineWidth', 1.5);
% hold on;
% plot(time_data, T_results(2, :), 'g', 'LineWidth', 1.5);
% plot(time_data, T_results(3, :), 'b', 'LineWidth', 1.5);
% plot(time_data, T_results(4, :), 'y', 'LineWidth', 1.5);
% xlabel('Time (hours)');
% ylabel('Temperature (K)');
% legend('T1', 'T2', 'T3', 'T4');
% grid on;

%%
% 
% % T_bottom_tank1 = (interp1(time_data_1min,T_upper_water_tank1_130423(60*7+1:60*19+1),time_data_8min)'+T_upper_water_tank2_130423((15*3+8)+1:(9*15+8)+1))/2;
% % T_middle_tank2 = (interp1(time_data_1min,T_buttom_water_tank2_130423(60*7+1:60*19+1),time_data_8min)'+T_upper_water_tank2_130423((15*3+8)+1:(9*15+8)+1))/2;
% 
% T_bottom_tank1 = T_upper_water_tank1_130423(60*7+1:60*19+1) - d_T1;
% T_middle_tank2 = T_upper_water_tank2_130423((15*3+8)+1:(9*15+8)+1) - d_T3;
% 
% 
% figure;
% subplot(5,1,1)
% plot(time_data, T_results(:,1),time_data_1min, T_upper_water_tank1_130423(60*7+1:60*19+1));
% legend("layer 1","upper water of tank 1")
% xlabel("Time (0=07:00; 12=19:00)")
% ylabel("Temp (°C)")
% grid on
% subplot(5,1,2)
% plot(time_data, T_results(:,2),time_data_1min, T_bottom_tank1);
% legend("layer 2","buttom water of tank 1")
% xlabel("Time (0=07:00; 12=19:00)")
% ylabel("Temp (°C)")
% grid on
% subplot(5,1,3)
% plot(time_data, T_results(:,3),time_data_8min, T_upper_water_tank2_130423((15*3+8)+1:(9*15+8)+1));
% legend("layer 3","upper water of tank 2")
% xlabel("Time (0=07:00; 12=19:00)")
% ylabel("Temp (°C)")
% grid on
% subplot(5,1,4)
% plot(time_data, T_results(:,4),time_data_8min, T_middle_tank2);
% legend("layer 4","buttom water of tank 2")
% xlabel("Time (0=07:00, 12=19:00)")
% ylabel("Temp (°C)")
% grid on
% subplot(5,1,5)
% plot(time_data, T_results(:,5),time_data_1min, T_buttom_water_tank2_130423(60*7+1:60*19+1));
% legend("layer 5","buttom water of tank 2")
% xlabel("Time (0=07:00, 12=19:00)")
% ylabel("Temp (°C)")
% grid on
% set(gcf,'position',[400,300,1000,600])
% sgtitle('Comparison of Tempertures (2023-04-13) ') 

%%
figure;
subplot(3,1,1)
plot(time_data, T_results(:,1),time_data_1min, T_upper_water_tank1_130423(60*7+1:60*19+1));
legend("layer 1","upper water of tank 1")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
grid on
subplot(3,1,2)
plot(time_data, T_results(:,3),time_data_8min, T_upper_water_tank2_130423((15*3+8)+1:(9*15+8)+1));
legend("layer 3","upper water of tank 2")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
grid on
subplot(3,1,3)
yyaxis left
plot(time_data, T_results(:,6),time_data_1min, T_buttom_water_tank2_130423(60*7+1:60*19+1));
xlabel("Time (0=07:00, 12=19:00)")
ylabel("Temp (°C)")
yyaxis right
plot(time_data_1min,m_p_dot_data,time_data,m_s_dot_interp)
legend("layer 5","buttom water of tank 2")
grid on
set(gcf,'position',[400,300,1000,600])
sgtitle('Comparison of Tempertures (2023-04-13) ')

%%
figure;
yyaxis left
plot(time_data_1min, T_buttom_water_tank2_130423(60*7+1:60*19+1));
yyaxis right
plot(time_data_1min,FR_water_hp_130423(60*7+1:60*19+1)*1000,0:12,V_water_130423(8:20)*1000,time_data,m_s_dot_interp)
grid on
set(gcf,'position',[400,300,1000,600])