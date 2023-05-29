%% the Euler method 
% mutil-layers
clear all
close all

%m cp dT/dt = m_dot * cp * (T1 - T2) = Q_dot 
%[kg] * [J/kg·K] [K/s] = [kg/s] * [J/kg·K] * [K] = [J/s] = [W]
%[kg] * [J/kg·K] [K/h] = [kg/h] * [J/kg·K] * [K] = [J/h] = 1/(60*60) * [J/s] = 1/(60*60) * [W] = 1/(60*60) * [J/s]

%R (T1-T2)
%[W/K] * [K] = [W] = [J/s] = (60*60) * [J/h] 


num_layer_tank1 = 2;
num_layer_tank2 = 4;


% Parameters and initial conditions
cp = 4186; % specific heat  [J/kg·K]        cp = 4217 - 0.3 * T
k= 0.6; % the thermal conductivity of water  [W/m·K]
r=0.65/2; % radius of tank  [m]
A= pi*r^2; % cross-sectional area between layers    [m^2]
H=1.64; % height of tank    [m]

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


R12=0;

m_c_dot= 1*1000; % flow rate of hot water circulation    [m^3/h]=1000*[L/h]=1000*[kg/h]
diff_T_c = 2; % temperature difference  [K]
T_s = 12; % cold water temperature      [°C]


%% Read data 2023-03-27
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



T_initial_layer_tank1 = zeros(1,num_layer_tank1);
T_initial_layer_tank2 = zeros(1,num_layer_tank2);

% measurement values
T_initial_layer_tank1(1)=T_upper_water_tank1(421);
T_initial_layer_tank2(1)=(T_upper_water_tank2(53)+T_upper_water_tank2(54))/2;
T_initial_layer_tank2(end)=T_buttom_water_tank2(421);


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



% %7:00 %1min:60*7+1=421  8min:60/8*7+1=53.5
% T1_initial=T_upper_water_tank1(421); % layer 1 initial temperature
% T3_initial=(T_upper_water_tank2(53)+T_upper_water_tank2(54))/2; % layer 3 initial temperature
% T6_initial=T_buttom_water_tank2(421); % layer 5 initial temperature
% 
% 
% d_T1 = (T1_initial - T3_initial) * (R1 / R1_tot);
% T2_initial = T1_initial - d_T1;     % layer 2 initial temperature
% d_T3 = (T3_initial - T6_initial) * (R3 / R2_tot);
% T4_initial= T3_initial - d_T3; % layer 4 initial temperature
% d_T4 = (T3_initial - T6_initial) * (R4 / R2_tot);
% T5_initial= T4_initial - d_T4; % layer 4 initial temperature



%% Known input data
%2023-03-27 7:00 ~ 19:00
%1min
m_p_dot_data=FR_water_hp(60*7+1:60*19+1)*1000;  %[L/h] 60*7+1~60*19
T_in_data=T_out_water_he(60*7+1:60*19+1);
% T_outdoor_data=T_outdoor(60*7+1:60*19+1);

%1h
m_s_dot_data=V_water(1*8:1*20)*1000;    %[L/h] 07:00~19:00
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

peak_shift = 0.5; 
% time_data_1h_interp = time_data_1h-0.5; %make the value in the middle of each hour
time_data_1h_interp = time_data_1h-peak_shift; %make the value in the middle of each hour
m_s_dot_interp = interp1(time_data_1h_interp, m_s_dot_data, time_span,'spline');
% m_s_dot_interp = interp1(time_data_1h, m_s_dot_data, time_span,'spline');
% m_s_dot_interp = interp1(time_data_1h, m_s_dot_data, time_span,'nearest');

% Adjust the minute-level data to make the total sums equal
total_hourly = sum(m_s_dot_data(2:end));
total_minute = sum(m_s_dot_interp(2:end))*(sample_time/60); % Don't divide by 60

difference = total_hourly - total_minute;
m_s_dot_interp_adjusted = m_s_dot_interp;
m_s_dot_interp_adjusted(2:end) = m_s_dot_interp(2:end) + (difference/(sample_time/60))/ (length(m_s_dot_interp)-1);
% Calculate the new total for the minute-level data
total_minute_adjusted = sum(m_s_dot_interp_adjusted(2:end))*(sample_time/60); % Don't divide by 60
disp(['Total water consumption (Hourly data): ', num2str(total_hourly), ' liters']);
disp(['Total water consumption (Minute data, adjusted): ', num2str(total_minute_adjusted), ' liters']);


%%
m_s_dot_interp=m_s_dot_interp_adjusted;



% Define ground true data
T_upper_water_tank1_interp = interp1(time_data_1min,T_upper_water_tank1(60*7+1:60*19+1),time_data);
T_upper_water_tank2_interp = interp1(time_data_8min,T_upper_water_tank2(round(60/8*7)+1:round(60/8*19)+1),time_data);
T_buttom_water_tank2_interp = interp1(time_data_1min,T_buttom_water_tank2(60*7+1:60*19+1),time_data);
data = [T_upper_water_tank1_interp',T_upper_water_tank2_interp',T_buttom_water_tank2_interp']; % Your historical data here


% Initialization
T=[];
for i=1:num_layer_tank1
    T=[T;T_initial_layer_tank1(i)];
end
for i=1:num_layer_tank2
    T=[T;T_initial_layer_tank2(i)];
end

% T1 = T1_initial; % initial temperature T1
% T2 = T2_initial; % initial temperature T2
% T3 = T3_initial; % initial temperature T3
% T4 = T4_initial; % initial temperature T4
% T5 = T5_initial; % initial temperature T5
% T6 = T6_initial; % initial temperature T6
% T = [T1; T2; T3; T4; T5; T6];



%% Euler method
% T_results = system_of_equations_Euler_6layers_5R(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, m1, m2, m3, m4, m5, m6, cp, R12, R23, R34, R45, R56, m_c_dot, diff_T_c, T_s);
% T_results = system_of_equations_Euler_6layers_5R_singlePipe(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, m1, m2, m3, m4, m5, m6, cp, R12, R23, R34, R45, R56, m_c_dot, diff_T_c, T_s);

T_results = system_of_equations_Euler_Nlayers_singlePipe(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, num_layer_tank1, num_layer_tank2, m_layer_tank1,m_layer_tank2, R_thermal_tank1, R_thermal_tank2, R12, cp, m_c_dot, diff_T_c, T_s);



%% Identify optimal parameters
% 1R  diff_T_c T_s m_c_dot
% Initial guess for R1 and R2
R0 = [R12];
P0 = [R0,diff_T_c,T_s,m_c_dot];

% Define the objective function
objective_function = @(params) objFun_Euler_Nlayers_1R_diffT_Ts_flowrate(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, num_layer_tank1, num_layer_tank2, m_layer_tank1,m_layer_tank2, R_thermal_tank1, R_thermal_tank2, cp, H, A, k, params, data);

% Call fminunc
Aeq = [0, 0,0,0];
beq = [500];
lb = [0, 0,0,800];
ub = [0.1, 5,30,1200];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
R_optimal = fmincon(objective_function, P0, [], [], Aeq, beq, lb, ub, [], options)

R12_optimal = R_optimal(1);
diff_T_c_optimal = R_optimal(2);
T_s_optimal = R_optimal(3);
m_c_dot_optimal = R_optimal(4);


%%
T_results = system_of_equations_Euler_Nlayers_singlePipe(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, num_layer_tank1, num_layer_tank2, m_layer_tank1,m_layer_tank2, R_thermal_tank1, R_thermal_tank2, R12_optimal, cp, m_c_dot_optimal, diff_T_c_optimal, T_s_optimal);



%%


err1 = sqrt(mean((T_results(:,1) - data(:,1)).^2));               %Root Mean Squared Error
err2 = sqrt(mean((T_results(:,num_layer_tank1+1) - data(:,2)).^2));               %Root Mean Squared Error
err3 = sqrt(mean((T_results(:,num_layer_tank1+num_layer_tank2) - data(:,3)).^2));               %Root Mean Squared Error
figure;
subplot(3,1,1)
plot(time_data, T_results(:,1),time_data_1min, T_upper_water_tank1(60*7+1:60*19+1));
legend("Corresponding layer","upper water of tank 1")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
title("Root Mean Squared Error: ", err1)
grid on
subplot(3,1,2)
plot(time_data, T_results(:,num_layer_tank1+1),time_data_8min, T_upper_water_tank2((15*3+8)+1:(9*15+8)+1));
legend("Corresponding layer","upper water of tank 2")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
title("Root Mean Squared Error: ", err2)
grid on
subplot(3,1,3)
plot(time_data, T_results(:,num_layer_tank1+num_layer_tank2),time_data_1min, T_buttom_water_tank2(60*7+1:60*19+1));
xlabel("Time (0=07:00, 12=19:00)")
ylabel("Temp (°C)")
legend("Corresponding layer","buttom water of tank 2")
title("Root Mean Squared Error: ", err3)
grid on
set(gcf,'position',[400,300,1000,600])
til = sprintf('Comparison of Tempertures (2023-05-01) (layers: %d + %d)',num_layer_tank1 , num_layer_tank2);
sgtitle(til)
% subplot(3,1,3)
% yyaxis left
% plot(time_data, T_results(:,6),time_data_1min, T_buttom_water_tank2(60*7+1:60*19+1));
% xlabel("Time (0=07:00, 12=19:00)")
% ylabel("Temp (°C)")
% yyaxis right
% plot(time_data_1min,m_p_dot_data,time_data,m_s_dot_interp)
% legend("layer 5","buttom water of tank 2")
% grid on
% set(gcf,'position',[400,300,1000,600])
% sgtitle('Comparison of Tempertures ')




%% flow rate interpolation visulization
% m_s_dot_interp_next = interp1(time_data_1h, m_s_dot_data, time_span,'next');
% figure;
% yyaxis left
% plot(time_data_1min, T_buttom_water_tank2(60*7+1:60*19+1),'LineWidth',2);
% xlabel("Time (0=07:00, 12=19:00)")
% ylabel("Temp (°C)")
% yyaxis right
% plot(time_data,m_s_dot_interp_next,"LineWidth",2)
% hold on
% plot(time_data,m_s_dot_interp,time_data_1min,FR_water_hp(60*7+1:60*19+1)*1000,'LineWidth',2)
% xlabel("Time (0=07:00, 12=19:00)")
% ylabel("Flow rate (L/h)")
% legend("Temp of buttom water of tank 2","Hour-level flow rate of supplied water","Interpolated flow rate of supplied water","Flow rate of heat pump")
% grid on
% hold off
% set(gcf,'position',[400,300,1000,600])
% title('Flow rate of supplied water (2023-05-01)')
% 
% %%
% figure;
% bar(0:12,V_water(8:20)*1000)
% hold on
% plot(time_data,m_s_dot_interp,"LineWidth",2)
% legend("Hour-level flow rate of supplied water","Interpolated flow rate of supplied water")
% grid on
% hold off
% xlabel("Time (0=07:00, 12=19:00)")
% ylabel("Flow rate (L/h)")
% set(gcf,'position',[400,300,1000,600])
% title('Flow rate of supplied water (2023-05-01)')
% 
% %%
% figure;
% % plot(0:12,V_water(8:20)*1000)
% plot(time_data,m_s_dot_interp_next,"LineWidth",2)
% hold on
% plot(time_data,m_s_dot_interp,"LineWidth",2)
% legend("Hour-level flow rate of supplied water","Interpolated flow rate of supplied water")
% grid on
% hold off
% xlabel("Time (0=07:00, 12=19:00)")
% ylabel("Flow rate (L/h)")
% set(gcf,'position',[400,300,1000,600])
% title('Flow rate of supplied water (2023-05-01)')
