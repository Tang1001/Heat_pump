function T_results=system_simulation_Euler_period(period_time, time_data,sample_time,num_layer_tank1,num_layer_tank2,T_out_water_he,FR_water_hp,T_buttom_water_tank2,T_upper_water_tank1,T_upper_water_tank2,V_water)

%4hour period
sample_num = period_time/(sample_time/60);

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



%% Interpolation of input data
%2023-03-27 7:00 ~ 19:00
%1min
m_p_dot_data=FR_water_hp(60*7+1:60*19+1)*1000;  %[L/h] 60*7+1~60*19
T_in_data=T_out_water_he(60*7+1:60*19+1);

%1h
m_s_dot_data=V_water(1*8:1*20)*1000;    %[L/h] 07:00~19:00

time_data_1min=0:1/60:12;
time_data_1h=0:12;
time_data_8min=0:8/60:12;


% Interpolation
m_p_dot_interp = interp1(time_data_1min, m_p_dot_data, time_data);
T_in_interp = interp1(time_data_1min, T_in_data, time_data);


peak_shift = 0.5; 
time_data_1h_interp = time_data_1h-peak_shift; %make the value in the middle of each hour
m_s_dot_interp = interp1(time_data_1h_interp, m_s_dot_data, time_data,'spline');
% m_s_dot_interp = interp1(time_data_1h, m_s_dot_data, time_span,'spline');
% m_s_dot_interp = interp1(time_data_1h, m_s_dot_data, time_span,'nearest');

% Adjust the minute-level data to make the total sums equal
total_hourly = sum(m_s_dot_data(2:end));
total_minute = sum(m_s_dot_interp(2:end))*(sample_time/60); % Don't divide by 60

difference = total_hourly - total_minute;
m_s_dot_interp_adjusted = m_s_dot_interp;
m_s_dot_interp_adjusted(2:end) = m_s_dot_interp(2:end) + (difference/(sample_time/60))/ (length(m_s_dot_interp)-1);
m_s_dot_interp=m_s_dot_interp_adjusted;


%% split data

time_data_period = zeros(12/period_time,(sample_num+1));
m_p_dot_interp_period = zeros(12/period_time,sample_num);
T_in_interp_period = zeros(12/period_time,sample_num);
m_s_dot_interp_period = zeros(12/period_time,sample_num);

for i=1:12/period_time
    time_data_period(i,:) = [period_time*(i-1):sample_time/60:period_time*i];
    m_p_dot_interp_period(i,:) = m_p_dot_interp(1+(i-1)*sample_num:i*sample_num);
    T_in_interp_period(i,:) = T_in_interp(1+(i-1)*sample_num:i*sample_num);
    m_s_dot_interp_period(i,:) = m_s_dot_interp(1+(i-1)*sample_num:i*sample_num);
end




%%
%7:00 %1min:60*7+1=421  8min:60/8*7+1=53.5
T_results_period = zeros(sample_num+1,6,12/period_time);


for i=1:12/period_time

T1_0=T_upper_water_tank1(60*(7+(i-1)*period_time)+1); % layer 1 initial temperature
T3_0=(T_upper_water_tank2(floor(60/8*(7+(i-1)*period_time)+1))+T_upper_water_tank2(ceil((60/8*(7+(i-1)*period_time)+1))))/2; % layer 3 initial temperature
T6_0=T_buttom_water_tank2(60*(7+(i-1)*period_time)+1); % layer 6 initial temperature

d_T1 = (T1_0 - T3_0) * (R1 / R1_tot);
T2_0 = T1_0 - d_T1;     % layer 2 initial temperature
d_T3 = (T3_0 - T6_0) * (R3 / R2_tot);
T4_0= T3_0 - d_T3; % layer 4 initial temperature
d_T4 = (T3_0 - T6_0) * (R4 / R2_tot);
T5_0= T4_0 - d_T4; % layer 5 initial temperature


T = [T1_0; T2_0; T3_0; T4_0; T5_0; T6_0];
T_results_period(:,:,i) = system_of_equations_Euler_6layers_5R_singlePipe(time_data_period(i,:), T, m_p_dot_interp_period(i,:), T_in_interp_period(i,:), m_s_dot_interp_period(i,:), m1, m2, m3_optimal, m4_optimal, m5_optimal, m6_optimal, cp, R12, R23_optimal, R34_optimal, R45_optimal, R56_optimal, m_c_dot_optimal, diff_T_c_optimal, T_s_optimal);

end


%%

T_results_0_12 = [];
for i=1:12/period_time
% reshapeT_results_period(:,end-1,i)
% T_results_period(:,:,i)
if i~=12/period_time
    T_results_0_12=[T_results_0_12;T_results_period(1:end-1,:,i)];
else
    T_results_0_12=[T_results_0_12;T_results_period(:,:,i)];
end
end


% T_results_0_12=[T_results_0_4(1:end-1,:);T_results_4_8(1:end-1,:);T_results_8_12];
err1 = sqrt(mean((T_results_0_12(:,1) - data(:,1)).^2));               %Root Mean Squared Error
err2 = sqrt(mean((T_results_0_12(:,3) - data(:,2)).^2));               %Root Mean Squared Error
err3 = sqrt(mean((T_results_0_12(:,6) - data(:,3)).^2));  
figure;
subplot(3,1,1)
plot(time_data, T_results_0_12(:,1),time_data_1min, T_upper_water_tank1(60*7+1:60*19+1));
for i=1:12/period_time-1
line([i*period_time,i*period_time],[0,100],'color','b','LineWidth',2,'linestyle','--')
end
legend("layer 1","upper water of tank 1")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
title("Root Mean Squared Error: ", err1)
grid on
subplot(3,1,2)
plot(time_data, T_results_0_12(:,3),time_data_8min, T_upper_water_tank2((15*3+8)+1:(9*15+8)+1));
for i=1:12/period_time-1
line([i*period_time,i*period_time],[0,100],'color','b','LineWidth',2,'linestyle','--')
end
legend("layer 3","upper water of tank 2")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
title("Root Mean Squared Error: ", err2)
grid on
subplot(3,1,3)
plot(time_data, T_results_0_12(:,6),time_data_1min, T_buttom_water_tank2(60*7+1:60*19+1));
for i=1:12/period_time-1
line([i*period_time,i*period_time],[0,100],'color','b','LineWidth',2,'linestyle','--')
end
xlabel("Time (0=07:00, 12=19:00)")
ylabel("Temp (°C)")
legend("layer 6","buttom water of tank 2")
title("Root Mean Squared Error: ", err3)
grid on
set(gcf,'position',[400,300,1000,600])
sgtitle('Comparison of Tempertures (2023-05-01)')


