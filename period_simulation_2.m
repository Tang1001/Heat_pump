%% 7:00~19:00 0~12


%%4hour period
period_time = 3;    %[h]
sample_num = period_time/(sample_time/60);

% time_data_period
% m_p_dot_interp_period
% T_in_interp_period
% m_s_dot_interp_period(i)
time_data_period = zeros(12/period_time,(sample_num+1));
m_p_dot_interp_period = zeros(12/period_time,sample_num);
T_in_interp_period = zeros(12/period_time,sample_num);
m_s_dot_interp_period = zeros(12/period_time,sample_num);

for i=1:12/period_time
    time_data_period(i,:) = [t_start+period_time*(i-1):sample_time/60:t_start+period_time*i];
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


T = [T1_0; T2_0; T3_0; T4_0; T5_0; T6_0]
% T_results_period(:,:,i) = system_of_equations_Euler_6layers_5R_singlePipe(time_data_period(i,:), T, m_p_dot_interp_period(i,:), T_in_interp_period(i,:), m_s_dot_interp_period(i,:), m1, m2, m3_optimal, m4_optimal, m5_optimal, m6_optimal, cp, R12, R23_optimal, R34_optimal, R45_optimal, R56_optimal, m_c_dot_optimal, diff_T_c_optimal, T_s_optimal);
T_results_period(:,:,i) = system_of_equations_Euler_6layers_5R_singlePipe(time_data_period(i,:), T, m_p_dot_interp_period(i,:), T_in_interp_period(i,:), m_s_dot_interp_period(i,:), m1, m2, m3, m4, m5, m6, cp, R12, R23, R34, R45, R56, m_c_dot, diff_T_c, T_s);

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
plot(time_data, T_results_0_12(:,3),time_data_8min, T_upper_water_tank2((15*3+8):(9*15+8)+1));
for i=1:12/period_time-1
line([i*period_time,i*period_time],[0,100],'color','b','LineWidth',2,'linestyle','--')
end
legend("layer 3","upper water of tank 2")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
title("Root Mean Squared Error: ", err2)
xlim([0 12])
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

