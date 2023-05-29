%% 7:00~19:00 0~12

%%4hour period
period_time = 4;    %[h]


time_data_0_4 = [t_start:sample_time/60:t_start+period_time];
time_data_4_8 = [t_start+period_time:sample_time/60:t_start+period_time*2];
time_data_8_12 = [t_start+period_time*2:sample_time/60:t_start+period_time*3];


m_p_dot_interp_0_4 = m_p_dot_interp(1:4*60/sample_time);
m_p_dot_interp_4_8 = m_p_dot_interp(4*60/sample_time+1:8*60/sample_time);
m_p_dot_interp_8_12 = m_p_dot_interp(8*60/sample_time+1:12*60/sample_time);

T_in_interp_0_4 = T_in_interp(1:4*60/sample_time);
T_in_interp_4_8 = T_in_interp(4*60/sample_time+1:8*60/sample_time);
T_in_interp_8_12 = T_in_interp(8*60/sample_time+1:12*60/sample_time);

m_s_dot_interp_0_4 = m_s_dot_interp(1:4*60/sample_time);
m_s_dot_interp_4_8 = m_s_dot_interp(4*60/sample_time+1:8*60/sample_time);
m_s_dot_interp_8_12 = m_s_dot_interp(8*60/sample_time+1:12*60/sample_time);



%%
%7:00 %1min:60*7+1=421  8min:60/8*7+1=53.5
T1_0=T_upper_water_tank1(421); % layer 1 initial temperature
T3_0=(T_upper_water_tank2(53)+T_upper_water_tank2(54))/2; % layer 3 initial temperature
T6_0=T_buttom_water_tank2(421); % layer 6 initial temperature


d_T1 = (T1_0 - T3_0) * (R1 / R1_tot);
T2_0 = T1_0 - d_T1;     % layer 2 initial temperature
d_T3 = (T3_0 - T6_0) * (R3 / R2_tot);
T4_0= T3_0 - d_T3; % layer 4 initial temperature
d_T4 = (T3_0 - T6_0) * (R4 / R2_tot);
T5_0= T4_0 - d_T4; % layer 5 initial temperature


T = [T1_0; T2_0; T3_0; T4_0; T5_0; T6_0];
T_results_0_4 = system_of_equations_Euler_6layers_5R(time_data_0_4, T, m_p_dot_interp_0_4, T_in_interp_0_4, m_s_dot_interp_0_4, m1, m2, m3_optimal, m4_optimal, m5_optimal, m6_optimal, cp, R12, R23_optimal, R34_optimal, R45_optimal, R56_optimal, m_c_dot_optimal, diff_T_c_optimal, T_s_optimal);

%%
%(7+4):00 %1min:60*(7+4)+1=421  8min:60/8*(7+4)+1=83.5
T1_4=T_upper_water_tank1(661); % layer 1 initial temperature
T3_4=(T_upper_water_tank2(83)+T_upper_water_tank2(84))/2; % layer 3 initial temperature
T6_4=T_buttom_water_tank2(661); % layer 6 initial temperature


d_T1 = (T1_4 - T3_4) * (R1 / R1_tot);
T2_4 = T1_4 - d_T1;     % layer 2 initial temperature
d_T3 = (T3_4 - T6_4) * (R3 / R2_tot);
T4_4= T3_4 - d_T3; % layer 4 initial temperature
d_T4 = (T3_4 - T6_4) * (R4 / R2_tot);
T5_4= T4_4 - d_T4; % layer 5 initial temperature

T = [T1_4; T2_4; T3_4; T4_4; T5_4; T6_4];
T_results_4_8 = system_of_equations_Euler_6layers_5R(time_data_4_8, T, m_p_dot_interp_4_8, T_in_interp_4_8, m_s_dot_interp_4_8, m1, m2, m3_optimal, m4_optimal, m5_optimal, m6_optimal, cp, R12, R23_optimal, R34_optimal, R45_optimal, R56_optimal, m_c_dot_optimal, diff_T_c_optimal, T_s_optimal);


%%
%(7+8):00 %1min:60*(7+8)+1=901  8min:60/8*(7+8)+1=113.5
T1_8=T_upper_water_tank1(901); % layer 1 initial temperature
T3_8=(T_upper_water_tank2(113)+T_upper_water_tank2(114))/2; % layer 3 initial temperature
T6_8=T_buttom_water_tank2(901); % layer 5 initial temperature


d_T1 = (T1_8 - T3_8) * (R1 / R1_tot);
T2_8 = T1_8 - d_T1;     % layer 2 initial temperature
d_T3 = (T3_8 - T6_8) * (R3 / R2_tot);
T4_8= T3_8 - d_T3; % layer 4 initial temperature
d_T4 = (T3_8 - T6_8) * (R4 / R2_tot);
T5_8= T4_8 - d_T4; % layer 4 initial temperature

T = [T1_8; T2_8; T3_8; T4_8; T5_8; T6_8];
T_results_8_12 = system_of_equations_Euler_6layers_5R(time_data_8_12, T, m_p_dot_interp_8_12, T_in_interp_8_12, m_s_dot_interp_8_12, m1, m2, m3_optimal, m4_optimal, m5_optimal, m6_optimal, cp, R12, R23_optimal, R34_optimal, R45_optimal, R56_optimal, m_c_dot_optimal, diff_T_c_optimal, T_s_optimal);


%%
T_results_0_12=[T_results_0_4(1:end-1,:);T_results_4_8(1:end-1,:);T_results_8_12];
err1 = sqrt(mean((T_results_0_12(:,1) - data(:,1)).^2));               %Root Mean Squared Error
err2 = sqrt(mean((T_results_0_12(:,3) - data(:,2)).^2));               %Root Mean Squared Error
err3 = sqrt(mean((T_results_0_12(:,6) - data(:,3)).^2));  
figure;
subplot(3,1,1)
plot(time_data, T_results_0_12(:,1),time_data_1min, T_upper_water_tank1(60*7+1:60*19+1));
line([4,4],[0,100],'color','b','LineWidth',2,'linestyle','--')
line([8,8],[0,100],'color','b','LineWidth',2,'linestyle','--')
legend("layer 1","upper water of tank 1")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
title("Root Mean Squared Error: ", err1)
grid on
subplot(3,1,2)
plot(time_data, T_results_0_12(:,3),time_data_8min, T_upper_water_tank2((15*3+8)+1:(9*15+8)+1));
line([4,4],[0,80],'color','b','LineWidth',2,'linestyle','--')
line([8,8],[0,80],'color','b','LineWidth',2,'linestyle','--')
legend("layer 3","upper water of tank 2")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
title("Root Mean Squared Error: ", err2)
grid on
subplot(3,1,3)
plot(time_data, T_results_0_12(:,6),time_data_1min, T_buttom_water_tank2(60*7+1:60*19+1));
line([4,4],[0,80],'color','b','LineWidth',2,'linestyle','--')
line([8,8],[0,80],'color','b','LineWidth',2,'linestyle','--')
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

