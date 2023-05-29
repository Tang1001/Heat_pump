%% Plot
%% Comparison of Temperatures
err1 = sqrt(mean((T_results(:,1) - data(:,1)).^2));               %Root Mean Squared Error
err2 = sqrt(mean((T_results(:,num_layer_tank1+1) - data(:,2)).^2));               %Root Mean Squared Error
err3 = sqrt(mean((T_results(:,num_layer_tank1+num_layer_tank2) - data(:,3)).^2));               %Root Mean Squared Error
figure;
subplot(3,1,1)
plot(time_data, T_results(:,1),time_data_1min, T_upper_water_tank1(60*7+1:60*19+1));
legend("layer 1","upper water of tank 1")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
title("Root Mean Squared Error: ", err1)
grid on
subplot(3,1,2)
plot(time_data, T_results(:,num_layer_tank1+1),time_data_8min, T_upper_water_tank2((15*3+8):(9*15+8)+1));
legend("layer 3","upper water of tank 2")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
xlim([0 12])
title("Root Mean Squared Error: ", err2)
grid on
subplot(3,1,3)
plot(time_data, T_results(:,num_layer_tank1+num_layer_tank2),time_data_1min, T_buttom_water_tank2(60*7+1:60*19+1));
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