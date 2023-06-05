%% Plot
%% Comparison of Temperatures
function plot_comparison_T_heatloss(time_data,T_results,time_data_1min,time_data_8min,num_layer_tank1,num_layer_tank2,T_upper_tank1_data,T_upper_tank2_data,T_bottom_tank2_data)



T_upper_tank1_interp = interp1(time_data_1min,T_upper_tank1_data,time_data);
T_upper_tank2_interp = interp1(time_data_8min,T_upper_tank2_data,time_data,"spline");
T_bottom_tank2_interp = interp1(time_data_1min,T_bottom_tank2_data,time_data);
T_ground_truth = [T_upper_tank1_interp',T_upper_tank2_interp',T_bottom_tank2_interp']; % Your historical data here


err1 = sqrt(mean((T_results(:,1) - T_ground_truth(:,1)).^2));               %Root Mean Squared Error
err2 = sqrt(mean((T_results(:,num_layer_tank1+1) - T_ground_truth(:,2)).^2));               %Root Mean Squared Error
err3 = sqrt(mean((T_results(:,num_layer_tank1+num_layer_tank2) - T_ground_truth(:,3)).^2));               %Root Mean Squared Error
figure;
subplot(3,1,1)
plot(time_data, T_results(:,1),time_data_1min, T_upper_tank1_data);
legend("layer 1","upper water of tank 1")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
xlim([time_data(1),time_data(end)])
title("Root Mean Squared Error: ", err1)
grid on
subplot(3,1,2)
plot(time_data, T_results(:,num_layer_tank1+1),time_data_8min, T_upper_tank2_data);
legend("layer 3","upper water of tank 2")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
xlim([time_data(1),time_data(end)])
title("Root Mean Squared Error: ", err2)
grid on
subplot(3,1,3)
plot(time_data, T_results(:,num_layer_tank1+num_layer_tank2),time_data_1min, T_bottom_tank2_data);
xlabel("Time (0=07:00, 12=19:00)")
ylabel("Temp (°C)")
xlim([time_data(1),time_data(end)])
legend("layer 6","buttom water of tank 2")
title("Root Mean Squared Error: ", err3)
grid on
set(gcf,'position',[400,300,1000,600])
sgtitle('Comparison of Tempertures (2023-05-01)')
end