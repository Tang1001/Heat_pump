%% Plot
%% Comparison of Temperatures
function plot_comparison_T_period(time_data,T_results,period_time,num_layer_tank1,num_layer_tank2,T_upper_tank1_data,T_upper_tank2_data,T_bottom_tank2_data)

% time_data_1min=0:1/60:12;
% time_data_8min=0-4/60:8/60:12+4/60;
time_data_1min=time_data(1):1/60:time_data(end);
% time_data_1h=floor(time_data(1)):ceil(time_data(end));
time_data_8min=(round((time_data(1)*60+60)/8)*8-60)/60:8/60:11+round((time_data(end)-11)*60/8)*8/60;  %7:04 ~ 18:32


T_upper_tank1_interp = interp1(time_data_1min,T_upper_tank1_data,time_data);
T_upper_tank2_interp = interp1(time_data_8min,T_upper_tank2_data,time_data);
T_bottom_tank2_interp = interp1(time_data_1min,T_bottom_tank2_data,time_data);
T_ground_truth = [T_upper_tank1_interp',T_upper_tank2_interp',T_bottom_tank2_interp']; % Your historical data here

% size(T_results(:,1))
% size(T_ground_truth(:,1))
err1 = sqrt(mean((T_results(:,1) - T_ground_truth(:,1)).^2));               %Root Mean Squared Error
err2 = sqrt(mean((T_results(:,num_layer_tank1+1) - T_ground_truth(:,2)).^2));               %Root Mean Squared Error
err3 = sqrt(mean((T_results(:,num_layer_tank1+num_layer_tank2) - T_ground_truth(:,3)).^2));  
figure;
subplot(3,1,1)
plot(time_data, T_results(:,1),time_data_1min, T_upper_tank1_data);
for i=1:ceil(12/period_time)+1
line([(i-1)*period_time,(i-1)*period_time],[0,100],'color','b','LineWidth',2,'linestyle','--')
end
legend("layer 1","upper water of tank 1")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
xlim([time_data(1) time_data(end)])
title("Root Mean Squared Error: ", err1)
grid on
subplot(3,1,2)
plot(time_data, T_results(:,num_layer_tank1+1),time_data_8min, T_upper_tank2_data);
for i=1:ceil(12/period_time)+1
line([(i-1)*period_time,(i-1)*period_time],[0,100],'color','b','LineWidth',2,'linestyle','--')
end
legend("layer 3","upper water of tank 2")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
xlim([time_data(1) time_data(end)])
title("Root Mean Squared Error: ", err2)
grid on
subplot(3,1,3)
plot(time_data, T_results(:,num_layer_tank1+num_layer_tank2),time_data_1min, T_bottom_tank2_data);
for i=1:ceil(12/period_time)+1
line([(i-1)*period_time,(i-1)*period_time],[0,100],'color','b','LineWidth',2,'linestyle','--')
end
xlabel("Time (0=07:00, 12=19:00)")
ylabel("Temp (°C)")
xlim([time_data(1) time_data(end)])
legend("layer 6","buttom water of tank 2")
title("Root Mean Squared Error: ", err3)
grid on
set(gcf,'position',[400,300,1000,600])
sgtitle('Comparison of Tempertures (2023-05-01)')
end