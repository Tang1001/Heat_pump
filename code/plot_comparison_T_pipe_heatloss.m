%% Plot
%% Comparison of Temperatures
function plot_comparison_T_pipe_heatloss(time_data,T_results,time_data_1min,T_pipe)

T_pipe_interp = interp1(time_data_1min,T_pipe,time_data)';

err1 = sqrt(mean((T_results - T_pipe_interp).^2))              %Root Mean Squared Error
figure;
plot(time_data, T_results,time_data_1min, T_pipe);
legend("predict","measurement")
xlabel("Time (0=07:00; 12=19:00)")
ylabel("Temp (°C)")
xlim([time_data(1),time_data(end)])
title("Root Mean Squared Error: ", err1)
grid on
% title('Comparison of Tempertures (2023-05-01)')
end