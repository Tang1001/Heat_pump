function [m_p_dot_interp,T_in_interp,m_s_dot_interp,T_upper_water_tank1_interp,T_upper_water_tank2_interp,...
    T_buttom_water_tank2_interp]=interpolation_data(time_data,sample_time,m_p_dot_data,T_in_data,m_s_dot_data,T_upper_water_tank1_data,T_upper_water_tank2_data,T_buttom_water_tank2_data) 
% Interpolation of data


%% Known input data
%2023-03-27 7:00 ~ 19:00
time_data_1min=0:1/60:12;
time_data_1h=0:12;
time_data_8min=0-4/60:8/60:12+4/60;


% Interpolation
m_p_dot_interp = interp1(time_data_1min, m_p_dot_data, time_data);
T_in_interp = interp1(time_data_1min, T_in_data, time_data);

% estimate min-level flow rate of supplying water
% m_s_dot_interp = interp1(time_data_1h, m_s_dot_data, time_span,'next');

peak_shift = 0.5; 
time_data_1h_interp = time_data_1h-peak_shift; %make the value in the middle of each hour
m_s_dot_interp = interp1(time_data_1h_interp, m_s_dot_data, time_data,'spline');

% Adjust the minute-level data to make the total sums equal
total_hourly = sum(m_s_dot_data(2:end));
total_minute = sum(m_s_dot_interp(2:end))*(sample_time/60); % Don't divide by 60

difference = total_hourly - total_minute;
m_s_dot_interp_adjusted = m_s_dot_interp;
m_s_dot_interp_adjusted(2:end) = m_s_dot_interp(2:end) + (difference/(sample_time/60))/ (length(m_s_dot_interp)-1);
% Calculate the new total for the minute-level data
% total_minute_adjusted = sum(m_s_dot_interp_adjusted(2:end))*(sample_time/60); % Don't divide by 60
% disp(['Total water consumption (Hourly data): ', num2str(total_hourly), ' liters']);
% disp(['Total water consumption (Minute data, adjusted): ', num2str(total_minute_adjusted), ' liters']);
m_s_dot_interp=m_s_dot_interp_adjusted;



T_upper_water_tank1_interp = interp1(time_data_1min,T_upper_water_tank1_data,time_data);
T_upper_water_tank2_interp = interp1(time_data_8min,T_upper_water_tank2_data,time_data);
T_buttom_water_tank2_interp = interp1(time_data_1min,T_buttom_water_tank2_data,time_data);



end


