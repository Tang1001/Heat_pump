function T_results_0_12=simulation_Euler_Nperiod(period_time,time_data,sample_time,num_layer_tank1,num_layer_tank2,m_layer_tank1, m_layer_tank2, R_tank12, m_c_dot, diff_T_c, T_s, T_in_data,m_p_dot_data,T_bottom_tank2_data,T_upper_tank1_data,T_upper_tank2_data,m_s_dot_data)

sample_num = period_time/(sample_time/60);

[m_p_dot_interp,T_in_interp,m_s_dot_interp,T_upper_tank1_interp,T_upper_tank2_interp,...
    T_bottom_tank2_interp]=interpolation_data(time_data,sample_time,m_p_dot_data,T_in_data,m_s_dot_data,T_upper_tank1_data,T_upper_tank2_data,T_bottom_tank2_data);


time_data_period = zeros((sample_num+1),12/period_time);
m_p_dot_interp_period = zeros(sample_num,12/period_time);
T_in_interp_period = zeros(sample_num,12/period_time);
m_s_dot_interp_period = zeros(sample_num,12/period_time);

T_upper_tank1_interp_period = zeros(sample_num,12/period_time);
T_upper_tank2_interp_period = zeros(sample_num,12/period_time);
T_bottom_tank2_interp_period = zeros(sample_num,12/period_time);


for i=1:12/period_time
    time_data_period(:,i) = [period_time*(i-1):sample_time/60:period_time*i];
    m_p_dot_interp_period(:,i) = m_p_dot_interp(1+(i-1)*sample_num:i*sample_num);
    T_in_interp_period(:,i) = T_in_interp(1+(i-1)*sample_num:i*sample_num);
    m_s_dot_interp_period(:,i) = m_s_dot_interp(1+(i-1)*sample_num:i*sample_num);

    T_upper_tank1_interp_period(:,i) =  T_upper_tank1_interp(1+(i-1)*sample_num:i*sample_num);
T_upper_tank2_interp_period(:,i) =  T_upper_tank2_interp(1+(i-1)*sample_num:i*sample_num);
T_bottom_tank2_interp_period(:,i) =  T_bottom_tank2_interp(1+(i-1)*sample_num:i*sample_num);


end




T_results_period = zeros(sample_num+1,num_layer_tank1+num_layer_tank2,12/period_time);
for i=1:12/period_time
T_initial=initial_T_setting(num_layer_tank1,num_layer_tank2,m_layer_tank1,m_layer_tank2,T_upper_tank1_interp_period(1,i),T_upper_tank2_interp_period(1,i),T_bottom_tank2_interp_period(1,i));
T_results_period(:,:,i) = system_of_equations_Euler_Nlayers(time_data_period(:,i), T_initial, m_p_dot_interp_period(:,i), T_in_interp_period(:,i), m_s_dot_interp_period(:,i), num_layer_tank1,num_layer_tank2, m_layer_tank1, m_layer_tank2, R_tank12, m_c_dot, diff_T_c, T_s);


end


T_results_0_12 = [];
for i=1:12/period_time
if i~=12/period_time
    T_results_0_12=[T_results_0_12;T_results_period(1:end-1,:,i)];
else
    T_results_0_12=[T_results_0_12;T_results_period(:,:,i)];
end
end


end
