function T_results=simulation_Euler_Nlayer(time_data,sample_time,num_layer_tank1,num_layer_tank2,m_layer_tank1, m_layer_tank2, R_tank12, m_c_dot, diff_T_c, T_s, T_in_data,m_p_dot_data,T_bottom_tank2_data,T_upper_tank1_data,T_upper_tank2_data,m_s_dot_data)


%% Interpolation
[m_p_dot_interp,T_in_interp,m_s_dot_interp,T_upper_tank1_interp,T_upper_tank2_interp,...
    T_bottom_tank2_interp]=interpolation_data(time_data,sample_time,m_p_dot_data,T_in_data,m_s_dot_data,T_upper_tank1_data,T_upper_tank2_data,T_bottom_tank2_data);


% Initial Temperature Setting

T_initial=initial_T_setting(num_layer_tank1,num_layer_tank2,m_layer_tank1,m_layer_tank2,T_upper_tank1_interp(1),T_upper_tank2_interp(1),T_bottom_tank2_interp(1));







%% Euler method
% T_results = system_of_equations_Euler_6layers_5R_singlePipe(time_data, T_initial, m_p_dot_interp, T_in_interp, m_s_dot_interp, m1, m2, m3, m4, m5, m6, cp, R12, R_tank12, R34, R45, R56, m_c_dot, diff_T_c, T_s);
T_results = system_of_equations_Euler_Nlayers(time_data, T_initial, m_p_dot_interp, T_in_interp, m_s_dot_interp, num_layer_tank1,num_layer_tank2, m_layer_tank1, m_layer_tank2, R_tank12, m_c_dot, diff_T_c, T_s);


end