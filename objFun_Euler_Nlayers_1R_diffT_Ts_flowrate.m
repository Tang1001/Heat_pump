function err = objFun_Euler_Nlayers_1R_diffT_Ts_flowrate(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data,  num_layer_tank1, num_layer_tank2, m_layer_tank1,m_layer_tank2, R_thermal_tank1, R_thermal_tank2, cp, H, A, k, params, data)
    R12 = params(1);
    diff_T_c = params(2);
    T_s = params(3);
    m_c_dot = params(4);
    
    
    % Solve the system of equations with the current R1 and R2
%     T_results = system_of_equations_Euler_6layers_5R(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, m1, m2, m3, m4, m5, m6, cp, R12, R23, R34, R45, R56, m_c_dot, diff_T_c, T_s);
T_results = system_of_equations_Euler_Nlayers_singlePipe(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, num_layer_tank1, num_layer_tank2, m_layer_tank1,m_layer_tank2, R_thermal_tank1, R_thermal_tank2, R12, cp, m_c_dot, diff_T_c, T_s);


    % Compare the model output with the actual data
%     err = mean((T(:,1) - data(:,1)).^2);    %Mean Squared Error
%     err = sqrt(mean((T_results(:,1) - data(:,1)).^2));    %Root Mean Squared Error
%     err = sqrt(mean((T_results(:,6) - data(:,3)).^2));    %Root Mean Squared Error
%     err = sum(sqrt(mean((T_results(:,[1,6]) - data(:,[1,3])).^2)));               %Root Mean Squared Error
    err = sum(sqrt(mean((T_results(:,[1,num_layer_tank1+1,num_layer_tank1+num_layer_tank2]) - data).^2)));               %Root Mean Squared Error
end