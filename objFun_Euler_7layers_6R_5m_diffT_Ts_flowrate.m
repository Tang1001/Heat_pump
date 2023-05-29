function err = objFun_Euler_7layers_6R_5m_diffT_Ts_flowrate(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, m1, m2, cp, params, data)
    R12 = params(1);
    R23 = params(2);
    R34 = params(3);
    R45 = params(4);
    R56 = params(5);
    R67 = params(6);
    
    m3 = params(7);
    m4 = params(8);
    m5 = params(9);
    m6 = params(10);
    m7 = params(11);
    diff_T_c = params(12);
    T_s = params(13);
    m_c_dot = params(14);
    
    
    % Solve the system of equations with the current R1 and R2
    T_results = system_of_equations_Euler_7layers_6R(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, m1, m2, m3, m4, m5, m6, m7, cp, R12, R23, R34, R45, R56, R67, m_c_dot, diff_T_c, T_s);


    % Compare the model output with the actual data
%     err = mean((T(:,1) - data(:,1)).^2);    %Mean Squared Error
%     err = sqrt(mean((T_results(:,1) - data(:,1)).^2));    %Root Mean Squared Error
%     err = sqrt(mean((T_results(:,6) - data(:,3)).^2));    %Root Mean Squared Error
    err = sum(sqrt(mean((T_results(:,[1,3,6]) - data).^2)));               %Root Mean Squared Error
end