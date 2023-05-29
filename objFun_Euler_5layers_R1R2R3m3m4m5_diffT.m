function err = objFun_Euler_5layers_R1R2R3m3m4m5_diffT(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, m1, m2, cp, m_c_dot, T_s, params, data)
    R12 = params(1);
    R34 = params(2);
    R45 = params(3);
    
    m3 = params(4);
    m4 = params(5);
    m5 = params(6);
    diff_T_c = params(7);
    
    
    % Solve the system of equations with the current R1 and R2
    T_results = system_of_equations_Euler_5layers(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, m1, m2, m3, m4, m5, cp, R12, R34, R45, m_c_dot, diff_T_c, T_s);


    % Compare the model output with the actual data
%     err = mean((T(:,1) - data(:,1)).^2);    %Mean Squared Error
%     err = sqrt(mean((T_results(:,1) - data(:,1)).^2));    %Root Mean Squared Error
    err = sum(sqrt(mean((T_results(:,[1,3,5]) - data).^2)));               %Root Mean Squared Error
end