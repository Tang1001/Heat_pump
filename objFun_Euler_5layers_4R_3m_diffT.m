function err = objFun_Euler_5layers_4R_3m_diffT(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, m1, m2, cp, m_c_dot, T_s, params, data)
    R12 = params(1);
    R23 = params(2);
    R34 = params(3);
    R45 = params(4);
    
    m3 = params(5);
    m4 = params(6);
    m5 = params(7);
    diff_T_c = params(8);
    
    
    % Solve the system of equations with the current R1 and R2
    T_results = system_of_equations_Euler_5layers_4R(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, m1, m2, m3, m4, m5, cp, R12, R23, R34, R45, m_c_dot, diff_T_c, T_s);


    % Compare the model output with the actual data
%     err = mean((T(:,1) - data(:,1)).^2);    %Mean Squared Error
%     err = sqrt(mean((T_results(:,1) - data(:,1)).^2));    %Root Mean Squared Error
    err = sum(sqrt(mean((T_results(:,[1,3,5]) - data).^2)));               %Root Mean Squared Error
end