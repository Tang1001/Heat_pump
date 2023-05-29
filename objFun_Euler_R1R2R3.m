function err = objFun_Euler_R1R2R3(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, m1, m2, m3, m4, cp, m_c_dot, diff_T_c, T_s, params, data)
    R1 = params(1);
    R2 = params(2);
    R3 = params(3);
    
    % Solve the system of equations with the current R1 and R2
    T_results = system_of_equations_Euler_4layers_R1R2R3(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, m1, m2, m3, m4, cp, R1, R2, R3, m_c_dot, diff_T_c, T_s);


    % Compare the model output with the actual data
%     err = mean((T(:,1) - data(:,1)).^2);    %Mean Squared Error
%     err = sqrt(mean((T_results(:,1) - data(:,1)).^2));    %Root Mean Squared Error
    err = sum(sqrt(mean((T_results(:,[1,3,4]) - data).^2)));               %Root Mean Squared Error
end