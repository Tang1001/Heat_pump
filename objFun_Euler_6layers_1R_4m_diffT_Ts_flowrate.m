function err = objFun_Euler_6layers_1R_4m_diffT_Ts_flowrate(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, m1, m2, cp, R12, H, A, k, params, data)
    R23 = params(1);
    
    m3 = params(2);
    m4 = params(3);
    m5 = params(4);
    m6 = params(5);
    diff_T_c = params(6);
    T_s = params(7);
    m_c_dot = params(8);
    
    

    % height of layer    [m]
    h3=H*(m3/500);
    h4=H*(m4/500);
    h5=H*(m5/500);
    h6=H*(m6/500);

    % thermal resistances    R_i = h_i / (k_i * A)
    R3=h3/(k*A);
    R4=h4/(k*A);
    R5=h5/(k*A);
    R6=h6/(k*A);


    % initial thermal efficiency between layers
    R34 = 1/(2 * (R3 * R4) / (R3 + R4));
    R45 = 1/(2 * (R4 * R5) / (R4 + R5));
    R56 = 1/(2 * (R5 * R6) / (R5 + R6));

    
    % Solve the system of equations with the current R1 and R2
%     T_results = system_of_equations_Euler_6layers_5R(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, m1, m2, m3, m4, m5, m6, cp, R12, R23, R34, R45, R56, m_c_dot, diff_T_c, T_s);
    T_results = system_of_equations_Euler_6layers_5R_singlePipe(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, m1, m2, m3, m4, m5, m6, cp, R12, R23, R34, R45, R56, m_c_dot, diff_T_c, T_s);


    % Compare the model output with the actual data
%     err = mean((T(:,1) - data(:,1)).^2);    %Mean Squared Error
%     err = sqrt(mean((T_results(:,1) - data(:,1)).^2));    %Root Mean Squared Error
%     err = sqrt(mean((T_results(:,6) - data(:,3)).^2));    %Root Mean Squared Error
%     err = sum(sqrt(mean((T_results(:,[1,6]) - data(:,[1,3])).^2)));               %Root Mean Squared Error
    err = sum(sqrt(mean((T_results(:,[1,3,6]) - data).^2)));               %Root Mean Squared Error
end