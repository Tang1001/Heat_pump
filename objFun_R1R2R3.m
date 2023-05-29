function err = objFun_R1R2R3(m_p_dot_data, T_in_data, m_s_dot_data, time_data_1min, time_data_1h, m1, m2, m3, m4, cp, m_c_dot, diff_T_c, T_s, params, time_span, initial_conditions, data)

    % Solve the ODE with the current parameters
%     [T,Y] = ode45(@(t,y) myODE(t, y, params), t, y0);
    [t, T] = ode45(@(t, T) system_of_equations_4layers_IDR1R2R3(t, T, m_p_dot_data, T_in_data, m_s_dot_data, time_data_1min, time_data_1h, m1, m2, m3, m4, cp, m_c_dot, diff_T_c, T_s, params), time_span, initial_conditions);


    % Compare the model output with the actual data
%     err = sum(sum((T(:,[1,3,4]) - data).^2));
%     err = mean((T(:,1) - data(:,1)).^2);    %Mean Squared Error
    err = sqrt(mean((T(:,4) - data(:,3)).^2));    %Root Mean Squared Error
end