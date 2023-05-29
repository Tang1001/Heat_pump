function err = objFun_heatloss_outdoor(m_p_dot_data, T_in_data, T_outdoor_data, m_s_dot_data, time_data_1min, time_data_1h, m1, m2, m3, m4, cp, m_c_dot, diff_T_c, T_s, params, time_span, initial_conditions, data)

    % Solve the ODE with the current parameters
%     [T,Y] = ode45(@(t,y) myODE(t, y, params), t, y0);
    [t, T] = ode45(@(t, T) system_of_equations_4layers_IDR1R2Const_outdoor(t, T, m_p_dot_data, T_in_data, T_outdoor_data, m_s_dot_data,  time_data_1min, time_data_1h, m1, m2, m3, m4, cp, m_c_dot, diff_T_c, T_s, params), time_span, initial_conditions);


    % Compare the model output with the actual data
%     err = sum(sum((T(:,[1,3,4]) - data).^2));
    err = sum((T(:,1) - data(:,1)).^2);
end