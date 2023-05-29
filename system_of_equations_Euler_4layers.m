function T_results = system_of_equations_Euler_4layers(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, m1, m2, m3, m4, cp, R1, R2, m_c_dot, diff_T_c, T_s)


T_results=zeros(length(time_data),4);
T=T_initial;
T_results(1,:)=T_initial;

dt = time_data(2) - time_data(1); % time step
for k = 1:length(time_data) - 1
    m_p_dot = m_p_dot_data(k);
    T_in = T_in_data(k);
    m_s_dot = m_s_dot_data(k);

    dT1dt = (m_p_dot * cp * (T_in - T(1)) - R1 * (T(1) - T(2))*(60*60) - (m_c_dot - m_s_dot) * cp * diff_T_c + m_s_dot * cp * (T(2) - T(1))) / (m1 * cp);
    dT2dt = (m_p_dot * cp * (T(1) - T(2)) + R1 * (T(1) - T(2))*(60*60) + m_s_dot * cp * (T(3) - T(2))) / (m2 * cp);
    dT3dt = (m_p_dot * cp * (T(2) - T(3)) - R2 * (T(3) - T(4))*(60*60) + m_s_dot * cp * (T(4) - T(3))) / (m3 * cp);
    dT4dt = (m_p_dot * cp * (T(3) - T(4)) + R2 * (T(3) - T(4))*(60*60) + m_s_dot * cp * (T_s - T(4))) / (m4 * cp);

    T = T + dt * [dT1dt; dT2dt; dT3dt; dT4dt]; % update the state

    % Store the results
    T_results(k + 1,:) = T;
end



end