function T_results = system_of_equations_Euler_7layers_6R(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, m1, m2, m3, m4, m5, m6, m7, cp, R12, R23, R34, R45, R56, R67, m_c_dot, diff_T_c, T_s)


T_results=zeros(length(time_data),7);
T=T_initial;
T_results(1,:)=T_initial;

dt = time_data(2) - time_data(1); % time step
for k = 1:length(time_data) - 1
    m_p_dot = m_p_dot_data(k);
    T_in = T_in_data(k);
    m_s_dot = m_s_dot_data(k);

    dT1dt = (m_p_dot * cp * (T_in - T(1)) - R12 * (T(1) - T(2))*(60*60) - (m_c_dot - m_s_dot) * cp * diff_T_c + m_s_dot * cp * (T(2) - T(1))) / (m1 * cp);
    dT2dt = (m_p_dot * cp * (T(1) - T(2)) + R12 * (T(1) - T(2))*(60*60) - R23 * (T(2) - T(3)) + m_s_dot * cp * (T(3) - T(2))) / (m2 * cp);
    dT3dt = (m_p_dot * cp * (T(2) - T(3)) - R34 * (T(3) - T(4))*(60*60) + R23 * (T(2) - T(3)) + m_s_dot * cp * (T(4) - T(3))) / (m3 * cp);
    dT4dt = (m_p_dot * cp * (T(3) - T(4)) + R34 * (T(3) - T(4))*(60*60) - R45 * (T(4) - T(5))*(60*60) + m_s_dot * cp * (T(5) - T(4))) / (m4 * cp);
    dT5dt = (m_p_dot * cp * (T(4) - T(5)) + R45 * (T(4) - T(5))*(60*60) - R56 * (T(5) - T(6))*(60*60) + m_s_dot * cp * (T(6) - T(5))) / (m5 * cp);
    dT6dt = (m_p_dot * cp * (T(5) - T(6)) + R56 * (T(5) - T(6))*(60*60) - R67 * (T(6) - T(7))*(60*60) + m_s_dot * cp * (T(7) - T(6))) / (m6 * cp);
    dT7dt = (m_p_dot * cp * (T(6) - T(7)) + R67 * (T(6) - T(7))*(60*60) + m_s_dot * cp * (T_s - T(7))) / (m7 * cp);


    T = T + dt * [dT1dt; dT2dt; dT3dt; dT4dt; dT5dt; dT6dt; dT7dt]; % update the state

    % Store the results
    T_results(k + 1,:) = T;
end



end