function dTdt = system_of_equations_4layers_IDR1R2R3(t, T, m_p_dot_data, T_in_data, m_s_dot_data, time_data_1min, time_data_1h, m1, m2, m3, m4, cp, m_c_dot, diff_T_c, T_s, params)
    T1 = T(1);
    T2 = T(2);
    T3 = T(3);
    T4 = T(4);
    
    R1 = params(1);
    R2 = params(2);
    R3 = params(3);
    
    m_p_dot = interp1(time_data_1min, m_p_dot_data, t);
    T_in = interp1(time_data_1min, T_in_data, t);
    m_s_dot = interp1(time_data_1h, m_s_dot_data, t,'previous');
    
    dT1dt = (m_p_dot * cp * (T_in - T1) - R1 * (T1 - T2)*(60*60) - (m_c_dot - m_s_dot) * cp * diff_T_c + m_s_dot * cp * (T2 - T1)) / (m1 * cp);
    dT2dt = (m_p_dot * cp * (T1 - T2) + R1 * (T1 - T2)*(60*60) - R3 * (T2 - T3)*(60*60) + m_s_dot * cp * (T3 - T2)) / (m2 * cp);
    dT3dt = (m_p_dot * cp * (T2 - T3) - R2 * (T3 - T4)*(60*60) + R3 * (T2 - T3)*(60*60) + m_s_dot * cp * (T4 - T3)) / (m3 * cp);
    dT4dt = (m_p_dot * cp * (T3 - T4) + R2 * (T3 - T4)*(60*60) + m_s_dot * cp * (T_s - T4)) / (m4 * cp);
    
    dTdt = [dT1dt; dT2dt; dT3dt; dT4dt];
end