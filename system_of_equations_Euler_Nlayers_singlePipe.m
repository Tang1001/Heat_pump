function T_results = system_of_equations_Euler_Nlayers_singlePipe(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, num_layer_tank1, num_layer_tank2, m_layer_tank1,m_layer_tank2, R_thermal_tank1, R_thermal_tank2, R12, cp, m_c_dot, diff_T_c, T_s)


T_results=zeros(length(time_data),num_layer_tank1+num_layer_tank2);
T=T_initial;
T_results(1,:)=T_initial;


heatpump_flag_data = (m_p_dot_data>0);

num_layer_tot = num_layer_tank1+num_layer_tank2;

dt = time_data(2) - time_data(1); % time step
for k = 1:length(time_data) - 1
    m_p_dot = m_p_dot_data(k);
    T_in = T_in_data(k);
    m_s_dot = m_s_dot_data(k);
    heatpump_flag = heatpump_flag_data(k);

    dTndt=zeros(num_layer_tank1+num_layer_tank2,1);

    %tank1
    dTndt(1) = (m_p_dot * cp * (T_in - T(1)) - R12 * (T(1) - T(2))*(60*60) - (m_c_dot - m_s_dot) * cp * diff_T_c + m_s_dot * (1 - heatpump_flag) * cp * (T(2) - T(1))) / (m_layer_tank1 * cp);
    for i=2:num_layer_tank1
        if i==num_layer_tank1
            dTndt(i) = ((m_p_dot-m_s_dot) * (heatpump_flag) * cp * (T(i-1) - T(i)) + R_thermal_tank1 * (T(i-1) - T(i))*(60*60) - R12 * (T(i) - T(i+1)) + m_s_dot * (1 - heatpump_flag) * cp * (T(i+1) - T(i))) / (m_layer_tank1 * cp);
        else
            dTndt(i) = ((m_p_dot-m_s_dot) * (heatpump_flag) * cp * (T(i-1) - T(i)) + R_thermal_tank1 * (T(i-1) - T(i))*(60*60) - R_thermal_tank1 * (T(i) - T(i+1)) + m_s_dot * (1 - heatpump_flag) * cp * (T(i+1) - T(i))) / (m_layer_tank1 * cp);
        end

    end
    %tank2
    for i=num_layer_tank1+1:num_layer_tot
        if i==num_layer_tank1+1
            dTndt(i) = ((m_p_dot-m_s_dot) * (heatpump_flag) * cp * (T(i-1) - T(i)) - R_thermal_tank2 * (T(i) - T(i-1))*(60*60) + R12 * (T(i-1) - T(i)) + m_s_dot * (1 - heatpump_flag) * cp * (T(i+1) - T(i))) / (m_layer_tank2 * cp);
        elseif i==num_layer_tot
            dTndt(i) = ((m_p_dot-m_s_dot) * (heatpump_flag) * cp * (T(i-1) - T(i)) + R_thermal_tank2 * (T(i-1) - T(i))*(60*60) + m_s_dot * (1 - heatpump_flag) * cp * (T_s - T(i))) / (m_layer_tank2 * cp);
        else
            dTndt(i) = ((m_p_dot-m_s_dot) * (heatpump_flag) * cp * (T(i-1) - T(i)) + R_thermal_tank2 * (T(i-1) - T(i))*(60*60) - R_thermal_tank2 * (T(i) - T(i+1))*(60*60) + m_s_dot * (1 - heatpump_flag) * cp * (T(i+1) - T(i))) / (m_layer_tank2 * cp);
        end

    end


    T = T + dt * dTndt; % update the state


    % Store the results
    T_results(k + 1,:) = T;
end



end