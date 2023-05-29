%27-5-2023
function T_results = system_of_equations_Euler_Nlayers(time_data, T_initial, m_p_dot_interp, T_in_interp, m_s_dot_interp, num_layer_tank1,num_layer_tank2, m_layer_tank1, m_layer_tank2, R_tank12, m_c_dot, diff_T_c, T_s)

num_layer=num_layer_tank1+num_layer_tank2;

% Parameters
cp = 4186; % specific heat  [J/kg·K]        cp = 4217 - 0.3 * T


[h_layer_tank1,h_layer_tank2,R_layer_tank1,R_layer_tank2,R_thermal_tank1,R_thermal_tank2]=parameters_tank(num_layer_tank1,num_layer_tank2,m_layer_tank1,m_layer_tank2);



T_results=zeros(length(time_data),num_layer);
 

T=T_initial;
T_results(1,:)=T_initial;


heatpump_flag_data = (m_p_dot_interp>0);


dt = time_data(2) - time_data(1); % time step
for k = 1:length(time_data) - 1
    m_p_dot = m_p_dot_interp(k);
    T_in = T_in_interp(k);
    m_s_dot = m_s_dot_interp(k);
    heatpump_flag = heatpump_flag_data(k);

    dTndt=zeros(num_layer_tank1+num_layer_tank2,1);

    %tank1
    dTndt(1) = (m_p_dot * cp * (T_in - T(1)) - R_thermal_tank1(1) * (T(1) - T(2))*(60*60) - (m_c_dot - m_s_dot) * cp * diff_T_c + m_s_dot * (1 - heatpump_flag) * cp * (T(2) - T(1))) / (m_layer_tank1(1) * cp);
    for i=2:num_layer_tank1
        if i==num_layer_tank1
            dTndt(i) = ((m_p_dot-m_s_dot) * (heatpump_flag) * cp * (T(i-1) - T(i)) + R_thermal_tank1(i-1) * (T(i-1) - T(i))*(60*60) - R_tank12 * (T(i) - T(i+1)) + m_s_dot * (1 - heatpump_flag) * cp * (T(i+1) - T(i))) / (m_layer_tank1(i) * cp);
        else
            dTndt(i) = ((m_p_dot-m_s_dot) * (heatpump_flag) * cp * (T(i-1) - T(i)) + R_thermal_tank1(i-1) * (T(i-1) - T(i))*(60*60) - R_thermal_tank1(i) * (T(i) - T(i+1)) + m_s_dot * (1 - heatpump_flag) * cp * (T(i+1) - T(i))) / (m_layer_tank1(i) * cp);
        end

    end
    %tank2
    for i=num_layer_tank1+1:num_layer
        if i==num_layer_tank1+1
            dTndt(i) = ((m_p_dot-m_s_dot) * (heatpump_flag) * cp * (T(i-1) - T(i)) - R_thermal_tank2(i-num_layer_tank1) * (T(i) - T(i-1))*(60*60) + R_tank12 * (T(i-1) - T(i)) + m_s_dot * (1 - heatpump_flag) * cp * (T(i+1) - T(i))) / (m_layer_tank2(i-num_layer_tank1) * cp);
        elseif i==num_layer
            dTndt(i) = ((m_p_dot-m_s_dot) * (heatpump_flag) * cp * (T(i-1) - T(i)) + R_thermal_tank2(i-num_layer_tank1-1) * (T(i-1) - T(i))*(60*60) + m_s_dot * (1 - heatpump_flag) * cp * (T_s - T(i))) / (m_layer_tank2((i-num_layer_tank1)) * cp);
        else
            dTndt(i) = ((m_p_dot-m_s_dot) * (heatpump_flag) * cp * (T(i-1) - T(i)) + R_thermal_tank2(i-num_layer_tank1-1) * (T(i-1) - T(i))*(60*60) - R_thermal_tank2(i-num_layer_tank1) * (T(i) - T(i+1))*(60*60) + m_s_dot * (1 - heatpump_flag) * cp * (T(i+1) - T(i))) / (m_layer_tank2((i-num_layer_tank1)) * cp);
        end

    end


%     dT1dt = (m_p_dot * cp * (T_in - T(1)) - R12 * (T(1) - T(2))*(60*60) - (m_c_dot - m_s_dot) * cp * diff_T_c + m_s_dot * (1 - heatpump_flag) * cp * (T(2) - T(1))) / (m1 * cp);
%     dT2dt = ((m_p_dot-m_s_dot) * (heatpump_flag) * cp * (T(1) - T(2)) + R12 * (T(1) - T(2))*(60*60) - R23 * (T(2) - T(3)) + m_s_dot * (1 - heatpump_flag) * cp * (T(3) - T(2))) / (m2 * cp);
%     dT3dt = ((m_p_dot-m_s_dot) * (heatpump_flag) * cp * (T(2) - T(3)) - R34 * (T(3) - T(4))*(60*60) + R23 * (T(2) - T(3)) + m_s_dot * (1 - heatpump_flag) * cp * (T(4) - T(3))) / (m3 * cp);
%     dT4dt = ((m_p_dot-m_s_dot) * (heatpump_flag) * cp * (T(3) - T(4)) + R34 * (T(3) - T(4))*(60*60) - R45 * (T(4) - T(5))*(60*60) + m_s_dot * (1 - heatpump_flag) * cp * (T(5) - T(4))) / (m4 * cp);
%     dT5dt = ((m_p_dot-m_s_dot) * (heatpump_flag) * cp * (T(4) - T(5)) + R45 * (T(4) - T(5))*(60*60) - R56 * (T(5) - T(6))*(60*60) + m_s_dot * (1 - heatpump_flag) * cp * (T(6) - T(5))) / (m5 * cp);
%     dT6dt = ((m_p_dot-m_s_dot) * (heatpump_flag) * cp * (T(5) - T(6)) + R56 * (T(5) - T(6))*(60*60) + m_s_dot * (1 - heatpump_flag) * cp * (T_s - T(6))) / (m6 * cp);
% 
%     T = T + dt * [dT1dt; dT2dt; dT3dt; dT4dt; dT5dt; dT6dt]; % update the state


    T = T + dt * dTndt; % update the state


    % Store the results
    T_results(k + 1,:) = T;
end



end