function err = objFun_Euler_NandNlayers(time_data, T_initial, m_p_dot_interp, T_in_interp, m_s_dot_interp,num_layer_tank1,num_layer_tank2, params, T_ground_truth)
    



R_tank12 = params(1);
m_layer_tank1=zeros(num_layer_tank1,1);
m_layer_tank2=zeros(num_layer_tank2,1);
for i=1:num_layer_tank1
    m_layer_tank1(i) = params(i+1);
end
for i=1:num_layer_tank2
    m_layer_tank2(i) = params(i+1+num_layer_tank1);
end

diff_T_c = params((1+num_layer_tank1+num_layer_tank2)+1);
T_s = params((1+num_layer_tank1+num_layer_tank2)+2);
m_c_dot = params((1+num_layer_tank1+num_layer_tank2)+3);


% 
% R_tank12 = params(1);
% m_layer_tank2=zeros(num_layer_tank2,1);
%     for i=1:num_layer_tank2
%         m_layer_tank2(i) = params(i+1);
%     end
%     
%     diff_T_c = params((1+num_layer_tank2)+1);
%     T_s = params((1+num_layer_tank2)+2);
%     m_c_dot = params((1+num_layer_tank2)+3);
%     
    
% [h_layer_tank1,h_layer_tank2,R_layer_tank1,R_layer_tank2,R_thermal_tank1,R_thermal_tank2]=parameters_tank(num_layer_tank1,num_layer_tank2,m_layer_tank1,m_layer_tank2);



    
    % Solve the system of equations with the current R1 and R2
%     T_results = system_of_equations_Euler_6layers_5R(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, m1, m2, m3, m4, m5, m6, cp, R12, R23, R34, R45, R56, m_c_dot, diff_T_c, T_s);
%     T_results = system_of_equations_Euler_6layers_5R_singlePipe(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, m1, m2, m3, m4, m5, m6, cp, R12, R23, R34, R45, R56, m_c_dot, diff_T_c, T_s);

    T_results = system_of_equations_Euler_Nlayers(time_data, T_initial, m_p_dot_interp, T_in_interp, m_s_dot_interp, num_layer_tank1,num_layer_tank2, m_layer_tank1, m_layer_tank2, R_tank12, m_c_dot, diff_T_c, T_s);


    % Compare the model output with the actual data
    err = sum(sqrt(mean((T_results(:,[1,num_layer_tank1+1,num_layer_tank1+num_layer_tank2]) - T_ground_truth).^2)));               %Root Mean Squared Error
end