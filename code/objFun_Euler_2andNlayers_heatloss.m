function err = objFun_Euler_2andNlayers_heatloss(time_data, T_initial, T_outdoor, m_s_dot_interp,num_layer_tank2, m_layer_tank1,m_layer_tank2_optimal, R_tank12_optimal, m_c_dot_optimal, diff_T_c_optimal, T_s_optimal, params, T_ground_truth)
    
num_layer_tank1=2;

% load("optimal_params.mat")
% R_tank12_optimal = P_optimal(1);
% m_layer_tank2_optimal=zeros(num_layer_tank2,1);
% for i=1:num_layer_tank2
%     m_layer_tank2_optimal(i) = P_optimal(i+1);
% end
% 
% diff_T_c_optimal = P_optimal((1+num_layer_tank2)+1);
% T_s_optimal = P_optimal((1+num_layer_tank2)+2);
% m_c_dot_optimal = P_optimal((1+num_layer_tank2)+3);



R_loss = params;
    
    
% [h_layer_tank1,h_layer_tank2,R_layer_tank1,R_layer_tank2,R_thermal_tank1,R_thermal_tank2]=parameters_tank(num_layer_tank1,num_layer_tank2,m_layer_tank1,m_layer_tank2);



    
    % Solve the system of equations with the current R1 and R2
%     T_results = system_of_equations_Euler_6layers_5R(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, m1, m2, m3, m4, m5, m6, cp, R12, R23, R34, R45, R56, m_c_dot, diff_T_c, T_s);
%     T_results = system_of_equations_Euler_6layers_5R_singlePipe(time_data, T_initial, m_p_dot_data, T_in_data, m_s_dot_data, m1, m2, m3, m4, m5, m6, cp, R12, R23, R34, R45, R56, m_c_dot, diff_T_c, T_s);

    T_results = system_of_equations_Euler_Nlayers_heatloss(time_data, T_initial, T_outdoor, m_s_dot_interp, num_layer_tank1,num_layer_tank2, m_layer_tank1, m_layer_tank2_optimal, R_tank12_optimal,R_loss, m_c_dot_optimal, diff_T_c_optimal, T_s_optimal);

    % Compare the model output with the actual data
%     err = sum(sqrt(mean((T_results(:,[1,num_layer_tank1+1,num_layer_tank1+num_layer_tank2]) - T_ground_truth).^2)));               %Root Mean Squared Error
err = sum(sqrt(mean((T_results(:,1) - T_ground_truth(:,1)).^2)));               %Root Mean Squared Error

end