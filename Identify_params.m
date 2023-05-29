function params_optimal=Identify_params(time_data,sample_time,num_layer_tank1,num_layer_tank2,T_in_interp,m_p_dot_interp,T_buttom_water_tank2_interp,T_upper_water_tank1_interp,T_upper_water_tank2_interp,m_s_dot_interp,params0,data)

% Define the objective function
objective_function = @(params0) objFun_Euler_Nlayers(time_data,sample_time,num_layer_tank1,num_layer_tank2,T_in_interp,m_p_dot_interp,T_buttom_water_tank2_interp,T_upper_water_tank1_interp,T_upper_water_tank2_interp,m_s_dot_interp, params0, data)
% Call fminunc
% Aeq = [0, 0,0,0];
% beq = [500];
lb = [0, 800,1,6];
ub = [0.1, 1200, 5, 18];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
params_optimal = fmincon(objective_function, params0, [], [], [], [], lb, ub, [], options);
% params_optimal = fmincon(objective_function, params0, [], [], Aeq, beq, lb, ub, [], options);

% R12_optimal = R_optimal(1);
% diff_T_c_optimal = R_optimal(2);
% T_s_optimal = R_optimal(3);
% m_c_dot_optimal = R_optimal(4);


end