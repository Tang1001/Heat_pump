function params_optimal=Identify_params_2and4(time_data,sample_time,T_in_interp,m_p_dot_interp,T_buttom_water_tank2_interp,T_upper_water_tank1_interp,T_upper_water_tank2_interp,m_s_dot_interp,params0,data)


% Define the objective function
objective_function = @(params0) objFun_Euler_6layers_2and4(time_data,sample_time,T_in_interp,m_p_dot_interp,T_buttom_water_tank2_interp,T_upper_water_tank1_interp,T_upper_water_tank2_interp,m_s_dot_interp, params0, data)
% Call fminunc
Aeq = [0, 0,0,0, 1,1,1,1];
beq = [500];
lb = [0, 800,1,6, 0,0,0,0];
ub = [0.1, 1200, 5, 18, 500,500,500,500];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
params_optimal = fmincon(objective_function, params0, [], [], Aeq, beq, lb, ub, [], options);






% 
% 
% 
% R23_optimal = R_optimal(1);
% m3_optimal = R_optimal(2);
% m4_optimal = R_optimal(3);
% m5_optimal = R_optimal(4);
% m6_optimal = R_optimal(5);
% diff_T_c_optimal = R_optimal(6);
% T_s_optimal = R_optimal(7);
% m_c_dot_optimal = R_optimal(8);
% 
% 
% 
% % height of layer    [m]
% h3_optimal=H*(m3_optimal/500);
% h4_optimal=H*(m4_optimal/500);
% h5_optimal=H*(m5_optimal/500);
% h6_optimal=H*(m6_optimal/500);
% 
% 
% % thermal resistances    R_i = h_i / (k_i * A)
% R3_optimal=h3_optimal/(k*A);
% R4_optimal=h4_optimal/(k*A);
% R5_optimal=h5_optimal/(k*A);
% R6_optimal=h6_optimal/(k*A);
% 
% 
% % Thermal efficiency between layers
% R34_optimal = 1/(2 * (R3_optimal * R4_optimal) / (R3_optimal + R4_optimal));
% R45_optimal = 1/(2 * (R4_optimal * R5_optimal) / (R4_optimal + R5_optimal));
% R56_optimal = 1/(2 * (R5_optimal * R6_optimal) / (R5_optimal + R6_optimal));
% 

%%
% T_results = system_of_equations_Euler_6layers_5R(time_data, T, m_p_dot_interp, T_in_interp, m_s_dot_interp, m1, m2, m3_optimal, m4_optimal, m5_optimal, m6_optimal, cp, R12, R23_optimal, R34_optimal, R45_optimal, R56_optimal, m_c_dot_optimal, diff_T_c_optimal, T_s_optimal);


end