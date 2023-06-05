function err = objFun_Euler_pipe_heatloss(time_data, T_initial, T_ambient,  params, T_ground_truth)
    



R_loss_pipe = params(1);
m_pipe = params(2);

T_results = system_of_equations_pipe_heatloss(time_data, T_initial, T_ambient, m_pipe, R_loss_pipe);

% Compare the model output with the actual data
err = sqrt(mean((T_results - T_ground_truth).^2));               %Root Mean Squared Error

end