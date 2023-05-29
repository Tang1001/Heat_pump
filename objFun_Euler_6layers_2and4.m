function err = objFun_Euler_6layers_2and4(time_data,sample_time,T_in_interp,m_p_dot_interp,T_buttom_water_tank2_interp,T_upper_water_tank1_interp,T_upper_water_tank2_interp,m_s_dot_interp, params, data)

%params: 1R  diff_T_c T_s m_c_dot   4m
% params

T_results=system_simulation_Euler_2and4(time_data,sample_time,T_in_interp,m_p_dot_interp,T_buttom_water_tank2_interp,T_upper_water_tank1_interp,T_upper_water_tank2_interp,m_s_dot_interp,params);



% Compare the model output with the actual data
err = sum(sqrt(mean((T_results(:,[1,3,6]) - data).^2)));               %Root Mean Squared Error
end