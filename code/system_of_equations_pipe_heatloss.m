function T_results = system_of_equations_pipe_heatloss(time_data, T_initial, T_ambient, m_pipe, R_loss_pipe)


% Properties
cp = 4186; % Specific heat capacity of water in J/kg/K


dt = time_data(2) - time_data(1); % time step [h]



T_results=zeros(length(time_data),1);


% Set initial conditions
T = T_initial;
T_results(1)=T_initial;

% Euler's method
for i = 1:(length(time_data) - 1)
    dT_dt = -R_loss_pipe*(T - T_ambient(i))/(m_pipe*cp)*(60*60);    %[W/K] * [K] = [W] = [J/s] = (60*60) * [J/h] 
    T = T + dt*dT_dt;
    T_results(i+1)=T;
end


end