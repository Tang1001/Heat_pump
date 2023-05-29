%% Initial Temperature Setting
function T_initial=initial_T_setting(num_layer_tank1,num_layer_tank2,m_layer_tank1,m_layer_tank2,T0_upper_layer_tank1,T0_upper_layer_tank2,T0_bottom_layer_tank2)

[h_layer_tank1,h_layer_tank2,R_layer_tank1,R_layer_tank2,R_thermal_tank1,R_thermal_tank2]=parameters_tank(num_layer_tank1,num_layer_tank2,m_layer_tank1,m_layer_tank2);
R_tank1_tot=sum(R_layer_tank1);
R_tank2_tot=sum(R_layer_tank2);



T_initial_layer_tank1 = zeros(num_layer_tank1,1);
T_initial_layer_tank2 = zeros(num_layer_tank2,1);
% measurement values
T_initial_layer_tank1(1)=T0_upper_layer_tank1;
T_initial_layer_tank2(1)=T0_upper_layer_tank2;
T_initial_layer_tank2(end)=T0_bottom_layer_tank2;

%Tank1
for i=1:num_layer_tank1-1
    d_T_tank1=(T_initial_layer_tank1(1) - T_initial_layer_tank2(1)) * (R_layer_tank1(i) / R_tank1_tot);
    T_initial_layer_tank1(i+1) = T_initial_layer_tank1(i) - d_T_tank1; 
end

%Tank2
for i=1:num_layer_tank2-1-1
%     (R_layer_tank2(i) / R_tank2_tot)
    d_T_tank2=(T_initial_layer_tank2(1) - T_initial_layer_tank2(end)) * (R_layer_tank2(i) / R_tank2_tot);
    T_initial_layer_tank2(i+1) = T_initial_layer_tank2(i) - d_T_tank2;
end

T_initial = [T_initial_layer_tank1;T_initial_layer_tank2];



% %7:00 %1min:60*7+1=421  8min:60/8*7+1=53.5
% T1_initial=T_upper_tank1_interp(1); % layer 1 initial temperature
% T3_initial=T_upper_tank2_interp(1); % layer 3 initial temperature
% T6_initial=T_bottom_tank2_interp(1); % layer 5 initial temperature
% 
% d_T1 = (T1_initial - T3_initial) * (R1 / R1_tot);
% T2_initial = T1_initial - d_T1;     % layer 2 initial temperature
% d_T3 = (T3_initial - T6_initial) * (R3 / R2_tot);
% T4_initial= T3_initial - d_T3; % layer 4 initial temperature
% d_T4 = (T3_initial - T6_initial) * (R4 / R2_tot);
% T5_initial= T4_initial - d_T4; % layer 4 initial temperature
% 
% T_initial = [T1_initial; T2_initial; T3_initial; T4_initial; T5_initial; T6_initial];

end