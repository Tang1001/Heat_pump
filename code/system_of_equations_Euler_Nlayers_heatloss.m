%27-5-2023
function T_results = system_of_equations_Euler_Nlayers_heatloss(time_data, T_initial,T_outdoor, m_s_dot_interp, num_layer_tank1,num_layer_tank2, m_layer_tank1, m_layer_tank2, R_tank12,R_loss, m_c_dot, diff_T_c, T_s)

num_layer=num_layer_tank1+num_layer_tank2;

% Parameters
cp = 4186; % specific heat  [J/kg·K]        cp = 4217 - 0.3 * T


[h_layer_tank1,h_layer_tank2,R_layer_tank1,R_layer_tank2,R_thermal_tank1,R_thermal_tank2,A_heatloss_tank1,A_heatloss_tank2]=parameters_tank_heatloss(num_layer_tank1,num_layer_tank2,m_layer_tank1,m_layer_tank2);


A_heatloss_tank1


T_results=zeros(length(time_data),num_layer);
 

T=T_initial;
T_results(1,:)=T_initial;




dt = time_data(2) - time_data(1); % time step
for k = 1:length(time_data) - 1
    m_s_dot = m_s_dot_interp(k);

    dTndt=zeros(num_layer_tank1+num_layer_tank2,1);


    %tank1
    dTndt(1) = (R_loss*A_heatloss_tank1(1)*(T_outdoor(k)-T(1))*(60*60) - R_thermal_tank1(1) * (T(1) - T(2))*(60*60) - (m_c_dot - m_s_dot) * cp * diff_T_c + m_s_dot * cp * (T(2) - T(1))) / (m_layer_tank1(1) * cp);
    for i=2:num_layer_tank1
        if i==num_layer_tank1
            dTndt(i) = (R_loss*A_heatloss_tank1(i)*(T_outdoor(k)-T(i))*(60*60) + R_thermal_tank1(i-1) * (T(i-1) - T(i))*(60*60) - R_tank12 * (T(i) - T(i+1)) + m_s_dot * cp * (T(i+1) - T(i))) / (m_layer_tank1(i) * cp);
        else
            dTndt(i) = (R_loss*A_heatloss_tank1(i)*(T_outdoor(k)-T(i))*(60*60) + R_thermal_tank1(i-1) * (T(i-1) - T(i))*(60*60) - R_thermal_tank1(i) * (T(i) - T(i+1))  + m_s_dot * cp * (T(i+1) - T(i)) ) / (m_layer_tank1(i) * cp);
        end

    end
    %tank2
    for i=num_layer_tank1+1:num_layer
        if i==num_layer_tank1+1
            dTndt(i) = (R_loss*A_heatloss_tank2(i-num_layer_tank1)*(T_outdoor(k)-T(i))*(60*60) - R_thermal_tank2(i-num_layer_tank1) * (T(i) - T(i-1))*(60*60) + R_tank12 * (T(i-1) - T(i))  + m_s_dot * cp * (T(i+1) - T(i)) ) / (m_layer_tank2(i-num_layer_tank1) * cp);
        elseif i==num_layer
            dTndt(i) = (R_loss*A_heatloss_tank2(i-num_layer_tank1)*(T_outdoor(k)-T(i))*(60*60) + R_thermal_tank2(i-num_layer_tank1-1) * (T(i-1) - T(i))*(60*60)  + m_s_dot * cp * (T_s - T(i)) ) / (m_layer_tank2((i-num_layer_tank1)) * cp);
        else
            dTndt(i) = (R_loss*A_heatloss_tank2(i-num_layer_tank1)*(T_outdoor(k)-T(i))*(60*60) + R_thermal_tank2(i-num_layer_tank1-1) * (T(i-1) - T(i))*(60*60) - R_thermal_tank2(i-num_layer_tank1) * (T(i) - T(i+1))*(60*60)  + m_s_dot * cp * (T(i+1) - T(i))) / (m_layer_tank2((i-num_layer_tank1)) * cp);
        end

    end

%     %tank1
%     dTndt(1) = (R_loss*A_heatloss_tank1(1)*(T_outdoor(k)-T(1)) - R_thermal_tank1(1) * (T(1) - T(2))*(60*60) - (m_c_dot - m_s_dot) * cp * diff_T_c ) / (m_layer_tank1(1) * cp);
%     for i=2:num_layer_tank1
%         if i==num_layer_tank1
%             dTndt(i) = (R_loss*A_heatloss_tank1(i)*(T_outdoor(k)-T(i)) + R_thermal_tank1(i-1) * (T(i-1) - T(i))*(60*60) - R_tank12 * (T(i) - T(i+1)) ) / (m_layer_tank1(i) * cp);
%         else
%             dTndt(i) = (R_loss*A_heatloss_tank1(i)*(T_outdoor(k)-T(i)) + R_thermal_tank1(i-1) * (T(i-1) - T(i))*(60*60) - R_thermal_tank1(i) * (T(i) - T(i+1)) ) / (m_layer_tank1(i) * cp);
%         end
% 
%     end
%     %tank2
%     for i=num_layer_tank1+1:num_layer
%         if i==num_layer_tank1+1
%             dTndt(i) = (R_loss*A_heatloss_tank2(i-num_layer_tank1)*(T_outdoor(k)-T(i)) - R_thermal_tank2(i-num_layer_tank1) * (T(i) - T(i-1))*(60*60) + R_tank12 * (T(i-1) - T(i))) / (m_layer_tank2(i-num_layer_tank1) * cp);
%         elseif i==num_layer
%             dTndt(i) = (R_loss*A_heatloss_tank2(i-num_layer_tank1)*(T_outdoor(k)-T(i)) + R_thermal_tank2(i-num_layer_tank1-1) * (T(i-1) - T(i))*(60*60) ) / (m_layer_tank2((i-num_layer_tank1)) * cp);
%         else
%             dTndt(i) = (R_loss*A_heatloss_tank2(i-num_layer_tank1)*(T_outdoor(k)-T(i)) + R_thermal_tank2(i-num_layer_tank1-1) * (T(i-1) - T(i))*(60*60) - R_thermal_tank2(i-num_layer_tank1) * (T(i) - T(i+1))*(60*60)) / (m_layer_tank2((i-num_layer_tank1)) * cp);
%         end
% 
%     end



    T = T + dt * dTndt; % update the state


    % Store the results
    T_results(k + 1,:) = T;
end



end