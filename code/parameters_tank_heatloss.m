function [h_layer_tank1,h_layer_tank2,R_layer_tank1,R_layer_tank2,R_thermal_tank1,R_thermal_tank2,A_heatloss_tank1,A_heatloss_tank2]=parameters_tank_heatloss(num_layer_tank1,num_layer_tank2,m_layer_tank1,m_layer_tank2)


% Parameters
k= 0.6; % the thermal conductivity of water  [W/m·K]
r=0.65/2; % radius of tank  [m]
A= pi*r^2; % cross-sectional area between layers    [m^2]
H=1.64; % height of tank    [m]

h_layer_tank1 = H*(m_layer_tank1/500);
h_layer_tank2 = H*(m_layer_tank2/500);

% surface area
A_heatloss_tank1 = A*h_layer_tank1;
A_heatloss_tank2 = A*h_layer_tank2;


% thermal resistances
%R_i = h_i / (k_i * A)
R_layer_tank1 = h_layer_tank1/(k*A);
R_layer_tank2 = h_layer_tank2/(k*A);


% heat transfer coefficient
R_thermal_tank1 = zeros(num_layer_tank1-1,1);
R_thermal_tank2 = zeros(num_layer_tank2-1,1);


for i=1:num_layer_tank1-1
    R_thermal_tank1(i)=1/(2 * (R_layer_tank1(i) * R_layer_tank1(i+1)) / (R_layer_tank1(i) + R_layer_tank1(i+1)));
end
for i=1:num_layer_tank2-1
    R_thermal_tank2(i)=1/(2 * (R_layer_tank2(i) * R_layer_tank2(i+1)) / (R_layer_tank2(i) + R_layer_tank2(i+1)));
end

end







