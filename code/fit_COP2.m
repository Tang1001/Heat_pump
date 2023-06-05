clear all
close all

% fit the function COP=a1+a2*T1+a3*T1^2+a4*T2+a5*T2^2+a6*T1*T2
%%
% 60 - 65 °C
cop1 = xlsread('../data/2. COP different outdoor temperatures.xlsx', 'T buiten 18,5', 'G3:G43');
T_in1 = xlsread('../data/2. COP different outdoor temperatures.xlsx', 'T buiten 18,5', 'I3:I43');
T_amb1 = 18.5;

% 65 - 70 °C
cop2 = xlsread('../data/2. COP different outdoor temperatures.xlsx', 'T buiten 12,5', 'G3:G138');
T_in2 = xlsread('../data/2. COP different outdoor temperatures.xlsx', 'T buiten 12,5', 'I3:I138');
T_amb2 = 12.5;

%%
T_IN = [T_in1; T_in2]; % Ambient air temperature data
T_AMB = [ones(length(T_in1),1)*T_amb1;ones(length(T_in2),1)*T_amb2];
COP = [cop1; cop2]; % Coefficient of Performance data


data = [T_IN,T_AMB];
model = [0 0; 1 0; 2 0; 0 1; 1 1];

fit_result = polyfitn(data, COP, model);
a1 = fit_result.Coefficients(1);
a2 = fit_result.Coefficients(2);
a3 = fit_result.Coefficients(3);
a4 = fit_result.Coefficients(4);
a5 = fit_result.Coefficients(5);
cop_P2=[a1;a2;a3;a4;a5]

save("cop_P2.mat",'cop_P2')

fprintf('COP = %.2f + %.2f * T_in + %.2f * T_in^2 + %.2f * T_amb + %.2f * T_in * T_amb\n', a1, a2, a3, a4, a5);


[T_in_grid, T_amb_grid] = meshgrid(20:1:65, 5:1:25);
COP_grid = a1 + a2 * T_in_grid + a3 * T_in_grid.^2 + a4 * T_amb_grid + a5 * T_in_grid .* T_amb_grid;
figure;
surf(T_in_grid, T_amb_grid, COP_grid);
xlabel('Temperature of inlet water (°C)');
ylabel('Temperature of ambient water (°C)');
zlabel('COP');
title('COP Function');
hold on;
scatter3(T_IN, T_AMB, COP, 'ro', 'filled', 'MarkerEdgeColor', 'k');
legend('COP Function', 'Data points');
hold off;
