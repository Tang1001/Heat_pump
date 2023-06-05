clear all
close all

% fit the function COP=a1+a2*T1+a3*T1^2+a4*T2+a5*T2^2+a6*T1*T2
%%
% 60 - 65 °C
cop1 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '60 - 65 °C (WEINIG GEGEVENS)', 'G3:G15');
T_out1 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '60 - 65 °C (WEINIG GEGEVENS)', 'H3:H15');
T_in1 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '60 - 65 °C (WEINIG GEGEVENS)', 'I3:I15');

% 65 - 70 °C
cop2 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '65 - 70 °C', 'G3:G130');
T_out2 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '65 - 70 °C', 'H3:H130');
T_in2 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '65 - 70 °C', 'I3:I130');

% 70 - 75 °C
cop3 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '70 - 75 °C', 'G4:G494');
T_out3 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '70 - 75 °C', 'H4:H494');
T_in3 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '70 - 75 °C', 'I4:I494');

% 75 - 80 °C
cop4 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '75 - 80 °C', 'G3:G1078');
T_out4 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '75 - 80 °C', 'H3:H1078');
T_in4 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '75 - 80 °C', 'I3:I1078');

% 80 - 85 °C
cop5 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '80 - 85 °C', 'G3:G232');
T_out5 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '80 - 85 °C', 'H3:H232');
T_in5 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '80 - 85 °C', 'I3:I232');

% 85 - 90 °C
cop6 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '85 - 90 °C', 'G3:G159');
T_out6 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '85 - 90 °C', 'H3:H159');
T_in6 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '85 - 90 °C', 'I3:I159');

% 90 - 95 °C
cop7 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '90 - 95 °C', 'G3:G57');
T_out7 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '90 - 95 °C', 'H3:H57');
T_in7 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '90 - 95 °C', 'I3:I57');

% > 95 °C (WEINIG GEGEVENS)
cop8 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '> 95 °C (WEINIG GEGEVENS)', 'G3:G15');
T_out8 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '> 95 °C (WEINIG GEGEVENS)', 'H3:H15');
T_in8 = xlsread('../data/1. COP different temperatures heat pump in and out.xlsx', '> 95 °C (WEINIG GEGEVENS)', 'I3:I15');

%%
T1 = [T_out1; T_out2; T_out3; T_out4; T_out5; T_out6; T_out7; T_out8]; % Temperature of inlet water data
T2 = [T_in1; T_in2; T_in3; T_in4; T_in5; T_in6; T_in7; T_in8]; % Ambient air temperature data
COP = [cop1; cop2; cop3; cop4; cop5; cop6; cop7; cop8]; % Coefficient of Performance data

% 
% data = [T1, T2];
% model = [0 0; 1 0; 2 0; 0 1; 0 2; 1 1];
% 
% fit_result = polyfitn(data, COP, model);
% a1 = fit_result.Coefficients(1);
% a2 = fit_result.Coefficients(2);
% a3 = fit_result.Coefficients(3);
% a4 = fit_result.Coefficients(4);
% a5 = fit_result.Coefficients(5);
% a6 = fit_result.Coefficients(6);
% 
% fprintf('COP = %.2f + %.2f * T1 + %.2f * T1^2 + %.2f * T2 + %.2f * T2^2 + %.2f * T1 * T2\n', a1, a2, a3, a4, a5, a6);
% 
% 
% [T1_grid, T2_grid] = meshgrid(65:1:95, 20:1:65);
% COP_grid = a1 + a2 * T1_grid + a3 * T1_grid.^2 + a4 * T2_grid + a5 * T2_grid.^2 + a6 * T1_grid .* T2_grid;
% figure;
% surf(T1_grid, T2_grid, COP_grid);
% xlabel('Temperature of outlet water (°C)');
% ylabel('Temperature of inlet water (°C)');
% zlabel('COP');
% title('COP Function');
% hold on;
% scatter3(T1, T2, COP, 'ro', 'filled', 'MarkerEdgeColor', 'k');
% legend('COP Function', 'Data points');
% hold off;


%%

data2 = [T2];
model2 = [0; 1; 2];

fit_result = polyfitn(data2, COP, model2)
a1 = fit_result.Coefficients(1);
a2 = fit_result.Coefficients(2);
a3 = fit_result.Coefficients(3);

cop_P1=[a1;a2;a3]

save("cop_P1.mat",'cop_P1')

fprintf('COP = %.2f + %.2f * T2 + %.2f * T2^2\n', a1, a2, a3);



[T1_grid, T2_grid] = meshgrid(65:1:95, 20:1:65);
COP_grid = a1 + a2 * T2_grid + a3 * T2_grid.^2;
figure;
plot(T2_grid, COP_grid);
xlabel('Temperature of inlet water (°C)');
ylabel('COP');
title('COP Function');
hold on;
scatter(T2, COP, 'ro', 'filled', 'MarkerEdgeColor', 'k');
% legend('COP Function', 'Data points');
% hold off;
