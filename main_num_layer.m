close all
clear all

%% Read data 2023-03-27
%1min
Temp_out_water_he = readcell('Data/Mon 1 May 2023 outlet water T in heat exchanger.csv');
Temp_out_water_he=cell2table(Temp_out_water_he(2:end,:));
[date_out_water_he,T_out_water_he]=read_csv(Temp_out_water_he);

FlowRate_water_hp = readcell('Data/Mon 1 May 2023 water flow rate in heat pump.csv');
FlowRate_water_hp=cell2table(FlowRate_water_hp(2:end,:));
[date_FR_water_hp,FR_water_hp]=read_csv(FlowRate_water_hp);

Temp_buttom_water_tank2 = readcell('Data/Mon 1 May 2023 water T in bottom 2 tank.csv');
Temp_buttom_water_tank2=cell2table(Temp_buttom_water_tank2(2:end,:));
[date_buttom_water_tank2,T_buttom_water_tank2]=read_csv(Temp_buttom_water_tank2);

Temp_upper_water_tank1 = readcell('Data/Mon 1 May 2023 water T in upper 1 tank.csv');
Temp_upper_water_tank1=cell2table(Temp_upper_water_tank1(2:end,:));
[date_upper_water_tank1,T_upper_water_tank1]=read_csv(Temp_upper_water_tank1);


% Temp_outdoor = readcell('Thu 13 April 2023 outdoor T.csv');
% Temp_outdoor=cell2table(Temp_outdoor(2:end,:));
% [date_outdoor,T_outdoor]=read_csv(Temp_outdoor);

%8min
Temp_upper_water_tank2 = readcell('Data/Mon 1 May 2023 water T in upper 2 tank.csv');
Temp_upper_water_tank2=cell2table(Temp_upper_water_tank2(2:end,:));
[date_upper_water_tank2,T_upper_water_tank2]=read_csv(Temp_upper_water_tank2);

%1h
Consp_water = readcell('Data/Mon 1 May 2023 water V consumption.csv');
Consp_water=cell2table(Consp_water(2:end,:));
[date_V_water,V_water]=read_csv(Consp_water);  %[m^3]



%%
%2023-03-27 7:00 ~ 19:00
%7:00 ~ 19:00 -- 0~12
% Time span for the simulation
t_start = 0;
t_end = 12;     % end time [h]
sample_time = 2;    % [min]
time_span = [t_start:sample_time/60:t_end];
time_data = time_span;


%1min
m_p_dot_data=FR_water_hp(60*7+1:60*19+1)*1000;  %[L/h] 60*7+1~60*19
T_in_data=T_out_water_he(60*7+1:60*19+1);

%1h
m_s_dot_data=V_water(1*8:1*20)*1000;    %[L/h] 07:00~19:00
% m_s_dot_data=V_water_130423(1*9:1*20)*1000;    %[L/h] 08:00~19:00

time_data_1min=0:1/60:12;
time_data_1h=0:12;
time_data_8min=0:8/60:12;


T_upper_water_tank1_interp = interp1(time_data_1min,T_upper_water_tank1(60*7+1:60*19+1),time_data);
T_upper_water_tank2_interp = interp1(time_data_8min,T_upper_water_tank2(round(60/8*7)+1:round(60/8*19)+1),time_data);
T_buttom_water_tank2_interp = interp1(time_data_1min,T_buttom_water_tank2(60*7+1:60*19+1),time_data);
data = [T_upper_water_tank1_interp',T_upper_water_tank2_interp',T_buttom_water_tank2_interp']; % Your historical data here


%%


num_layer_tank1 = [2:10];
num_layer_tank2 = [2:20];

err_upper_tank1 = zeros(length(num_layer_tank1),length(num_layer_tank2));
err_upper_tank2 = zeros(length(num_layer_tank1),length(num_layer_tank2));
err_bottom_tank2 = zeros(length(num_layer_tank1),length(num_layer_tank2));

for i=1:length(num_layer_tank1)
    for j=1:length(num_layer_tank2)
        T_results=system_simulation_Euler(time_data,sample_time,num_layer_tank1(i),num_layer_tank2(j),T_out_water_he,FR_water_hp,T_buttom_water_tank2,T_upper_water_tank1,T_upper_water_tank2,V_water);

        err_upper_tank1(i,j) = sqrt(mean((T_results(:,1) - data(:,1)).^2));               %Root Mean Squared Error
        err_upper_tank2(i,j) = sqrt(mean((T_results(:,num_layer_tank1(i)+1) - data(:,2)).^2));               %Root Mean Squared Error
        err_bottom_tank2(i,j) = sqrt(mean((T_results(:,num_layer_tank1(i)+num_layer_tank2(j)) - data(:,3)).^2));               %Root Mean Squared Error

    end

end


%%
figure;
subplot(3,1,1)
h=heatmap(num_layer_tank2,num_layer_tank1,err_upper_tank1)
h.Colormap = winter;
xlabel('The number of layers in tank2')
ylabel('The number of layers in tank1')
title('Error of the upper temperature in tank1')
subplot(3,1,2)
h=heatmap(num_layer_tank2,num_layer_tank1,err_upper_tank2)
h.Colormap = winter;
xlabel('The number of layers in tank2')
ylabel('The number of layers in tank1')
title('Error of the upper temperature in tank2')
subplot(3,1,3)
h=heatmap(num_layer_tank2,num_layer_tank1,err_bottom_tank2)
h.Colormap = winter;
xlabel('The number of layers in tank2')
ylabel('The number of layers in tank1')
title('Error of the bottom temperature in tank1')
set(gcf,'position',[300,100,1200,800])

figure;
h=heatmap(num_layer_tank2,num_layer_tank1,err_upper_tank1+err_upper_tank2+err_bottom_tank2)
h.Colormap = winter;
xlabel('The number of layers in tank2')
ylabel('The number of layers in tank1')
title('Total error')
set(gcf,'position',[500,200,800,500])



% the optimal number of layers: 2+4

%%
figure;
subplot(3,1,1)
surf(num_layer_tank2, num_layer_tank1, err_upper_tank1);
colorbar;  % Adds a colorbar to the figure
xlabel('The number of layers in tank2')
ylabel('The number of layers in tank1')
zlabel('Error of the upper temperature in tank1')
subplot(3,1,2)
surf(num_layer_tank2, num_layer_tank1, err_upper_tank2);
colorbar;  % Adds a colorbar to the figure
xlabel('The number of layers in tank2')
ylabel('The number of layers in tank1')
zlabel('Error of the upper temperature in tank2')
subplot(3,1,3)
surf(num_layer_tank2, num_layer_tank1, err_bottom_tank2);
colorbar;  % Adds a colorbar to the figure
xlabel('The number of layers in tank2')
ylabel('The number of layers in tank1')
zlabel('Error of the bottom temperature in tank1')
set(gcf,'position',[500,100,700,800])

figure;
surf(num_layer_tank2, num_layer_tank1, err_upper_tank1+err_upper_tank2+err_bottom_tank2);
colorbar;  % Adds a colorbar to the figure
xlabel('The number of layers in tank2')
ylabel('The number of layers in tank1')
zlabel('Total error')
set(gcf,'position',[500,200,700,500])




%%
%% 7:00~19:00 0~12
period_time = 6;    %[h]
T_results=system_simulation_Euler_period(period_time, time_data,sample_time,num_layer_tank1,num_layer_tank2,T_out_water_he,FR_water_hp,T_buttom_water_tank2,T_upper_water_tank1,T_upper_water_tank2,V_water);





