close all
clear all

%% Read data
%1min
Temp_out_water_he = readcell('../Data/Mon 1 May 2023 outlet water T in heat exchanger.csv');
Temp_out_water_he=cell2table(Temp_out_water_he(2:end,:));
[date_out_he,T_out_he]=read_csv(Temp_out_water_he);

FlowRate_water_hp = readcell('../Data/Mon 1 May 2023 water flow rate in heat pump.csv');
FlowRate_water_hp=cell2table(FlowRate_water_hp(2:end,:));
[date_FR_hp,FR_hp]=read_csv(FlowRate_water_hp);

Temp_bottom_water_tank2 = readcell('../Data/Mon 1 May 2023 water T in bottom 2 tank.csv');
Temp_bottom_water_tank2=cell2table(Temp_bottom_water_tank2(2:end,:));
[date_bottom_tank2,T_bottom_tank2]=read_csv(Temp_bottom_water_tank2);

Temp_upper_water_tank1 = readcell('../Data/Mon 1 May 2023 water T in upper 1 tank.csv');
Temp_upper_water_tank1=cell2table(Temp_upper_water_tank1(2:end,:));
[date_upper_tank1,T_upper_tank1]=read_csv(Temp_upper_water_tank1);


%8min
Temp_upper_water_tank2 = readcell('../Data/Mon 1 May 2023 water T in upper 2 tank.csv');
Temp_upper_water_tank2=cell2table(Temp_upper_water_tank2(2:end,:));
[date_upper_tank2,T_upper_tank2]=read_csv(Temp_upper_water_tank2);

%1h
Consp_water = readcell('../Data/Mon 1 May 2023 water V consumption.csv');
Consp_water=cell2table(Consp_water(2:end,:));
[date_V_water,V_water]=read_csv(Consp_water);  %[m^3]



%%
% set the number of layers

num_layer_tank1=2;
num_layer_tank2=4;

% initial masses of layers   L~kg
m_layer_tank1=[250;250];
m_layer_tank2=[100;150;150;100];

R_tank12=0;
m_c_dot= 1*1000; % flow rate of hot water circulation    [m^3/h]=1000*[L/h]=1000*[kg/h]
diff_T_c = 2; % temperature difference  [K]
T_s = 12; % cold water temperature      [°C]


%% Choose data between 7:05 ~ 18:30
%7:05 ~ 18:30
%0+5/60~11.5
% Time span for the simulation
t_start = 0+5/60;
t_end = 11.5;     % end time [h]
sample_time = 5;    % [min]
time_span = [t_start:sample_time/60:t_end];
time_data = time_span;

%%
%1min
m_p_dot_data=FR_hp(60*7+1+5:60*18.5+1)*1000;  %[L/h] 60*7+1+5~60*(7+11.5)+1
T_in_data=T_out_he(60*7+1+5:60*18.5+1);
T_upper_tank1_data = T_upper_tank1(60*7+1+5:60*18.5+1);
T_bottom_tank2_data = T_bottom_tank2(60*7+1+5:60*18.5+1);

%8min
T_upper_tank2_data = T_upper_tank2(round(60/8*7)+1:round(60/8*18.5)+1); %07:04~18:32

%1h
m_s_dot_data=V_water(1*8:1*20)*1000;    %[L/h] 07:00~19:00


[m_p_dot_interp,T_in_interp,m_s_dot_interp,T_upper_tank1_interp,T_upper_tank2_interp,...
    T_bottom_tank2_interp]=interpolation_data(time_data,sample_time,m_p_dot_data,T_in_data,m_s_dot_data,T_upper_tank1_data,T_upper_tank2_data,T_bottom_tank2_data);

% Define ground truth data
% T_ground_truth = [T_upper_tank1_interp',T_upper_tank2_interp',T_bottom_tank2_interp']; % Your historical data here


%% 
% period time
period_time = 3;    %[h]
sample_num = period_time/(sample_time/60);

% time_data_period = zeros((sample_num+1),12/period_time);
% m_p_dot_interp_period = zeros(sample_num,12/period_time);
% T_in_interp_period = zeros(sample_num,12/period_time);
% m_s_dot_interp_period = zeros(sample_num,12/period_time);
% 
% T_upper_tank1_interp_period = zeros(sample_num,12/period_time);
% T_upper_tank2_interp_period = zeros(sample_num,12/period_time);
% T_bottom_tank2_interp_period = zeros(sample_num,12/period_time);
% 
% 
% for i=1:12/period_time
%     time_data_period(:,i) = [t_start+period_time*(i-1):sample_time/60:t_start+period_time*i];
%     m_p_dot_interp_period(:,i) = m_p_dot_interp(1+(i-1)*sample_num:i*sample_num);
%     T_in_interp_period(:,i) = T_in_interp(1+(i-1)*sample_num:i*sample_num);
%     m_s_dot_interp_period(:,i) = m_s_dot_interp(1+(i-1)*sample_num:i*sample_num);
% 
%     T_upper_tank1_interp_period(:,i) =  T_upper_tank1_interp(1+(i-1)*sample_num:i*sample_num);
% T_upper_tank2_interp_period(:,i) =  T_upper_tank2_interp(1+(i-1)*sample_num:i*sample_num);
% T_bottom_tank2_interp_period(:,i) =  T_bottom_tank2_interp(1+(i-1)*sample_num:i*sample_num);
% 
% 
% end


%%
load("optimal_params.mat")
R_tank12_optimal = P_optimal(1);
m_layer_tank2_optimal=zeros(num_layer_tank2,1);
for i=1:num_layer_tank2
    m_layer_tank2_optimal(i) = P_optimal(i+1);
end

diff_T_c_optimal = P_optimal((1+num_layer_tank2)+1);
T_s_optimal = P_optimal((1+num_layer_tank2)+2);
m_c_dot_optimal = P_optimal((1+num_layer_tank2)+3);



%%
T_results=[];
initial_index=1;
for i=1:length(time_data)
    if (rem(i,sample_num)==1)
        initial_index=i;
        T_initial=initial_T_setting(num_layer_tank1,num_layer_tank2,m_layer_tank1,m_layer_tank2,T_upper_tank1_interp(initial_index),T_upper_tank2_interp(initial_index),T_bottom_tank2_interp(initial_index));
    end

    if (rem(i-initial_index+1,sample_num)==0)
        T_results_period = system_of_equations_Euler_Nlayers(time_data(initial_index:i), T_initial, m_p_dot_interp(initial_index:i), T_in_interp(initial_index:i), m_s_dot_interp(initial_index:i), num_layer_tank1,num_layer_tank2, m_layer_tank1, m_layer_tank2_optimal, R_tank12_optimal, m_c_dot_optimal, diff_T_c_optimal, T_s_optimal);
        T_results = [T_results;T_results_period];
    end

    if (i==length(time_data))
        T_results_period = system_of_equations_Euler_Nlayers(time_data(initial_index:i), T_initial, m_p_dot_interp(initial_index:i), T_in_interp(initial_index:i), m_s_dot_interp(initial_index:i), num_layer_tank1,num_layer_tank2, m_layer_tank1, m_layer_tank2_optimal, R_tank12_optimal, m_c_dot_optimal, diff_T_c_optimal, T_s_optimal);
        T_results = [T_results;T_results_period];
    end

        
end



plot_comparison_T_period(time_data,T_results,period_time,num_layer_tank1,num_layer_tank2,T_upper_tank1_data,T_upper_tank2_data,T_bottom_tank2_data)

%%

save("model_T_010523","T_results")