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


Temp_out_water_hp = readcell('../Data/Mon 1 May 2023 outlet water T in heat pump.csv');
Temp_out_water_hp=cell2table(Temp_out_water_hp(2:end,:));
[date_out_hp,T_out_hp]=read_csv(Temp_out_water_hp);

Temp_in_water_hp = readcell('../Data/Mon 1 May 2023 inlet water T in heat pump.csv');
Temp_in_water_hp=cell2table(Temp_in_water_hp(2:end,:));
[date_in_hp,T_in_hp]=read_csv(Temp_in_water_hp);



%8min
Temp_upper_water_tank2 = readcell('../Data/Mon 1 May 2023 water T in upper 2 tank.csv');
Temp_upper_water_tank2=cell2table(Temp_upper_water_tank2(2:end,:));
[date_upper_tank2,T_upper_tank2]=read_csv(Temp_upper_water_tank2);

%1h
Consp_water = readcell('../Data/Mon 1 May 2023 water V consumption.csv');
Consp_water=cell2table(Consp_water(2:end,:));
[date_V_water,V_water]=read_csv(Consp_water);  %[m^3]

Heat_thermal = readcell('../Data/Mon 1 May 2023 thermal energy.csv');
Heat_thermal=cell2table(Heat_thermal(2:end,:));
[date_thermal,Q_thermal]=read_csv(Heat_thermal);


%% Choose data between 7:00 ~ 19:00
%7:00 ~ 19:00
%0~12
% Time span for the simulation
t_start = 0;
t_end = 12;     % end time [h]
sample_time = 2;    % [min]
time_span = [t_start:sample_time/60:t_end];
time_data = time_span;

time_data_1min=time_data(1):1/60:time_data(end);
time_data_1h=time_data(1):time_data(end);
time_data_8min=time_data(1)-4/60:8/60:time_data(end)+4/60;

%%
%1min
m_p_dot_data=FR_hp(60*7+1:60*19+1)*1000;  %[L/h] 60*7+1~60*19
T_in_data=T_out_he(60*7+1:60*19+1);
T_upper_tank1_data = T_upper_tank1(60*7+1:60*19+1);
T_bottom_tank2_data = T_bottom_tank2(60*7+1:60*19+1);

T_out_hp_data = T_out_hp(60*7+1:60*19+1);
T_in_hp_data = T_in_hp(60*7+1:60*19+1);


%8min
T_upper_tank2_data = T_upper_tank2(round(60/8*7):round(60/8*19)+1); %06:56~19:04

%1h
m_s_dot_data=V_water(1*8:1*20)*1000;    %[L/h] 07:00~19:00
Q_hp_data=Q_thermal(1*8:1*20);    %[kWh] 07:00~19:00



[m_p_dot_interp,T_in_interp,m_s_dot_interp,T_upper_tank1_interp,T_upper_tank2_interp,...
    T_bottom_tank2_interp]=interpolation_data(time_data,sample_time,m_p_dot_data,T_in_data,m_s_dot_data,T_upper_tank1_data,T_upper_tank2_data,T_bottom_tank2_data);



T_out_hp_interp = interp1(time_data_1min,T_out_hp_data,time_data);
T_in_hp_interp = interp1(time_data_1min,T_in_hp_data,time_data);


load("optimal_params.mat")  % P_optimal: R_tank12_optimal  m_layer_tank2_optimal  diff_T_c_optimal T_s_optimal m_c_dot_optimal[L/h]
load("model_T_010523.mat")  % T_results_0_12 (i,6)

%%
Q_he_estimate_min=zeros(length(m_p_dot_interp),1);

diff_T_he = 2.5; % difference between inlet water on heat pump side and outlet water buffer side
T_out_he_es=zeros(length(m_p_dot_interp),1);
T_in_he_es=zeros(length(m_p_dot_interp),1);

for i=1:length(m_p_dot_interp)
    if m_p_dot_interp(i)~=0
        
        % predictive temperature of inlet water in the heat exchanger
        T_in_he = T_results_0_12(i,6)*(1-m_s_dot_interp(i)/m_p_dot_interp(i))+P_optimal(1+4+2)*m_s_dot_interp(i)/m_p_dot_interp(i);
        % predictive temperature of outlet water in the heat exchanger
        T_out_he =T_out_hp_interp(i)-diff_T_he;
        T_out_he_es(i) = T_out_he;
        T_in_he_es(i) = T_in_he;
%         T_out_he = T_in_interp(i);

        if T_in_he>T_out_he
            continue;
        end


        %[kg/h]*[J/kg·K]*[K]=[J/h]=1/(60*60) * [J/s]=1/(60*60) * [W] = =1/(60*60) * 1/1000 *[kW] 
        %[kW]*sample/60 hour = [kWh]
        Q_he_estimate_min(i)=m_p_dot_interp(i)*4184*(T_out_he-T_in_he)*1/(60*60)* 1/1000 * sample_time/60;
    end
end

Q_he_estimate_hour=zeros(12,1);
Q=0;
j=1;
for i=1:length(Q_he_estimate_min)
    if j<30
        Q=Q+Q_he_estimate_min(i);
        j=j+1;
    else
        Q=Q+Q_he_estimate_min(i);
        Q_he_estimate_hour(i/30)=Q;
        Q=0;
        j=1;
    end
end

%%
thermal_hp_13_min=zeros(length(m_p_dot_data),1);
for i=1:length(m_p_dot_data)
    if m_p_dot_data(i)~=0
        thermal_hp_13_min(i)=m_p_dot_data(i)*4184*(T_out_hp_data(i)-T_in_hp_data(i))*1/(60*60)* 1/1000*1/60;
    end
end
thermal_hp_13_hour=zeros((length(m_p_dot_data)-1)/60,1);
Q=0;
j=1;
for i=1:length(thermal_hp_13_min)-1
    if j<60
        Q=Q+thermal_hp_13_min(i);
        j=j+1;
    else
        Q=Q+thermal_hp_13_min(i);
        thermal_hp_13_hour(i/60)=Q;
        Q=0;
        j=1;
    end

end

%%
err1=sqrt(mean((Q_hp_data(2:end) - Q_he_estimate_hour).^2))
err2=sqrt(mean((Q_hp_data(2:end) - thermal_hp_13_hour).^2)) 
figure;
plot(time_data_1h,Q_hp_data,'-o',time_data_1h,[0;Q_he_estimate_hour],'-o',time_data_1h,[0;thermal_hp_13_hour],'-o')
legend('Thermal Meter','Calculation based on temperatures')
xlabel('Time (0=07:00, 12=19:00)')
ylabel('Heat (kWh)')
title("Root Mean Squared Error: ", err1)
set(gcf,'position',[500,500,1000,400])
grid on



%%

T_in_hp = zeros(length(m_p_dot_interp),1);
T_in_hp_diff = zeros(length(m_p_dot_interp)-1,1);
for i=1:length(m_p_dot_interp)
    if m_p_dot_interp(i)~=0
        
        % predictive temperature of inlet water in the heat exchanger
        T_in_he = T_results_0_12(i,6)*(1-m_s_dot_interp(i)/m_p_dot_interp(i))+P_optimal(1+4+2)*m_s_dot_interp(i)/m_p_dot_interp(i);
        % predictive temperature of outlet water in the heat exchanger
        T_out_he =T_out_hp_interp(i)-diff_T_he;
%         T_out_he = T_in_interp(i);
%         T_in_hp(i) = T_out_hp_interp(i) - (T_out_he-T_in_he);
%         T_in_hp(i) = T_in_he-diff_T_he;
        T_in_hp_diff(i) = T_in_hp(i-1)-T_in_hp_interp(i);
        
    end
end

% T_in_hp_diff = T_in_hp(2:end) - T_in_hp_interp(1:end-1)';

figure;
plot(time_data,T_in_hp_interp,'-o',time_data,T_in_hp,'-o')
legend('Measurement','Calculation based on temperatures')
xlabel('Time (0=07:00, 12=19:00)')
ylabel('Temp (°C)')
set(gcf,'position',[500,500,1000,400])
grid on

figure;
plot(time_data,[T_in_hp_diff;0])
xlabel('Time (0=07:00, 12=19:00)')
ylabel('Temp (°C)')
set(gcf,'position',[500,500,1000,400])
grid on

%%
T_in_hp = zeros(length(m_p_dot_interp),1);
T_in_hp_diff = zeros(length(m_p_dot_interp)-1,1);
T_in_he_es2=zeros(length(m_p_dot_interp),1);
for i=1:length(m_p_dot_interp)
    if m_p_dot_interp(i)~=0
        
        
        T_out_he =T_out_hp_interp(i)-diff_T_he;
        T_in_he_es2(i) = T_out_he - (T_out_hp_interp(i) - T_in_hp_interp(i));
%         T_out_he = T_in_interp(i);
%         T_in_hp(i) = T_out_hp_interp(i) - (T_out_he-T_in_he);
%         T_in_hp(i) = T_in_he-diff_T_he;
        
    end
end

figure;
plot(time_data,T_in_he_es2,time_data,T_in_he_es)
xlabel('Time (0=07:00, 12=19:00)')
ylabel('Temp (°C)')
set(gcf,'position',[500,500,1000,400])
grid on

%%
figure;
plot(time_data,T_out_he_es,'-o',time_data,T_in_interp,'-o')
xlabel('Time (0=07:00, 12=19:00)')
ylabel('Temp (°C)')
set(gcf,'position',[500,500,1000,400])
grid on

