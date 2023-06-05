clear all
close all
%% Read data
%1min
Temp_out_he = readcell('../Data/Mon 1 May 2023 outlet water T in heat exchanger.csv');
Temp_out_he=cell2table(Temp_out_he(2:end,:));
[date_out_he,T_out_he]=read_csv(Temp_out_he);

FlowRate_hp = readcell('../Data/Mon 1 May 2023 water flow rate in heat pump.csv');
FlowRate_hp=cell2table(FlowRate_hp(2:end,:));
[date_FR_hp,FR_hp]=read_csv(FlowRate_hp);

Temp_in_hp = readcell('../Data/Mon 1 May 2023 inlet water T in heat pump.csv');
Temp_in_hp=cell2table(Temp_in_hp(2:end,:));
[date_in_hp,T_in_hp]=read_csv(Temp_in_hp);

Temp_out_hp = readcell('../Data/Mon 1 May 2023 outlet water T in heat pump.csv');
Temp_out_hp=cell2table(Temp_out_hp(2:end,:));
[date_out_hp,T_out_hp]=read_csv(Temp_out_hp);



Temp_buttom_tank2 = readcell('../Data/Mon 1 May 2023 water T in bottom 2 tank.csv');
Temp_buttom_tank2=cell2table(Temp_buttom_tank2(2:end,:));
[date_buttom_tank2,T_buttom_tank2]=read_csv(Temp_buttom_tank2);

Temp_upper_tank1 = readcell('../Data/Mon 1 May 2023 water T in upper 1 tank.csv');
Temp_upper_tank1=cell2table(Temp_upper_tank1(2:end,:));
[date_upper_tank1,T_upper_tank1]=read_csv(Temp_upper_tank1);


Temp_ambient = readcell('../Data/Mon 1 May 2023 ambient air T.csv');
Temp_ambient=cell2table(Temp_ambient(2:end,:));
[date_mbient,T_amb]=read_csv(Temp_ambient);

%8min
Temp_upper_tank2 = readcell('../Data/Mon 1 May 2023 water T in upper 2 tank.csv');
Temp_upper_tank2=cell2table(Temp_upper_tank2(2:end,:));
[date_upper_tank2,T_upper_tank2]=read_csv(Temp_upper_tank2);

%1h
Consp_water = readcell('../Data/Mon 1 May 2023 water V consumption.csv');
Consp_water=cell2table(Consp_water(2:end,:));
[date_V_water,V_water]=read_csv(Consp_water);  %[m^3]


%% Choose data between 7:00~19:00
%1min
m_p_dot_data=FR_hp(60*7+1:60*19+1)*1000;  %[L/h] 60*7+1~60*19
T_in_data=T_out_he(60*7+1:60*19+1);
T_amb_data=T_amb(60*7+1:60*19+1);
T_upper_water_tank1_data = T_upper_tank1(60*7+1:60*19+1);
T_buttom_water_tank2_data = T_buttom_tank2(60*7+1:60*19+1);


T_in_hp_data = T_in_hp(60*7+1:60*19+1);
T_out_hp_data = T_out_hp(60*7+1:60*19+1);


%8min
T_upper_water_tank2_data = T_upper_tank2(round(60/8*7):round(60/8*19)+1); %06:56~19:04


%1h
m_s_dot_data=V_water(1*8:1*20)*1000;    %[L/h] 07:00~19:00


%% Time span for the simulation
t_start = 0;
t_end = 12;     % end time [h]
sample_time = 2;    % [min]
time_span = [t_start:sample_time/60:t_end];
time_data = time_span;

time_data_1min=0:1/60:12;
time_data_1h=0:12;
time_data_8min=0-4/60:8/60:12+4/60;

%%
%m cp dT/dt = m_dot * cp * (T1 - T2) = Q_dot 
%[kg] * [J/kg·K] [K/s] = [kg/s] * [J/kg·K] * [K] = [J/s] = [W]
%[kg] * [J/kg·K] [K/h] = [kg/h] * [J/kg·K] * [K] = [J/h] = 1/(60*60) * [J/s] = 1/(60*60) * [W] = 1/(60*60) * [J/s]


%P
P_hp = 10*1000;  %[kW]

T_out_hp_diff1 = zeros(length(m_p_dot_data),1);
T_out_hp_es1=zeros(length(m_p_dot_data),1);
T_out_hp_diff2 = zeros(length(m_p_dot_data),1);
T_out_hp_es2=zeros(length(m_p_dot_data),1);

for i=1:length(m_p_dot_data)
    if m_p_dot_data(i)~=0
        n1=COP1(T_in_hp_data(i));
        Q1=n1*P_hp;
        T_out_hp_es1(i)=Q1/(0.8*1000/3600*4184)+T_in_hp_data(i);
        T_out_hp_diff1(i) = T_out_hp_data(i)-T_out_hp_es1(i);

        n2=COP2(T_in_hp_data(i),T_amb_data(i));
        Q2=n2*P_hp;
        T_out_hp_es2(i)=Q2/(0.8*1000/3600*4184)+T_in_hp_data(i);
        T_out_hp_diff2(i) = T_out_hp_data(i)-T_out_hp_es2(i);

%         if T_out_hp_es2(i)>500
%             Q
%             FR_hp(i)
%             T_in_hp(i)
%         end
    end
end


figure()
plot(time_data_1min,T_out_hp_data,'-o',time_data_1min,T_out_hp_es1,'-o',time_data_1min,T_out_hp_es2,'-o')
legend('Measurement','Estimate1','Estimate2')
xlabel('Time')
ylabel('Temp (°C)')
set(gcf,'position',[500,500,1000,400])
grid on


figure()
plot(time_data_1min,T_out_hp_diff1,time_data_1min,T_out_hp_diff2)
legend('Estimate1','Estimate2')
xlabel('Time')
ylabel('Temp (°C)')
set(gcf,'position',[500,500,1000,400])
grid on

