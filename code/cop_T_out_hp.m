%%COP thermal
% thermal
clear all

%Read data
%1min
Temp_out_he_13 = readcell('../Data/week 13 outlet water T in heat exchanger.csv');
Temp_out_he_13=cell2table(Temp_out_he_13(2:end,:));
Temp_in_hp_13 = readcell('../Data/week 13 inlet water T in heat pump.csv');
Temp_in_hp_13=cell2table(Temp_in_hp_13(2:end,:));
Temp_out_hp_13 = readcell('../Data/week 13 outlet water T in heat pump.csv');
Temp_out_hp_13=cell2table(Temp_out_hp_13(2:end,:));
FlowRate_hp_13 = readcell('../Data/week 13 water flow rate in heat pump.csv');
FlowRate_hp_13=cell2table(FlowRate_hp_13(2:end,:));

[date_out_he_13,T_out_he_13]=read_csv(Temp_out_he_13);
[date_in_hp_13,T_in_hp_13]=read_csv(Temp_in_hp_13);
[date_out_hp_13,T_out_hp_13]=read_csv(Temp_out_hp_13);
[date_FR_hp_13,FR_hp_13]=read_csv(FlowRate_hp_13);

%1h
Prod_heat_13 = readcell('../Data/week 13 thermal energy production.csv');
Prod_heat_13=cell2table(Prod_heat_13(2:end,:));
[date_heat_13,Q_heat_13]=read_csv(Prod_heat_13);


%%
%m cp dT/dt = m_dot * cp * (T1 - T2) = Q_dot 
%[kg] * [J/kg·K] [K/s] = [kg/s] * [J/kg·K] * [K] = [J/s] = [W]
%[kg] * [J/kg·K] [K/h] = [kg/h] * [J/kg·K] * [K] = [J/h] = 1/(60*60) * [J/s] = 1/(60*60) * [W] = 1/(60*60) * [J/s]


%P
P_hp = 10*1000;  %[kW]

T_out_hp_13_diff = zeros(length(T_out_hp_13),1);
T_out_hp_13_es=zeros(length(T_out_hp_13),1);
for i=1:length(date_FR_hp_13)
    if FR_hp_13(i)~=0
        COP=COP1(T_in_hp_13(i));
        Q=COP*P_hp;
        T_out_hp_13_es(i)=Q/(0.8*1000/3600*4184)+T_in_hp_13(i);
        T_out_hp_13_diff(i) = T_out_hp_13(i)-T_out_hp_13_es(i);
        if T_out_hp_13_es(i)>500
            Q
            FR_hp_13(i)
            T_in_hp_13(i)
        end
    end
end

figure()
plot(date_FR_hp_13(60*7+1:60*19+1),T_out_hp_13(60*7+1:60*19+1),'-o',date_FR_hp_13(60*7+1:60*19+1),T_out_hp_13_es(60*7+1:60*19+1),'-o')
legend('Measurement','Estimate')
xlabel('Time')
ylabel('Temp (°C)')
set(gcf,'position',[500,500,1000,400])
grid on


figure()
plot(date_FR_hp_13(60*7+1:60*19+1),T_out_hp_13_diff(60*7+1:60*19+1))
legend('Measurement','Estimate')
xlabel('Time')
ylabel('Temp (°C)')
set(gcf,'position',[500,500,1000,400])
grid on

%%
figure()
plot(date_FR_hp_13(60*7+1:60*19+1),FR_hp_13(60*7+1:60*19+1))
xlabel('Time')
ylabel('Temp (°C)')
set(gcf,'position',[500,500,1000,400])
grid on
