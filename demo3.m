%% week 13 2023 
% thermal
clear all

%Read data
%1min
Temp_out_water_he_13 = readcell('week 13 outlet water T in heat exchanger.csv');
Temp_out_water_he_13=cell2table(Temp_out_water_he_13(2:end,:));
Temp_in_water_hp_13 = readcell('week 13 inlet water T in heat pump.csv');
Temp_in_water_hp_13=cell2table(Temp_in_water_hp_13(2:end,:));
Temp_out_water_hp_13 = readcell('week 13 outlet water T in heat pump.csv');
Temp_out_water_hp_13=cell2table(Temp_out_water_hp_13(2:end,:));
FlowRate_water_hp_13 = readcell('week 13 water flow rate in heat pump.csv');
FlowRate_water_hp_13=cell2table(FlowRate_water_hp_13(2:end,:));

[date_out_water_he_13,T_out_water_he_13]=read_csv(Temp_out_water_he_13);
[date_in_water_hp_13,T_in_water_hp_13]=read_csv(Temp_in_water_hp_13);
[date_out_water_hp_13,T_out_water_hp_13]=read_csv(Temp_out_water_hp_13);
[date_FR_water_hp_13,FR_water_hp_13]=read_csv(FlowRate_water_hp_13);

%1h
Prod_heat_13 = readcell('week 13 thermal energy production.csv');
Prod_heat_13=cell2table(Prod_heat_13(2:end,:));
[date_heat_13,Q_heat_13]=read_csv(Prod_heat_13);


%% water T in heat pump
figure;
plot(date_out_water_hp_13((1440+1):2*1440),T_out_water_hp_13((1440+1):2*1440),date_in_water_hp_13((1440+1):2*1440),T_in_water_hp_13((1440+1):2*1440))
hold on
plot(date_FR_water_hp_13((1440+1):2*1440),FR_water_hp_13((1440+1):2*1440)*100);
legend("outlet water","inlet water");
set(gcf,'position',[500,500,1000,400])
grid on
hold off

%%
figure;
plot(date_out_water_he_13((1440+1):2*1440),T_out_water_he_13((1440+1):2*1440));
legend("outlet water","inlet water");
set(gcf,'position',[500,500,1000,400])
grid on


%%
thermal_hp_13_min=zeros(length(date_FR_water_hp_13),1);
for i=1:length(date_FR_water_hp_13)
    if FR_water_hp_13(i)~=0
        thermal_hp_13_min(i)=(FR_water_hp_13(i)*1000/60)*(4184/3600000)*(T_out_water_hp_13(i)-T_in_water_hp_13(i));
    end
end
figure()
plot(date_FR_water_hp_13,thermal_hp_13_min)
xlabel('Time')
ylabel('Heat (kWh)')
title('Heat energy (min level)')
set(gcf,'position',[500,500,1000,400])
grid on

thermal_hp_13_hour=zeros(length(date_FR_water_hp_13)/60,1);
Q=0;
j=1;
for i=1:length(thermal_hp_13_min)
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
figure()
plot(date_heat_13,thermal_hp_13_hour,'-o',date_heat_13,Q_heat_13,'-o')
legend('Calculation based on temperatures','Thermal Meter')
xlabel('Time')
ylabel('Heat (kWh)')
title('Heat energy')
set(gcf,'position',[500,500,1000,400])
grid on
