%% 2023-03-27 week 19
clear all
close all

%Read data
%1min
Temp_out_water_he_19 = readcell('week 19 outlet water T in heat exchanger.csv');
Temp_out_water_he_19=cell2table(Temp_out_water_he_19(2:end,:));
[date_out_water_he_19,T_out_water_he_19]=read_csv(Temp_out_water_he_19);
Temp_out_water_hp_19 = readcell('week 19 outlet water T in heat pump.csv');
Temp_out_water_hp_19=cell2table(Temp_out_water_hp_19(2:end,:));
[date_out_water_hp_19,T_out_water_hp_19]=read_csv(Temp_out_water_hp_19);

FlowRate_water_hp_19 = readcell('week 19 water flow rate in heat pump.csv');
FlowRate_water_hp_19=cell2table(FlowRate_water_hp_19(2:end,:));
[date_FR_water_hp_19,FR_water_hp_19]=read_csv(FlowRate_water_hp_19);


%%

difference_he_19=zeros(length(date_FR_water_hp_19),1);
for i=1:length(date_FR_water_hp_19)
    if FR_water_hp_19(i)~=0
        difference_he_19(i)=T_out_water_hp_19(i)-T_out_water_he_19(i);
    end
end
figure()
plot(date_FR_water_hp_19(60*7+1:60*19+1),difference_he_19(60*7+1:60*19+1))
xlabel('Time')
ylabel('Temperature Difference (°C)')
yyaxis right
plot(date_FR_water_hp_19(60*7+1:60*19+1),FR_water_hp_19(60*7+1:60*19+1));
ylabel('water flow rate (m^3/h)')
title('Temperature Difference between outlet water of the heat pump and outlet water of the heat exchanger')
set(gcf,'position',[500,500,1000,400])
grid on

