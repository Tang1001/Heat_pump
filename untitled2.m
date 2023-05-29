%% week 13 2023 
% temperature of inlet water in heat exchanger
clear all

%Read data
%1min
Temp_in_water_hp_13 = readcell('week 13 inlet water T in heat pump.csv');
Temp_in_water_hp_13=cell2table(Temp_in_water_hp_13(2:end,:));
Temp_out_water_hp_13 = readcell('week 13 outlet water T in heat pump.csv');
Temp_out_water_hp_13=cell2table(Temp_out_water_hp_13(2:end,:));
FlowRate_water_hp_13 = readcell('week 13 water flow rate in heat pump.csv');
FlowRate_water_hp_13=cell2table(FlowRate_water_hp_13(2:end,:));
Temp_out_water_he_13 = readcell('week 13 outlet water T in heat exchanger.csv');
Temp_out_water_he_13=cell2table(Temp_out_water_he_13(2:end,:));

Temp_b_water_2tank_13 = readcell('week 13 water T in bottom 2 tank.csv');
Temp_b_water_2tank_13=cell2table(Temp_b_water_2tank_13(2:end,:));

Consp_water_13 = readcell('week 13 water V consumption.csv');
Consp_water_13=cell2table(Consp_water_13(2:end,:));


[date_in_water_hp_13,T_in_water_hp_13]=read_csv(Temp_in_water_hp_13);
[date_out_water_hp_13,T_out_water_hp_13]=read_csv(Temp_out_water_hp_13);
[date_FR_water_hp_13,FR_water_hp_13]=read_csv(FlowRate_water_hp_13);
[date_out_water_he_13,T_out_water_he_13]=read_csv(Temp_out_water_he_13);

[date_b_water_2tank_13,T_b_water_2tank_13]=read_csv(Temp_b_water_2tank_13);

[date_water_13,V_water_13]=read_csv(Consp_water_13);

%%
%Q=dot(m)*Cp*T


T_in_water_he_13_min=zeros(length(date_FR_water_hp_13),1);

for i=1:length(date_FR_water_hp_13)
    if FR_water_hp_13(i)~=0 && T_in_water_hp_13(i)<62 && (T_out_water_hp_13(i)-T_in_water_hp_13(i))>0
%         T_in_water_he_13_min(i)=(T_out_water_hp_13(i)-T_in_water_hp_13(i));
        T_in_water_he_13_min(i)=T_out_water_he_13(i)-0.9*(T_out_water_hp_13(i)-T_in_water_hp_13(i));
        if T_in_water_he_13_min(i)<0
            T_in_water_he_13_min(i)=0;
        end
    end
end

figure()
% yyaxis left
% plot(date_out_water_he_13,T_in_water_he_13_min,date_out_water_he_13,T_out_water_he_13)
% plot(date_out_water_he_13,T_out_water_he_13,date_out_water_he_13,(T_out_water_hp_13-T_in_water_hp_13))
plot(date_out_water_he_13,T_in_water_he_13_min,date_b_water_2tank_13,T_b_water_2tank_13)
% yyaxis right
% plot(date_water_13,V_water_13)
% xlabel('Time')
% ylabel('Heat (kWh)')
% title('Heat energy (min level)')
set(gcf,'position',[500,500,1800,400])
grid on


%% calculate the real time of supplying water

FR_water_supply_13=zeros(length(FR_water_hp_13),1);

for i=1:length(FR_water_supply_13)
    if FR_water_hp_13(i)~=0
        FR_water_supply_13(i)=FR_water_hp_13(i)*(T_b_water_2tank_13(i)-T_in_water_he_13_min(i))/(T_in_water_he_13_min(i)-12);
    end
end
  

df=zeros(length(T_in_water_he_13_min),1);
for i=1:length(T_in_water_he_13_min)
    if FR_water_hp_13(i)~=0
        df(i)=T_in_water_he_13_min(i)-12;
    end
end

figure()
% plot(date_out_water_he_13,T_in_water_he_13_min,date_out_water_he_13,T_b_water_2tank_13)
plot(T_in_water_he_13_min)
% plot(T_b_water_2tank_13-T_in_water_he_13_min)
set(gcf,'position',[500,500,1800,400])
grid on