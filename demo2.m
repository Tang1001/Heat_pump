%% week 13 2023 and Mon 17 April 2023
% COP

clear all

%Read data
%1h
Prod_heat_13 = readcell('week 13 thermal energy production.csv');
Prod_heat_13=cell2table(Prod_heat_13(2:end,:));
Prod_elec_13 = readcell('week 13 electrical energy production.csv');
Prod_elec_13=cell2table(Prod_elec_13(2:end,:));
[date_heat_13,Q_heat_13]=read_csv(Prod_heat_13);
[date_elec_13,W_elec_13]=read_csv(Prod_elec_13);


Prod_heat_elec_23_4_17 = readcell('Mon 17 April 2023 heat and electricity.csv');
Prod_heat_23_4_17=cell2table(Prod_heat_elec_23_4_17(2:end,[1,2]));
Prod_elec_23_4_17=cell2table(Prod_heat_elec_23_4_17(2:end,[1,3]));
[date_heat_23_4_17,Q_heat_23_4_17]=read_csv(Prod_heat_23_4_17);
[date_elec_23_4_17,W_elec_23_4_17]=read_csv(Prod_elec_23_4_17);

%%
figure()
plot(date_heat_13,Q_heat_13/10,'-o')
% plot(date_water_12(1:24),Consp_water_12(1:24))  %one day
xlabel('Time')
ylabel('Heat (kWh)')
title('Week 13 heat production')
set(gcf,'position',[500,500,1000,400])
grid on

figure()
plot(date_elec_13,W_elec_13,'-o')
% plot(date_water_12(1:24),Consp_water_12(1:24))  %one day
xlabel('Time')
ylabel('Heat (kWh)')
title('Week 13 electricity production')
set(gcf,'position',[500,500,1000,400])
grid on

figure()
plot(date_elec_23_4_17,W_elec_23_4_17,'-o',date_heat_23_4_17,Q_heat_23_4_17,'-o')
% plot(date_water_12(1:24),Consp_water_12(1:24))  %one day
xlabel('Time')
ylabel('Heat (kWh)')
title('17/4/23 energy production')
set(gcf,'position',[500,500,1000,400])
grid on

%%
figure()
plot(date_elec_23_4_17,Q_heat_23_4_17./W_elec_23_4_17,'-o')
% plot(date_water_12(1:24),Consp_water_12(1:24))  %one day
xlabel('Time')
ylabel('Heat (kWh)')
title('COP')
set(gcf,'position',[500,500,1000,400])
grid on


%还是要用之前实验的数据去拟合COP方程