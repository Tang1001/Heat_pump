%% week 12 2023 and week 50 2020
% water consumption and heat consumption

clear all

%Read data
%1h
V_water_12 = readcell('week 12 water consumption.csv');
V_water_12=cell2table(V_water_12(2:end,:));
V_water_50 = readcell('week 50 water consumption.csv');
V_water_50=cell2table(V_water_50(2:end,:));
[date_water_12,Consp_water_12]=read_csv(V_water_12);
[date_water_50,Consp_water_50]=read_csv(V_water_50);


% process data  
date_water_12=date_water_50;
Consp_water_12=[Consp_water_12(1:146);0;Consp_water_12(147:end)];

%8min
T_water_12 = readcell('week 12 hot water T.csv');
T_water_12 = cell2table(T_water_12(2:end,:));
[date_water_12_Temp,Temp_water_12]=read_csv(T_water_12);
T_water_50 = readcell('week 50 hot water T.csv');
T_water_50 = cell2table(T_water_50(2:end,:));
[date_water_50_Temp,Temp_water_50]=read_csv(T_water_50);



% 8min to 1hour
[aver_date_water_12_Temp,aver_Temp_water_12]=aver_min2hour(date_water_12_Temp,Temp_water_12);
[aver_date_water_50_Temp,aver_Temp_water_50]=aver_min2hour(date_water_50_Temp,Temp_water_50);
% aver_date_water_12_Temp=aver_date_water_12_Temp';
% aver_date_water_50_Temp=aver_date_water_50_Temp';

% process data  
aver_date_water_12_Temp=aver_date_water_50_Temp;
aver_Temp_water_12=[aver_Temp_water_12(1:146),aver_Temp_water_12(146),aver_Temp_water_12(147:end)];


%%
% figure(1)
% plot(date_water_12,Consp_water_12)
% % plot(date_water_12(1:24),Consp_water_12(1:24))  %one day
% xlabel('Time')
% ylabel('Water Total')
% set(gcf,'position',[500,500,1000,400])
% grid on

close all
Consp_heat_12=Consp_water_12*1000*4184.*(aver_Temp_water_12'-12); %J
Consp_heat_12=Consp_heat_12/3600000;    %kWh
figure()
plot(date_water_12,Consp_heat_12,'-*')
% plot(date_water_12(1:24),Consp_water_12(1:24))  %one day
xlabel('Time')
ylabel('Heat (kWh)')
title('Week 12 heat consumption')
set(gcf,'position',[500,500,1000,400])
grid on


Consp_heat_12_heatmap=reshape(Consp_heat_12,[24,7])';
figure()
heatmap(Consp_heat_12_heatmap)
xlabel('hour')
ylabel('day')
title('Week 12 heat consumption')
set(gcf,'position',[500,500,1000,400])



Consp_heat_50=Consp_water_50*1000*4184.*(aver_Temp_water_50'-12); %J
Consp_heat_50=Consp_heat_50/3600000;    %kWh
figure()
plot(date_water_50,Consp_heat_50,'-*')
% plot(date_water_12(1:24),Consp_water_12(1:24))  %one day
xlabel('Time')
ylabel('Heat (kWh)')
title('Week 50 heat consumption')
set(gcf,'position',[500,500,1000,400])
grid on


Consp_heat_50_heatmap=reshape(Consp_heat_50,[24,7])';
figure()
heatmap(Consp_heat_50_heatmap)
xlabel('hour')
ylabel('day')
title('Week 50 heat consumption')
set(gcf,'position',[500,500,1000,400])

%%
figure()
subplot(2,1,1)
plot(date_water_12,Consp_water_12,'-*')
xlabel('Time')
ylabel('Water (m^3)')
title('Week 12 water consumption')
grid on
subplot(2,1,2)
plot(date_water_12,Consp_heat_12,'-*')
% plot(date_water_12(1:24),Consp_water_12(1:24))  %one day
xlabel('Time')
ylabel('Heat (kWh)')
title('Week 12 heat consumption')
grid on
set(gcf,'position',[500,500,1000,400])


figure()
subplot(2,1,1)
xvalues = {'0:00','1:00','2:00','3:00','4:00','5:00','6:00','7:00','8:00','9:00','10:00','11:00','12:00','13:00','14:00','15:00','16:00','17:00','18:00','19:00','20:00','21:00','22:00','23:00'};
yvalues = {'Monday','Tuesday','Wednesday','Thursday','Firday','Saturday','Sunday'};
heatmap(xvalues,yvalues,Consp_heat_12_heatmap)
% xlabel('hour')
% ylabel('day')
title('Week 12 heat consumption (kWh)')

subplot(2,1,2)
xvalues = {'0:00','1:00','2:00','3:00','4:00','5:00','6:00','7:00','8:00','9:00','10:00','11:00','12:00','13:00','14:00','15:00','16:00','17:00','18:00','19:00','20:00','21:00','22:00','23:00'};
yvalues = {'Monday','Tuesday','Wednesday','Thursday','Firday','Saturday','Sunday'};
heatmap(xvalues,yvalues,Consp_heat_50_heatmap)
% xlabel('hour')
% ylabel('day')
title('Week 50 heat consumption (kWh)')
set(gcf,'position',[500,500,1000,400])