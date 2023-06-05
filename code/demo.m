figure;
plot(time_data,T_out_he_es,'-o',time_data,T_in_interp,'-o')
legend('Estimate','Measurement')
xlabel('Time (0=07:00, 12=19:00)')
ylabel('Temp (°C)')
set(gcf,'position',[500,500,1000,400])
grid on

figure;
plot(time_data,T_out_he_es-T_in_interp','-o')
legend('Estimate','Measurement')
xlabel('Time (0=07:00, 12=19:00)')
ylabel('Temp (°C)')
set(gcf,'position',[500,500,1000,400])
grid on

%%
figure;
plot(time_data,T_in_he_es,'-o',time_data,T_in_interp,'-o',time_data,T_out_hp_interp,'-o',time_data,T_in_hp_interp,'-o')
legend('T in he','T in buffer','T out hp')
xlabel('Time (0=07:00, 12=19:00)')
ylabel('Temp (°C)')
set(gcf,'position',[500,500,1000,400])
grid on

%%

figure;