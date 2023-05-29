function [date,value]=read_csv(tab)

date=table2array(tab(:,1));
value=string(table2array(tab(:,2)));
value = double(replace(value,",","."));
