% Compare alpha values from different sources

clear all
close all

t=-60:30;

A_w=1./(4.792-3.63e-2*t-1.897e-4*t.^2);
A_k=1./(1.041-2.606e-2*t+5.313e-4*t.^2);
A_x=1./(7.919e-2-2.226e-3*t+8.201e-5*t.^2);

alpha=1./(0.0003859.*t.^2-0.02718.*t+1.034);

figure('renderer','painters');
hold on

plot(t,A_w,'linewidth',1.5);
plot(t,A_k,'linewidth',1.5);
%plot(t,A_x,'linewidth',1.5);
plot(t,alpha,'linewidth',1.5);

legend('W patent','Ka patent','Ka paper');

xlabel('Temperature')
ylabel('alpha')