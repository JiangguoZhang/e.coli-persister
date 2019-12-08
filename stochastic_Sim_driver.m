%% This is the script to generate the model from 
%  Edo Kussell, etal. 2005
%  It is adapted from the codes in HW6
%
clc; close; clear all;
steps = 40000000;
[xx1,tt1]=wildtype(steps);
[xx2,tt2]=hipQ(steps);
figure(1)
sgtitle('Bacterial Persistence')
%subplot(2,2,1)
%%
plot(tt1,xx1)
hold on;
plot(tt2,xx2)
xlabel('time/hours')
ylabel('population size')
legend('wt','wt persisters','hipQ','hipQ persisters','Location','best');
xlim([0 200])
ax = gca;
ax.LineWidth = 1.5;
set(gca,'yscale','log')