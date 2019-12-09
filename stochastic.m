clear; clc; close all;
%% This is the script to generate the model from 
%  Edo Kussell, etal. 2005
%
steps = 10000000;
[xx1,tt1]=reactionmodel(steps);
figure(1)
%sgtitle('Bacterial Persistence')
color = [[repelem(linspace(1,0,5), 5),0];zeros(1,26);[repmat(linspace(1,0,5),1,5),0]];
for i = 1:25
    plot(tt1,xx1(i,:) + xx1(i+25,:), 'color', color(:,i));
    hold on;
end
plot(tt1,xx1(51,:), 'color', 'green');
title("red=large N-->P, blue=large P-->N, green=toxin");
xlabel('time/hours')
ylabel('population size')