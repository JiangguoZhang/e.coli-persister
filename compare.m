%Simulate growth and death of different cells under various antibiotic
%conditions
%levels: faster, fast, medium, slow, slower (FfmsS)
clear, clc
close all
labels = {'Faster', 'Fast', 'Medium', 'Slow', 'Slower'};
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'NonNegative', 1:31);
%options = odeset('NonNegative', 1:51);
tspan = [0 10000];
y0 = zeros(51,1);
for i = 1:2:50
    y0(i) = 1;
    y0(i+1) = 0;
end
y0(51) = 1000;
[t,y] = ode45(@ddt, tspan, y0, options);
figure, plot(t,y(:,1:50))
title(sprintf('Population Growth for Toxin = %d (AU), Winner: ', y0(51)));
xlabel('Time')
ylabel('Population Size (Percent of Carrying Capacity)')
res = zeros(25, 1);
for i = 1:25
    res(i) = y(end, 2*i-1) + y(end, 2*i);
end
[~, I] = sort(res)
sort(res)
y(end,51)
title(sprintf('Population Growth for Toxin = %d (AU)', y0(51)));
ylim([0, max(res)*1.3])
%labels(I);
figure, plot(t,y(:,51))
function dt = ddt(t,y)
K = 100;
F = 0.1;
f = 0.01;
m = 0.001;
s = 0.0001;
S = 0;
ksmall = [F, f, m, s, S];
k = zeros(1, 50);
for i = 1:5
    for j = 1:5
        k(10*(i-1)+2*j-1) = ksmall(i);
        k(10*(i-1)+2*j) = ksmall(j);
    end
end
mu_Nmax = 0.1;
mu_Pmax = 0;
if t > 0
    tox = y(51);
else
    tox = 0;
end
tot = sum(y(1:50));
mu_N = mu_Nmax *(K - tot)/K;
mu_P = mu_Pmax *(K - tot)/K;
dt = zeros(51,1);
for i = 1:2:50
    dt(i) = mu_N*y(i) - k(i)*y(i) + k(i+1)*y(i+1) - tox * y(i) - 0.01 * y(i);
    dt(i+1) = mu_P*y(i+1) - k(i+1)*y(i+1) + k(i)*y(i) - 0.01 * y(i+1);
    dt(51) = dt(51) - tox*y(i);
end
dt(51) = dt(51) - tox*0.01;
end