%% This is the script to simulate the deterministic model from 
%  Edo Kussell, etal. 2005
clc; clear; close all
% Set parameters
global a b mu_n mu_p

%% Simulate wt condition
n = 50000;
p = 0;
vars = [n, p];
options = odeset('RelTol',1e-8,'AbsTol',1e-10, 'Events', @capacityFcn);
figure(1);
t_end = 400;
t_g = 20;
t_s = 2.5;
t_pos = [];
case_pos = [];
t_last = 0;
a = 1.2E-6;
b = 0.1;
while t_last < t_end
    % growth
    t_span = [t_last, t_last + t_g];    
    mu_n = 2;
    mu_p = 0;
    [t_temp, case_temp] = ode45(@dvardt, t_span, vars, options);
    t_pos = [t_pos; t_temp];
    case_pos = [case_pos; case_temp];
    t_last = t_temp(end);
    while t_last < t_span(end)
        t_span2 = [t_last, t_span(end)];
        vars = case_temp(end, :)*10/11;
        [t_temp, case_temp] = ode45(@dvardt, t_span2, vars, options);
        t_pos = [t_pos; t_temp];
        case_pos = [case_pos; case_temp];
        t_last = t_temp(end);
    end
    vars = case_temp(end, :);
    % antibiotic
    t_span = [t_last, t_last + t_s];
    mu_n = -4;
    mu_p = -0.4;
    [t_temp, case_temp] = ode45(@dvardt, t_span, vars, options);
    t_pos = [t_pos; t_temp];
    case_pos = [case_pos; case_temp];
    t_last = t_temp(end);
    while t_last < t_span(end)
        t_span2 = [t_last, t_span(end)];
        vars = case_temp(end, :)*10/11;
        [t_temp, case_temp] = ode45(@dvardt, t_span2, vars, options);
        t_pos = [t_pos; t_temp];
        case_pos = [case_pos; case_temp];
        t_last = t_temp(end);
    end
    vars = case_temp(end, :);
end

%% Simulate hipQ condition
n = 50000;
p = 0;
vars = [n, p];
options = odeset('RelTol',1e-8,'AbsTol',1e-10, 'Events', @capacityFcn);
t_end = 400;
t_g = 20;
t_s = 2.5;
t_pos2 = [];
case_pos2 = [];
t_last = 0;
a = 0.001;
b = 1E-4;
while t_last < t_end
    % growth
    t_span = [t_last, t_last + t_g];    
    mu_n = 2;
    mu_p = 0.2;
    [t_temp, case_temp] = ode45(@dvardt, t_span, vars, options);
    t_pos2 = [t_pos2; t_temp];
    case_pos2 = [case_pos2; case_temp];
    t_last = t_temp(end);
    while t_last < t_span(end)
        t_span2 = [t_last, t_span(end)];
        vars = case_temp(end, :)*10/11;
        [t_temp, case_temp] = ode45(@dvardt, t_span2, vars, options);
        t_pos2 = [t_pos2; t_temp];
        case_pos2 = [case_pos2; case_temp];
        t_last = t_temp(end);
    end
    vars = case_temp(end, :);
    % antibiotic
    t_span = [t_last, t_last + t_s];
    mu_n = -4;
    mu_p = -0.4;
    [t_temp, case_temp] = ode45(@dvardt, t_span, vars, options);
    t_pos2 = [t_pos2; t_temp];
    case_pos2 = [case_pos2; case_temp];
    t_last = t_temp(end);
    while t_last < t_span(end)
        t_span2 = [t_last, t_span(end)];
        vars = case_temp(end, :)*10/11;
        [t_temp, case_temp] = ode45(@dvardt, t_span2, vars, options);
        t_pos2 = [t_pos2; t_temp];
        case_pos2 = [case_pos2; case_temp];
        t_last = t_temp(end);
    end
    vars = case_temp(end, :);
end
plot(t_pos,case_pos);
hold on;
plot(t_pos2,case_pos2);
xlabel('time/hours')
ylabel('population size')
legend('wt','wt persisters','hipQ','hipQ persisters','Location','best');
xlim([0 200])
ylim([1, 1E6])
ax = gca;
ax.LineWidth = 2;
set(gca,'yscale','log')


function derivative = dvardt(t, var)
    global a b mu_n mu_p
    n = var(1);
    p = var(2);
    dndt = mu_n * n - a * n + b * p;
    dpdt = mu_p * p - b * p + a * n;
    derivative = [dndt; dpdt];
end

function [position, isternimal, direction] = capacityFcn(t, var)
    Nmax = 110000;
    position = Nmax - var(1) - var(2);
    isternimal = 1;
    direction = -1;
end