function [xx,tt] = hipQ(stepN)
    % initial condition
    n = 50000; 
    p = 0; 
    % parameters
    Nmax = 110000;
    Ni = 100000;
    a = 0.001;
    b = 1E-4;
    t_g = 20;
    t_s = 2.5;
    mu_n = 2;
    mu_p = 0.2;
    %initial values 
    t = 0;  %time 
    %storage 
    xx = zeros(2,stepN); % normal and persister
    tt = zeros(1,stepN); 
    l = zeros(1, 4);
    % simulate for stepN
    counter = 1;
    while counter <= stepN 
        if mod(t,(t_g + t_s)) < t_g
            mu_n = 2;
            mu_p = 0.2;
        else
            mu_n = -4;
            mu_p = -0.4;
        end
        %update propensities 
        l(1) = abs(mu_n) * n;  % growth of normal cells
        l(2) = abs(mu_p) * p;  % growth of persister cells
        l(3) = a * n;  % switch from normal to persister
        l(4) = b * p; % switch from persister to normal
        l_total = sum(sum(l)); 
        % now throw the dart
        r1 = rand;
        tau = (1/l_total) * log(1/r1); 
        t = t + tau;
        r2 = rand; 
        comparison = l_total * r2; 
        sum_ls = 0; 
        q = 0; 
        for i = 1:4
            sum_ls = sum_ls + l(i);
            if sum_ls > comparison  
                q = i; break;  
            end
        end
        if q == 1 % growth of normal cells
            n = n + sign(mu_n); 
        elseif q == 2 % growth of persister cells
            p = p + sign(mu_p); 
        elseif q == 3 % switch from normal to persister 
            n = n - 1;
            p = p + 1;
        elseif q == 4 % switch from persister to normal
            p = p - 1;
            n = n + 1;
        end
        Ntot = p + n;
        if Ntot > Nmax
            mean_n = n*Ni/Ntot;
            n = poissrnd(mean_n);
            mean_p = p*Ni/Ntot;
            p = poissrnd(mean_p);
        end
        % store the updated states and time
        x = [n; p];
        tt(1, counter) = t; 
        xx(:, counter) = x; 
        counter = counter + 1;
    end;
end