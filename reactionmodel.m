function [xx,tt] = reactionmodel(stepN)
    % initial condition
    n = ones(25,1) * 10000;
    p = zeros(25,1);
    tox = 250000;
    % growth rate
    mu_n = 0.1;
    mu_p = 0;
    % death rate
    d_n = 0.01;
    d_p = 0.01;
    % toxin to death rate
    tox_n = 1;
    tox_p = 0;
    % capacity
    K = 1000000;
    %initial values 
    t = 0;  %time 
    %storage 
    xx = zeros(51,stepN); % normal and persister
    tt = zeros(1,stepN); 
    l = zeros(8, 25);
    % simulate for stepN
    counter = 1;
    while counter <= stepN 
        %update propensities 
        tot = sum(n) + sum(p);
        if (tot == 0)
            tt = tt(:, 1:counter-1);
            xx = xx(:, 1:counter-1);
            break;
        end
        for type = 1:25
            l(1, type) = mu_n*(K-tot)/K*n(type);        % growth N
            l(2, type) = mu_p*(K-tot)/K*p(type);        % growth P
            l(3, type) = 10^-fix((type-1)/5)*n(type);   % N --> P
            l(4, type) = 10^-mod((type-1),5)*p(type);   % P --> N
            l(5, type) = d_n*n(type);                   % death N
            l(6, type) = d_p*p(type);                   % death P
            l(7, type) = tox_n * tox * n(type)/tot; % toxin N
            l(8, type) = tox_p * tox * p(type)/tot; % toxin P
        end
        l_total = sum(sum(l));
        % now throw the dart
        r1 = rand;
        tau = (1/l_total) * log(1/r1); 
        t = t + tau;
        r2 = rand; 
        comparison = l_total * r2;
        sum_ls = 0; 
        q = 0; 
        for i = 1:200
            sum_ls = sum_ls + l(i);
            if sum_ls > comparison  
                q = i;
                break;  
            end
        end
        et = ceil(q/8); % effective type
        switch mod(q, 8)
            case 1  % growth N
                n(et) = n(et) + 1;
            case 2  % growth P
                p(et) = p(et) + 1;
            case 3  % N --> P
                n(et) = n(et) - 1;
                p(et) = p(et) + 1;
            case 4  % P --> N
                p(et) = p(et) - 1;
                n(et) = n(et) + 1;
            case 5  % death N
                n(et) = n(et) - 1;
            case 6  % death P
                p(et) = p(et) - 1;
            case 7  % toxin N
                n(et) = n(et) - 1;
                tox = tox - 1;
            case 0  % toxin P
                p(et) = p(et) - 1;
                tox = tox - 1;
        end
        % store the updated states and time
        x = [n ; p; tox];
        tt(1, counter) = t; 
        xx(:, counter) = x; 
        counter = counter + 1;
    end
end