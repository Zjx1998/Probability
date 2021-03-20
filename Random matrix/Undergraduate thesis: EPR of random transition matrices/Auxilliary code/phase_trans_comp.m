%--------------------------------------------------------------------------
% Comparison of the noise's effect on epr, whether the noise-epr curve is
% uniquely determined by initial value of epr, i.e. ground truth
%--------------------------------------------------------------------------


global P1 P2


% transition kernel of the Markov chain, the invariant measure is given by
% [0.2 0.3 0.5]^T, P_{ij} = P(X_{n + 1} = j|X_n = i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE #1-1 naive detail balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P1      = ones(3)/3;
pi      = ones(3, 1)/3;
n       = size(P1, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE #1-2 detail balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P2      = [ 0.5     0.4     0.1;
            0.3     0.4     0.3;
            0.1     0.4     0.5];
pi      = [ 0.3;    0.4;    0.3];
n       = size(P2, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE #2 not detail balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P2       = [ 0.1     0.4     0.12;
%             0.1     0.1     0.5;
%             0.8     0.5     0.38]';
% pi      = [ 0.2;    0.3;    0.5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE #3 not detail balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P       = [ 0.6     0.3     0.1;
%             0.75    0.05    0.2;
%             0.45    0.05    0.5];
% pi      = [ 0.6;    0.2;    0.2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE #4 high dimension model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n       = 100;
% P       = rand(n, n);
% P       = P ./ sum(P, 2);
% [V, e]  = eig(P');
% pi      = V(:, 1) / sum(V(:, 1));


Q           = pi .* P2;
epr     	= sum(sum(Q .* log(Q ./ Q')));
eps_        = 1e-2;
num         = 1e4;                                                      % number of samples
n_num       = 30;                                                          % number of noise size
n           = size(P2, 1);
noise_sz1   = linspace(0.01, 0.2, n_num);                                   % size of random noise
noise_sz2   = linspace(0.02, 1, n_num);                                    % size of random noise
epr_arr1    = zeros(num, 1);                                               % array of epr samples
epr_arr2    = zeros(num, 1);                                               % array of epr samples
epr_mean1   = zeros(n_num, 1);                                             % array of mean of epr 
epr_mean2   = zeros(n_num, 1);                                             % array of mean of epr
epr_var1    = zeros(n_num, 1);                                             % array of var of epr
epr_var2    = zeros(n_num, 1);                                             % array of var of epr


% Calculating the epr after a random noise
%profile on
for j = 1:n_num 

    
    for i = 1:num


        xi          = rand(n, n) * noise_sz1(j);
        xi          = xi - mean(xi);
        P_          = xi' + P1;
        [V, e]      = eig(P_');
        for k = 1:n
            if (abs(e(k, k) - 1) < eps_)
                break
            end
        end
        if (abs(sum(V(:, k))) < eps_)
            i = i - 1;
            continue
        end
        pi_         = V(:, k) / sum(V(:, k));
        Q_          = pi_ .* P_;
        epr_arr1(i) = sum(sum(Q_ .* log(Q_ ./ Q_')));
        
        
        xi          = rand(n, n) * noise_sz1(j);
        xi          = xi - mean(xi);
        P_          = xi' + P2;
        [V, e]      = eig(P_');
        for k = 1:n
            if (abs(e(k, k) - 1) < eps_)
                break
            end
        end
        if (abs(sum(V(:, k))) < eps_)
            i = i - 1;
            continue
        end
        pi_         = V(:, k) / sum(V(:, k));
        Q_          = pi_ .* P_;
        epr_arr2(i) = sum(sum(Q_ .* log(Q_ ./ Q_')));


    end
    epr_mean1(j)    = mean(epr_arr1);
    epr_var1(j)     = var(epr_arr1);
    epr_mean2(j)    = mean(epr_arr2);
    epr_var2(j)     = var(epr_arr2);
    
    
end
%profile viewer


hold on;
plot(noise_sz1, epr_mean1, '-pg', 'LineWidth', 2, 'markersize',10);
plot(noise_sz1, epr_mean2, '-dr', 'LineWidth', 2, 'markersize',10);
legend('epr mean','epr var');
xlabel("Magnitude of noise");
ylabel("Statistics");
set(gcf,'outerposition',get(0,'screensize'));
set(gca,'fontsize',40,'fontname','Times');