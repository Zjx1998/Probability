%--------------------------------------------------------------------------
% Variance and other statistics of entropy production rate under random
% noise of transition matrix
%--------------------------------------------------------------------------


global P1 P2


% transition kernel of the Markov chain, the invariant measure is given by
% [0.2 0.3 0.5]^T, P_{ij} = P(X_{n + 1} = j|X_n = i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE #1 detail balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P1       = [ 0.5     0.4     0.1;
%             0.3     0.4     0.3;
%             0.1     0.4     0.5];
% pi      = [ 0.3;    0.4;    0.3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE #2 not detail balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P1       = [ 0.1     0.4     0.12;
            0.1     0.1     0.5;
            0.8     0.5     0.38]';
pi      = [ 0.2;    0.3;    0.5];
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


Q           = pi .* P1;
epr     	= sum(sum(Q .* log(Q ./ Q')));
eps_        = 1e-2;
num         = 100000;                                                      % number of samples
n_num       = 30;                                                          % number of noise size
n           = size(P1, 1);
noise_sz1   = linspace(0.01, 0.15, n_num);                                    % size of random noise
noise_sz2   = linspace(0.02, 10, n_num);                                    % size of random noise
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
        
        
%         P_          = rand(n, n) * noise_sz2(j) + P1;
%         P_          = P_ ./ sum(P_, 2);
%         [V, e]      = eig(P_');
%         for k = 1:n
%             if (abs(e(k, k) - 1) < eps_)
%                 break
%             end
%         end
%         if (abs(sum(V(:, k))) < eps_)
%             i = i - 1;
%             continue
%         end
%         pi_         = V(:, k) / sum(V(:, k));
%         Q_          = pi_ .* P_;
%         epr_arr2(i) = sum(sum(Q_ .* log(Q_ ./ Q_')));
        
        
%         P_          = rand(n, n) * noise_sz2(j) + P;
%         P_          = P_ ./ sum(P_, 2);
%         [V, e]      = eig(P_');
%         for k = 1:n
%             if (abs(e(k, k) - 1) < eps)
%                 break
%             end
%         end
%         pi_         = V(:, k) / sum(V(:, k));
%         Q_          = pi_ .* P_;
%         epr_arr2(i) = sum(sum(Q_ .* log(Q_ ./ Q_')));


    end
    epr_mean1(j)    = mean(epr_arr1);
    epr_var1(j)     = var(epr_arr1);
    epr_mean2(j)    = mean(epr_arr2);
    epr_var2(j)     = var(epr_arr2);
    
    
end
%profile viewer


%subplot(1, 2, 1)
hold on;
plot(noise_sz1, epr_mean1, '-pg', 'LineWidth', 2, 'markersize',10);
plot(noise_sz1, epr_var1, '-dr','LineWidth', 2, 'markersize',10);
plot(noise_sz1, epr * ones(n_num, 1));
legend('epr mean','epr var','true epr');
xlabel("Magnitude of noise");
ylabel("Statistics");
set(gcf,'outerposition',get(0,'screensize'));
set(gca,'fontsize',40,'fontname','Times');

% subplot(1, 2, 2)
% hold on;
% plot(noise_sz2, epr_mean2, '-pg', 'LineWidth', 2, 'markersize',10);
% plot(noise_sz2, epr_var2, '-dr','LineWidth', 2, 'markersize',10);
% %plot(noise_sz2, epr * ones(n_num, 1));
% legend('epr mean','epr var','true epr');
% xlabel("Magnitude of noise");
% ylabel("Statistics");
% set(gcf,'outerposition',get(0,'screensize'));
% set(gca,'fontsize',40,'fontname','Times');
% axis([0.02 1 0 0.35]);
%print(figure(1), "Stat_epr.jpg" , '-djpeg', '-r500');



%opts        = [];
%output      = Epr_MC(opts);