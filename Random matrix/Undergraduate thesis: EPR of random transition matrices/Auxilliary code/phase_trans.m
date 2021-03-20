%--------------------------------------------------------------------------
% Phase transition of the noisy model, trying to establish the phase
% transition of the effect of noise on estimation of epr.
% The fact is that there does not exist a phase transition phenomenon, 
% however for positive epr, small noise will always decrease the epr
%--------------------------------------------------------------------------


global P


% transition kernel of the Markov chain, the invariant measure is given by
% [0.2 0.3 0.5]^T, P_{ij} = P(X_{n + 1} = j|X_n = i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE #1-1 naive detail balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P       = ones(3)/3;
pi      = ones(3, 1)/3;
n       = size(P, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE #1-2 detail balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P       = [ 0.5     0.4     0.1;
%             0.3     0.4     0.3;
%             0.1     0.4     0.5];
% pi      = [ 0.3;    0.4;    0.3];
% n       = size(P, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE #2 not detail balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P       = [ 0.1     0.4     0.12;
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


n           = 3;
p_num       = 10;
n_num       = 20;
num         = 1e3;
eps_        = 1e-2;
prob        = linspace(0.01, P(1, 1) + P(1, 2) - 0.01, p_num);             % sequence of choice for transition matrix
noise_sz1   = linspace(0.01, 1, n_num);                                    % size of random noise
noise_sz2   = linspace(0.01, 0.1, n_num);                                    % size of random noise
epr_arr1    = zeros(num, 1);                                               % array of epr samples
epr_arr2    = zeros(num, 1);                                               % array of epr samples
epr_mean1   = zeros(n_num, p_num);                                             % array of mean of epr 
epr_mean2   = zeros(n_num, p_num);                                             % array of mean of epr
epr_var1    = zeros(n_num, p_num);                                             % array of var of epr
epr_var2    = zeros(n_num, p_num);                                             % array of var of epr


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:p_num


    % Change the model for transition matrix
    P(1, 1) = prob(i);
    P(1, 2) = 2/3 - prob(i);
    % Calculating the epr using well-known formula
    %Q       = pi .* P;
    %epr     = sum(sum(Q .* log(Q ./ Q')));
    %disp(epr);
    
    
    % Calculating the epr after a random noise
    %profile on
    for j = 1:n_num 


        for k = 1:num

            
            % Model of choice #1
            P_          = rand(n, n) * noise_sz1(j) + P;
            P_          = P_ ./ sum(P_, 2);
            % Model of choice #2
%             noise       = rand(n, n) * noise_sz1(j);
%             noise       = noise - 
%             P_          = P + noise;
%             xi          = rand(n, n) * noise_sz1(j);
%             xi          = xi - mean(xi);
%             P_          = xi' + P;
            [V, e]      = eig(P_');
            for l = 1:n
                if (abs(e(l, l) - 1) < eps_)
                    break
                end
            end
            if (abs(sum(V(:, l))) < eps_)
                k = k - 1;
                continue
            end
            pi_         = V(:, l) / sum(V(:, l));
            Q_          = pi_ .* P_;
            epr_arr1(k)  = sum(sum(Q_ .* log(Q_ ./ Q_')));
            
            
            % Model of choice #1
            P_          = rand(n, n) * noise_sz2(j) + P;
            P_          = P_ ./ sum(P_, 2);
            % Model of choice #2
%             noise       = rand(n, n) * noise_sz1(j);
%             noise       = noise - 
%             P_          = P + noise;
%             xi          = rand(n, n) * noise_sz1(j);
%             xi          = xi - mean(xi);
%             P_          = xi' + P;
            [V, e]      = eig(P_');
            for l = 1:n
                if (abs(e(l, l) - 1) < eps_)
                    break
                end
            end
            if (abs(sum(V(:, l))) < eps_)
                k = k - 1;
                continue
            end
            pi_         = V(:, l) / sum(V(:, l));
            Q_          = pi_ .* P_;
            epr_arr2(k) = sum(sum(Q_ .* log(Q_ ./ Q_')));


        end
        epr_mean1(j, i)    = mean(epr_arr1);
        epr_var1(j, i)     = var(epr_arr1);
        epr_mean2(j, i)    = mean(epr_arr2);
        epr_var2(j, i)     = var(epr_arr2);


    end
    %profile viewer

    
end


subplot(1, 2, 1)
hold on;
for i = 1:p_num
    plot(noise_sz1, epr_mean1(:, i), '-pg', 'LineWidth', 2, 'markersize',10);
    plot(noise_sz1, epr_var1(:, i), '-dr','LineWidth', 2, 'markersize',10);
end
%legend('epr mean','epr var');
xlabel("Magnitude of noise");
ylabel("Statistics");
set(gcf,'outerposition',get(0,'screensize'));
set(gca,'fontsize',40,'fontname','Times');


subplot(1, 2, 2)
hold on;
for i = 1:p_num
    plot(noise_sz2, epr_mean2(:, i), '-pg', 'LineWidth', 2, 'markersize',10);
    plot(noise_sz2, epr_var2(:, i), '-dr','LineWidth', 2, 'markersize',10);
end
%legend('epr mean','epr var');
xlabel("Magnitude of noise");
ylabel("Statistics");
set(gcf,'outerposition',get(0,'screensize'));
set(gca,'fontsize',40,'fontname','Times');
%print(figure(1), "phase_trans1.jpg" , '-djpeg', '-r500');