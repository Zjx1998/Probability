%--------------------------------------------------------------------------
% Expectation of statistics of RMT and relation with dimension n 
% i.e. E[log\kappa_n] = log n + c, more details could be found in 
% EIGENVALUES AND CONDITION NUMBERS OF RANDOM MATRICES
%--------------------------------------------------------------------------


n_num           = 8;                                                               
n_arr           = 2 .^ linspace(1, n_num, n_num);                          % dimension of random matrices
%n_max           = n_arr(n_num);
num             = 10000;                                                    % number of sample
%eigen_H         = zeros(n_max, num);                                       % eigenvalue for Harr ensemble
%eigen_G         = zeros(n_max, num);                                       % eigenvalue for GOE
max_sing_G      = zeros(num, 1);                                           % max singular value for GOE
min_sing_G      = zeros(num, 1);                                           % min singular value for GOE
cond_G          = zeros(num, 1);                                           % conditional number for GOE
E_logcond_G     = zeros(n_num, 1);                                         % expectation of log conditional number
E_logmaxsing_G  = zeros(n_num, 1);                                         % expectation of log max singular value
E_logminsing_G  = zeros(n_num, 1);                                         % expectation of log min singular value
gen             = 1;                                                       % indicator for generator
sym_c           = 2;                                                       % choice of symmetrization


tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP: iterate RM with different dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:8
    
    
    n = n_arr(j);
    for i = 1:num


        % this check the semicircle law of GOE
        switch gen
            case 1
                % entries given by Gaussian
                A_1             = randn(n, n);
            case 2
                % entries given by exponential
                A_1             = -log(rand(n, n));
            case 3
                % entries given by Cauchy
                A_1             = tan(pi * (rand(n, n) - 0.5 ));
            otherwise
                print('Unknown generator');
        end


        % symmetrization
        switch sym_c
            case 1
                % symmetrize using (A + A^T)/2, with diagonal elements ~ N(0,\sqrt(2)), other ~ N(0,1)
                A               = (A_1 + A_1')/ sqrt(2 * n);
            case 2
                % symmetrize using AA^T
                A               = A_1 * A_1';
            otherwise
                print('Unknown symmetrization method');
        end
        [V, D]          = eig(A);
        %eigen_G(:, i)   = diag(D);
        max_sing_G(i)   = max(diag(D));
        min_sing_G(i)   = min(diag(D));
        cond_G(i)       = cond(A);


        % this check the semicircle law of correspondent incorrect sampling function of
        % Harr measure
    %     A               = Harr_O(n);
    %     eigen_H(:, i)   = eig(A);


    end
    E_logcond_G(j)      = mean(log(cond_G));                                   
    E_logmaxsing_G(j)   = mean(log(max_sing_G));                                    
    E_logminsing_G(j)   = mean(log(min_sing_G));    
    
    
end
toc;


hold on;
grid on;
plot(log(n_arr), E_logcond_G,'LineWidth',2);
plot(log(n_arr), E_logmaxsing_G,'LineWidth',2);
plot(log(n_arr), E_logminsing_G,'LineWidth',2);
legend('expectation of log conditional number',...
       'expectation of log max singular value',...
       'expectation of log min singular value');
set(gca,'fontsize',20,'fontname','Times');
xlabel('log n');
ylabel('log E');