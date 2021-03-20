%--------------------------------------------------------------------------
% Convergence rate of important quantities such as invariant distribution 
% for random transition matrix as the dimension n tends to infinity. This
% script can be viewed as a substep of 'RM_epr.m'
%--------------------------------------------------------------------------



eps_    = 1e-2;
num     = 10000;
num_b   = 5;
low_n   = 4;
high_n  = 40;
num_n   = high_n - low_n + 1;
inv_arr = zeros(num, 1);
epr_avr = zeros(num_b, num_n);
epr_var = zeros(num_b, num_n);
bias    = linspace(0, 1, num_b);
n       = linspace(low_n, high_n, num_n);

tic;
for k = 1:num_n
    for j = 1:num_b
        for i = 1:num
            
            
            P           = rand(n(k), n(k)) + bias(j);
            P           = P ./ sum(P, 2);
            [V, e]      = eig(P');
            for l = 1:n(k)
                if (abs(e(l, l) - 1) < eps_)
                    break
                end
            end
            if (abs(sum(V(:, l))) < eps_)
                i = i - 1;
                continue
            end
            pi          = V(:, l) / sum(V(:, l));
            inv_arr(i)  = pi(1);
            
            
        end
        inv_avr(j, k)   = mean(inv_arr);
        inv_var(j, k)   = var(inv_arr);
        
        
    end
end
toc;


subplot(1, 2, 1);
hold on; 
for i = 1:num_b
    plot(n,  (abs(inv_avr(i, :) - 1 ./ n) .* n), 'LineWidth', 2);
end
xlabel("Dimension of matrix n");
ylabel("1/(epr bias)");
set(gcf,'outerposition',get(0,'screensize'));
set(gca,'fontsize',40,'fontname','Times');
legend( 'epr mean, bias = 0', ...
        'epr mean, bias = 0.25', ...
        'epr mean, bias = 0.5', ...
        'epr mean, bias = 0.75', ...
        'epr mean, bias = 1');
 

subplot(1, 2, 2);
hold on;
for i = 1:num_b
    plot(n, n .* n .* (inv_var(i, :)), 'LineWidth', 2);
end
xlabel("Dimension of matrix n");
ylabel("1/(epr var)");

set(gcf,'outerposition',get(0,'screensize'));
set(gca,'fontsize',40,'fontname','Times');
legend( 'epr var, bias = 0', ...
        'epr var, bias = 0.25', ...
        'epr var, bias = 0.5', ...
        'epr var, bias = 0.75', ...
        'epr var, bias = 1');
%print(figure(1), "Rm_epr2.jpg" , '-djpeg', '-r500');