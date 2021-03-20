%--------------------------------------------------------------------------
% Convergence rate of entropy production rate for random transition matrix
% as the dimension n tends to infinity
%--------------------------------------------------------------------------



eps_    = 1e-2;
num     = 10000;
num_b   = 5;
low_n   = 4;
high_n  = 40;
num_n   = high_n - low_n + 1;
epr_arr = zeros(num, 1);
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
            Q           = pi .* P;
            epr_arr(i)  = sum(sum(Q .* log(Q ./ Q')));
            %epr_avr(j)  = mean(epr_arr);
            
            
        end
        epr_avr(j, k)   = mean(epr_arr);
        epr_var(j, k)   = var(epr_arr);
        
        
    end
end
toc;


subplot(1, 2, 1);
hold on; 
for i = 1:num_b
    plot(n,  1 ./ (abs(epr_avr(i, :) - 0.5 + f(bias(i)))), 'LineWidth', 2);
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
    plot(n, 1 ./ (epr_var(i, :)), 'LineWidth', 2);
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

    
function [ a ] = f(b)
    % function calculate the limit of entropy production rate under biased transition matrix 
        if b == 0
            a   = 0;
        else
            a   = b .* (b + 1) ./ (2 * b + 1) .* log((b + 1)./b);
        end
end