function [ T_cov ] = Coupon_Cov( n, iter )
% Ce script calcule une sequence de temp couverture en une torus annulus 
% T_N par MCMC

    T_cov       = zeros(iter, 1);
    tic;
    for i = 1:iter
        n_V     = n;
        V       = zeros(n, 1);
        t       = 0;
        rnd     = 1;
        while n_V ~= 0
            t   = t + 1;
            if V(rnd) == 0
                n_V     = n_V - 1;
            end
            V(rnd)      = V(rnd) + 1;
            rnd = ceil(rand() * n);
        end
        T_cov(i)        = t;
        disp(t / n / log(n));
    end
    toc;


end