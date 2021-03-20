%--------------------------------------------------------------------------
% Calculate the derivative of the epr under random noise, i.e.
%               d/dt|_{t = 0} epr(P + t\xi).
%--------------------------------------------------------------------------
% Ground truth matrix
P       = [ 0.1     0.4     0.12;
            0.1     0.1     0.5;
            0.8     0.5     0.38];
pi      = [ 0.2;    0.3;    0.5];
eps_    = 1e-4;
      
        
num_t   = 10;
t       = linspace(1e-4, 1e-2, num_t);
num     = 1000000;
Epr_    = zeros(num, 1); 
epr_    = zeros(num_t, 1);


for j = 1:num_t
    for i = 1:num


        % random noise
        xi      = rand(3) - 0.5;
        xi      = xi - mean(xi);
        P_      = (P + xi * t(j))';     
        [V, e]  = eig(P_');
        for l   = 1:3
            if (abs(e(l, l) - 1) < eps_)
                break
            end
        end
        if (abs(sum(V(:, l))) < eps_)
            i   = i - 1;
            continue
        end
        pi_     = V(:, l) / sum(V(:, l));
        Q_      = pi_ .* P_;
        Epr_(i) = sum(sum(Q_ .* log(Q_ ./ Q_')));


    end
    epr_(j)     = mean(Epr_);
    
    
end


Q       = pi .* P';
epr     = sum(sum(Q .* log(Q ./ Q')));
disp(epr);
plot(t, epr_);
