%--------------------------------------------------------------------------
% Distribution of the invariant distribution under different ensembles of
% transition matrices
%--------------------------------------------------------------------------


global P


% ensemble of transition kernel of the Markov chain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE #1 detail balance
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
num                 = 1000000;
n                   = 3;
epr_arr             = zeros(num, 1);
inv_dist            = zeros(n, num);
eps_                = 1e-3;


for i = 1:num
    
    
    P               = rand(n, n);
    P               = P ./ sum(P, 2);
    [V, e]          = eig(P');
    for l = 1:n
        if (abs(e(l, l) - 1) < eps_)
            break
        end
    end
    if (abs(sum(V(:, l))) < eps_)
        i = i - 1;
        continue;
    end
    pi_             = V(:, l) / sum(V(:, l));
    Q_              = pi_ .* P;
    epr_arr(i)      = sum(sum(Q_ .* log(Q_ ./ Q_')));
    inv_dist(:, i)  = pi_;

    
end


% histogram for epr
% [y, x]              = hist(epr_arr,500);
% y                   = y / length(epr_arr) / mean(diff(x));
% bar(x,y,1);


% histogram for invariant distribution
inv_dist            = reshape(inv_dist, n * num, 1);
[y, x]              = hist(inv_dist,500);
y                   = y / length(inv_dist) / mean(diff(x));
bar(x,y,1);