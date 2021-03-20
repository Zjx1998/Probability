%--------------------------------------------------------------------------
% Estimation of entropy production rate using thermodynamic uncertainty
% relation
%--------------------------------------------------------------------------


global P


% transition kernel of the Markov chain, the invariant measure is given by
% [0.2 0.3 0.5]^T, P_{ij} = P(X_{n + 1} = j|X_n = i)
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
P       = [ 0.1     0.4     0.12;
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



% Calculating the epr using well-known formula
n       = size(P, 1);
Q       = pi .* P;
epr     = sum(sum(Q .* log(Q ./ Q')));
disp(epr);


%profile on;
% Estimating epr using thermodynamic uncertainty relation, the relation
% between the estimation using different edges and different epochs,
% batchsizes
opts    = [];
e_num   = 10;
b_num   = 1;
epoch   = linspace(100000, 1000000, 10);
batch   = [100; 1000; 10000];
epr_arr = zeros(n, n, e_num, b_num);
t       = 1;

tic;
output  = repmat(struct('curr_mean',{0}, 'curr_var',{0}, 'est_epr',{0} ,'time',{0}), n ,n, e_num, b_num);
for i = 1:n
    for j = i + 1:n
        
        
            for e = 1:e_num
                for b = 1:b_num
                    opts.source         = i;
                    opts.target         = j;
                    opts.epoch          = epoch(e);
                    opts.batch          = batch(b + 1);
                    output(i, j, e, b)  = Epr_MC(opts);
                    epr_arr(i, j, e, b) = output(i, j, e, b).est_epr;
                end
            end
            
            
    end
end
toc;
%profile viewer;

%opts        = [];
%output      = Epr_MC(opts);