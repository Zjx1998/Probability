% num     = 20;
% n_arr   = linspace(100, 2000, num);
% error   = zeros(num, 1);
% 
% 
% tic;
% for i = 1:20
%     n           = n_arr(i);
%     P           = rand(n, n);
%     P           = P ./ sum(P, 2);
%     [V, e]      = eig(P');
%     pi          = V(:, 1) / sum(V(:, 1));
%     pi_         = ones(n, 1) / n;
%     error(i)    = sum(abs(pi - pi_));
% end
% toc;
% 
% 
% Q       = pi .* P;
% epr     = sum(sum(Q .* log(Q ./ Q')));


