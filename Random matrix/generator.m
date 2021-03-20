%--------------------------------------------------------------------------
% A collection of generating code for RMT
%--------------------------------------------------------------------------


% \gamma kurtosis, c.f. BEYOND UNIVERSALITY IN RANDOM MATRIX THEORY
% Normal, \gamma = 0
A           = randn(n);
A           = randn(n + v, n) + im * randn(n + v, n);
% Uniform, \gamma = - 1.2
A           = (rand-0.5)*sqrt(12);
A           = ((rand(n + v, n) - 0.5) + im * rand(n + v, n) - 0.5)) * sqrt(12);
% Bernoulli, \gamma = - 2
A           = sign(randn(n));
A           = sign(randn(n + v, n)) + im * sign(randn(n + v, n));
% Gamma, \gamma = 6
A           = rand(Gamma(), n) - 1;
A           = (rand(Gamma(),n + v, n) - 1) + im * (rand(Gamma(),n + v, n) - 1);