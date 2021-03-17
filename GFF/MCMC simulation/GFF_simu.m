% Ce script simule le champ libre de Gaussian. La taille du champ est n +
% 2 avec la condition de limite du Dirichlet type. 

%clc;
step            = 1e-2;
n               = 50;
d               = 2;
iter_num        = 1e5;
sample_num      = 1e4;
%d               = 1/n;
n2              = n * n;
% using Langevin dynamics or spectral method for sampling
choice          = true;
%choice          = false;

% initialization of the 2D Laplacian restricted to the function class which
% vanishes outside the domain
T               = -2 * eye(n);
for i           = 2 : n - 1
	T(i, i - 1) = 1;
	T(i, i + 1) = 1;
end
T               = T / (2 * d);
L1_Res          = T;
% setting boundary condition
L1_Res(1, 2)    = 1 / (2 * d);
L1_Res(n, n-1)  = 1 / (2 * d);
L1_Res          = sparse(L1_Res);
%disp(L1_Res);
    
I               = eye(n);
L2_Res          = kron(I, L1_Res) + kron(L1_Res, I);
%disp(L2_Res);
%G               = -inv(L2_Res);
%disp(L2_Res);

tic;
if choice
    % spectral method using diagonalization
%     [V,D]           = eig(L2_Res);
%     %disp(diag(V' * V))
%     d               = sqrt(-diag(D));
%     x               = randn(n2, 1);
%     %disp(size(x));
%     %disp(size(d));
%     %disp(reshape(d, n, n));
%     h               = V * x .* d;
    
    % spectral method using Cholesky factorization
    R               = chol(-L2_Res);
    x               = randn(n2, 1);
    h               = R \ x;
else
    % Convergence phase
    h               = zeros(n2, 1);
    for i = 1:iter_num
        h           = h + L2_Res * h * step + randn(n2, 1) * sqrt(step);
    end
end    
toc;

h_ext       = zeros(n + 2, n + 2);
h_ext(2:n + 1, 2: n + 1) = reshape(h, n, n);
s           = surf(h_ext, 'FaceAlpha',0.9);
set(gcf,'outerposition',get(0,'screensize'));
colorbar
colormap hot
s.EdgeColor = 'none';


% %Sampling phase, using Langevin dynamics to sampling
% sample_0        = zeros(sample_num, 2);
% x               = floor(n/2) * n + floor(n/2);                             % point 1
% y               = floor(n/3) * n + floor(n/3);                             % point 2
% for i = 1:sample_num
%     h           = L2_Res * h * step + randn(n2, 1) * sqrt(2 * step);
%     sample_0(i, 1) = h(x);
%     sample_0(i, 2) = h(y);
%     subplot(2,2,i);

    % plot the result
%     h_ext       = zeros(n + 2, n + 2);
%     h_ext(2:n + 1, 2: n + 1) = reshape(h, n, n);
%     s           = surf(h_ext, 'FaceAlpha',0.9);
%     set(gcf,'outerposition',get(0,'screensize'));
%     colorbar
%     colormap hot
%     s.EdgeColor = 'none';


% end
% 
% %disp(cov(sample_0))
% %disp(mean(sample_0(:, 1)))
% %disp([   Green_func(x, x, n) Green_func(x, y, n); ...
% %         Green_func(y, x, n) Green_func(y, y, n)]);
% %disp([G(x, x) G(x, y); G(y, x) G(y, y)]);