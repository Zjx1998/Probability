% Ce script simule le champ de 4D membrane. La taille du champ est n +
% 2 avec la condition de limite du Dirichlet type. 

clc;
step            = 1e-4;
n               = 15;
d               = 2;
iter_num        = 1e5;
%d               = 1/n;
n2              = n * n;
n4              = n2 * n2;
choice          = true;


T               = -2 * eye(n);
for i           = 2 : n - 1
	T(i, i - 1) = 1;
	T(i, i + 1) = 1;
end
T               = T / (d * d);

% Laplacian restricted to the function class which vanishes outside the domain 
L1_Res          = T;
% setting boundary condition
L1_Res(1, 2)    = 1 / (d * d);
L1_Res(n-1, n)  = 1 / (d * d);
L1_Res          = sparse(L1_Res);

I               = eye(n);
L2_Res          = kron(I, L1_Res) + kron(L1_Res, I);
I               = eye(n2);
L4_Res          = kron(I, L2_Res) + kron(L2_Res, I);
disp(size(L4_Res))

tic;
if choice
    % spectral method using diagonalization
%     [V,D]           = eig(L4_Res);
%     %disp(diag(V' * V))
%     d               = sqrt(-diag(D));
%     x               = randn(n4, 1);
%     %disp(size(x));
%     %disp(size(d));
%     %disp(reshape(d, n, n));
%     h               = V * x .* d;
    
    % spectral method using Cholesky factorization
    R               = chol(-L4_Res);
    x               = randn(n4, 1);
    h               = R \ x;
else
    % Convergence phase
    h               = zeros(n2, 1);
    for i = 1:iter_num
        h           = h + L2_Res * h * step + randn(n2, 1) * sqrt(step);
    end
end    
toc;


% plot of the whole field
% h_ext       = zeros(n2 + 2, n2 + 2);
% h_ext(2:n2 + 1, 2: n2 + 1) = reshape(real(h), n2, n2);

% plot of a slice of the field
h_ext       = zeros(n + 2, n + 2);
h_ext(2:n + 1, 2: n + 1) = reshape(real(h()), n, n);


s           = surf(h_ext, 'FaceAlpha',0.9);
set(gcf,'outerposition',get(0,'screensize'));
colorbar
colormap hot
s.EdgeColor = 'none';