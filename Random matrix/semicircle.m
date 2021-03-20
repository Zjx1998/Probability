%--------------------------------------------------------------------------
% Semicircle law for the GOE and beyond:
% 1. Semicircle law for GOE
% 2. Min singular value law
% 3. Conditional number law
%-------------------------------------------------------------------------- 


n       = 100;                                                             % dimension of random matrices
num     = 10000;                                                          % number of sample
eigen_H = zeros(n, num);                                                   % eigenvalue for Harr ensemble
eigen_G = zeros(n, num);                                                   % eigenvalue for GOE
sing_G  = zeros(num, 1);                                                   % min singular value for GOE
cond_G  = zeros(num, 1);                                                   % conditional number for GOE
gen     = 1;                                                               % indicator for generator
sym_c   = 2;                                                               % choice of symmetrization


tic;
for i = 1:num
    
    
    % this check the semicircle law of GOE
    switch gen
        case 1
            % entries given by Gaussian
            A_1             = randn(n, n);
        case 2
            % entries given by exponential
            A_1             = -log(rand(n, n));
        case 3
            % entries given by Cauchy
            A_1             = tan(pi * (rand(n, n) - 0.5 ));
        otherwise
            print('Unknown generator');
    end
    
    
    % symmetrization
    switch sym_c
        case 1
            % symmetrize using (A + A^T)/2, with diagonal elements ~ N(0,\sqrt(2)), other ~ N(0,1)
            A               = (A_1 + A_1')/ sqrt(2 * n);
        case 2
            % symmetrize using AA^T
            A               = A_1 * A_1';
        otherwise
            print('Unknown symmetrization method');
    end
    [V, D]          = eig(A);
    eigen_G(:, i)   = diag(D);
    sing_G(i)    = max(diag(D));
    cond_G(i)    = cond(A);
    
  
    % this check the semicircle law of correspondent incorrect sampling function of
    % Harr measure
%     A               = Harr_O(n);
%     eigen_H(:, i)   = eig(A);
    
    
end
toc;


%eigen_H = reshape(eigen_H, num * n, 1);
eigen_G = reshape(eigen_G, num * n, 1);
eigen_G = sort(eigen_G);
cond_G  = sort(cond_G);
subplot(3,1,1)
hist(eigen_G, 500);
subplot(3,1,2)
hist(sing_G, 500);
subplot(3,1,3)
hist(cond_G(1:60000), 500);