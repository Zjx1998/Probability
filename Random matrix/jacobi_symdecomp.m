%--------------------------------------------------------------------------
% Jacobian calculation of the symmetric eigenvalue decomposition:
% A = QDQ', A symmetric, D diagonal, Q orthogonal
% dD * Q'dQ / dA
%-------------------------------------------------------------------------- 


n           = 30;                                                          % Size of the matrix
X           = randn(n);
A           = (X + X')/n;                                                  % Generate a symmetric matrix
[Q,L]       = eig(A);                                                      % Compute its eigenvalues/eigenvectors
L           = diag(L);
JacMatrix   = zeros(n * (n + 1)/2);                                        % Initialize Jacobian matrix
epsilon     = 1e-7; 
idx         = 1;
mask        = triu(ones(n),1); 
mask        = logical(mask(:));                                            % Upper triangular mask
factor      = 1;                                                           % Factor of varianble transformation


for i = 1:n
    for j = i:n
        
        
        %%% Perturbation Matrix
        % Initialize perturbation
        Eij                         = zeros(n); 
        Eij(i,j)                    = 1;
        Eij(j,i)                    = 1;
        Ap                          = A + epsilon * Eij; 
        
        
        %%% Eigenvalues and Eigenvectors
        [Qp,Lp]                     = eig(Ap);
        % Eigenvalue perturbation
        dL                          = (diag(Lp) - L)/epsilon; 
        % Eigenvector perturbation
        QdQ                         = Q' * (Qp-Q)/epsilon; 
        
        
        %%% The Jacobian Matrix
        % Eigenvalue part of Jacobian
        JacMatrix(1:n,idx)          = dL; 
        % Eigenvector part of Jacobian
        JacMatrix((n+1):end,idx)    = QdQ(mask); 
        % Increment column counter
        idx                         = idx + 1; 
        
        
        % Calculate the factor
        if i ~= j
            factor                      = factor * abs(L(i) - L(j));
        end
        
        
    end
end


format long
disp(factor * det(JacMatrix))