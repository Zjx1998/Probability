function [  ] = Laplacian2D( n )
% Laplacian in 1D with different boundary condition
% n is the dimension of Laplacian

    % There are two different kinds of scaling, the first one is more
    % common in the probability and graph-theoretic setting while the
    % second one is more common in scientific computing, i.e. d equals to
    % 1/n, leaving the domain invariant under the change of n.
    d               = 2;
    %d               = 1/n;
    n2              = n * n;

    T               = -2 * eye(n);
    for i           = 2 : n - 1
        T(i, i - 1) = 1;
        T(i, i + 1) = 1;
    end
    T               = T / (d * d);

    % Laplacian restricted to the function class which vanishes outside the domain 
    L1_Res          = T;
    % setting boundary condition
    L1_Res(1, 2)    = -1 / (d * d);
    L1_Res(n, n-1)  = -1 / (d * d);
    L1_Res          = sparse(L1_Res);
    
    I               = eye(n);
    L2_Res          = kron(I, L1_Res) + kron(L1_Res, I);
    
    
    % Dirichlet boundary condition Laplacian
    L1_Diri         = T;
    % setting boundary condition
    L1_Diri(1, 1)   = 1;
    L1_Diri(n, n)   = 1;
    L1_Diri         = sparse(L1_Diri);
    
    I               = eye(n);
    L2_Diri         = kron(I, L1_Diri) + kron(L1_Diri, I);
    % setting boundary condition
    L2_Diri(1:n, 1:n)               = eye(n);
    L2_Diri(1:n, 1:n)               = eye(n);
    for i = 1:n - 2
        L2_Diri(i*n + 1, :)         = zeros(1, n2);
        L2_Diri(i*n+1,i*n+1)        = 1;
        L2_Diri(i*n+n-1, :)         = zeros(1, n2);
        L2_Diri(i*n+n-1,i*n+n-1)    = 1;
    end
    
    format short e
    fprintf('n = %d\n', n);
    fprintf('condition number for restricted Laplacian: %g\n', ...
            condest(L2_Res));
    fprintf('condition number for Laplacian in 2D with Dirichlet BC: %g\n', ...
            condest(L2_Diri));
  
    
    
    
    
      % solving Laplacian
%     % parallelized version for computing laplacian matrix
%     laplacian = - (kron(I,T) + kron(T,I));
%     [U, S, V] = svds(laplacian, 10);
%     
%     % An interesting thing is that if you truly do the following things
%     % which make the laplacian the same as its discretized version, it will
%     % be singular thus we cannot just do qr decomposition to it.
%         %for i = 1:n*n
%         %    laplacian(i,i) = - sum(laplacian(i,:)) + laplacian(i,i);   
%         %end
%     
%     [L, R] = qr(laplacian/(d^2));    
%     u = R^(-1) * L^(-1) * f;
%     
%     %u_m = reshape(u,n,n);
%     %mesh(u_m);

end