function [ ] = Laplacian1D( n )
% Laplacian in 1D with different boundary condition
% n is the dimension of Laplacian

    % There are two different kinds of scaling, the first one is more
    % common in the probability and graph-theoretic setting while the
    % second one is more common in scientific computing, i.e. d equals to
    % 1/n, leaving the domain invariant under the change of n.
    %d               = 1;
    d               = 1/n;

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
    
    % Dirichlet boundary condition Laplacian
    L1_Diri         = T;
    % setting boundary condition
    L1_Diri(1, 1)   = 1;
    L1_Diri(n, n)   = 1;
    L1_Diri         = sparse(L1_Diri);
    
    % Neuman boundary condition Laplacian, where the directional
    % derivatives are both outwards at two boundaries
    L1_Neu          = T;
    % setting boundary condition
    L1_Neu(1, 1)    = 1 / d;
    L1_Neu(1, 2)    = -1 / d;
    L1_Neu(n, n)    = 1 / d;
    L1_Neu(n, n-1)  = -1 / d;
    L1_Neu          = sparse(L1_Neu);
    
    % Mixed boundary condition Laplacian
    L1_Mix          = T;
    % setting boundary condition
    L1_Mix(1, 1)    = 1 / d;
    L1_Mix(1, 2)    = -1 / d;
    L1_Mix(n, n)    = 1;
    L1_Mix          = sparse(L1_Mix);
    
    format short e
    fprintf('n = %d\n', n);
    fprintf('condition number for restricted Laplacian: %g\n', ...
            condest(L1_Res));
    fprintf('condition number for Laplacian in 1D with Dirichlet BC: %g\n', ...
            condest(L1_Diri));
    fprintf('condition number for Laplacian in 1D with Dirichlet BC: %g\n', ...
            condest(L1_Neu));
    fprintf('condition number for Laplacian in 1D with Dirichlet BC: %g\n', ...
            condest(L1_Mix));
    
end