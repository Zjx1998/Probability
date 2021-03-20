%--------------------------------------------------------------------------
% Distribution of the spacing between eigenvalues of random matrix
%-------------------------------------------------------------------------- 


n       = 100;                                                             % dimension of random matrices
num     = 100000;                                                           % number of sample
space   = zeros(n - 1, num);                                               % spacing between eigenvaluea for GOE
space_1 = zeros(num, 1);                                                   % spacing between first two eigenvalueas for GOE


tic;
%profile on
for i = 1:num
    
    
    % this check the semicircle law of GOE
    A_1             = randn(n, n);
    A               = (A_1 + A_1') / sqrt(2 * n); 
    eig_A           = eig(A); 
    space(:, i)     = Sp(sort(eig_A));
    space_1(i, 1)   = Sp(sort(eig_A(1:2)));
    
    
end
%profile viewer
toc;


subplot(2, 1, 1);
space   = reshape(space, num * (n - 1), 1);
hist(space, 500);
subplot(2, 1, 2);
hist(space_1, 50);