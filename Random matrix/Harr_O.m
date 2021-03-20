function [ Q ] = Harr_O( n )
% This function generates samples from the Harr measure of the cpt Lie
% group O(n) by applying qr decomposition to random Gaussian matrix, this
% is known as the Circular Orthogonal Ensemble (COE)


    A       = randn(n,n);
    [Q, R]  = qr(A);
    
    
end