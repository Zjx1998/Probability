n   = 100;
num = 1e5;
A   = rand(n, n);
disp(eig(A));
for i = 1:num
    [Q, R]  = qr(A);
    A       = R * Q;
end
disp(eig(A));
disp(A);