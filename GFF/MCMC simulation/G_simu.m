n               = 5;
step            = 3e-2;
iter_num        = 1e7;
sample_num      = 1e4;
%d               = 1/n;
n2              = n * n;
h               = rands(n, 1);

for i = 1:iter_num
    h           = h - h * step + randn(n, 1) * sqrt(2 * step);
end

%Sampling phase
sample_0        = zeros(sample_num, 2);
x               = 2;
y               = 4;
for i = 1:sample_num
    h           = h - h * step + randn(n, 1) * sqrt(2 * step);
    sample_0(i, 1) = h(x);
    sample_0(i, 2) = h(y);
    %subplot(2,2,i);
    %h_ext       = zeros(n + 2, n + 2);
    %h_ext(2:n + 1, 2: n + 1) = reshape(h, n, n);
    %s           = surf(h_ext, 'FaceAlpha',0.9);
    %s.EdgeColor = 'none';
end

disp(cov(sample_0))
disp(mean(sample_0(:, 1)))