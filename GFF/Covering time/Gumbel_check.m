% Deux normalisations de temp couverture 
% n1T_cov = \sqrt{T_cov} / n - \log(n)
% n2T_cov = T_cov / n^2 - \log^2(n)


%load('covertime2.mat')
%n1T_cov = sqrt(T_cov) / n - log(n);
n2T_cov = T_cov / n / n - power(log(n), 2) * 4 / pi;
% This is for Coupon collector problem
%n2T_cov = T_cov / n - log(n);

nT_cov  = n2T_cov;
mu_     = mean(nT_cov);
sigma   = var(nT_cov);

beta = pi / sqrt(sigma * 6);                    % inverse scale parameter
mu      = mu_ - vpa(eulergamma / beta);

x_      = linspace(-10, 10, 1000);
y_      = beta * exp(-beta * (x_ - mu) - exp(-beta * (x_ - mu)));
%y_      = evpdf(x_, mu, 1/beta);


[y, x]  = hist(nT_cov,300);
y       = y / length(T_cov) / mean(diff(x));
bar(x,y,1);
ylabel('Density');
xlabel('Cover time');
title('Empirical distribution of normalized cover time for T_{50}');
set(gca,'fontsize',20,'fontname','Times');
set(gca,'fontsize',30,'fontname','Times');
hold on
plot(x_, y_, 'LineWidth',5);