% This script reproduce the figure of "CT30-1.jpg, CT30-2.jpg, CT50-1.jpg, CT50-2.jpg"
% Choose one data to load in and then use a specific choice of
% normalization

% Deux normalisation
% n1T_cov = \sqrt{T_cov} / n - \log(n)
% n2T_cov = T_cov / n^2 - \log^2(n)


load('Tcovn30.mat')
%load('Tcovn50.mat')
choice = 2;

% two choices of normalization
if choice == 1
    nT_cov = sqrt(T_cov) / n - log(n);
else
    nT_cov = T_cov / n / n - power(log(n), 2);% * 4 / pi;
end
% This is for Coupon collector problem
%n2T_cov = T_cov / n - log(n);

mu_     = mean(nT_cov);
sigma   = var(nT_cov);

beta = pi / sqrt(sigma * 6);                    % inverse scale parameter
mu      = mu_ - vpa(eulergamma / beta);

d       = max(nT_cov) + 0.1;
x_      = linspace(-d, d, 1000);
y_      = beta * exp(-beta * (x_ - mu) - exp(-beta * (x_ - mu)));
%y_      = evpdf(x_, mu, 1/beta);



[y, x]  = hist(nT_cov,300);
y       = y / length(T_cov) / mean(diff(x));
bar(x,y,1);
ylabel('Density');
xlabel('Cover time');

if choice == 1
    title("Empirical distribution of the normalized cover time of first type with n = " + n);
    axis([-1.5 1.5 0 3]);
else
    title("Empirical distribution of the normalized cover time of second type with n = " + n);
    axis([-10 10 0 0.4]);
end

set(gca,'fontsize',40,'fontname','Times');
hold on
plot(x_, y_, 'LineWidth', 5);
legend('Empirical distribution', 'Gumbel fitting')