% Parameters for random walk
iter    = 1e5;
X_t     = 0;
Y_t     = 0;
n       = 80;
L_t     = ones(n, n) * -1;                                                     % Local time matrix
L_t(1,1)= 0;
n_V     = n * n - 1;                                                       % # of not-v node
t       = 0;


% Generating the Torus
aminor  = 1.;                                                              % Torus minor radius
Rmajor  = 3.;                                                              % Torus major radius

theta   = linspace(-pi, pi, n)   ;                                         % Poloidal angle
phi     = linspace(0., 2.*pi, n) ;                                         % Toroidal angle

[q, p]  = meshgrid(phi, theta);
x       = (Rmajor + aminor.*cos(p)) .* cos(q);
y       = (Rmajor + aminor.*cos(p)) .* sin(q);
z       = aminor.*sin(p);


set(gcf,'outerposition',get(0,'screensize'));
surf(x, y, z, L_t);
title([t, n_V]);
axis equal
colorbar
colormap hot
while n_V ~= 0
    %disp(t)
    t   = t + 1;
	rnd = rand();
	if rnd > 0.75
        X_t = mod(X_t + 1, n);
    elseif rand > 0.5
        X_t = mod(X_t - 1, n);
    elseif rand > 0.25
        Y_t = mod(Y_t + 1, n);
    else
        Y_t = mod(Y_t - 1, n);
    end
    if(L_t(X_t + 1, Y_t + 1) == -1)
        n_V = n_V - 1;
    end
    L_t(X_t + 1, Y_t + 1)   = L_t(X_t + 1, Y_t + 1) / 1.1;
    
    if t < 300 || mod(t, 10) == 0
        %surf(x, y, z, - 1 ./ (L_t+1));
        surf(x, y, z, L_t);
        title([t, n_V]);
        axis equal
        colorbar
        %colormap hot
    end
    pause(1e-10)
    
end 