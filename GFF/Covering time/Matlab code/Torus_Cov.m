function [ T_cov ] = Torus_Cov( n, iter )
% Ce script calcule une sequence de temp couverture en une torus annulus 
% T_N par MCMC, 

% Something wrong with this prpgram since the covering number sampled by
% this algorithm is different from the python and C code

    % Parameters for random walk
    X_t     = 0;
    Y_t     = 0;
    T_cov   = zeros(iter, 1);


    % Generating the Torus
    aminor  = 1.;                           % Torus minor radius
    Rmajor  = 3.;                           % Torus major radius

    theta   = linspace(-pi, pi, n)   ;      % Poloidal angle
    phi     = linspace(0., 2.*pi, n) ;      % Toroidal angle

    [q, p]  = meshgrid(phi, theta);
    x       = (Rmajor + aminor.*cos(p)) .* cos(q);
    y       = (Rmajor + aminor.*cos(p)) .* sin(q);
    z       = aminor.*sin(p);

    tic;
    for i = 1:iter
        L_t     = zeros(n, n);              % Local time matrix
        n_V     = n * n;                    % # of not-v node
        t       = 0;
        while n_V ~= 0
            %disp(t)
            if(L_t(X_t + 1, Y_t + 1) == 0)
                n_V = n_V - 1;
            end
            L_t(X_t + 1, Y_t + 1)   = L_t(X_t + 1, Y_t + 1) + 1;
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

        %     if t < 300 || mod(t, 10) == 0
        %         surf(x, y, z, L_t)
        %         title([t, n_V])
        %         axis equal
        %     end
        %     pause(1e-10)

        end
        T_cov(i)    = t;
    end
    toc;

end

