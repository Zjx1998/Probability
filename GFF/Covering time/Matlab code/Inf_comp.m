function [ ] = Inf_comp( eps )
% 

    iter    = 1e7;
    sam_num = 1e5;
    N       = zeros(sam_num, 1);
    for j = 1:sam_num
        S       = 0;                                                           % partial sum
        for i = 1:iter
            rnd = rand();
            if rnd > 0.5
                S   = S + 2;
            end
            S   = S - 1;
            if abs(S/i) < eps
                break;
            end
        end
        N(j)    = i;
    end
    N       = N * eps * eps;
    [y, x]  = hist(N,100);
    y       = y / length(N) / mean(diff(x));
    bar(x,y,1); 
    %histogram(N,100,'FaceColor','yellow');
    %histogram(N,100,'DisplayStyle','stairs');
        

end