function [ cov ] = Green_func(x ,y , n)
% This function numerically calculate the value of 2D Green function by the
% defintion via simple random walk.

    sample_num      = 1e5;
    x1              = floor((x + n - 1)/n);
    x2              = mod(x - 1, n) + 1;
    y1              = floor((y + n - 1)/n);
    y2              = mod(y - 1, n) + 1;
    count_          = 0;
    for i = 1:sample_num
        pt          = [x1, x2];                                            % starting from x
        while ~BC(pt, n)                                                   % not hitting boundary
            if sum(abs(pt - [y1, y2])) == 0                                % hitting y
                count_  = count_ + 1;
            end
            pt      = next(pt);                                            % sample next walk
        end
    end
    cov             = count_ / sample_num;
      
    function [ state ] = BC(pt, n)
        state       = false;
        if abs((pt(1) - n - 1) * (pt(2) - n - 1) * pt(1) * pt(2)) < eps
            state   = true;
        end
    end

    function [ pt_ ] = next(pt)
        p       = rand();
        pt_     = pt;
        if(p < 0.25)
            pt_(1)  = pt_(1) + 1;
        elseif(p < 0.5)
            pt_(1)  = pt_(1) - 1;
        elseif(p < 0.75)
            pt_(2)  = pt_(2) + 1;
        else
            pt_(2)  = pt_(2) - 1;
        end
    end

end