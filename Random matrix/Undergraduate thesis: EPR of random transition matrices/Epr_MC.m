function [ output ] = Epr_MC( opts )
%--------------------------------------------------------------------------
% Estimation of entropy production rate using thermodynamic uncertainty
% relation
%--------------------------------------------------------------------------


    global P

    
    if ~isfield(opts,'epoch');      opts.epoch      = 100000;     end
	if ~isfield(opts,'batch');      opts.batch      = 10000;     end
	if ~isfield(opts,'source');     opts.source     = 2;        end
    if ~isfield(opts,'target');     opts.target     = 3;        end

    
	% copy paramter
    n           = size(P, 1);
	epoch       = opts.epoch;                                              % simulation epoch
    batch       = opts.batch;                                              % parallel sample path
    source      = opts.source;                                             % target state for current
    target      = opts.target;                                             % target state for current

    
    % auxiliary parameter in optimization
    tr          = @trans;                                                  % function handle for transition
    c           = @curr;                                                   % function handle for current
    state       = ones(batch, 1);                                         % parameters to optimize                           
    current     = zeros(batch, 1);                                         % approximate Hessian
    
    
    % set up print format
% 	if itPrint > 0
% 	    if ispc; str1 = '  %10s'; str2 = '  %7s';
% 	    else     str1 = '  %10s'; str2 = '  %7s'; end
% 	    stra = ['%5s', str2, str2, str2, str2, str2 '\n'];
% 	    str_head = sprintf(stra, 'iter', 'obj', 'obj2', 'mu', 's', 't');
% 	    str_num = ['%4d %2.1e %+2.1e %+2.1e %2.1e %2.1e \n'];    
% 	end


    %tic;
    prob    = rand(batch, epoch);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MAIN LOOP: first stage and second stage
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for t = 1:epoch
        
        
        old_state   = state;
        state       = tr(state, prob(:, t));
        current     = current + c(old_state, state);

                
%         % save information for graphic output
%         if opts.plott
%             time_temp   = toc;
%             obj_full    = f_handle.get_f(x);
%             res_full    = norm(F_handle.get_F(x,f_handle.get_grad(x),1));
%             eval_time   = eval_time + toc - time_temp;
% 
%             plott.epoch(iter)   = nr_epoch;
%             plott.time(iter)    = toc - eval_time;
%             plott.obj(iter)     = obj_full;
%             plott.res(iter)     = res_full;
%             plott.tau(iter)     = tau;
%         end
        
        
    end
    %toc;
    

    output.curr_mean    = mean(current);
    output.curr_var     = var(current);
    output.est_epr      = output.curr_var * epoch/(output.curr_mean)^2;
    output.time         = toc;
    
    
    function [ new_s ] = trans(s, p)
    % parallel version of the function which calculate the transition Markov chain        
        new_s       = ones(batch, 1);
        for iter_s  = 1:n
            p       = p - P(s, iter_s);
            new_s   = new_s + (p > 0);
        end        
    end


    function [ cc ] = curr(o_s, n_s)
    % parallel version of the function which calculate the transition Markov chain        
        cc = (o_s == source) .* (n_s == target) - (n_s == source) .* (o_s == target);        
    end
    

end