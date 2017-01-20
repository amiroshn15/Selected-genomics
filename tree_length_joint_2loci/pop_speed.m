function f = pop_speed(t, type, option)

% t = vector or matrix of times

% type = type of speed 
% type = 0 returns lambda(t)
% type = 1 returns Lambda(t)=int_0^t lambda(s) ds
% type = 2 returns Lambda^{-1}(t)


% option = the choice of a particular speed function
% option = 1 lambda(t) = piecewise-constant function
% option = 2 continuous smooth function
% option = "else" lambda(t)=1, Lambda(t)=t, Lambda_inv(a)=a

% lambda(t) = 1/2 from 0 to 0.03
% lambda(t) = 5 from 0.03 to 0.12
% lambda(t) = 1/2 from 0.12 to \infty

% option is either integer number 0,1,2,3,...
% or it is a structure option.opt, option.param; the latter contains the
% parameters.


if isstruct(option)
    
    opt_val = option.val;            
    
else
    
    opt_val = option;   
    
end

    
if (opt_val==1) % piecewise-continuous lambda(t)

        % intervals of continuity, 0.0 must be included:
            
        if isstruct(option)
           
            if isfield(option,'param')
                
                if isfield(option.param,'times')
            
                    times = option.param.times;
        
                end
    
                if isfield(option.param,'lambda_vals')

                    lambda_vals = option.param.lambda_vals;

                end

                 
            end
            
        
        else
            
            times  = [0.0, 0.05, 0.2 ];
            
            lambda_vals = [0.5, 4, 0.5]; 
            
        end           
    
        % compute cumulutive speed at the nodes times(i)
        
        
        Lambda_vals = pcont_cum_speed_node_vals(times,lambda_vals);
        
        
    if (type==1)        

        % this function requires the values at the nodes Lambda(times(i))
        % and the last slope on the interval [times(end),\infty]:

        f = pcont_cum_speed(t, times, Lambda_vals, lambda_vals(end));

    elseif (type==2)       
        
        % the same as the comment above:
        
        f = pcont_cum_speed_inv(t, times, Lambda_vals, lambda_vals(end));
        
    else
        
        f = pconst_speed(t,times,lambda_vals);          
        
    end
    
    
elseif (opt_val==2) % smooth continuous lambda(t) = 2 - 1/(1+t)^2
    
    if (type==1)
    
        f = 2 .*t - t ./ (1+t);

    elseif (type==2)
        
        f = (1/4) * ( -(1-t) + ((1-t).^2 + 8*t ).^(0.5));
        
    else
        
        f = 2  - (1+t).^(-2);
        
    end
               
elseif (opt_val==3) % peice-wise continuous lambda(t) given by pcexp_speed(t)
    
    
    if isstruct(option)
        
        if isfield(option,'param')
            
            prm = option.param;
            
            % by accident in this option I implmented c(x)=1/lambda(x)
            % to handle this mistake one simply need to change parameters
            % to N_A -> 1/N_A, N_B=1/N_B and g-> -g
            % these are parameters for the speed function
            
            prm.N_A = 1/prm.N_A;
            
            prm.N_B = 1/prm.N_B;
            
            prm.g  = -prm.g;
            
        else
            assert(false,'field param is missing in the structure option');
        end
        
        
    else  % these are parameters for the speed function
    
        prm.N_A = 1/2;  
        
        prm.N_B = 1/0.25;  
        
        prm.T_B = 0.15;
        
        prm.T_G = 0.025;
        
        prm.g = (-1)*200;
        
    end
    
    if (prm.g==0)
        
        p=[];
        
        p.times = [0.0, prm.T_B];  p.lambda_vals = [ prm.N_B, prm.N_A]; 
        
        lambda_opt.val = 1;      % piecewise-constant lambda(t)
        
        lambda_opt.param = p;
        
        f = pop_speed(t,type,lambda_opt);

    else
        
    
        if (type==1)

            % this function requires the values at the nodes Lambda(times(i))
            % and the last slope on the interval [times(end),\infty]:

            f = arrayfun(@(t) pcexp_cum_speed(t,prm),t);               

        elseif (type==2)       

            % it is not implemented yet:

            f = pcexp_cum_speed_inv(t,prm);

        else

            f = arrayfun(@(t) pcexp_speed(t,prm),t); 

        end
    
    end
        
else    % trivial case lambda(t)=1
    
        if ( (type==1) || (type==2) )
            
            f = t;
            
        else            
             f=ones(size(t));
        end        
    
end

end