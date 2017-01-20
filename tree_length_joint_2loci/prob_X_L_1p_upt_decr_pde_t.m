function [F,x,t,P_X,N_t_int,vv_inv_M,I,mdata] = prob_X_L_1p_upt_decr_pde_t(rho,pop_speed_option,Qc,Qr,vv, dom_bnd, solver_options, display_flag)


if isstruct(solver_options)
    
    if isfield(solver_options,'ode_options') % first we check if solver_options has a subfield ode_options
        
        ode_options = solver_options.ode_options;
        
    else % if not then we assume that solver_options is the structure containing only fields of odeset()
        
        ode_options = solver_options;
        
        solver_options = [];
        
        solver_options.options = ode_options;
        
    end
            
else
    assert(false,'pde options must be structure');
end


if ~exist('display_flag')
    
    display_flag=0;
    
end

if ( numel(dom_bnd) ~= 2 )
    assert(false,'Incorrect input for dom_bnd=[x_max, t_max]!');
end

x_max = dom_bnd(1);

t_max = dom_bnd(2);


K = size(Qc,1); 

if ( x_max<=0 )
    assert(false,'x_max must be > 0');
end

if ( t_max<=0 )
    assert(false,'t_max must be > 0');
end

if (numel(vv)~=K)
   assert(false,'state space vector vv must be of size Kn');
end

% check if maximal value of vv is positive

max_vv = max(vv);

if (max_vv<=0)
    assert(false,'maximal value of vv must be positive');
end

t_int = x_max/max_vv;

% y_0 = zeros(1,K); y_0(1)=1; % initial condition 

y_0 = eye(1,K); % initial condition 

% y_0 = ones(1,K)/K; % initial condition 

if (t_int>=t_max)    
    [t,P_X] = prob_X_ode(rho, pop_speed_option, Qc, Qr, [0.0, t_max], y_0, ode_options );     
    N_t_int = length(t);    
else   
   [t1, P_X1] = prob_X_ode(rho, pop_speed_option, Qc, Qr, [0.0, t_int], y_0, ode_options );
    N_t_int = numel(t1);
    y_0 = P_X1(end,:);
    [t2, P_X2] = prob_X_ode(rho, pop_speed_option, Qc, Qr, [t_int,t_max], y_0, ode_options );  
    t = [t1(1:end);t2(2:end)];       
    P_X = [P_X1(1:end,:);P_X2(2:end,:)];      
end

% TODO Check this slight improvement:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(solver_options,'ode_options_fine') % recompute on the mesh t values of P_X with better accuracy

    [t,P_X] = prob_X_ode( rho, pop_speed_option, Qc, Qr, t, eye(1,K) , solver_options.ode_options_fine );
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_t = numel(t)

[F,x,vv_inv_M,I,mdata] = prob_X_L_1p_upt_decr_pde_ti(rho,pop_speed_option,Qc,Qr,vv,t,P_X, N_t_int, solver_options, display_flag);

end