
function [P_X_L,P_X,SE_X_L] = prob_X_L_1p_upt_sim( rho, pop_speed_option, Qc, Qr, vv, x , t , N_traj_arr, silent_display )


if (nargout>=1)
    expect.P_X_L = true;
else
    expect.P_X_L = false;
end

if (nargout>=2)
    expect.P_X = true;
else
    expect.P_X = false;
end

if (nargout>=3)
    expect.SE_X_L = true;
else
    expect.SE_X_L = false;
end


% This functions returns the empirically computed probabilities:
% 
%  { P( X(t(i)) = state )}_{i=1..length(t)} 
%  { P( X(t(i)) = state, L(t)<x )}_{i=1..length(t)} 


% Here:
%
% rho = recombination rate
%
% pop_speed_option = the choice of the speed function lambda(t):
%
%  0 = trivial, lambda(t) = 1;  
%  1 = piecewise constant;
%  2 = smooth.
%
% t = array of times   >= 0
% x = array of lengths >= 0

    time_start = cputime;

    % check if matrices Qc and Qr have appropiate sizez
    
    K = size(Qc,1); 

    if ( (K < 2) ) 
        assert(false,'Incorrect input!');
    end

    if (size(Qc,1) ~= K) || (size(Qc,2) ~= K) || (size(Qr,1) ~= K) || (size(Qr,2) ~= K)
        assert(false,'Generators Qc, Qr have incorrect sizes');
    end

    if (numel(vv)~=K)
       assert(false,'state space vector vv must be of size K');
    end

    % TODO: No need to check if vv is decreasing and nonnegative. 
    
    if ~isempty(find(vv<0,1))
        assert(false,'vv must be nonnegative');
    end

    [ii,jj] = find(Qc);

    decreasing=1;
    up_triang=1;

    for k=1:numel(ii)
        if (ii(k)~=jj(k))
            if (vv(ii(k))<vv(jj(k)))
                decreasing=0;
            end
        end
        if (ii(k)>jj(k)), up_triang=0; end
    end

    [ii,jj] = find(Qr);

    for k=1:numel(ii)
        if (ii(k)~=jj(k))
            if (vv(ii(k))<vv(jj(k)))
                decreasing=0;
            end
        end    
        if (ii(k)>jj(k)), up_triang=0; end
    end

    if (decreasing==0)
        assert(false,'v(process(t)) must be non-increasing');
    end

    if (up_triang==0)
        assert(false,'each component of Q must be upper-triangular');
    end
    
    % TODO: Take care of the situation rho=0
    
    % get the diagonal
    
    % the routine works with Q=lambda(t)*Qc + Qr*rho. To make things fast
    % we will redefine Qr by Qr * rho.
    
    Qr  = rho*Qr;       
    q_c = -diag(Qc);    
    q_r = -diag(Qr);

    qq_r = q_r;
    
%     qq_r = rho*q_r;

    % detemine sets Delta_c and Delta_r
    
    Delta_c = cell(1,K);            
    Delta_r = cell(1,K);
    
    P_c = cell(1,K);
    P_r = cell(1,K);            
    
    abs_state = zeros(1,K);
    
    for i=1:K
        
        tmp = Qc(i,:); tmp(i)=0;  Delta_c{i}=find(tmp);        
        
        tmp = Qr(i,:); tmp(i)=0;  Delta_r{i}=find(tmp);
        
        if ~isempty(Delta_c{i})
            P_c{i}=Qc(i,Delta_c{i});
        end

        if ~isempty(Delta_r{i})
            P_r{i}=Qr(i,Delta_r{i});
        end

        % check that Delta_c and Delta_r are non-intersecting     
        
        if ~isempty(intersect(Delta_c{i},Delta_r{i}))
            assert(false,'Connectivity matrices for Qc and Qr must have no coommon entries except diagonals');
        end
    
        % determine is 'i' is absorbing
        if isempty(Delta_c{i}) && isempty(Delta_r{i})
           abs_state(i)=1;
        end
    end
           
        
    % check that the last state is absorbing    
    if ~(isempty(Delta_c{K}) && isempty(Delta_r{K}))
        assert(false,'last state must be absorbing');
    end
    
    % %


    if (sum(N_traj_arr<1)>0) || (sum(diff(N_traj_arr)<1)>0)
     assert (false, 'N_{traj} must be array of positive numbers ');
    end

    if (isempty(t))
        assert (false, 't must be nonempty array ');
    end

    if (sum(t<0) > 0 )
        assert (false, 'times t must be nonnegative');
    end  

    if (sum(x<0) > 0 )
        assert (false, 'lengths L_a must be nonnegative');
    end

    
    if ( nargin < 9 )
        display_flag = 1;
    else
        display_flag = ~(silent_display);
    end

    
    % some anonimous useful functions:

    Lambda     = @(t) pop_speed(t, 1, pop_speed_option);

    Lambda_inv = @(t) pop_speed(t, 2, pop_speed_option);
    
    
    num_N_traj_arr = numel(N_traj_arr);
    
    N_t   = length(t);
   
    N_x   = length(x);

    hits_X = zeros(N_t,K);
   
    hits_X_L = zeros(N_x, N_t, K);
    
    P_X = zeros(N_t,K,num_N_traj_arr);

    P_X_L = zeros(N_x, N_t, K, num_N_traj_arr);

    
  % precomputing coalescence rates (for speeding up)
  
 
for r=1:num_N_traj_arr

    if (display_flag)
        disp(' ')
    end
    
    if (r==1)
        N_traj = N_traj_arr(1);
        if (display_flag)
            disp(strcat('N_{traj_interval} = [','1,',num2str(N_traj_arr(r)),']:'))
        end     
        
    else
        N_traj = N_traj_arr(r) - N_traj_arr(r-1);
        if (display_flag)
            disp(strcat('N_{traj_interval} = [',num2str(N_traj_arr(r-1)+1),',',num2str(N_traj_arr(r)),']'))        
        end
    end

    if (display_flag)
       disp(strcat('Precomputing i.i.d. Random Variables' ))
    end
    
   % Generate exponential random variables:
   % Maximal number of jumps cannot exceed K-1 and hence the chain visits
   % no more than K states.

    E = zeros(N_traj, K, 2);
    E(:,:,1) = exprnd( 1.0, N_traj, K );    
    E(:,:,2) = exprnd( 1.0, N_traj, K );    
   

    rnd_Delta_c = cell(1,K);
    rnd_Delta_r = cell(1,K);
        
    
    for i=1:K
        
        if ~isempty(Delta_c{i})
            rnd_Delta_c{i} = Delta_c{i}(rand_choose2(P_c{i},N_traj));
        end
        
        if ~isempty(Delta_r{i})
            rnd_Delta_r{i} = Delta_r{i}(rand_choose2(P_r{i},N_traj));
        end
        
    end
    
    % Simulate trajectories:

    if (display_flag)
        disp(strcat('Generating Trajectories: N_traj= ', num2str(N_traj)));
    end

    if (display_flag)
        disp( sprintf('progress: 0 %%'));    
    end
    
    dprogress = 10;
    
    progress = dprogress;
    
    
    for j=1:N_traj

        if (display_flag)
            if (ceil(100*(j/N_traj))>progress)
              disp( sprintf('progress: %d %%',progress));
              progress = progress+dprogress;
            end
        end

        % preallocate n+2 elements for embedded chain. 

        % note: n+2 is maximally possible length (including the initial state)

        T  = zeros(1,K); %  contains jump times. T(end)=T_{Ult. MRCA}        
        XX = zeros(1,K); %  jump chain for t >= 0 on the state space 1,2,...,K
        LL = zeros(1,K);   %  Path-integral value up to N_state-th state: L(i)=int_0^{T{i-1}} vv(X(t)) dt
        
        % TODO: Determine the maximal number of jumps in the trajectory:
                
        % initiate variables:

        m = 1;          % Counter of the number of states the chain jumped through 
                              % N_state-1 is the number of jumps up to now :
                              
        XX(1)   = 1;          % initial state, start always from state 1.
        T(1)    = 0.0;        % initial time of the initial state
        LL(1)   = 0.0;        % the value of the integral up to the current state is zero
        
        expect_absorbing_state = 1;
        
% % This Part Computes Trajectory
% %[[

        if ~abs_state(XX(m))  
        
            while(expect_absorbing_state==1)           

                 X_m = XX(m);
                                                                                                
                % step 1
                
                if isempty(Delta_c{X_m}) 
                    tau_1=inf;
                else
                    tau_1 = Lambda_inv( (E(j, m, 1)/q_c(X_m)) + Lambda(T(m)) ) - T(m);
                end

                
                % step 2
                
                if isempty(Delta_r{X_m}) 
                    tau_2=inf;
                else
                    tau_2 = E(j, m, 2)/(qq_r(X_m));
                end
                
                                
                % step 3
                
                [tau, idx] = min([tau_1,tau_2]); % the inter-event time        
                
                % step 4  
                                    
                 if (idx==1) % chains jumps to Delta_c - set
                    
                     XX(m+1)=rnd_Delta_c{X_m}(j);
                    
                 else % chains jumps to Delta_r - set
                    
                     XX(m+1)=rnd_Delta_r{X_m}(j);
                     
                 end
                
                % Compute the path-integral up to XX(m+1): 
                                
                LL(m+1) = LL(m) + vv(XX(m))*tau;
                        
                % Record the number of jumps and the time of the jump
                 
                T(m+1) =  T(m) + tau;  % time of entering into new state
                
                if abs_state(XX(m+1)), expect_absorbing_state = 0;  end   
            
                m = m+1; % the number of the states the chain visited up to now
                               
            end
                                    
        end
                 
% % ]]
        
        N_states = m; % N_state now contains the number of states in the jump chain
        
        % XX(N_states) contains an absorbing state where chain stops. T(N_states)=T_{abs}.
                       
        idx_t_m_1 = 1;
        
        expect_quit = 1;
        
        m=1;
        
        while(expect_quit==1)                       
            
            if (m==N_states)
               idx_t_m_2 = N_t;               
            else
               idx_t_m_2 = (idx_t_m_1-1) + find((t(idx_t_m_1:end)<T(m+1)),1,'last');
            end                            
            
            if ~isempty(idx_t_m_2)

                X_m = XX(m);                

                hits_X(idx_t_m_1:idx_t_m_2, X_m) = hits_X(idx_t_m_1:idx_t_m_2,X_m) + 1;
                
                L_t_m = LL(m) +  vv(X_m) * (t(idx_t_m_1:idx_t_m_2)-T(m));

                jj = 1;

                for ii=idx_t_m_1:idx_t_m_2                             

                 hits_X_L(:,ii,X_m) = hits_X_L(:,ii,X_m) + (L_t_m(jj)<x)';                          

                 jj = jj+1;
                end

                idx_t_m_1 = idx_t_m_2+1;

            end                         
            
            if (idx_t_m_1<=N_t)                
                m=m+1;                
            else                
                expect_quit=0;
            end

        end
        
    end       

if (display_flag)
    disp( sprintf('progress: %d %%',progress));        
end


P_X(:,:,r) = hits_X/N_traj_arr(r);

P_X_L(:,:,:,r) = hits_X_L/N_traj_arr(r);

end


if (expect.SE_X_L)
    
    SE_X_L = (P_X_L - P_X_L.^2);

    for r=1:numel(N_traj_arr)
        SE_X_L(:,:,:,r)=sqrt(SE_X_L(:,:,:,r)/(N_traj_arr(r)-1));
    end
    
end


if (display_flag)
    time_elapsed = cputime - time_start
end

end

