function [F,x,vv_inv_M,I,mdata] = prob_X_L_1p_upt_decr_pde_ti(rho,pop_speed_option,Qc,Qr,vv,t,P_X, N_t_int, solver_options, display_flag)


if ~isstruct(solver_options)

    assert(false,'solver_options must be structure.');
    
end

%TODO: ignore saving I and mdata if not requested
% max_num_arg_ouput = nargout('prob_X_L_1p_upt_decr_pde_ti');
% num_arg_ouput=nargout;


if ~exist('display_flag')
    display_flag = 0;
end


% This function computes P( L_n(t) < s, A_n(t) = j )  for  j = 1 ,..., n
% where L_a(t) = int_0^t N_a(u)*1_{N_a(u)>=2} du
% Here n = number of individuals
% (xx,tt)-mesh which is a rectangular mesh determined by 
% x-interval [0,x_max] and t-interval [0,t_max] as [0,x_max] x [0,t_max]

% Here the generator has representation:
% Q = lambda * Qc + rho * Qr !!!


Kn = size(Qc,1); 

if ( (Kn < 2) ) 
    assert(false,'Incorrect input!');
end


if (size(Qc,1) ~= Kn) || (size(Qc,2) ~= Kn) || (size(Qr,1) ~= Kn) || (size(Qr,2) ~= Kn)
    assert(false,'Generators Qc, Qr have incorrect sizes');
end

if (numel(vv)~=Kn)
   assert(false,'state space vector vv must be of size Kn');
end


% check if v is decreasing and nonnegative

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


max_vv = max(vv);

speed_set_vv = union(sort(vv),[]);


% assign the number of the time slice where cone stops

N_t = numel(t);

if (N_t_int>N_t)
    N_t_int = N_t;
end


% Anonimous functions

lambda = @(t) pop_speed(t, 0, pop_speed_option); % population speed lambda

Lambda = @(t) pop_speed(t, 1, pop_speed_option); % cumulative speed Lambda


Qr = rho * Qr;

Qc_hat_T   = (Qc - diag(diag(Qc)))';

Qr_hat_T = (Qr - diag(diag(Qr)))';

q_c   = -diag(Qc);

q_rho = -diag(Qr);

genact_trans = @(m,t,Y) ( lambda(t) .* (Qc_hat_T(m,:) * Y) + Qr_hat_T(m,:) * Y );

% set variables
    
disp(' ')

disp('PDE computation for P( X(t)=s, L_{v}(t) < x ).')

time_start = cputime;

interp_option = 'linear';


% Mesh t = 0 = t_1 < t_2 < ... t_{max} is dictated by the ODE solver


% default value of x-mesh multiplier

mult_x   = max_vv;  % This number allows to add more points on x-mesh at each time t_i, mustbe >= 1 


if isfield(solver_options,'x_mesh_mult')
    
    if ~isempty(solver_options.x_mesh_mult)
        
        if (solver_options.x_mesh_mult<=1)
            
            mult_x = 1.0;
            
        else
            
            mult_x = solver_options.x_mesh_mult;
            
        end
        
    end

end

mult_x_e = mult_x; % (default) number allows to add more points on x-mesh at each extracted time 

if isfield(solver_options,'x_mesh_mult_extr')
    
    if ~isempty(solver_options.x_mesh_mult_extr)
        
        if (solver_options.x_mesh_mult_extr<=1)
            
            mult_x_e = 1.0;
            
        else
            
            mult_x_e = solver_options.x_mesh_mult_extr;
            
        end
        
    end
    
else
    
    solver_options.x_mesh_mult_extr = [];
    
end

if ~isfield(solver_options,'idx_times_e')
    solver_options.idx_times_e = [];
end


mult_x

mult_x_e


% Create fixed mesh in the triangular region 0 < s < nt
% and preallocate memory for u-values, g-values, s-values
% Here u=F_k, g=F_{k+1} is known


x = cell(1,N_t);
F = cell(1,N_t);
I = zeros(Kn,N_t,3);


x{1} = 0;
F{1} = zeros(Kn,1);    % Note that the "trace" for values s in Psi_d(n) is                        
                       % actually zero. Thus, in calculations we take these 
                       % values as zero. 

mdata = cell(Kn,N_t); % here for each i=2..N_t and m=1..Kn we will store the values of F
                      % at the points x{i}(j) with j=N_12(i)+1...N_34(i)-1
                      % and also at the corresponding points at earlier times 
                      % obtained by projection of x{i} backwards in time.


% creating a trapezoidal time-space mesh


x_max = max_vv * t(N_t_int);

mesh_option = 1;

for i=2:N_t

  for m=2:Kn
    mdata{m,i} = struct('xx',[],'tt',[],'uu',[],'gg',[],'length_1','length_2');
  end
  
 if (mesh_option==1)
               
     if ~isempty(find(i==solver_options.idx_times_e,1))
         
        mult_x_i = mult_x_e;
        
     else
         
        mult_x_i = mult_x;
        
     end
     
 
     if (i<=N_t_int)
      
       x{i} = linspace(0, max_vv * t(i), ceil(mult_x_i*(i-1)+1) );

     else

       x{i} = linspace(0,x_max, ceil( x_max*(mult_x_i*(i-1)+1)/(max_vv*t(i))));

     end
     
 else
     
     if (i<=N_t_int)
      
       x{i} = t(1:i)*max_vv;

     else

       x{i} = t(1:N_t_int)*(x_max/t(N_t_int));
       
       x{i}(end)=x_max;

     end
     
               
 end
 
 
 % TODO: Improve the choice of mesh step. For now a temporary solution:
  
 
    x_add = t(i)*(1:1:max_vv)/max_vv;
    
    idx  = find(x_add<x{i}(end));
   
    if ~isempty(idx)
        x{i} = union(x{i},x_add(idx));
    end
    
%     np=4;
%     
%     x_monitor=linspace(0,x{i}(end),np*10+1);
%       
%     x{i}=union(x{i},x_monitor);

    
%   if (numel(speed_set_vv)>1)
%          x_vv = speed_set_vv(1:(end-1))*t(i);
%          idx  = find(x_vv<x{i}(end));
%          if ~isempty(idx)
%              x{i} = union(x{i},x_vv(idx))';
%          end
%   end
   
 

  F{i} = zeros(Kn,length(x{i}));
          
end

% Precomputing e^(-(H_m(t_i)-H_m(t_{i-1}))) = 
%              = e^( - q_c(m) *(Lambda(t_i)-Lambda(t_{i-1})) - q_rho(m)*(t_i-t_{i-1})):

E = exp( - (  diff( Lambda(t) ) * q_c' + diff(t) * q_rho' ) ); 

% pre-computing cumulative speed values Lambda(t_{i})

Lambda_vals = Lambda(t);


% Input zone indices for elements in vv_inv_M(n)

vv_inv_M = find(max_vv==vv); % elements for which probability is 0 for x<Mt and P(A_{rho}(t)=*) for x > Mt.


for m=1:Kn
    
    I(m,1,:) = [1; 1; 1];

    % Setting up indices for s_m in Psi_d(n)
    
    if (ismember(m, vv_inv_M))
        
      for i=2:N_t
        I(m,i,:) = [length(x{i}); length(x{i}); length(x{i})];       
      end

    end
    
    % Setting up the boundary conditions (trace) for s_m not in Psi_d(n)
    
    if (~ismember(m,vv_inv_M))
        
        % Set up initial data for u = F_m on the line { x = n*t } up to t=t(N_t_int)
        
        for i=1:N_t_int
          F{i}(m,end) = P_X(i,m); % setting up the ODE boundary condition
        end
        
    end
    
end

        
% CDF for s_m in Psi_d(n) does not contribute anything to computation of other
% components as they all are zero in the region of computations. 

dprogress = 10;
    
progress = dprogress;

if (display_flag)
    disp(' '); disp('Computing 1-D time slices:');
end

for i=2:N_t
    
    if (display_flag)
        if (ceil(100*(i/N_t))>progress)
          disp( sprintf('progress: %d %%',progress));
          progress = progress+dprogress;
        end
    end    
    
    for m=1:Kn % Loop over all elements
    
     if ~ismember(m, vv_inv_M)   
        
         u = zeros(1,length(x{i}));
        
         % slope of characteristics
         
         v = vv(m);
         
         % Setting up the boundary indeces (that depend on k) 

        % Boundary indeces:

        dt = t(i)-t(i-1);
        
        %TODO: check if the method with float_eps variable incorporated in
        %inequalities works
                
        if (i<=N_t_int)
            
            I1 = find(x{i}<=v*t(i),1,'last');
            
%             I2 = find((x{i}<=(max_vv*t(i-1)+v*dt)),1,'last');

            I2 = find(x{i}<=(v*t(i)+(max_vv-v)*t(i-1)),1,'last');
            
            I3 = length(x{i});
        
            if (I1>I2)
                assert(false,'indices I1 and I2 for zones 1,2 are incorrectly computed!')
            end
            
            if (I2>I3)
                assert(false,'indices I2 and I3 for zones 2,3 are incorrectly computed!')
            end
            
        else
            
            I1 = find((x{i}<=v*t(i)),1,'last');

            I2 = length(x{i});

            I3 = length(x{i});
        
        end
               
        % Boundary indices for an s_m element 
        
        I(m,i,:) = [I1; I2; I3];
             
     
        % Compute values of u=F_m at each t(i), i=2,...N_t.

        % Zone-0: u = 0 in Zone-1

        gg_down = [];            

        gg_up   = [];          

        uu_down = [];

        % First, we get values of s-points from the time t(i-1)

        
        XI_up_1  = (I(m,i,1)+1):(I(m,i,2)); % for index I(m,i,1) u-value=0

        if ~isempty(XI_up_1)

             XI_down = I(m,i-1,1):I(m,i-1,3); 

             dt = t(i)-t(i-1);

             % Trace back x-points at t=t(i-1)

             x_up   = x{i}(XI_up_1);

             x_down = x_up - v * dt;

             % for interpolation purpose
             if (x{i-1}(XI_down(end))<x_down(end))                
                x_down(end)=x{i-1}(XI_down(end));
             end

             % Interpolate u, F, and g values at t = t(i-1) 
             
             % TO DO: NO NEED TO INTERPOLATE THE WHOLE VECTOR %

             F_down = (interp1(x{i-1}(XI_down), (F{i-1}(:,XI_down))', x_down, interp_option))';

             u_down = F_down(m,:);

             g_up   = genact_trans(m, t(i), F{i}(:,XI_up_1) );

             g_down = genact_trans(m, t(i-1), F_down );

             % Compute u_up in Zone-2:
             % Note in E(i-1,m) there is a shift in the index. In the
             % notes I denote this as E_{m,i).

             u(XI_up_1) = u_down * E(i-1,m) + (dt/2) * ( g_up + E(i-1,m) * g_down );                 

             gg_up   = g_up;

             gg_down = g_down;                                  

             uu_down = u_down;  

        end


         % Zone-2: Zone-2 appears only for t<=t_intermediate

         XI_up_2=[];
         
         if (i<=N_t_int)
             
             XI_up_2 = (I(m,i,2)+1):(I(m,i,3)-1);

             if   ~isempty(XI_up_2) 

                % Compute times on the line x = n*t

                x_up   = x{i}(XI_up_2);

                t_down = ( x_up - v * t(i) ) / ( max_vv - v );

                % Interpolate u and g values at t_down

                % for interpolation purpose
                if (t_down(1)<t(i-1))                
                 t_down(1)=t(i-1);                                  
                end

                F_down = interp1([t(i-1),t(i)], [F{i-1}(:,end),F{i}(:,end)]', t_down, interp_option)';

                u_down = F_down(m,:);

                g_down = genact_trans(m, t_down, F_down );

                g_up   = genact_trans(m, t(i), F{i}(:,XI_up_2) );

                % Compute u_up in Zone-3

                dt = t(i) - t_down; % vector of dt's

                % Precompute exp(-(H_m(t_i)-H_m(t_{i-1}))) 
                % = exp( - q_c(m) * (Lambda(t_i)-Lambda(t_{i-1})) - q_rho(m)*(t_i-t_{i-1})):

                E_bar_m = exp( -( (Lambda_vals(i)-Lambda(t_down)) * q_c(m) + dt * q_rho(m)) ); 

                u(XI_up_2) = u_down .* E_bar_m + (dt/2) .* ( g_up + E_bar_m .* g_down );            

                gg_down = [gg_down,g_down];                        

                gg_up = [gg_up,g_up];

                uu_down = [uu_down,u_down];  

             end
             
         end

         % update F_m term by entering u-values :

         bnd_i = (i<=N_t_int);
         
         F{i}(m,1:(end-bnd_i)) = u(1:(end-bnd_i));
         
         
        % record the values of u at x{i} for
        % i=I(m,i,1)+1 ... I(m,i,3)-1 as well as the values at
        % corresponding points at an earlier times t=t(i-1) 

        
        mdata_length_1 = numel(XI_up_1);
         
        mdata_length_2 = numel(XI_up_2);
         
        XI_up = [ XI_up_1, XI_up_2 ];

        mdata_length = mdata_length_1 + mdata_length_2;
        
        if (mdata_length>0) % otherwise there is no new data

             mdata{m,i}.xx = zeros(2,mdata_length);
             mdata{m,i}.tt = zeros(2,mdata_length);
             mdata{m,i}.uu = zeros(2,mdata_length);                                  
             mdata{m,i}.gg = zeros(2,mdata_length);

             % data at time t=t(i-1)

             if (mdata_length_1>0)
                mdata{m,i}.xx(1,1:mdata_length_1) = x_down; % these values are from Zone 2
                mdata{m,i}.tt(1,1:mdata_length_1) = t(i-1) * ones(1,mdata_length_1);
             end

             if (mdata_length_2>0)                    
                mdata{m,i}.xx(1,(mdata_length_1+1):mdata_length) = x{i}(XI_up_2) - v * dt;                
                mdata{m,i}.tt(1,(mdata_length_1+1):mdata_length) = t_down;                                
             end

             mdata{m,i}.uu(1,:) = uu_down; 
 
             mdata{m,i}.gg(1,:) = gg_down;

             % data at time t=t(i)

             mdata{m,i}.tt(2,:) = t(i) * ones(1,mdata_length);                                  

             mdata{m,i}.gg(2,:) = gg_up;

             mdata{m,i}.xx(2,:) = x{i}(XI_up);

             mdata{m,i}.uu(2,:) = u(XI_up);

        end

      end
       
    end
    
end

    if (display_flag)
          disp( sprintf('progress: 100 %%'));
          time_elapsed = cputime - time_start
    end


% setting up the ODE boundary condition to s_m in Psi_d(n)
% this is the boundary where jump happens. Though the trace for the PDE
% is zero for these components we set up these values anyway:

end