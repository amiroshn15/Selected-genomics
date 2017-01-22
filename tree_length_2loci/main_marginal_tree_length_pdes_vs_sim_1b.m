%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code has been downloaded from :
% https://github.com/amiroshn15/Selected-genomics/tree_length_2loci
%
% Authors:
% Alexey Miroshnikov <amiroshn@gmail.com>
% Matthias Steinrücken <steinruecken@schoolph.umass.edu>
%
% License: 	GPL-2 | GPL-3 [expanded from: GPL (>=2)]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Set 1b] In this script we compare computations of
%
%          F(t,x):= { F_k=P( A(t)=k, L(t) < x ) }_{k=1,2,...,n}                 
%                   with t in [0,t_max] and x in [0, x_max]
%
%    by use of PDEs and stochastic simulations.
%
%    A(t,w)- ancestral (marginal) process for one locus on a chromosose 
%            for the case of a variable size population with the coalescent 
%            rate function lambda(t) which is set in 'pop_speed_opt.val'.
%
%    L(t,w)= int_0^t A(s,w) ds - accumulated length of the tree.
%
%    Here: (A(t),L(t)) is an inhomogenous time-continuous Markov process.
%
%   The description of the process (A(t),L(t)) and PDE equations for CDF
%   described in the journal article:
%
%    Alexey Miroshnikov, Matthias Steinrücken
%    "The marginal and joint distributions of the total
%       tree lengths across loci in populations with variable size".
%
%    Preprint is available at https://arxiv.org/pdf/1609.08880v3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1b.1 PDE Computations of F = P( A(t)=k , L(t)<x )

orig_path = path();    % save original path

addpath(genpath(pwd)); % add folders where the script is to search path

n = 5; % number of individuals

state_space = n:(-1):1; % return the state space S={n,n-1,...1}

v = state_space; v(end)=0; % define a state space function v(s)

Q=A_n_gen(n); % get an ancestral generator (responsible for coalescence)

% The state number m for which we compare F_m PDE vs Simulations
% Note that if m=n then the state S_m=S_n=1 is an absorbing state.

m = n; 

% set the coalescent rate function (piecewise-exponential lambda(t))

pop_speed_opt = [];

pop_speed_opt.val = 3; % piecewise-exponential lambda(t)

% coalescent rate function parameters 

prm.N_A = 2;  
prm.N_B = 0.25;  
prm.T_B = 0.15;  
prm.T_G = 0.025;  
prm.g = 200;
pop_speed_opt.param = prm;

% preliminary time boundary (later is changed according to prob_A)

t_max_prelim = 5.0; % maximal t-interval

% Probability lists for (x,t)-mesh. For each value prob_A(i) in the list 
% prob_A we find the time t_i such that P(A(t)=1)=prob_A(i). We then 
% compare PDE solutions to simulated ancestral process for each time t_i at
% multiple x-values.

prob_A = [0.1, 0.2, 0.3, 0.4, 0.5]; 

ode_options = odeset('NormControl','on','RelTol', 3e-13,...
                     'AbsTol', 3e-12, 'MaxStep', 0.0005);


% first check if there is time for which P_A_m(t) > max(prob_A)

rho = 0.0; % recombination rate is zero when computing marginal.

[T,P_A] = prob_X_ode( rho, pop_speed_opt, Q, Q, ...
                      [0,t_max_prelim], eye(1,n), ode_options );

P_A_m = sum(P_A(:,m),2);

% Find times t for which P(A(t)=1)<= prob_A(end)

t = zeros(size(prob_A));

for r=1:numel(prob_A)
    
    idx = find(P_A_m>=prob_A(r),1);
    
    if ~isempty(idx)

        t(r) = round(T(idx),1); % maximal t-value
        
        err = P_A_m(idx)-prob_A(r);
                
    else
        assert(false,'choose time larger!!!')
    end

end
 
display(t,'times used for comparison')

t_max = t(end);

x_max = t(end) * n;

ode_options = odeset('NormControl','on','RelTol', 3e-11,'AbsTol', 3e-10);

[F,X,T,P_A,N_t_int,v_inv_max] = ...
    prob_X_L_1p_upt_decr_pde_t(rho,pop_speed_opt,...
                               Q,Q,v,[x_max,t_max],ode_options,1);

path(orig_path);

%% Plot F_m(x,t) computed via PDE

orig_path = path();    % save original path

addpath(genpath(pwd)); % add folders where the script is to search path

% make a finer (xx,tt)-mesh
 
tt = linspace(t(1), t(end), 401);

xx = linspace(0.0, x_max, 401);

% Interpolate PDE on (xx,tt)-mesh

P_A_L_pde_fine = ...
prob_X_L_1p_upt_decr_interp_t(F,X,T,P_A,N_t_int,v_inv_max,xx,tt,'linear');

% Plot the surface on (xx,tt)-mesh

h_fig = figure('units','normalized','outerposition',[0.2 0.1 0.6 0.9]); 

set(gca,'FontSize',15);  % set the fontsize of axis

surface(xx,tt,P_A_L_pde_fine(:,:,n),...
        'FaceColor','interp','EdgeColor','none'); 
    
colormap(jet(256)); view(3);

xlabel('x'); ylabel('t');

title(['P_{PDE}( A(t) = ',num2str(state_space(m)),', L(t) < x ),',...
        ' n=',num2str(n), ', \lambda=piecewise-exponential']);

path(orig_path);

%% Create (x,t)-mesh used in simulations. Interpolate and plot F_m-pde.

orig_path = path();    % save original path

addpath(genpath(pwd)); % add folders where the script is to search path

% The mesh (x,t) has t-points which are computed above.

N_t = numel(t); 

N_x = 10;  % number of points in x-direction

x_domain_portion = 0.5; % prob. values are very flat for most of the domain

x = linspace(1.0, x_max*x_domain_portion, N_x);

% Interpolate PDE on (x,t)-mesh

P_A_L_pde = ...
   prob_X_L_1p_upt_decr_interp_t(F,X,T,P_A,N_t_int,v_inv_max,x,t,'linear');

Fm_pde = P_A_L_pde(:,:,m)';

% Plot the surface at (x,t)-mesh 

h_fig = figure('units','normalized','outerposition',outpos); hold on; 

set(gca,'FontSize',15);  % set the fontsize of axis

surface(x,t,Fm_pde'); colormap(jet(256)); view(3);

xlabel('x'); ylabel('t');

title(['P_{PDE}( A(t) = ',num2str(state_space(m)),', L(t) < x ),',...
        ' n=',num2str(n), ', \lambda=piecewise-exponential']);

path(orig_path);

%% Perform Simulations. Simulates trajectories and Computes Length

orig_path = path();    % save original path

addpath(genpath(pwd)); % add folders where the script is to search path

rand('seed',100);

% array that contains the number of trajectories

N_traj = 4.^[5,6,7,8,9];

[ P_A_L_sim, P_A_sim, SE_A_L] = ...
    prob_X_L_1p_upt_sim( 0.0, pop_speed_opt, Q, Q, v, x, t, N_traj,0);

Fm_sim = squeeze(P_A_L_sim(:,:,m,:));   % solution for m-th component

SEm = squeeze(SE_A_L(:,:,m,:));         % standard error

% Set a desired subset for the comparison table

N_traj_set = N_traj(1:end);

idx_traj = zeros(size(N_traj_set));

for j=1:numel(N_traj_set)

    idx_traj(j)=find(N_traj==N_traj_set(j));
    
end

path(orig_path);

%% Plot log-difference of PDE and Simulation version of F_m 

for jj=1:numel(idx_traj)
    
    j=idx_traj(jj);
    
    surf(x,t,log10(abs(Fm_pde(:,:)'-Fm_sim(:,:,j)')),'FaceColor','flat');    
     
    colormap(jet(256)); hold on;   
end

set(gca,'FontSize',10);  % set the fontsize of axis

xlabel('x'); ylabel('t');

title(['log|P_{PDE}-P_{SIM}|, A(t)=', num2str(state_space(m)),...
        ', n=',num2str(n), ', \lambda=piecewise-exponential']);

%% Make an Comparison Table PDEs versus Simulations

% Chose a subset of x-points for each time slice for comparison table.

prop_x = 0.5; % proportion of subset x-points

idx_x =1:floor(1/prop_x):numel(x);

num_subset_x = numel(idx_x);

Ridx = zeros(N_t,num_subset_x);

for i=1:N_t
    Ridx(i,:)=idx_x;
end

TABLE = [];

for i=1:N_t
        
    for rr=1:numel(Ridx(i,:))        
        r = Ridx(i,rr);
        row = [t(i), x(r), Fm_pde(r,i)];
        for jj=1:numel(idx_traj)
            j=idx_traj(jj);
            row = [row, Fm_sim(r,i,j), N_traj(j), SEm(r,i,j) ];
        end            
        TABLE = [ TABLE; row ];
    end    
end

% save table

filename = 'set_1b_TABLE_P_A_L.txt';

fID = fopen(filename,'w');

header = ['t    ','|x    ','|P-PDE     '];

for l=1:numel(idx_traj)
    header = [header,'|P-SIM     ','|NTR       ', '|SER       '];
end

header=[header,'|'];

 fprintf(fID,header);       
 fprintf(fID,'\n');       
 
for l=1:size(TABLE,1)        
 fprintf(fID,' %1.1f |',TABLE(l,1:2));
 fprintf(fID,' %1.6f |',TABLE(l,3));
 for jj=1:numel(idx_traj)
      fprintf(fID,' %1.6f |',TABLE(l,1+jj*3));
      fprintf(fID,' %8.0u |',TABLE(l,2+jj*3));
      fprintf(fID,' %1.6f |',TABLE(l,3+jj*3));
 end     
 fprintf(fID,'\n');       
end

fclose(fID);

% ouput TABLE on display

fID = fopen(filename,'r');

tline = fgetl(fID);

while ischar(tline)

    disp(tline)
    tline = fgetl(fID);
end

fclose(fID);

%% Plot standard error as a function of N_traj for subset (x,t)-point
        
np=1:(N_t*N_x);

figure(); hold on; j=1;

for i=1:N_t
    for rr=1:numel(Ridx(i,:))
        r=Ridx(i,rr);

        plot(log10(N_traj),log10(squeeze(SEm(r,i,:))),'-*',...
            'Color',[double(np(j)/np(end)),0.6*double(np(j)/np(end)),...
            0.6*double(np(j)/np(end))]); j=j+1;
    end
end

title('Log-SERR( A(t) = 1, L(t) < x )');

xlabel(strcat('log10(N_{traj})',', n=',num2str(n),', state s=1'));  

hold off;
