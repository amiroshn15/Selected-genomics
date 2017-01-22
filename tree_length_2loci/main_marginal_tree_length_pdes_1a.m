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
% [Set 1a] In this script we compute and then plot 
%
%          F(t,x):= { P( A(t)=k, L(t) < x ) }_{k=1,2,...,n}                 
%                   with t in [0,t_max] and x in [0, x_max]
%    where
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
%   Alexey Miroshnikov, Matthias Steinrücken
%   "The marginal and joint distributions of the total
%      tree lengths across loci in populations with variable size".
%
%   Preprint is available at https://arxiv.org/pdf/1609.08880v3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1a.1 Compute P( A(t)=k, L(t) < x ), k=1,2,...,n                 

orig_path = path();    % save original path

addpath(genpath(pwd)); % add folders where the script is to search path

n = 5; % number of individuals

state_space = n:(-1):1; % return the set of state space enumerated as 
                        % S={ n, n-1,...,2,1}
                        
v = state_space; v(end)=0;

Q=A_n_gen(n); % get the ancestral generator for n individuals

% parameters for population rate function lambda(t) of option 3:

prm.N_A = 2; prm.N_B = 0.25; prm.T_B = 0.15; prm.T_G = 0.025; prm.g = 200;

pop_speed_opt = [];

pop_speed_opt.val = 3; % = piecewise-exponential lambda(t) with bottleneck.

pop_speed_opt.param = prm;

% % time ans space domain boundaries

t_max = 0.5;      % maximal t-interval
 
x_max = n*t_max;  % maximal x-interval

ode_options = odeset('NormControl','on','RelTol', 3e-11,...
                     'AbsTol', 3e-11,'MaxStep', 0.25*10^(-2));

rho = 0.0; % recombination rate is zero when computing marginal.

[F,X,T,P_A,N_t_int,v_inv_max] = ...
    prob_X_L_1p_upt_decr_pde_t(rho,pop_speed_opt,...
                               Q,Q,v,[x_max,t_max],ode_options,1);

% interpolation onto a rectangular (x,t)-mesh

N_t = 1001;
 
N_x = (N_t-1)+1;
 
t = linspace(0.0, t_max, N_t);

x = linspace(0.0, x_max, N_x);

P_A_L_pde = prob_X_L_1p_upt_decr_interp_t(F, X, T, P_A,N_t_int,...
                                          v_inv_max,x,t,'linear');
                                                                        
path(orig_path); % return original path

%% 1a.2 Plot P(A(t)=k,L(t)<x)

fntsz = 20;

lnwidth = 1;

for m=1:n

    k = n - m + 1; 
    
    % plot surface:
    
    h_fig=figure('units','normalized','outerposition',[0.2 0.1 0.6 0.9]);
    
    hold on;
    
    xlim([x(1), x(end)]);
    
    surface(x,t,P_A_L_pde(:,:,m),'FaceColor','interp','EdgeColor','none');
    colormap(jet(256));    
    
    set(gca,'FontSize',fntsz);  % set the fontsize of axis
    
    colorbar('location','eastoutside','FontSize',fntsz); % make colorbar
                
    % plot lines x=(k*1_{k>1})*t and x=n*t:
    
    c = max(max(P_A_L_pde(:,:,m)));
    
    lin1 = plot3([0,n*t_max], [0,t_max],[c,c], 'Color', [0.0,0.0,0.0], ...
                 'LineWidth', lnwidth); hold on;
    
    if (k==1)
        lin2 = plot3([0,2*t_max], [0,t_max],[c,c], '--','Color', ...
                     [0.0,0.0,0.0], 'LineWidth', lnwidth); hold on;
    else
        if (k~=n)
            lin2 = plot3([0,v(m)*t_max], [0,t_max],[c,c], 'Color',...
                         [0.8,0.0,0.0], 'LineWidth', lnwidth); hold on;
        end
    end
    
    % make legends:
        
    title(['P( A(t) = ', num2str(k), ', L(t) < x ), ','n=',num2str(n)]);
    
    lgtext2 = 'x = 5\cdott';
    
    if (k~=1)
        lgtext3 = ['x = ',num2str(k),'\cdott'];
    else
        lgtext3 = ['x = 2\cdott'];
    end
    
    if (k~=n)
        h_leg=legend([lin1,lin2],lgtext2,lgtext3,'location','northwest');
    else
        h_leg=legend([lin1],lgtext2, 'location','northwest');
    end
    
    set(h_leg,'FontSize',fntsz);
    
    legend('boxon')   
        
    xlabel('x');
    
    ylabel('t');
    
    pause(0.2);
      
end

%% 1a.3 Plot coalescent rate lambda(t)

fntsz = 20;

orig_path = path();    % save original path

addpath(genpath(pwd)); % add folders to search path

h_fig=figure('units','normalized','outerposition',[0.2 0.1 0.6 0.9]);

hold on;

set(gca,'FontSize',fntsz);  % set the fontsize of axis

tt=0:0.001:0.2;
    
yy=pop_speed(tt,0,pop_speed_opt);
    
plot(tt,yy, 'LineWidth',2); hold on;

ttxt = ['coalescent rate \lambda(t)'];

title(ttxt);

xlabel('t');

path(old_path); % return original path

%% 1a.4 Plot relative population size c(t)=lambda^{-1}(t)

fntsz = 20;

orig_path = path();    % save original path

addpath(genpath(pwd)); % add folders to search path

h_fig=figure('units','normalized','outerposition',[0.2 0.1 0.6 0.9]);

hold on;

set(gca,'FontSize',fntsz);  % set the fontsize of axis

tt=0:0.001:0.2;
    
yy=pop_speed(tt,0,pop_speed_opt);
    
plot(tt,1./yy, 'LineWidth',2); hold on;

ttxt = ['relative population size c(t)=1/{\lambda(t)}'];

title(ttxt);

xlabel('t');

path(old_path); % return original path