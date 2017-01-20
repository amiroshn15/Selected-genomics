function FF = prob_X_L_1p_upt_decr_interp_t( F, x, t, P_X, N_t_int, vv_inv_M, xx,tt, interp_option)

if ~exist('interp_option')
    interp_option = 'linear';
end

xx_max = xx(end);

tt_max = tt(end);

t_max = t(end);

x_max = x{N_t_int}(end);

N_t = numel(t);

N_tt = numel(tt);

N_xx = numel(xx);


if (N_t_int<N_t) && (xx_max>x_max)
    assert(false,strcat('pde solution is unknown for x>',num2str(x_max),'. Change values of xx-mesh.'))
end

if (tt_max>t_max)
    assert(false,strcat('pde solution is unknown for t>',num2str(t_max),'. Change values of tt-mesh.'))
end

K=size(P_X,2);

FF = zeros( N_tt, N_xx, K );

vv_inv_M_complement = setdiff(1:K,vv_inv_M);

PP_X_v_inv_M = (interp1( t, P_X(:,vv_inv_M), tt, interp_option,'extrap'));

M = x{N_t_int}(end)/t(N_t_int);


for idx_m=1:numel(vv_inv_M)
    
    m=vv_inv_M(idx_m);
    
    for r=1:N_tt
        off_diag_idx = find(xx > M*tt(r));
        
        if ~isempty(off_diag_idx)

            if numel(vv_inv_M)==1
                FF(r,off_diag_idx,m) = PP_X_v_inv_M(r);
            else
                FF(r,off_diag_idx,m) = PP_X_v_inv_M(r,idx_m);
            end
        end
    end      
end

% %

prob_x_unif = zeros(N_t,N_xx,K);

for idx_m=1:numel(vv_inv_M_complement)
    
    m=vv_inv_M_complement(idx_m);
    
    for i=1:N_t
        if (xx(end)>x{i}(end))            
            prob_x_unif(i,:,m) = interp1([x{i}, xx(end)], [F{i}(m,:),F{i}(m,end)], xx, interp_option,'extrap');
        else            
            prob_x_unif(i,:,m) = ( interp1( x{i}, F{i}(m,:), xx , interp_option,'extrap') );
        end
    end
   
    for j=1:N_xx
        FF(:,j,m) = interp1(t, prob_x_unif(:,j,m), tt ,interp_option,'extrap');
    end
    
end

end