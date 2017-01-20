function Q = A_n_gen(n)

    if (n < 2)
        assert (false, 'n muste be >= 2');
    end
        
%     q = (0:1:(n-1))' .* (1:n)' * 0.5 ;

    q = ((n-1):(-1):0)' .* (n:(-1):1)' * 0.5; 

    Q=sparse(n,n);
    
    for i=1:(n-1)
        Q(i,i)=-q(i); Q(i,i+1)=q(i);
    end
    
%     Z1 = [ -q, [0; q(2:end)] ];
% 
%     FJ = spdiags(Z1,[0,1], n, n );
% 
%     Q=sparse(n,n);
%     
%     Q = FJ';

end