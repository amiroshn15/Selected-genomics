% This function teturns the generator 'Q' for the ancestral process A(t)
% under constant size population. Note, the time-dependent 
% generator for the process A(t; lambda) under variable population size
% that corresponds to the coalescent rate lambda(t) is given by Q*lambda(t)

function Q = A_n_gen(n)

    if (n < 2)
        assert (false, 'n muste be >= 2');
    end
    
    q = ((n-1):(-1):0)' .* (n:(-1):1)' * 0.5; 

    Q=sparse(n,n);
    
    for i=1:(n-1)
        Q(i,i)=-q(i); Q(i,i+1)=q(i);
    end   

end