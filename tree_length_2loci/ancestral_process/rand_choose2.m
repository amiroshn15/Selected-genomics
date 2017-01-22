function f=rand_choose2(W,n,m)

if ( nargin < 1 )
    assert (false, 'vector of weights is not provided');
end

if ( nargin == 1 )
    n=1;
    m=1;
end

if ( nargin == 2 )
    m=1;
end

if ( ( sum(W < 0) < 1 ) && (n>0) && (m>0) )

    f = zeros(n,m);
    
    cumW = cumsum(W);

    TotW = cumW(end);
    
    CumProbW = cumW/TotW;

    U = rand(n,m);

    for i=1:n
        for j=1:m
            f(i,j) =  find((CumProbW > U(i,j)), 1, 'first');
        end
    end
    
else
    assert (false, 'negative weight given');
end


end