function f=rand_choose(W)

if ( sum(W < 0) < 1 )

    cumW = cumsum(W);

    Tot_W = cumW(end);
    
    U = rand;
    
    f = find( ( cumW/Tot_W ) > U , 1,'first');    

else
    assert (false, 'negative weight given');
end

end