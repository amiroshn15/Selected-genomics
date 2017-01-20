% this is the cumulative speed function based on the speed function lambda_1 used in the article

function val = pcexp_cum_speed(t,p)
  
  if (t >= p.T_B)
      
    val =  (p.N_B/p.g) * (exp(p.T_G*p.g)-1) + p.N_B*(p.T_B-p.T_G)+p.N_A*(t-p.T_B);
    
  elseif (t >= p.T_G) 
      
    val =  (p.N_B/p.g) * (exp(p.T_G*p.g)-1) + p.N_B*(t-p.T_G);
      
  else 
      
    val =  (p.N_B/p.g) * exp(p.T_G*p.g)*(1-exp(-t*p.g));
      
  end
  
end