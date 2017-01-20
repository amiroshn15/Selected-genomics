% speed function called lambda_1 in the article

function val = pcexp_speed(t,p)
   
  if (t >= p.T_B)
      
    val =  (p.N_A);  
    
  elseif (t >= p.T_G) 
      
      val =  (p.N_B);    
      
  else 
      
      val =  p.N_B * exp((t-p.T_G)*(-p.g));    
      
  end
  
end
