% this is the cumulative speed function based on the speed function lambda_1 used in the article

function val = pcexp_cum_speed_inv(y,p)
   
  y_G = (p.N_B/p.g)*(exp(p.T_G*p.g)-1);

  y_B = (p.N_B/p.g)*(exp(p.T_G*p.g)-1) + p.N_B*(p.T_B-p.T_G);
  
  if (y >= y_B)
            
      val =  (y-y_B)/p.N_A + p.T_B;
   
    
  elseif (y >= y_G)       
      
      val =  (y-y_G)/p.N_B + p.T_G;
      
  else 
      
      val =  -(1/p.g)*log(1-y*p.g*exp(-p.T_G*p.g)/p.N_B);
      
  end
  
end
