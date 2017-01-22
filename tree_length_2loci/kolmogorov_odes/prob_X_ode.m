function [T,Y] = prob_X_ode(rho, pop_speed_option, Qc,Qr, t_range, y_0, ode_options )

    % some anonimous useful functions:
     
    if (rho < 0)
        assert (false, 'rho muste be >= 0.0');
    end

    % anonimous population speed function
    
    lambda = @(t) pop_speed(t, 0, pop_speed_option); % population speed 

    % use Qc, Qr   that decompose the generator Q = lambda(t) * Qc + rho*Qr

    QcT   = Qc';
    
    QrT = rho * Qr';
    
    % setup the ode solver 

    rhs = @(t,y) ( lambda(t) * QcT * y + QrT * y );       
    
    
    if ~isfield(ode_options,'solver')
        
        ode_options.solver = @ode45;
        
    end
    
    [T,Y] = ode_options.solver(rhs, t_range, y_0, ode_options);
    
end