function T = cohesive_law(delta_max, delta , delta_c, sigma_c, alpha)
    T_max = sigma_c*(1 - delta_max/delta_c);
    
    if delta == delta_max 
        T = T_max;
        
    elseif delta >= 0 && delta < delta_max
        T = T_max*delta/delta_max;
        
    elseif delta < 0
        T = alpha*delta;
        
    else 
        T = 0;
        
    end
end
