function e_rev = reversible_energy(delta_max, sigma_c, delta_c, delta, alpha)
    e_rev = 0;
    for i=1:length(delta_max)
        T = cohesive_law(delta_max(i), delta(i) , delta_c(i+1), sigma_c(i), alpha);
        e_rev = e_rev + 0.5*T*delta(i); 
    end
end
    