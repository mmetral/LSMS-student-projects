function e_dis = energy_dissipated(delta_max, sigma_c, A)
        e_dis = 0;
        for i=1:length(delta_max)
            e_dis = e_dis + 0.5*delta_max(i)*sigma_c(i+1)*A;
        end
end