function e_con = energy_contact(alpha, delta)
    e_con=0;
    for i=1:length(delta)
        if delta(i)<0
            e_con=e_con+alpha*delta(i);
        end
    end
end

            