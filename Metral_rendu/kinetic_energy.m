function e_kin = kinetic_energy(M, spe, E_kin_start)
    e_kin = 0.5*spe'*M*spe - E_kin_start(1);
end