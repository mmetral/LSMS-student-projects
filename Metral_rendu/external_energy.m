function e_ext = external_energy(f_ext, E_ext, U, t, n, M, Acc, k_initial)
    
    f_ext(1) = M(1,1) * Acc(1, t) + k_initial(1) * (U(1, t) - U(2, t));
    f_ext(n+1) = M(1,1) * Acc(n+1, t)+ k_initial(end) * (U(n+1, t) - U(n, t));
    
    e_ext = E_ext(end)- f_ext(1)*(U(1, t)-U(1, t-1)) - f_ext(n+1)*(U(n+1, t) - U(n+1, t-1));
end