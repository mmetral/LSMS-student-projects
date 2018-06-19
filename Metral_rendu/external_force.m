function f_ext = external_force(n, acc, u, M, k_initial)
    f_ext = zeros(n+1, 1);

    f_ext(1) = M(1,1) * acc(1) + k_initial(1) * (u(1) - u(2));
    f_ext(n+1) = M(n+1, n+1) * acc(n+1) + k_initial(end) * (u(n+1) - u(n));

end

    