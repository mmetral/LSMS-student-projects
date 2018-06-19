function [u, spe, acc] = dynamic_equation(u_in, spe_in, acc_in, M_vec, dt, spe_0, n, time, pos_tot, k_initial, delta_max, delta, delta_c, sigma_c, alpha, nodes_live) %u_in = displacement at time - 1
    %predicator
    u=zeros(n+1, 1);
    u = u_in + dt * spe_in + 0.5 * (dt^2) * acc_in;
    u(1,1) = -spe_0 * time;
    u(n+1,1) = spe_0 * time;
    
    spe_p = zeros(n+1, 1);
    spe_p = spe_in + dt*acc_in;
    spe_p(1,1)= -spe_0;
    spe_p(n+1,1)= spe_0;
    
    acc_p = zeros(n+1, 1);
    acc_p = acc_in;
    acc_p(1,1)=0;
    acc_p(n+1,1)=0;
    
    %solve
    f_int = intern_force(pos_tot, n, k_initial, u, delta_max, delta, delta_c, sigma_c, alpha, nodes_live);
    f_ext = external_force(n, acc_p, u, M_vec, k_initial);
    
    acc_term=zeros(n+1, 1);
    for i=1:n+1
        acc_term(i)=M_vec(i)*acc_p(i);
    end
    
    res = f_ext - f_int - acc_term;
    
    da=zeros(n+1, 1);
    for i=1:n+1
        da(i)=res(i)/M_vec(i);
    end
    
    %da = M \ res; %delta acceleration
    
    %corrector
    acc=zeros(n+1, 1);
    acc = acc_p + da;
    acc(1,1) = 0;
    acc(n+1,1) = 0;
    
    spe = zeros(n+1, 1);
    spe = spe_p + 0.5 * dt * da;
    spe(1,1) = -spe_0;
    spe(n+1,1) = spe_0;
end



       
       

    