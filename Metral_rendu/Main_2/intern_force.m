function f_int = intern_force(pos_tot, n, k_initial, u, delta_max, delta, delta_c, sigma_c, alpha, nodes_live) 
       f_int = zeros(n+1, 1);
       
       %sur les extrémités, on a que K*u (pas d'introduction avec la force 
       f_int(1) = k_initial(1)*(u(1) - u(2));
       f_int(n+1) = k_initial(end) * (u(n+1) - u(n));
       
       %tous les noeuds qui sont dans pos_tot ont leur force interne donnée
       %par la loi de cohésion + Ku
       if ~isempty(pos_tot)
           for i = pos_tot
               j = index_node(i, pos_tot);
               T = cohesive_law(delta_max(i-1), delta(i-1) , delta_c(i), sigma_c(i), alpha);
               f_int(j) = k_initial(i-1)*(u(j)-u(j-1)) - T;
               f_int(j+1) = k_initial(i)*(u(j+1)-u(j+2)) + T;
           end
       end
       
       %tous les noeuds en dehors des extrémités qui ne sont pas cassés :
       %Ku
       if ~isempty(nodes_live)
           for i = nodes_live
               j = index_node(i, pos_tot);
               f_int(j) = k_initial(i-1)*(u(j)-u(j-1)) + k_initial(i)*(u(j)-u(j+1));
           end
       end
end

       
       