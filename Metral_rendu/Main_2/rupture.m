function [pos_rup , pos]= rupture(sigma_c, nodes_live, pos_tot, A, k_initial, u)
    pos_rup = [];
    pos=[];
    
    if ~isempty(nodes_live)
        for i=nodes_live %exprimé dans la base initiale
            j = index_node(i, pos_tot);
            
            f_spring_av = 0.5*(k_initial(i-1)*(u(j)-u(j-1)) + k_initial(i)*(u(j+1)-u(j)));
            
            if f_spring_av >= sigma_c(i)*A
               pos_rup=[pos_rup j];
               
               pos=[pos i];
               pos=sort(pos);
            end  
        end
    end
end