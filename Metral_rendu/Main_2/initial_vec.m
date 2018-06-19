function [u_in, spe_in, acc_in]=initial_vec(n, pos_rup, u, acc, spe)
    u_in=u;
    spe_in=spe;
    acc_in=acc;
    
    n_adap=n-length(pos_rup);
    
    for j=1:length(pos_rup)
        i=pos_rup(j);
        
        u_in(1:i)=u_in(1:i); 
        u_in(i+2:n_adap+2)=u_in(i+1:end);
        u_in(i+1)=u_in(i);

        spe_in(1:i)=spe_in(1:i); 
        spe_in(i+2:n_adap+2)=spe_in(i+1:end);
        spe_in(i+1)=spe_in(i);
        
        acc_in(1:i)=acc_in(1:i); 
        acc_in(i+2:n_adap+2)=acc_in(i+1:end);
        acc_in(i+1)=acc_in(i);
        
        n_adap=n_adap+1;
        
    end
        
end