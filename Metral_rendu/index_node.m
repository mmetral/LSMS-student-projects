function j = index_node(i, pos_tot) 
% this fonction allows us to know the position of nodes j when we
% concider all nodes added, knowing the initial position i (without
% concidering the nodes added)
    if ~isempty(pos_tot)
        ind=find(pos_tot<i);
        if isempty(ind)
            j=i;
        else
            j=i+length(ind);
        end
    else
        j=i;
    end
end