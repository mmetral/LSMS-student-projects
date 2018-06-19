function j = index_node(i, pos_tot) 
% cette fonction sert � connaitre la position d'un noeud apr�s rupture
% connaissant sa position initiale
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