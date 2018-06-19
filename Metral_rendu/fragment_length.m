function l_frag = fragment_length(nodes_dead, s)
    l_frag=[];
    s_nodes_dead=sort(nodes_dead);
    for i=1:length(s_nodes_dead)
        if i==1 %premier fragment
            l_frag=[l_frag sum(s(1:s_nodes_dead(i)-1))];
        else
            l_frag=[l_frag sum(s(1:s_nodes_dead(i)-1))-sum(s(1:s_nodes_dead(i-1)-1))];
        end
    end
    l_frag=[l_frag sum(s(1:s_nodes_dead(length(s_nodes_dead))))-sum(s(1:s_nodes_dead(length(s_nodes_dead))-1))]; %on rajoute le dernier fragment pour bien en avoir n       
end