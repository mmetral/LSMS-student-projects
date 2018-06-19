%% this section allows us to calculate the average lenght and the fragment size distribution 

%% fragements length
l_frag = fragment_length(nodes_dead, s);

%% arithmetic average 
l_moy_arit = sum(l_frag)/length(l_frag);

%% weigth average
%on condière que les plus gros morceaux ont plus d'importance et donc on utilise des poids
l_max=max(l_frag);
poids=l_frag./l_max;

l_moy_pond=sum(poids.*l_frag)/sum(poids);

%% fragment size distribution
r_length=l_frag./l_moy_arit;
[counts,centers] = hist(r_length);
bar(centers,counts)

