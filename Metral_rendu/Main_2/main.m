%% this code is similar to the first one : the differences are we compute only dissipated energy, and we use a mass vector and no a mass matrix
clear all
close all
clc

%% geometric variables

l=0.5; %bar length [m]
A=0.05; %area [m2]

%% variables for finite elements

nb_frag=[1];
n=100; %number of springs

s = length_spring(n, l); %length of springs [m]
% s = zeros(n,1);
% s(:) = l/n;

%% material parameters

E=310e9; %elastic modulus [Pa]
nu=0.27; %poisson's ratio [-]
rho=3300; %mass density [kg/m3]

epsilon_rate_0=5*10^2; %[1/s]

alpha=0; %3.815e15; %[Pa/m] : penalty for contact
G_c=30; %dissipated energy [J/m2]

%sigma_c = 300e6*ones(n+1, 1); %same sigma_c for each spring
sigma_c=299*10^6+2*10^6.*rand(n+1, 1); %uniform distributed critical stress [Pa] 300+-1MPa
%sigma_c((n+2)/2)=300*10^6; %only the node in the middle can break

delta_c=2*G_c./sigma_c; %critical strain for each spring [m]

%% vector for finite element

k_initial=(E*A./s);

M_vec=zeros(n+1, 1);
for i=2:n
    M_vec(i)=rho*A*0.5*(s(i)+s(i-1));
end
M_vec(1)=rho*A*0.5*s(1);
M_vec(n+1)=rho*A*0.5*s(n);
%M=diag(M_vec); %mass vector

%% time

c = sqrt(E/rho); %celerity [m/s]
t_final = 5*10^-6; %[s]
c_dt= 0.1; %coefficient for the time step
dt = (c_dt/c)*min(s); %time step [s]
T = 0:dt:t_final; %time vector []

%% vector to caracterize nodes

nodes_all = 1:n+1;
nodes_live=2:n; %vector contains nodes that can be still broke
nodes_dead=[];
not_broken=nodes_live;

delta_max=zeros(length(nodes_live), 1); %maximum delta for node which can be broke
delta=zeros(length(nodes_live), 1);

pos=[]; %position of nodes broken (initial position)
pos_tot=[]; %vector with all the position where nodes were broken
pos_rup=[];

%% initial condition

spe_0=epsilon_rate_0*1; %[m/s] vo=epsilon_ref*r

u_in=zeros(n+1, 1);

grad_spe=@(i) (2*spe_0/l)*sum(s(1:i-1))-spe_0;

spe_in=zeros(n+1, 1);
for i=1:n+1
    if i==1
        spe_in(i)=-spe_0;
    else
        spe_in(i)=grad_spe(i);
    end
end

acc_in=zeros(n+1, 1);

f_ext = zeros(n+1, 1);

f_int = zeros(n+1, 1);

%% connectivity 

con = zeros(n, 2); %connectivity matrix
for e = 1:n
    con(e, 1) = e;
    con(e, 2) = e + 1;
end


%% dissipated energy
E_dis=[0];

%% main
for t=2:length(T)
    
    time=T(t);
    
    %displacement/speed/acceleration for a fixed time
    [u, spe, acc] = dynamic_equation(u_in, spe_in, acc_in, M_vec, dt, spe_0, n, time, pos_tot, k_initial, delta_max, delta, delta_c, sigma_c, alpha, nodes_live);
    
    %force at time
    f_ext = external_force(n, acc, u, M_vec, k_initial);
    f_int = intern_force(pos_tot, n, k_initial, u, delta_max, delta, delta_c, sigma_c, alpha, nodes_live);
   
    e_dis = energy_dissipated(delta_max, sigma_c, A);
    E_dis = [E_dis e_dis];
    
    
    %rupture ?
    [pos_rup , pos]= rupture(sigma_c, nodes_live, pos_tot, A, k_initial, u);
    
    %we add node in pos_rup : pos_rup is define with all nodes
    if length(pos_rup)>1
        for j=2:length(pos_rup)
            pos_rup(j)=pos_rup(j)+length(pos_rup(1:j-1));
        end
    end
        
    %define total rupture positions concidering the new ones
    Index0=[];
    
    if ~isempty(pos)
        if ~isempty(pos_tot)
                for i=1:length(pos_tot)
                    index=find(pos==pos_tot(i));
                    Index0=[Index0 index];
                end
        end
        add=pos;
        add(Index0)=[];    
        pos_tot=[pos_tot add];
        pos_tot=sort(pos_tot);
    end
    
    %define nodes position which can be still broken
    Index1=[];
    if ~isempty(pos)
            for i=1:length(pos)
                index=find(nodes_live==pos(i));
                Index1=[Index1 index];
            end
    end
    
    nodes_live(Index1)=[];
    
    if ~isempty(pos_rup)
            
            for j=1:length(pos_rup)
                
                i=pos_rup(j);
                n=n+1;
             
                %mass matrix
                M_vec(i+2:n+1)=M_vec(i+1:end);
                M_vec(i)=M_vec(i)/2;
                M_vec(i+1)=M_vec(i);
                %M=diag(M_vec);
                
                %new connectivity matrix
                con = zeros(n, 2); %connectivity matrix
                for e = 1:n
                    con(e, 1) = e;
                    con(e, 2) = e + 1;
                end   
            end
            
             %initial value
            [u_in, spe_in, acc_in]=initial_vec(n, pos_rup, u, acc, spe);
    else
            u_in = u;
            spe_in = spe;
            acc_in = acc;
    end
    
    %define the displacement for cohesive springs
    if ~isempty(pos_tot)
        for i = pos_tot
            j = index_node(i, pos_tot);
            delta(i-1)=u_in(j+1)-u_in(j);  
        end
    end
    
    %delta max
    if ~isempty(not_broken)
       for i=pos_tot
            j = index_node(i, pos_tot);
            delta_max(i-1) = max(delta_max(i-1), delta(i-1));
            delta_max(i-1) = min(delta_max(i-1), delta_c(i));
       end
    end
   
   
   %number of fragments
   nb = nb_frag(end);
   
   node_supp = [];
   if ~isempty(not_broken)
       for j = 1:length(not_broken)
           i = not_broken(j);
           if delta(i-1) >= delta_c(i)
               nb = nb+1;
               nodes_dead = [nodes_dead i];
               node_supp = [node_supp j];
           end
       end
   end
   not_broken(node_supp)=[];
   nb_frag=[nb_frag nb];
   
end

energy_dissipated = E_dis(end)

%% this section allows us to calculate the average lenght and the fragment size distribution 

%% fragements length
l_frag = fragment_length(nodes_dead, s);

%% arithmetic average 
l_moy_arit = sum(l_frag)/length(l_frag)

%% weigth average
%on condière que les plus gros morceaux ont plus d'importance et donc on utilise des poids
l_max=max(l_frag);
poids=l_frag./l_max;

l_moy_pond=sum(poids.*l_frag)/sum(poids)

%% fragment size distribution
r_length=l_frag./l_moy_arit;
[counts,centers] = hist(r_length);
bar(centers,counts)
 
   
   