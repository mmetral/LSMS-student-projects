%% this section define all parameters, and compute displacement U, speed Spe, acceleration Acc for each time (columns) and nodes (rows) 
% + Stress along the bar for each time + Energy terms

clear all
close all
clc

%% geometric variables

l=0.5; %bar length [m]
A=0.05; %area [m2]

%% variables for finite elements

nb_frag=[1]; %number of fragments
n=50; %number of springs

%s = length_spring(n, l); %springs length  : random mesh [m]
s = zeros(n,1);
s(:) = l/n; %springs length  : uniform mesh [m]

%% material parameters

E=310e9; %elastic modulus [Pa]
nu=0.27; %poisson's ratio [-]
rho=3300; %mass density [kg/m3]

epsilon_rate_0=5*10^2; % strain rate [1/s]

alpha=0; %[Pa/m] : penalty coefficient for contact
G_c=30; %dissipated energy [J/m2]

%sigma_c = 300e6*ones(n+1, 1); %critical stress : uniform
sigma_c=299*10^6+2*10^6.*rand(n+1, 1); %critical stress [Pa] 300+-1MPa : uniform distributed 
%sigma_c((n+2)/2)=300*10^6; %only the node in the middle can break

delta_c=2*G_c./sigma_c; %critical strain for each spring [m]

%% vector for finite element

%initial stiffness
k_initial=(E*A./s);

%mass matrix
M_vec=zeros(n+1, 1);
for i=2:n
    M_vec(i)=rho*A*0.5*(s(i)+s(i-1));
end
M_vec(1)=rho*A*0.5*s(1);
M_vec(n+1)=rho*A*0.5*s(n);
M=diag(M_vec); 

%% time

c = sqrt(E/rho); %celerity [m/s]
t_final = 1*10^-5; %[s] %final time
c_dt= 0.001; %coefficient for the time step
dt = (c_dt/c)*min(s); %time step [s]
T = 0:dt:t_final; %time vector []

%% vector to caracterize nodes

nodes_all = 1:n+1; 
nodes_live=2:n; %vector contains nodes that can be still broke
nodes_dead=[]; %vector contains nodes that is broken
not_broken=nodes_live; %vector contains nodes which have have crack but are not broken

delta_max=zeros(length(nodes_live), 1); %maximum delta for node which can be broke
delta=zeros(length(nodes_live), 1);

pos=[]; %position of nodes broken (initial position)
pos_tot=[]; %vector with all the position where nodes were broken
pos_rup=[]; %position of nodes broken (concidering nodes which I had)

%% initial condition

spe_0=epsilon_rate_0*1; %[m/s] speed applied

u_in=zeros(n+1, 1); %initial displacement

grad_spe=@(i) (2*spe_0/l)*sum(s(1:i-1))-spe_0; %initial gradient speed

spe_in=zeros(n+1, 1); %initial speed
for i=1:n+1
    if i==1
        spe_in(i)=-spe_0;
    else
        spe_in(i)=grad_spe(i);
    end
end

acc_in=zeros(n+1, 1); %initial acceleration

f_ext = zeros(n+1, 1); %initial external force
f_int = zeros(n+1, 1); %initial intern force

%% connectivity 

con = zeros(n, 2); %connectivity matrix
for e = 1:n
    con(e, 1) = e;
    con(e, 2) = e + 1;
end


%% vector for plot
%energy
E_pot=[0]; %potentiel
E_kin_start=[0.5*spe_in'*M*spe_in]; %initial kinetic energy 
E_kin=[0]; %variation of kinetic energy 
E_dis=[0]; %dissipated energy
E_ext=[0]; %external energy
E_con=[0]; %contact energy : throughtout this section, no contact is considered
E_rev=[0]; %reversible energy
E_tot=[0]; %total energy

%stress
Sigma=[0]; %stress

%plot total displacement/speed/acceleration
U=zeros(2*length(nodes_live)+2, length(T));
Spe=zeros(2*length(nodes_live)+2, length(T));
Acc=zeros(2*length(nodes_live)+2, length(T));

X=zeros(2*length(nodes_live)+2, length(T));
n_in=n;
X(1:n_in+1, 1)=1:n_in+1;
N=zeros(1, length(T));
N(1)=n_in;

%% main
for t=2:length(T)
    
    time=T(t);
    
    %displacement/speed/acceleration for a fixed time
    [u, spe, acc] = dynamic_equation(u_in, spe_in, acc_in, M, dt, spe_0, n, time, pos_tot, k_initial, delta_max, delta, delta_c, sigma_c, alpha, nodes_live);
    
    %force at time
    f_ext = external_force(n, acc, u, M, k_initial);
    f_int = intern_force(pos_tot, n, k_initial, u, delta_max, delta, delta_c, sigma_c, alpha, nodes_live);
    
    %sigma 
    Sigma=[Sigma f_ext(end)/A]; %1D problem
    
    % plot : displacement/speed/acceleration
    U(1:n+1, t) = u;
    Spe(1:n+1, t) = spe;
    Acc(1:n+1, t) = acc;
    
    if ~isempty(pos_tot)
        vec=sort([1:n_in+1 pos_tot]);
        X(1:n+1, t)=vec;
    else
        X(1:n+1, t)=1:n_in+1;
    end
    N(t)=n;
    
    
    %energy
    
    e_ext = external_energy(f_ext, E_ext, U, t, n, M, Acc, k_initial);
    E_ext=[E_ext e_ext];
    
    e_kin = kinetic_energy(M, spe, E_kin_start);
    E_kin = [E_kin e_kin];
    
    e_pot = potential_energy(f_int, u);
    E_pot = [E_pot e_pot];
    
    e_dis = energy_dissipated(delta_max, sigma_c, A);
    E_dis = [E_dis e_dis];
    
    e_con = energy_contact(alpha, delta);
    E_con = [E_con e_con];
    
    e_rev = reversible_energy(delta_max, sigma_c, delta_c, delta, alpha);
    E_rev=[E_rev e_rev];
    
    e_tot=e_ext+e_kin+e_pot+e_dis+e_con+e_rev;
    E_tot=[E_tot e_tot];
  
    
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
    
    %define new vectors when we add nodes
    if ~isempty(pos_rup)
            
            for j=1:length(pos_rup)
                
                i=pos_rup(j);
                n=n+1;
             
                %mass matrix
                M_vec(i+2:n+1)=M_vec(i+1:end);
                M_vec(i)=M(i, i)/2;
                M_vec(i+1)=M(i, i)/2;
                M=diag(M_vec);
                
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

 
   
   