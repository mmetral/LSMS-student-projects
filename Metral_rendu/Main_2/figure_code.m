%% figure : number of fragment
figure 
plot(T(1:12325), nb_frag, '-')
hold on
xlabel('time [s]')
ylabel('number of fragments')
grid minor
title('number of fragments over time')

%% figure : plot of stress
figure
plot(T, Sigma, '-')
hold on
xlabel({'time [s]'}, 'FontSize', 14)
ylabel({'stress [Pa]'}, 'FontSize', 14)
grid minor

%% figure : energy without contact

figure 
plot(T, E_ext, 'b')
hold on
plot(T, E_kin, 'r')
hold on
plot(T, E_pot, 'g')
hold on
plot(T, E_dis, 'c')
hold on
plot(T, E_rev, 'y')
hold on
plot(T, E_tot, 'k')
legend({'external energy', 'kinetic energy', 'potential energy (spring)', 'dissipated energy', 'reversible energy', 'total energy'}, 'FontSize', 10, 'Location', 'eastoutside'  )
grid minor
xlabel({'time [s]'}, 'FontSize', 14)
ylabel({'energy terms [J]'}, 'FontSize', 14)

%% displacement
figure
grid minor

for t=1:length(T)
    n=N(t);
    plot(X(1:n+1, t), U(1:n+1 , t))
    ylim([min(min(U)) max(max(U))]);
    drawnow
end

%% speed
figure
grid minor

for t=1:length(T)
    n=N(t);
    plot(X(1:n+1, t), Spe(1:n+1 , t))
    ylim([min(min(Spe)) max(max(Spe))]);
    drawnow
end

%% acceleration
figure
grid minor

for t=1:length(T)
    n=N(t);
    plot(X(1:n+1, t), Acc(1:n+1 , t))
    ylim([min(min(Acc)) max(max(Acc))]);
    drawnow
end

%% test
figure 
plot(T, E_ext)
grid minor