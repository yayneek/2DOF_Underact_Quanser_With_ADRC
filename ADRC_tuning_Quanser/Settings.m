clear; clc;
%Ustawienia początkowe dla symulacji:
dt = 0.001;
l = 0.127;
L = l*2;
m1 = 0.065;
m2 = m1*0.5;
g = 9.81;
J = 9.5e-5;
M = diag([m1,m1,J, m1,m1,J, m1,m1,J, m1,m1,J, m2,m2,m2]);
initial_cond.M = M;
initial_cond.dt = dt;
initial_cond.l = l;
initial_cond.L = L;
initial_cond.m1 = m1;
initial_cond.m2 = m2;
initial_cond.g = g;
initial_cond.J = J;
save("Initial_cond.mat","initial_cond");