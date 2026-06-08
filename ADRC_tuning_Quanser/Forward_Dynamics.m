clear; clc;
%Wczytanie wartości do symulacji
settings = load("Initial_cond.mat").initial_cond;
dt = settings.dt;
g = settings.g;
l = settings.l;
L = settings.L;
m1 = settings.m1;
m2 = settings.m2;
J = settings.J;
M = diag([m1,m1,J,m1,m1,J,m1,m1,J,m1,m1,J,m2,m2,m2]);
singals = load('output.mat').output;
u = singals.u;
traj = singals.traj;

%Definicja czasu:
t0 = 0;
tk = 10;
t = t0:dt:tk;

%Szukane wektory
q = zeros(15,length(t));
v = zeros(15,length(t));
lagrange = zeros(11, length(t));

%Warunki początkowe:
q0 = [0; -l/2; 3/2*pi; l/2; -l; 0; l; -l/2; pi/2; 3/2*l;0;0;-L; l;-l];
v0 = zeros(15,1);
P0 = [l;-l]; %Położenie początkowe ładunku
q(:,1) = q0;
v(:,1) = v0;

%Sprawdzenie, czy warunki początkowe spełniają więzy:
if norm(Phi(q0,l,L)) > 1e-10
    Phi(q0,l,L)
    error('Warunki początkowe nie spełniają równania więzów!')
end
initial_guess = [q0;v0];
options = odeset('RelTol',1e-12,'AbsTol',1e-16);

[t, sol] = ode45(@(t, y) ode_system(t, y, u, M, l, L, t0, tk), t, initial_guess, options);

%Wyniki:
q = sol(:, 1:15)';
v = sol(:, 16:30)';

figure()
plot(t, q(14,:), t, q(15,:))
title('xy(t)')
grid on
legend('x','y')

figure()
plot(t, q(14,:) - traj(1,:), t, q(15,:) - traj(2,:));
title('e(t)')
legend('ex','ey')
grid on

function dydt = ode_system(t, y, u, M, l, L, t0, tk)

   
    q = y(1:15);
    v = y(16:30);

    m2 = M(13,13);
    g = 9.81;

    %Interpolacja sterowania u(t)
    u1 = interp1(linspace(t0, tk, length(u(1,:))), u(1,:), t);
    u2 = interp1(linspace(t0, tk, length(u(2,:))), u(2,:), t);


    F = [0;0;u1; 0;0;0; 0;0;0; 0;0;u2; -m2*g;0;0];
    C = Phiq(q, l);
    wiezy = Phi(q, l, L);

    alpha = 100;
    beta = 100;

    A = [M C'; C zeros(11,11)];
    b = [F; Gamma(q, v, l) - 2*alpha*C*v - beta^2*wiezy];
    x = A\b;
    
    a = x(1:15);
    lagrange = x(16:end);
    dydt = [v; a];
end


%% Obliczenie mnożników lagrange'a i położenia punktu charakterystycznego
function lambda = compute_lagrange(q, v, t, u, M, l, L, t0, tk)
    % Tutaj powtarzamy (prawie) wszystko to, co robimy w ode_system:
    %
    %  1. Obliczenie u1, u2 na chwili t
    %  2. F = [...], C = ..., wiezy = ...
    %  3. Macierz A = [M C'; C 0]
    %  4. Wektor b = [F; Gamma - 2*alpha*C*v - beta^2*wiezy]
    %  5. Rozwiązanie x = A\b
    %  6. lambda = x(16:end);

    m2 = M(13,13);
    g = 9.81;
    F = [0;0;u(1); 0;0;0; 0;0;0; 0;0;u(2); -m2*g;0;0];

    C = Phiq(q, l);
    wiezy = Phi(q, l, L);
    
    alpha = 100;
    beta = 100;

    A = [M, C'; C, zeros(11, 11)];

    size([M, C'])

    b = [ F; Gamma(q, v, l) - 2*alpha*C*v - beta^2*wiezy ];
    size(F)
    x = A\b;
    % Ostatnie 11 składowych to mnożniki:
    lambda = x(16:end);
end


lambda_all = zeros(11, length(t)); % jeśli chcemy mieć 11 mnożników w czasie

for i = 1:length(t)
    lambda_all(:,i) = compute_lagrange(q(:,i), v(:,i), t, u(:,i), M, l, L, t0, tk);
end

Px = zeros(1, length(t));
Py = zeros(1, length(t));
for i = 1:length(t)
    Px(i) = (cos(q(3,i)) * l + cos(q(6,i)) * l);
    Py(i) = (sin(q(3,i)) * l + sin(q(6,i)) * l);
end

%%
folder = 'C:\Users\janek\Desktop\Publikacja\Kod\Wykresy';
Case = 'line';
close all

figure()
plot(t, u)
title('$u(t)$','Interpreter','latex', 'FontSize',16)
grid on
legend('$u_1$','$u_2$','Interpreter','latex', 'FontSize',16, 'location','northwest')
ylabel('$u[Nm]$','Interpreter','latex', 'FontSize',16)
xlabel('$t[s]$','Interpreter','latex', 'FontSize',16)

% filename = 'u(t)_';
% fullfile_pdf = fullfile(folder, [append(filename, Case) '.pdf']);
% exportgraphics(gcf, fullfile_pdf);



figure()
plot(t, q(14,:) - traj(1,:), t, q(15,:) - traj(2,:));
title('$e(t)$','Interpreter','latex', 'FontSize',16)
grid on
legend('$e_x$','$e_y$','Interpreter','latex', 'FontSize',16, 'location','northwest')
ylabel('$e[m]$','Interpreter','latex', 'FontSize',16)
xlabel('$t[s]$','Interpreter','latex', 'FontSize',16)

% filename = 'e(t)_';
% fullfile_pdf = fullfile(folder, [append(filename, Case) '.pdf']);
% exportgraphics(gcf, fullfile_pdf);



figure()
plot(t, lambda_all(11,:))
grid on
title('$\lambda (t) -$ Tension force in the rope','Interpreter','latex', 'FontSize',16)
ylabel('$\lambda[N]$','Interpreter','latex', 'FontSize',16)
xlabel('$t[s]$','Interpreter','latex', 'FontSize',16)

% filename = 'lag(t)_';
% fullfile_pdf = fullfile(folder, [append(filename, Case) '.pdf']);
% exportgraphics(gcf, fullfile_pdf);


figure()
plot(q(14,:), q(15,:), Px, Py)
grid on
title('Location of the load','Interpreter','latex', 'FontSize',16)
xlabel('$x[m]$','Interpreter','latex', 'FontSize',16)
ylabel('$y[m]$','Interpreter','latex', 'FontSize',16)
legend('Load','Characteristic point','Interpreter','latex', 'FontSize',16)

% filename = 'XY_';
% fullfile_pdf = fullfile(folder, [append(filename, Case) '.pdf']);
% exportgraphics(gcf, fullfile_pdf);


figure()
plot(t, q(13,:))
grid on
title('$z(t)$','Interpreter','latex', 'FontSize',16)
xlabel('$t[s]$','Interpreter','latex', 'FontSize',16)
ylabel('$z[m]$','Interpreter','latex', 'FontSize',16)

% filename = 'z(t)_';
% fullfile_pdf = fullfile(folder, [append(filename, Case) '.pdf']);
% exportgraphics(gcf, fullfile_pdf);

figure()
plot(traj(1,:), traj(2,:), q(14,:), q(15,:))
grid on
title('Desired and output trajectories','Interpreter','latex', 'FontSize',16)
xlabel('$x[m]$','Interpreter','latex', 'FontSize',16)
ylabel('$y[m]$','Interpreter','latex', 'FontSize',16)
legend('Desired trajectory','Output of the system','Interpreter','latex', 'FontSize',16)
% 
% filename = 'traj_';
% fullfile_pdf = fullfile(folder, [append(filename, Case) '.pdf']);
% exportgraphics(gcf, fullfile_pdf);

%% Animacja 3D

% Animacja3d(q,l,L,t,dt,traj);