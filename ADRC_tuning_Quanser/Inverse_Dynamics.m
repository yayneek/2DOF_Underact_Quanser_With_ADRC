clear; clc;
close all;
%Wczytanie wartości do symulacji
settings = load("Initial_cond.mat");
settings = settings.initial_cond;
dt = settings.dt;
l = settings.l;
L = settings.L;
m1 = settings.m1;
m2 = settings.m2;
J = settings.J;
M = settings.M;
g = 9.81;

%Definicja czasu:
t0 = 0;
tk = 10;
t = t0:dt:tk;

%Szukane wektory
q = zeros(15,length(t));
v = zeros(15,length(t));
u = zeros(2,length(t));
lagrange = zeros(11, length(t));
lagrange(11,1) = m2*g;

%Warunki początkowe:
q0 = [0; -l/2; 3/2*pi; l/2; -l; 0; l; -l/2; pi/2; 3/2*l;0;0;-L; l;-l];
v0 = zeros(15,1);
P0 = [l;-l];        %Położenie początkowe ładunku
Pk = [l/2; -l/2];   % Położenie końcowe ładunku 
q(:,1) = q0;
v(:,1) = v0;


%Sprawdzenie, czy warunki początkowe spełniają więzy:
if norm(Phi(q0,l,L)) > 1e-10
    Phi(q0,l,L)
    error('Warunki początkowe nie spełniają równania więzów!')
end

%Trajektoria ładunku:
[traj, traj_d, traj_dd] = Trajectory_Circle(t0,dt,tk,P0, l/5, 5);
X = Trajectory_Trapz(t0,dt,tk,P0(1),Pk(1),5);
Y = Trajectory(t0,dt,tk,P0(2),Pk(2),10);

%Rozwiązanie iteracyjne:
initial_guess = [q0;v0;u(:,1);lagrange(:,1)];
options = optimset('Display','none', 'TolFun', 1e-10, 'TolX', 1e-7, 'MaxIter', 10000, 'MaxFunEvals', 50000);

N = 1;

for i=1:length(t)-1
func = @(unknowns)ID(unknowns, q, v, i, dt, l, L, traj, M, N);
[solution, fval(:,i), exitflag(i)] = fsolve(func, initial_guess, options);
initial_guess = solution;

q(:,i+1) = solution(1:15);
v(:,i+1) = solution(16:30);
u(:,i+1) = solution(31:32);
lagrange(:,i+1) = solution(33:end);
fvalnorm(i) = norm(fval(:,i));
if exitflag(i) < 1
    break;
end
end


figure()
plot(t, q(14,:), t, q(15,:));
grid on
title('xy')
legend('x','y')

figure()
plot(t, u)
title('u(t)')
grid on

output.u = u;
output.traj = traj;
save('output.mat',"output");

u1 = u(1,:);
u2 = u(2,:);
x_traj = traj(1,:);
y_traj = traj(2,:);
b_hat = 41874;
omega_c = 1;

function solution = ID(unknowns, q, v, i, dt, l, L, traj, M, N)
Q = unknowns(1:15);
V = unknowns(16:30);
u = unknowns(31:32);
lagrange = unknowns(33:end);


%Macierz sygnałów sterujących
B = zeros(15,2);
B(3,1) = 1;
B(12,2) = 1;

%Siły przyłożone do mechanizmu:
g = 9.81;
f = zeros(15,1);
f(13,1) = -M(15,15)*g;

%Tłumienie w mechaniźmie:
d = zeros(15,1);

%Jakobian więzów geometrycznych:
C = Phiq(Q,l);
S = [
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]    
];
D = Orthogonal_Complement(q(:,i),l);

solution = [
Phi(Q,l,L)
traj(:,i+1) - [Q(14);Q(15)]
[D'; S*inv(M); C*inv(M)]*(M*((V-v(:,i))/dt) + d - f -B*u + C'*lagrange)
(Q - q(:,i))/dt - V
];

    if N == 2 && i >2
                solution = [
            (3*Q - 4 * q(:, i) + q(:,i-1))/(2*dt) - V;
            Phi(Q,l,L)
            traj(:,i+1) - [Q(14);Q(15)]
            [D'; S*inv(M); C*inv(M)]*(M*(3*V - 4*v(:,i) + v(:,i-1))/(2*dt) + d - f -B*u + C'*lagrange)
        ];
                
    elseif N == 3 && i > 3
                solution = [
            (11*Q - 18*q(:,i) + 9*q(:,i-1) -2*q(:,i-2))/(6*dt) - V;
            Phi(Q,l,L)
            traj(:,i+1) - [Q(14);Q(15)]
            [D'; S*inv(M); C*inv(M)]*(M*(11*V - 18*v(:,i) + 9*v(:,i-1) - 2*v(:,i-2))/(6*dt) + d - f -B*u + C'*lagrange)
        ];
    elseif N == 4 && i > 4
                solution = [
            Phi(Q,l,L)
            traj(:,i+1) - [Q(14);Q(15)]
            (25*Q - 48*q(:,i) + 36*q(:,i-1) - 16*q(:,i-2) + 3*q(:,i-3))/(12*dt) - V;
            [D'; S*inv(M); C*inv(M)]*(M*(25*V - 48*v(:,i) + 36*v(:,i-1) - 16*v(:,i-2) + 3*v(:,i-3))/(12*dt) + d - f -B*u + C'*lagrange)
        ];
    elseif N == 5 && i > 5
                solution = [
           Phi(Q,l,L)
           traj(:,i+1) - [Q(14);Q(15)]
           (137*Q - 300*q(:,i) + 300*q(:,i-1) - 200*q(:,i-2) + 75*q(:,i-3) - 12*q(:,i-4))/(60*dt) - V;
           [D'; S*inv(M); C*inv(M)]*(M*(137*V - 300*v(:,i) + 300*v(:,i-1) - 200*v(:,i-2) + 75*v(:,i-3) - 12*v(:,i-4))/(60*dt) + d - f -B*u + C'*lagrange)
        ];
    end

end
