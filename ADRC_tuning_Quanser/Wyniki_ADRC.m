Inverse_Dynamics;
%Wykresy3D;
b_hat = 10^5;
omega_c = 10^0.3;
close all;
out = sim('ADRC_Quanser_final.slx',10);

q = out.coordinates.data;
error_fb_ffw = out.error_fb_ffw.data';
error_fb = out.error_fb.data';
error_ffw = out.error_ffw.data';
u_ffw = u;
u_fb = out.u_feedback.data';

e_norm_fb_ffw = zeros(1, length(t));
e_norm_fb = zeros(1, length(t));
e_norm_ffw = zeros(1, length(t));

for i =1:length(t)
    e_norm_ffw(i) = norm(error_ffw(:,i));
    e_norm_fb_ffw(i) = norm(error_fb_ffw(:,i));
    e_norm_fb(i) = norm(error_fb(:,i));
end

% Dla układu FB+FFW:
ISE_fb_ffw = cumtrapz(t, e_norm_fb_ffw.^2);   % całka z e^2
IAE_fb_ffw = cumtrapz(t, abs(e_norm_fb_ffw)); % całka z |e|

% Dla układu FB:
ISE_fb = cumtrapz(t, e_norm_fb.^2);
IAE_fb = cumtrapz(t, abs(e_norm_fb));

% Dla układu FFW:
ISE_ffw = cumtrapz(t, e_norm_ffw.^2);
IAE_ffw = cumtrapz(t, abs(e_norm_ffw));

figure()
plot(t, ISE_fb_ffw, t, ISE_ffw, t, ISE_fb)
grid on
title('ISE')
legend('FFW+FB','FFW','FB')

figure()
plot(t, IAE_fb_ffw, t, IAE_ffw, t, IAE_fb)
grid on
title('IAE')
legend('FFW+FB','FFW','FB')

figure()
plot(t, error_fb_ffw, t, error_ffw)
legend('FFW+FB','FFW')
%% Sygnały sterujące
%u_ffw
folder = 'C:\Users\Jan Piliszczuk\Desktop\plots';
figure()
plot(t, u_ffw)
title('$u_{ff}(t)$', 'Interpreter', 'latex', 'FontSize',16);
ylabel('$u$[Nm]', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('$t$[s]', 'Interpreter', 'latex', 'FontSize', 16);
filename = 'uffw(t).pdf';
grid on
legend('$u_1$','$u_2$','Interpreter','latex', 'FontSize',16, 'location','southwest')
fullfile_pdf = fullfile(folder, filename);
exportgraphics(gcf,fullfile_pdf);

%u_fb
figure()
plot(t, u_fb)
title('$u_{fb}(t)$', 'Interpreter', 'latex', 'FontSize',16);
filename = 'ufb(t).pdf';
ylabel('$u$[Nm]', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('$t$[s]', 'Interpreter', 'latex', 'FontSize', 16);
grid on
legend('$u_1$','$u_2$','Interpreter','latex', 'FontSize',16, 'location','southwest')
fullfile_pdf = fullfile(folder, filename);
exportgraphics(gcf,fullfile_pdf);

%u_ffw_fb
figure()
plot(t, u_fb+u_ffw)
title('$u_{ff}(t) + u_{fb}$', 'Interpreter', 'latex', 'FontSize',16);
filename = 'uffwfb(t).pdf';
ylabel('$u$[Nm]', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('$t$[s]', 'Interpreter', 'latex', 'FontSize', 16);
grid on
legend('$u_1$','$u_2$','Interpreter','latex', 'FontSize',16, 'location','southwest')
fullfile_pdf = fullfile(folder, filename);
exportgraphics(gcf,fullfile_pdf);
%% Wskaźniki całkowe
%ISE 
figure()
plot(t, ISE_fb_ffw, t, ISE_ffw, t, ISE_fb)
grid on
ylabel('$ISE$[m]', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('$t$[s]', 'Interpreter', 'latex', 'FontSize', 16);
title('ISE$(t)$', 'Interpreter', 'latex', 'FontSize', 16);
legend('combined','feedforward', 'feedback', 'Interpreter', 'latex', 'FontSize', 16, 'location', 'northwest');

%IAE
figure()
plot(t, IAE_fb_ffw, t, IAE_ffw, t, IAE_fb)
grid on
ylabel('$IAE$[m]', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('$t$[s]', 'Interpreter', 'latex', 'FontSize', 16);
title('IAE$(t)$', 'Interpreter', 'latex', 'FontSize', 16);
legend('combined','feedforward', 'feedback', 'location','northwest','Interpreter', 'latex', 'FontSize', 16);
%% Współrzędne ładunku 
figure()
plot(traj(1,:), traj(2,:),traj(1,:)+error_fb_ffw(1,:), traj(2,:)+error_fb_ffw(2,:), traj(1,:)+error_ffw(1,:), traj(2,:)+error_ffw(2,:),traj(1,:)+error_fb(1,:), traj(2,:)+ error_fb(2,:))
title('','Interpreter', 'latex', 'FontSize', 16);
legend('reference trajectory','combined','ffw','fb','Interpreter', 'latex', 'FontSize', 16);
grid on

%% Uchyby


figure();
hold on
filename = 'e(t)_ffwfb.pdf';
plot(t, error_fb_ffw(1, :), '-',  'LineWidth', 1);  
plot(t, error_fb_ffw(2, :), '-', 'LineWidth', 1);
grid on
title('$e(t)$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$e$ [m]', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 16);
legend('$e_x$','$e_y$','Interpreter','latex', 'FontSize',16, 'location','southwest')
box on
fullfile_pdf = fullfile(folder, filename);
exportgraphics(gcf,fullfile_pdf);

%%

figure();
hold on;

% --- Właściwe krzywe ---
% 1) Błąd łączony (fb + ffw)
c1 = "#0072BD";  % kolor "combined"
p1 = plot(t, error_fb_ffw(1, :), '-',  'LineWidth', 1.5, 'Color', c1);  
p2 = plot(t, error_fb_ffw(2, :), '--', 'LineWidth', 1.5, 'Color', c1);

% 2) Błąd feedforward (ffw)
c2  = "#D95319"; % kolor "ffw"
p3 = plot(t, error_ffw(1, :), '-',  'LineWidth', 1.5, 'Color', c2);
p4 = plot(t, error_ffw(2, :), '--', 'LineWidth', 1.5, 'Color', c2);

% 3) Błąd feedback (fb)
c3 = "#EDB120";  % kolor "fb"
p5 = plot(t, error_fb(1, :), '-',  'LineWidth', 1.5, 'Color', c3);
p6 = plot(t, error_fb(2, :), '--', 'LineWidth', 1.5, 'Color', c3);

% --- Formatowanie wykresu ---
grid on;
title('$e(t)$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$e$ [m]', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 16);

% --- Pierwsza legenda opisująca schematy sterowania ---
h(1) = plot(NaN, NaN, 'LineWidth', 2, 'Color', c1);
h(2) = plot(NaN, NaN, 'LineWidth', 2, 'Color', c2);
h(3) = plot(NaN, NaN, 'LineWidth', 2, 'Color', c3);
l1 = legend(h, {'Combined', 'Feedforward', 'Feedback'}, ...
            'Location', 'northwest', ...
            'Interpreter', 'latex', ...
            'FontSize', 12);
box on
% --- Druga legenda opisująca styl linii ---
% Utwórz niewidoczne osie tylko do drugiej legendy
ah1 = axes('Position', get(gca, 'Position'), 'Color', 'none', 'Visible', 'off');
hold(ah1, 'on');
h2(1) = plot(ah1, NaN, NaN, '-',  'LineWidth', 1.5, 'Color', 'k');
h2(2) = plot(ah1, NaN, NaN, '--', 'LineWidth', 1.5, 'Color', 'k');
l2 = legend(ah1, h2, {'$e_x$', '$e_y$'}, ...
            'Location', 'northeast', ...
            'Interpreter', 'latex', ...
            'FontSize', 12);
set(l2, 'Color', 'white');
hold off;

filename = 'e(t)_combined.pdf';
fullfile_pdf = fullfile(folder, filename);
exportgraphics(gcf,fullfile_pdf);

%% Wykres ładunku oraz punktu charakterystycznego i porównanie z trajektorią:
q = out.coordinates.Data;
traj = out.trajectory.Data;
close all
r = [l;0];
R = @(theta)[cos(theta) -sin(theta); sin(theta) cos(theta)];
for i =1:length(t)
    char(:,i) = R(q(3,i))*r + R(q(6,i))*r; 
end

figure()
plot(q(14,:),q(15,:), char(1,:), char(2,:))
title('Location of the load', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('$x[m]$', 'Interpreter','latex' , 'FontSize', 16);
ylabel('$y[m]$', 'Interpreter','latex' , 'FontSize', 16);
legend('Load','Characteristic point', 'Interpreter','latex' , 'FontSize', 16)
grid on
filename = 'XY_char.pdf';
fullfile_pdf = fullfile(filename);
exportgraphics(gcf, fullfile_pdf);


figure()
plot(traj(1,:), traj(2,:),q(14,:),q(15,:))
title('Desired and output trajectories', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('$x[m]$', 'Interpreter','latex' , 'FontSize', 16);
ylabel('$y[m]$', 'Interpreter','latex' , 'FontSize', 16);
legend('Desired trajectory','Output of the load', 'Interpreter','latex' , 'FontSize', 16)
grid on
filename = 'XY_load.pdf';
fullfile_pdf = fullfile(filename);
exportgraphics(gcf, fullfile_pdf);