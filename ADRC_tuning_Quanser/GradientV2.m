clear; clc;
load("ProperVals.mat");

dW = 0.05;
dB = 0.05;
b_hat_range = -0.2:dB:5.2;
omega_c_range = -3.2:dW:0.7;
Nb = length(b_hat_range);
Nomega = length(omega_c_range);

ISE = [validSettings.ISE];
ISE = reshape(ISE, [Nomega, Nb]);

%% Minimum ISE only in selected range:
% b_hat in [0, 5], omega_c in [-3, 0.5]
b_mask = (b_hat_range >= 0) & (b_hat_range <= 5);
w_mask = (omega_c_range >= -3) & (omega_c_range <= 0.5);

ISE_sub = ISE(w_mask, b_mask);
[min_ISE, idx_min] = min(ISE_sub(:));
[row_min, col_min] = ind2sub(size(ISE_sub), idx_min);

b_sub = b_hat_range(b_mask);
w_sub = omega_c_range(w_mask);

b_min = b_sub(col_min);
omega_c_min = w_sub(row_min);

fprintf('Minimum ISE in selected range:\n');
fprintf('ISE_min = %.6e\n', min_ISE);
fprintf('b_hat   = %.4f\n', b_min);
fprintf('omega_c = %.4f\n', omega_c_min);

%% Basic surface plot:
figure;
surf(b_hat_range, omega_c_range, ISE);
title('ISE surface');
xlabel('b_{hat}');
ylabel('\omega_c');
zlabel('ISE');
set(gca, 'ZScale', 'log');
clim([1e-5, 1e-1]);
colorbar

%% Gradient calculation on full range
[dISEdW, dISEdB] = gradient(ISE, dW, dB);
MagGrad = hypot(dISEdB, dISEdW);

% Gradient magnitude:
figure;
imagesc(b_hat_range, omega_c_range, MagGrad);
title('Magnitude of Gradient of ISE');
xlabel('b_{hat}');
ylabel('\omega_c');
set(gca, 'YDir', 'normal');
set(gca, 'ColorScale', 'log');
clim([1e-5, 1e-1]);
colorbar

% Gradient ratios:
figure;
RatioGrad = abs(dISEdB ./ dISEdW);
imagesc(b_hat_range, omega_c_range, RatioGrad);
title('Ratio of Gradient of ISE');
xlabel('b_{hat}');
ylabel('\omega_c');
set(gca, 'YDir', 'normal');
set(gca, 'ColorScale', 'log');
clim([1e-2, 1e2]);
colorbar

%% Hessian on full range
Hessian = zeros(Nomega, Nb);
[d2ISEd2W, d2ISEdWdB] = gradient(dISEdW, dW, dB);
[d2ISEdBdW, d2ISEd2B] = gradient(dISEdB, dW, dB);

for i = 1:Nomega
    for j = 1:Nb
        matrix = [d2ISEd2B(i, j), d2ISEdBdW(i, j);
                  d2ISEdWdB(i, j), d2ISEd2W(i, j)];
        lambda = eig(matrix);
        LambdaMax = abs(max(lambda));
        LambdaMin = abs(min(lambda));
        Hessian(i, j) = LambdaMax / LambdaMin;
    end
end

figure;
imagesc(b_hat_range, omega_c_range, Hessian);
title('Hessian');
xlabel('b_{hat}');
ylabel('\omega_c');
set(gca, 'YDir', 'normal');
set(gca, 'ColorScale', 'log');
clim([1e-3 1e3])
colorbar

%% Plots in selected range: b in [0,5], omega_c in [-3,0.5]

b_mask = (b_hat_range >= 0) & (b_hat_range <= 5);
w_mask = (omega_c_range >= -3) & (omega_c_range <= 0.5);

b_plot = b_hat_range(b_mask);
w_plot = omega_c_range(w_mask);

ISE_plot       = ISE(w_mask, b_mask);
MagGrad_plot   = MagGrad(w_mask, b_mask);
RatioGrad_plot = RatioGrad(w_mask, b_mask);
Kappa_plot     = Hessian(w_mask, b_mask);

%% 1. Surface plot ISE(b, omega)
figure;
hfig = figure;
pictureWidth = 20;
hwRatio = 0.65;
set(findall(hfig, '-property','FontSize'), 'FontSize', 16)
set(findall(hfig, '-property','Interpreter'), 'Interpreter', 'Latex');
set(findall(hfig, '-property','TickLabelInterpreter'), 'TickLabelInterpreter', 'Latex')
set(hfig, 'Units', 'centimeters', 'Position',[3 3 pictureWidth hwRatio*pictureWidth])
surf(b_plot, w_plot, ISE_plot);
xlabel('$log_{10}(\hat{b})$', 'Interpreter','latex','FontSize',16);
ylabel('$log_{10}(\omega_c)$', 'Interpreter','latex','FontSize',16);
title('Surface plot of $ISE$ index', 'Interpreter','latex','FontSize',16);
zlabel('ISE');
set(gca, 'ColorScale','log');
set(gca, 'ZScale', 'log');
clim([1e-4, 1e-1]);
colorbar
hold on;
plot3(b_min, omega_c_min, min_ISE, 'rx', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'r', 'LineWidth',5);
exportgraphics(gcf, 'ISE_surface.pdf', 'ContentType', 'vector');
%% 2. MagGrad(b, omega)
figure;
hfig = figure;
pictureWidth = 20;
hwRatio = 0.65;
set(findall(hfig, '-property','FontSize'), 'FontSize', 16)
set(findall(hfig, '-property','Interpreter'), 'Interpreter', 'Latex');
set(findall(hfig, '-property','TickLabelInterpreter'), 'TickLabelInterpreter', 'Latex')
set(hfig, 'Units', 'centimeters', 'Position',[3 3 pictureWidth hwRatio*pictureWidth])
imagesc(b_plot, w_plot, MagGrad_plot);
xlabel('$log_{10}(\hat{b})$', 'Interpreter','latex','FontSize',16);
ylabel('$log_{10}(\omega_c)$', 'Interpreter','latex','FontSize',16);
title('$||\nabla ISE||$', 'Interpreter','latex','FontSize',16);
set(gca, 'YDir', 'normal');
set(gca, 'ColorScale', 'log');
clim([1e-5 1])
colorbar
hold on
plot3(b_min, omega_c_min, min_ISE, 's', ...
    'MarkerSize', 6, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

drawAnnotation(MagGrad_plot(row_min, col_min), col_min, row_min, ...
    b_plot, w_plot, b_min, omega_c_min, '$\left\| \nabla ISE \right\|$');
exportgraphics(gcf, 'MagGrad.pdf', 'ContentType', 'vector');
%% 3. RatioGrad(b, omega)
figure;
hfig = figure;
pictureWidth = 20;
hwRatio = 0.65;
set(findall(hfig, '-property','FontSize'), 'FontSize', 16)
set(findall(hfig, '-property','Interpreter'), 'Interpreter', 'Latex');
set(findall(hfig, '-property','TickLabelInterpreter'), 'TickLabelInterpreter', 'Latex')
set(hfig, 'Units', 'centimeters', 'Position',[3 3 pictureWidth hwRatio*pictureWidth])
imagesc(b_plot, w_plot, RatioGrad_plot);
title('$\frac{ISE_{\hat{b}}}{ISE_{\omega_c}}$', 'Interpreter','latex','FontSize',16);
xlabel('$log_{10}(\hat{b})$', 'Interpreter','latex','FontSize',16);
ylabel('$log_{10}(\omega_c)$', 'Interpreter','latex','FontSize',16);
set(gca, 'YDir', 'normal');
set(gca, 'ColorScale', 'log');
colorbar
clim([1e-1 1e2])
hold on
plot3(b_min, omega_c_min, min_ISE, 's', ...
    'MarkerSize', 6, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
drawAnnotation(RatioGrad_plot(row_min, col_min), col_min, row_min, ...
    b_plot, w_plot, b_min, omega_c_min, ...
    '$\left|\frac{\partial \mathrm{ISE}/\partial b}{\partial \mathrm{ISE}/\partial \omega_c}\right|$');
exportgraphics(gcf, 'RatioGrad.pdf', 'ContentType', 'vector');
%% 4. Kappa(b, omega)
figure;
hfig = figure;
pictureWidth = 20;
hwRatio = 0.65;
set(findall(hfig, '-property','FontSize'), 'FontSize', 16)
set(findall(hfig, '-property','Interpreter'), 'Interpreter', 'Latex');
set(findall(hfig, '-property','TickLabelInterpreter'), 'TickLabelInterpreter', 'Latex')
set(hfig, 'Units', 'centimeters', 'Position',[3 3 pictureWidth hwRatio*pictureWidth])
imagesc(b_plot, w_plot, Kappa_plot);
xlabel('$log_{10}(\hat{b})$', 'Interpreter','latex','FontSize',16);
ylabel('$log_{10}(\omega_c)$', 'Interpreter','latex','FontSize',16);
title('ISE: $\kappa = \frac{\lambda_{max}}{\lambda_{min}}$', 'Interpreter','latex','FontSize',16)
set(gca, 'YDir', 'normal');
set(gca, 'ColorScale', 'log');
clim([1e-3 1e3])
colorbar
hold on
plot3(b_min, omega_c_min, min_ISE, 's', ...
    'MarkerSize', 6, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

drawAnnotation(Kappa_plot(row_min, col_min), col_min, row_min, ...
    b_plot, w_plot, b_min, omega_c_min, '$\kappa$');
exportgraphics(gcf, 'Kappa.pdf', 'ContentType', 'vector');
%% Helper: mark optimal point with annotation box and arrowfunction drawAnnotation(val, ix, iy, b_hat, omega_c, xOpt, yOpt, labelText)
function drawAnnotation(val, ix, iy, b_hat, omega_c, xOpt, yOpt, labelText)

    % Stała pozycja textboxa względem osi:
    % [x y] w jednostkach normalized osi
    boxPos = [0.68 0.35];

    dB = b_hat(2) - b_hat(1);
    dW = omega_c(2) - omega_c(1);

    % Tekst jako cell array - dużo stabilniejsze dla LaTeX w MATLAB
    txt = { ...
        '\textbf{ISE optimal point}', ...
        ['$\log_{10}(\hat{b}) = ' sprintf('%.2f', b_hat(ix)) '$'], ...
        ['$\log_{10}(\omega_c) = ' sprintf('%.2f', omega_c(iy)) '$'], ...
        [labelText ' = ' sprintf('%.5g', val)] ...
    };

    % Textbox osadzony w osiach - stałe położenie względem wykresu
    text(boxPos(1), boxPos(2), txt, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'Interpreter', 'latex', ...
        'FontSize', 11, ...
        'BackgroundColor', 'w', ...
        'EdgeColor', 'k', ...
        'Margin', 6);

    % Marker punktu optymalnego - mała czerwona ramka
    hold on
    plot(xOpt, yOpt, 's', ...
        'MarkerSize', 4, ...
        'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', 'r', ...
        'LineWidth', 1.2);

    % Konwersja punktu optymalnego z danych osi do normalized figure units
    ax = gca;
    fig = gcf;
    axPos = ax.Position;
    xLim = ax.XLim;
    yLim = ax.YLim;

    % Delikatne przesunięcie końca strzałki względem markera
    x0 = xOpt + dB/4;
    y0 = yOpt + dW/4;

    xN = axPos(1) + (x0 - xLim(1)) / (xLim(2) - xLim(1)) * axPos(3);
    yN = axPos(2) + (y0 - yLim(1)) / (yLim(2) - yLim(1)) * axPos(4);

    % Początek strzałki: liczony od stałego boxPos osadzonego w osiach
    xBoxN = axPos(1) + boxPos(1) * axPos(3);
    yBoxN = axPos(2) + boxPos(2) * axPos(4);

    annotation(fig, 'arrow', [xBoxN xN], [yBoxN yN], ...
        'LineWidth', 1.5, ...
        'Color', 'k');
end