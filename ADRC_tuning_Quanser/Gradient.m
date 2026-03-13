clear; clc;
% Loading the data:
data32 = load("ExtendedOmegaRange.mat");
data50 = load("valid_settings.mat");
datab5Plus = load('b_hat5+.mat');
datab0minus = load("b_hat0-.mat");

datab0minus = datab0minus.validSettings;
datab5Plus = datab5Plus.validSettings;
data50 = data50.validSettings;
data32 = data32.validSettings;

ISE32 = [data32.ISE];
ISE50 = [data50.ISE];
ISEb5Plus = [datab5Plus.ISE];
ISEb05minus = [datab0minus.ISE];

nB = 100;
nW = 100;

% Optimal points:
ISEWopt = log10(1.9179);

ISEBopt = log10(100000);

b_hat = linspace(0,5,nB);
omega_c32 = linspace(-3, 2, nW);
omega_c50 = linspace(-5, 0, nW);

dB = b_hat(2) - b_hat(1);
dW = omega_c50(2) - omega_c50(1);
b_hat5PlusRange = [5+dB, 5+2*dB, 5+3*dB];
b_hat0minusRange = [-3*dB, -2*dB, -dB];

%Reshaping:
ISE32 = reshape(ISE32, [nW, nB]);
ISE50 = reshape(ISE50, [nW, nB]);
ISEb5Plus = reshape(ISEb5Plus, [nW, 3]);
ISEb05minus = reshape(ISEb05minus, [nW, 3]);


[dISEdB, dISEdW] = gradient(ISE32, dB, dW);

% Calculate the magnitude of the gradients
magnitudeGradientISE = hypot(dISEdB, dISEdW);

% Calcualte 
ISEgradientRatio = abs(dISEdB./dISEdW);
[~, ix] = min(abs(b_hat - ISEBopt));
[~, iy] = min(abs(omega_c32 - ISEWopt));

%% Surface plot of extended b range:
ISETotal = [ISE32, ISEb5Plus];
b = [b_hat, b_hat5PlusRange];
figure;
imagesc(b, omega_c32, ISETotal);
xlabel('$log_{10}(\hat{b})$', 'Interpreter','latex','FontSize',16);
ylabel('$log_{10}(\omega_c)$', 'Interpreter','latex','FontSize',16);
title('Heatmap of $ISE$ index', 'Interpreter','latex','FontSize',16);
set(gca, 'ColorScale','log');
set(gca,'YDir','normal');
clim([1e-5 1e1])
colorbar

%% Surface plot:
% Displaying surface of ISE index:

omega_cTot = [omega_c50, omega_c32(61:end)];
ISETot = [ISE50; ISE32(61:end, :)];

figure;
surf(10.^b_hat, 10.^omega_cTot, ISETot);
xlabel('$log_{10}(\hat{b})$', 'Interpreter','latex','FontSize',16);
ylabel('$log_{10}(\omega_c)$', 'Interpreter','latex','FontSize',16);
title('Heatmap of $ISE$ index', 'Interpreter','latex','FontSize',16);
set(gca, 'ColorScale','log');
colorbar;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'ZScale', 'log');
clim([1e-5 1e1])
colorbar
%% ISE Heatmap:
figure;
imagesc(b_hat, omega_cTot, ISETot);
xlabel('$log_{10}(\hat{b})$', 'Interpreter','latex','FontSize',16);
ylabel('$log_{10}(\omega_c)$', 'Interpreter','latex','FontSize',16);
title('Heatmap of $ISE$ index', 'Interpreter','latex','FontSize',16);
set(gca,'YDir','normal');  
set(gca, 'ColorScale','log');
clim([5e-4 1e-2])
colorbar;

%% Ratios:
% Visualizing ratios of the gradient components:
figure;
imagesc(b_hat, omega_cTot, ISEgradientRatio)
xlabel('$log_{10}(\hat{b})$', 'Interpreter','latex','FontSize',16);
ylabel('$log_{10}(\omega_c)$', 'Interpreter','latex','FontSize',16);
title('Heatmap of gradient component ratio for ISE', 'Interpreter','latex','FontSize',16);
set(gca,'YDir','normal');  
set(gca, 'ColorScale','log');
clim([1e-4 1e1])
colorbar;
hold on;
rectangle('Position',[ISEBopt, ISEWopt, dB/2, dW/2], ...
          'EdgeColor','r', 'LineWidth',1);

val = ISEgradientRatio(iy, ix);
drawAnnotation(val,ix,iy,b_hat,omega_c32,ISEBopt,ISEWopt);
grid on
%% Magnitude:
% Visualizing the magnitude of the gradients
figure;
imagesc(b_hat, omega_cTot, magnitudeGradientISE);
xlabel('$log_{10}(\hat{b})$', 'Interpreter','latex','FontSize',16);
ylabel('$log_{10}(\omega_c)$', 'Interpreter','latex','FontSize',16);
title('Heatmap of $\left\| \nabla \mathrm{ISE} \right\|$', 'Interpreter','latex','FontSize',16);
set(gca,'YDir','normal');   
set(gca, 'ColorScale','log');
clim([1e-4 1])
colorbar;

hold on
rectangle('Position',[ISEBopt, ISEWopt, dB/2, dW/2], ...
          'EdgeColor','r', 'LineWidth',2);

val = magnitudeGradientISE(iy, ix);
drawAnnotation(val,ix,iy,b_hat,omega_cTot,ISEBopt,ISEWopt);
grid on
%% Hessian analysis:
% ISE:

[d2ISEd2B, d2ISEdBdW] = gradient(dISEdB);
[d2ISEdWdB, d2ISEd2W] = gradient(dISEdW);
ISEHessian = zeros(100,100);

for i = 1:nB
    for j =1:nW
        matrix = [d2ISEd2B(i, j), d2ISEdBdW(i, j); d2ISEdWdB(i, j), d2ISEd2W(i, j)];
        lambda = eig(matrix);
        LambdaMax = abs(max(lambda));
        LambdaMin = abs(min(lambda));
        ISEHessian(i, j) = LambdaMax/LambdaMin;
    end
end



figure;
imagesc(b_hat, omega_cTot, ISEHessian);
set(gca,'YDir','normal');   % <- to "odbija" mapę w pionie
xlabel('$log_{10}(\hat b)$', 'Interpreter','latex','FontSize',16);
ylabel('$log_{10}(\omega_c)$', 'Interpreter','latex','FontSize',16);
title('ISE: $\kappa = \frac{\lambda_{max}}{\lambda_{min}}$', 'Interpreter','latex','FontSize',16)
set(gca, 'ColorScale','log');
clim([1e-3 1e3])
colorbar;
hold on
rectangle('Position',[ISEBopt ISEWopt, dB, dW], ...
          'EdgeColor','r', 'LineWidth',1);
grid on



function drawAnnotation(val, ix, iy, b_hat, omega_c,ISEBopt, ISEWopt,txt)
boxPos = [0.75 0.25];
dB = b_hat(2) - b_hat(1);
dW = omega_c(2) - omega_c(1);
    txt = sprintf(['\\bf Optimum\\rm\n' ...
               '$\\log_{10}(\\hat b)=%.3f$\n' ...
               '$\\log_{10}(\\omega_c)=%.3f$\n' ...
               '$\\left| \\nabla ISE \\right|=%.3g$'], ...
               b_hat(ix), omega_c(iy), val);

text(boxPos(1),boxPos(2),txt,'Units','normalized', ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'Interpreter','latex','FontSize',11, ...
    'BackgroundColor','w','EdgeColor','k','Margin',6);

ax = gca; fig = gcf;
axPos = ax.Position; xLim = ax.XLim; yLim = ax.YLim;

x0 = ISEBopt + dB/4;
y0 = ISEWopt + dW/4;

xN = axPos(1) + (x0 - xLim(1)) / (xLim(2) - xLim(1)) * axPos(3);
yN = axPos(2) + (y0 - yLim(1)) / (yLim(2) - yLim(1)) * axPos(4);

xBoxN = axPos(1) + boxPos(1)*axPos(3);
yBoxN = axPos(2) + boxPos(2)*axPos(4);

annotation(fig,'arrow',[xBoxN xN],[yBoxN yN], ...
    'LineWidth',1.5,'Color','k');
end