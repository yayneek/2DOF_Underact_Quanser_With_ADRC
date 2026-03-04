clear; clc;
% Loading the data:
data = load("valid_settings.mat");
data = data.validSettings;
IAE = [data.IAE];
ISE = [data.ISE];
nB = 100;
nW = 100;

% Optimal points:
ISEWopt = log10(1);
ISEBopt = log10(44306);
IAEWopt = log10(1);
IAEBopt = log10(39442);

b_hat = linspace(0,5,nB);
omega_c = linspace(-5,0,nW);

%Reshaping:
IAE = reshape(IAE, [nW, nB]);
ISE = reshape(ISE, [nW, nB]);

dB = b_hat(2) - b_hat(1);
dW = omega_c(2) - omega_c(1);
[dISEdB, dISEdW] = gradient(ISE, dB, dW);
[dIAEdB, dIAEdW] = gradient(IAE, dB, dW);


% Calculate the magnitude of the gradients
magnitudeGradientISE = hypot(dISEdB, dISEdW);
magnitudeGradientIAE = hypot(dIAEdB, dIAEdW);

% Calcualte 
ISEgradientRatio = abs(dISEdB./dISEdW);
IAEgradientRatio = abs(dIAEdB./dIAEdW);
%% Ratios:
% Visualizing ratios of the gradient components:
figure;
imagesc(b_hat, omega_c, ISEgradientRatio)
xlabel('$log_{10}(\hat{b})$', 'Interpreter','latex','FontSize',16);
ylabel('$log_{10}(\omega_c)$', 'Interpreter','latex','FontSize',16);
title('Heatmap of gradient component ratio for ISE', 'Interpreter','latex','FontSize',16);
set(gca,'YDir','normal');  
set(gca, 'ColorScale','log');
clim([1e-3 1])
colorbar;
hold on;
rectangle('Position',[ISEBopt, ISEWopt, dB/2, dW/2], ...
          'EdgeColor','r', 'LineWidth',1);
boxPos = [0.5 0.5];  % [x y] w Units='normalized' osi (0..1). Zmień numerycznie.

[~, ix] = min(abs(b_hat - ISEBopt));
[~, iy] = min(abs(omega_c - ISEWopt));

val = magnitudeGradientISE(iy, ix);

txt = sprintf(['\\bf Optimum\\rm\n' ...
               '$\\log_{10}(\\hat b)=%.3f$\n' ...
               '$\\log_{10}(\\omega_c)=%.3f$\n' ...
               '$\\left\\|\\nabla\\,\\mathrm{ISE}\\right\\|=%.3g$'], ...
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
grid on


figure;
imagesc(b_hat, omega_c, IAEgradientRatio)
xlabel('$log_{10}(\hat{b})$', 'Interpreter','latex','FontSize',16);
ylabel('$log_{10}(\omega_c)$', 'Interpreter','latex','FontSize',16);
title('Heatmap of gradient component ratio for IAE', 'Interpreter','latex','FontSize',16);
set(gca,'YDir','normal');  
set(gca, 'ColorScale','log');
clim([1e-3 1])
colorbar;
hold on;
rectangle('Position',[ISEBopt, ISEWopt, dB/2, dW/2], ...
          'EdgeColor','r', 'LineWidth',1);
grid on

%% Magnitude:
% Visualizing the magnitude of the gradients
figure;
imagesc(b_hat, omega_c, magnitudeGradientISE);
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




figure;
imagesc(b_hat, omega_c, magnitudeGradientIAE);
set(gca,'YDir','normal');   % <- to "odbija" mapę w pionie
xlabel('$log_{10}(\hat{b})$', 'Interpreter','latex','FontSize',16);
ylabel('$log_{10}(\omega_c)$', 'Interpreter','latex','FontSize',16);
title('Heatmap of $\left\| \nabla \mathrm{IAE} \right\|$', 'Interpreter','latex','FontSize',16);
set(gca, 'ColorScale','log');
clim([1e-4 1])
colorbar;

hold on
rectangle('Position',[IAEBopt, IAEWopt, dB/2, dW/2], ...
          'EdgeColor','r', 'LineWidth',2);

%% Hessian analysis:
% Compute the Hessian matrices for IAE and ISE
% ISE:

[d2ISEd2B, d2ISEdBdW] = gradient(dISEdB, dB, dW);
[d2ISEdWdB, d2ISEd2W] = gradient(dISEdW, dB, dW);
ISEHessian = zeros(100,100);

for i = 1:nB
    for j =1:nW
        matrix = [d2ISEd2B(i, j), d2ISEdBdW(i, j); d2ISEdWdB(i, j), d2ISEd2W(i, j)];
        lambda = eig(matrix);
        % if lambda(1) > 0 && lambda(2) > 0
        %     ISEHessian(i, j) = 1; % Mark positive definiteness
        % elseif lambda(1) < 0 && lambda(2) < 0
        %     ISEHessian(i, j) = -1; % Mark non-positive definiteness
        % else
        %     ISEHessian(i, j) = 0; %Saddle point
        % end
        LambdaMax = abs(max(lambda));
        LambdaMin = abs(min(lambda));
        ISEHessian(i, j) = LambdaMax/LambdaMin;
    end
end



figure;
imagesc(b_hat, omega_c, ISEHessian);
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


%IAE:
[d2IAEd2B, d2IAEdBdW] = gradient(dIAEdB, dB, dW);
[d2IAEdWdB, d2IAEd2W] = gradient(dIAEdW, dB, dW);
IAEHessian = zeros(100,100);

for i = 1:nB
    for j =1:nW
        matrix = [d2IAEd2B(i, j), d2IAEdBdW(i, j); d2IAEdWdB(i, j), d2IAEd2W(i, j)];
        lambda = eig(matrix);
        % if lambda(1) > 0 && lambda(2) > 0
        %     IAEHessian(i, j) = 1; % Mark positive definiteness
        % elseif lambda(1) < 0 && lambda(2) < 0
        %     IAEHessian(i, j) = -1; % Mark non-positive definiteness
        % else
        %     IAEHessian(i, j) = 0; %Saddle point
        % end
        LambdaMax = abs(max(lambda));
        LambdaMin = abs(min(lambda));
        IAEHessian(i,j) = LambdaMax/LambdaMin;
    end
end

figure;
imagesc(b_hat, omega_c, IAEHessian);
set(gca,'YDir','normal');   % <- to "odbija" mapę w pionie
xlabel('$log_{10}(\hat b)$', 'Interpreter','latex','FontSize',16);
ylabel('$log_{10}(\omega_c)$', 'Interpreter','latex','FontSize',16);
title('IAE: $\kappa = \frac{\lambda_{max}}{\lambda_{min}}$', 'Interpreter','latex','FontSize',16)
set(gca, 'ColorScale','log');
clim([1e-3 1e3]);
colorbar;
hold on
rectangle('Position',[IAEBopt-dB, IAEWopt-dW, dB, dW], ...
          'EdgeColor','r', 'LineWidth',1);
grid on

