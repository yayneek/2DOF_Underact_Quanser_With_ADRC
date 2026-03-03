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


% Visualizing the magnitude of the gradients
figure;
imagesc(b_hat, omega_c, magnitudeGradientISE);
xlabel('b_{hat}');
ylabel('omega_c');
zlabel('Magnitude of Gradient');
title('ISE Gradient Magnitude Surface');
set(gca,'YDir','normal');   % <- to "odbija" mapę w pionie
set(gca, 'ColorScale','log');
clim([1e-4 1])
colorbar;
hold on
rectangle('Position',[ISEBopt-dB, ISEWopt-dW, dB, dW], ...
          'EdgeColor','r', 'LineWidth',1);

figure;
imagesc(b_hat, omega_c, magnitudeGradientIAE);
set(gca,'YDir','normal');   % <- to "odbija" mapę w pionie
xlabel('b_{hat}');
ylabel('omega_c');
zlabel('Magnitude of Gradient');
title('IAE Gradient Magnitude Surface');
set(gca, 'ColorScale','log');
clim([1e-4 1])
colorbar;
hold on
rectangle('Position',[IAEBopt-dB, IAEWopt-dW, dB, dW], ...
          'EdgeColor','r', 'LineWidth',1);

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
xlabel('$log_{10}(\hat b)$', 'Interpreter','latex');
ylabel('$log_{10}(\omega_c)$', 'Interpreter','latex');
title('ISE: $\kappa = \frac{\lambda_{max}}{\lambda_{min}}$', 'Interpreter','latex')
set(gca, 'ColorScale','log');
clim([1e-3 1e3])
colorbar;
hold on
rectangle('Position',[ISEBopt-dB, ISEWopt-dW, dB, dW], ...
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
        LambdaMax = max(lambda);
        LambdaMin = min(lambda);
        IAEHessian(i,j) = LambdaMax/LambdaMin;
    end
end

figure;
imagesc(b_hat, omega_c, IAEHessian);
set(gca,'YDir','normal');   % <- to "odbija" mapę w pionie
xlabel('$log_{10}(\hat b)$', 'Interpreter','latex');
ylabel('$log_{10}(\omega_c)$', 'Interpreter','latex');
title('IAE: $\kappa = \frac{\lambda_{max}}{\lambda_{min}}$', 'Interpreter','latex')
set(gca, 'ColorScale','log');
clim([1e-3 1e3]);
colorbar;
hold on
rectangle('Position',[IAEBopt-dB, IAEWopt-dW, dB, dW], ...
          'EdgeColor','r', 'LineWidth',1);
grid on