clear; clc;

% Loading the data:
data32      = load("ExtendedOmegaRange.mat");      % omega_c in [-3, 2]
data50m     = load("valid_settings.mat");          % omega_c in [-5, 0]
data05p     = load("valid_settings.mat");      % omega_c in [0, 5]  <-- podmien na swoja nazwe pliku
datab5Plus  = load("b_hat5+.mat");
datab0minus = load("b_hat0-.mat");

data32      = data32.validSettings;
data50m     = data50m.validSettings;
data05p     = data05p.validSettings;
datab5Plus  = datab5Plus.validSettings;
datab0minus = datab0minus.validSettings;

ISE32       = [data32.ISE];
ISE50m      = [data50m.ISE];
ISE05p      = [data05p.ISE];
ISEb5Plus   = [datab5Plus.ISE];
ISEb05minus = [datab0minus.ISE];

nB = 100;
nW = 100;

% Optimal points:
ISEWopt = log10(1.9179);
ISEBopt = log10(100000);

b_hat      = linspace(0,5,nB);
omega_c32  = linspace(-3, 2, nW);
omega_c50m = linspace(-5, 0, nW);
omega_c05p = linspace(0, 5, nW);


dW = omega_c50m(2) - omega_c50m(1);
dB = b_hat(2) - b_hat(1);
b_hat5PlusRange  = [5+dB, 5+2*dB, 5+3*dB];
b_hat0minusRange = [-3*dB, -2*dB, -dB];

% Reshaping:
ISE32       = reshape(ISE32,       [nW, nB]);   % 100 x 100
ISE50m      = reshape(ISE50m,      [nW, nB]);   % 100 x 100
ISE05p      = reshape(ISE05p,      [nW, nB]);   % 100 x 100
ISEb5Plus   = reshape(ISEb5Plus,   [nW, 3]);    % 100 x 3
ISEb05minus = reshape(ISEb05minus, [nW, 3]);    % 100 x 3

%% Final omega range
idxLeft  = find(omega_c50m < omega_c32(1));
idxLeft  = idxLeft(end-3:end);                  % 4 probki tuz przed -3

idxRight = find(omega_c05p > omega_c32(end));   % probki powyzej 2

omega_cTotal = [omega_c50m(idxLeft), omega_c32, omega_c05p(idxRight)];
ISEmid       = [ISE50m(idxLeft,:); ISE32; ISE05p(idxRight,:)];
ISEleft  = [ISEb05minus(idxLeft,:); ISEb05minus; ISEb05minus(idxRight,:)];
ISEright = [ISEb5Plus(idxLeft,:);   ISEb5Plus;   ISEb5Plus(idxRight,:)];

b_hatTotal = [b_hat0minusRange, b_hat, b_hat5PlusRange];
ISETotal   = [ISEleft, ISEmid, ISEright];

[~, ix] = min(abs(b_hatTotal - ISEBopt));
[~, iy] = min(abs(omega_cTotal - ISEWopt));

% Surface plot
figure;
surf(b_hatTotal, omega_cTotal, ISETotal);

xlabel('$\log_{10}(\hat{b})$', 'Interpreter','latex','FontSize',16);
ylabel('$\log_{10}(\omega_c)$', 'Interpreter','latex','FontSize',16);
zlabel('$ISE$', 'Interpreter','latex','FontSize',16);
title('Heatmap of $ISE$ index', 'Interpreter','latex','FontSize',16);

set(gca, 'ColorScale','log');
set(gca, 'ZScale', 'log');
clim([1e-5 1e1]);
colorbar;

%% Gradient calculation:
[dISEdB, dISEdW] = gradient(ISETotal, b_hatTotal, omega_cTotal);
magnitudeGradientISE = hypot(dISEdB, dISEdW);
gradientRatioISE = abs(dISEdB ./ dISEdW);

%% Hessian:
[d2ISEd2B, d2ISEdBdW] = gradient(dISEdB, b_hatTotal, omega_cTotal);
[d2ISEdWdB, d2ISEd2W] = gradient(dISEdW, b_hatTotal, omega_cTotal);

Hessian = zeros(size(ISETotal));

for i = 1:size(ISETotal,1)
    for j = 1:size(ISETotal,2)
        matrix = [d2ISEd2B(i,j), d2ISEdWdB(i,j); ...
                  d2ISEdWdB(i,j), d2ISEd2W(i,j)];
        lambda = eig(matrix);
        LambdaMax = max(abs(lambda));
        LambdaMin = min(abs(lambda));

        if LambdaMin == 0
            Hessian(i,j) = NaN;
        else
            Hessian(i,j) = LambdaMax / LambdaMin;
        end
    end
end

%% Plotting:
idxW = omega_cTotal >= -3 & omega_cTotal <= 0.5;
idxB = b_hatTotal   >= 0  & b_hatTotal   <= 5;

omegaPlot = omega_cTotal(idxW);
bPlot     = b_hatTotal(idxB);

ISEplot   = ISETotal(idxW, idxB);
gradMag   = magnitudeGradientISE(idxW, idxB);
gradRatio = gradientRatioISE(idxW, idxB);
kappaPlot = Hessian(idxW, idxB);

% ISE surface plot:
figure;
surf(bPlot, omegaPlot, ISEplot);
xlabel('$\log_{10}(\hat{b})$', 'Interpreter','latex')
ylabel('$\log_{10}(\omega_c)$', 'Interpreter','latex')
zlabel('$ISE$', 'Interpreter','latex')
title('$ISE$', 'Interpreter','latex')
set(gca,'ColorScale','log')
clim([1e-4 1])
colorbar
set(gca,'ZScale','log')

% ISE gradient magnitude:
figure;
imagesc(bPlot, omegaPlot, gradMag);
axis xy
xlabel('$\log_{10}(\hat{b})$', 'Interpreter','latex')
ylabel('$\log_{10}(\omega_c)$', 'Interpreter','latex')
title('$|\nabla ISE|$', 'Interpreter','latex')
colorbar
set(gca,'ColorScale','log')
clim([1e-4 1])

figure;
imagesc(bPlot, omegaPlot, gradRatio);
axis xy
xlabel('$\log_{10}(\hat{b})$', 'Interpreter','latex')
ylabel('$\log_{10}(\omega_c)$', 'Interpreter','latex')
title('$|dISE/d\hat{b}| / |dISE/d\omega_c|$', 'Interpreter','latex')
colorbar
set(gca,'ColorScale','log')
clim([1e-4 1])

figure;
imagesc(bPlot, omegaPlot, kappaPlot);
axis xy
xlabel('$\log_{10}(\hat{b})$', 'Interpreter','latex')
ylabel('$\log_{10}(\omega_c)$', 'Interpreter','latex')
title('$\kappa = |\lambda_{\max}|/|\lambda_{\min}|$', 'Interpreter','latex')
colorbar
set(gca,'ColorScale','log')
clim([1e-3 1e3])
rectangle('Position',[ISEBopt, ISEWopt, dB/2, dW/2], ...
          'EdgeColor','r', 'LineWidth',2);