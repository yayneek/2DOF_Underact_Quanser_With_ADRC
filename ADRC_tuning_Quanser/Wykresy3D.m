
load('ExtendedOmegaRange.mat')

b_hat_range = linspace(0,5,100);
omega_c_range = linspace(-3,2,100);
x = zeros(1, 100);
y = zeros(1, 100);
z_ISE = zeros(1, 100000);

z_ISE_matrix = zeros(100, 100);



for i=1:10000
    z_ISE(i) = validSettings(i).ISE;
end


for i =1:100
    z_ISE_matrix(:,i) = z_ISE(1 + (i-1)*100: i*100);

    x(i) = 10^b_hat_range(i);
    y(i) = 10^omega_c_range(i);
end




[min_val, ind] = min(z_ISE_matrix(:));
min_val_ISE = min_val
[row_ISE, col_ISE] = ind2sub(size(z_ISE_matrix), ind); 

b_hat_ISE = x(col_ISE)
omega_c_ISE = y(row_ISE)

% Rysowanie wykresu
figure()
surf(x, y, z_ISE_matrix)
xlabel('$\hat{b}$', 'Interpreter', 'latex', 'FontSize',16); 
ylabel('$\omega_c$', 'Interpreter', 'latex', 'FontSize',16);
zlabel('ISE', 'Interpreter', 'latex', 'FontSize',16);
title('$ISE(\omega_c, \hat{b})$', 'Interpreter', 'latex', 'FontSize',16);
filename = 'X2Y2_ISE.pdf';
set(gca, 'ColorScale','log');
clim([1e-5 1e1])
colorbar;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'ZScale', 'log');
hold on;
plot3(x(col_ISE), y(row_ISE), min_val, 'rx', 'MarkerSize', 10, 'LineWidth', 5, 'MarkerFaceColor', 'r');
hold off;

