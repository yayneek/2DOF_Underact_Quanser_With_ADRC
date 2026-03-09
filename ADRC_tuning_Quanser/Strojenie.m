Inverse_Dynamics;
close all;
validSettings = struct('b_hat', {}, 'omega_c', {}, 'ISE', {});
b_hat_range = linspace(0,5,100);
omega_c_range = linspace(-3,2,100);
licznik = 0;

for i = 1:100
    for j = 1:100
        licznik = licznik +1;
            b_hat = 10^(b_hat_range(i));
            omega_c = 10^(omega_c_range(j));
            out = sim('ADRC_Quanser.slx',10);
            data = out.error_fb_ffw.data';
            for k = 1:length(t)
                uchyb(k) = norm(data(:,k));
            end
            ISE = trapz(t, uchyb.^2);
            fprintf('b_hat: %.4f \t omega_c: %.4f \t ISE: %d \t iteracja: %.2f%%\n', b_hat, omega_c, ISE, licznik*100/(length(b_hat_range)*length(omega_c_range)));
            validSettings(end+1).b_hat = b_hat;
            validSettings(end).omega_c = omega_c;
            validSettings(end).ISE = ISE;
    end
end

save("ExtendedOmegaRange", "validSettings");