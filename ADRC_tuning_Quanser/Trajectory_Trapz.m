function [s, sd, sdd] = Trajectory_Trapz(t0, dt, tk, starting_point, rest_point, setting_time)


t = t0:dt:tk;
tau  = setting_time/4;
T = 4*tau+dt:dt:tk;

t1 = t0:dt:tau;
t2 = tau+dt:dt:3*tau;
t3 = 3*tau+dt:dt:4*tau;

s1 = starting_point + ((rest_point- starting_point)/(setting_time - tau))*((-5*t1.^8)/(2*tau.^7) + (10*t1.^7)/tau^6 - (14*t1.^6)/tau^5 + (7*t1.^5)/tau^4);
s2 = starting_point + ((rest_point - starting_point)/(setting_time - tau))*(t2- tau/2);
s3 = rest_point + ((rest_point - starting_point)/(setting_time - tau))*((-5*(setting_time-t3).^8)/(2*tau.^7) + (10*(setting_time-t3).^7)/tau^6 - (14*(setting_time-t3).^6)/tau^5 + (7*(setting_time-t3).^5)/tau^4)*(-1);

s = [s1 s2 s3 ones(1, length(T))*rest_point];
sd = gradient(s)./gradient(t);
sdd = gradient(sd)./gradient(t);
% 
% figure()
% plot(t, s)
% grid on
% 
% figure()
% plot(t, sd)
% grid on
% 
% figure()
% plot(t, sdd)
% grid on
end