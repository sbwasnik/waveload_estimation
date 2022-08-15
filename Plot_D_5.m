clear all;
close all;
clc;

% 0 for extrapolated airy wave theory
% 1 for wheeler's stretching
% 2 for both
method = 2;

%%% Given Data
H = 18;
a = H/2;
d = 85;
T = 14;
D = 5;
L = get_wavelength(d,T);
k = 2*pi/L;
w = 2*pi/T;
x = 0;
C_d = 0.7;
C_m = 2;
t_max = 20;
dt = linspace(0,t_max,1000);
Force_a = zeros(1,numel(dt));
Force_w = Force_a;
numel_sec = 1000;
f_airy = zeros(numel_sec, numel(dt));
f_wheeler = f_airy;
Z = f_airy; 

for i = 1:numel(dt)
    t = dt(i);
    eta = a*cos(k*x - w*t);
    z = linspace(eta,-d,1000);
    Z(:,i) = z;
    f = zeros(1,numel(z));

    if method == 0  % extrapolated airy wave theory
        for j = 1:numel(z)
            if eta > 0 && z(j) > 0
                f(j) = get_force(C_d,C_m,D,H,L,0,d,x,T,t);
            else
                f(j) = get_force(C_d,C_m,D,H,L,z(j),d,x,T,t);
            end       
        end
        f_airy(:,i) = f;
        Force_a(i) = trapz(z,f); 

    elseif method == 1   % wheeler's stretching
        z_p = (z-eta)*(d/(d+eta));
        for j= 1:numel(z_p)
            f(j) = get_force(C_d,C_m,D,H,L,z_p(j),d,x,T,t);
        end
        f_wheeler(:,i) = f;
    
        Force_w(i) = trapz(z_p,f); 

    elseif method == 2  % plot both
        %%% for airy
        for j = 1:numel(z)
            if eta > 0 && z(j) > 0
                f(j) = get_force(C_d,C_m,D,H,L,0,d,x,T,t);
            else
                f(j) = get_force(C_d,C_m,D,H,L,z(j),d,x,T,t);
            end       
        end
        f_airy(:,i) = f;
        Force_a(i) = trapz(z,f);
        
        %%% for wheeler
        z_p = (z-eta)*(d/(d+eta));
        for j= 1:numel(z_p)
            f(j) = get_force(C_d,C_m,D,H,L,z_p(j),d,x,T,t);
        end
        f_wheeler(:,i) = f;
        Force_w(i) = trapz(z_p,f); 
    end
end

%%% Graphing
if method == 0
    figure;
    hold on;
    plot(dt, Force_a, '.','LineWidth',2);
    title('Estimation of waveloads using Extapolated Airy Wave Theory');
    xlabel('time (s)')
    ylabel('Force (N)')
    grid on;
    hold off;
elseif method == 1
    figure;
    hold on;
    plot(dt, Force_w, '.','LineWidth',2);
    title('Estimation of waveloads using Wheeler''s stretching');
    xlabel('time (s)')
    ylabel('Force (N)')
    grid on;
    hold off;
elseif method == 2
    figure;
    hold on;
    plot(dt, Force_a, '.','LineWidth',2);
    plot(dt, Force_w, '.','LineWidth',2);
    legend('Extrapolated Airy Wave Theory','Wheeler''s Stretching')
    title('Comparison of Airy and Wheeler Waves');
    xlabel('time (s)')
    ylabel('Force (N)')
    grid on;
    hold off;
end