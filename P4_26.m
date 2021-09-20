%% Problem 4.15 Adv. Orbital Mech. [curtis, 2010]
clc
clear all
close all
disp('A.Asgharpoor     email: A.Asgharpoor@ut.ac.ir')
disp('FNST')
disp('===================================================================================')
disp('Adv. Orbital Mech.')
disp('P4.26')
fprintf('\n')


%% Calculate r and v precisely 72 hours

syms E_n theta_n
    R   = [-2429.1;     4555.1;     4577.0];
    V   = [-4.7689;     -5.6113;    3.0535];


%% Initial conditions
    mu  = 398600;
    t_0 = 0;
    r_e = 6378;
    J_2 = 1.08263e-3;
    Dt  = 72*3600;
    r   = norm(R);
    v   = norm(V);
    vr  = dot(R,V) / r; 
    
%% Calculate the specific angular momegaentum
    
    H   = cross(R,V);
    h   = norm(H);
    
%% Calculate the inclination

    i   = acos(H(3) / h) * (180/pi);
    
%% Calculate node line vector

    k   = [0,       0,      1];
    N   = cross(k,H);
    n   = norm(N);

%% Calculate the right ascension of the ascending node
    
if N(2) >= 0
        omega = acos(N(1) / n) * (180/pi);
else
        omega = 360 - acos(N(1) / n) * (180/pi);
end
    
%% Calculate the eccentricity
    
    E   = (1 / mu)*((v^2 - (mu / r))*R - vr*V);
    e   = norm(E);
    
%% Calculate the argument of perigee
    
if E(3) >= 0
        w = acos(dot(N,E)/(n*e)) * (180/pi);
else
        w = 360 - acos(dot(N,E)/(n*e)) * (180/pi);
end
    
%% Calculate the true anomegaaly
    
if vr >= 0
        theta = acos(dot(E/e,R/r)) * (180/pi);
else
        theta = 360 - acos(dot(E/e,R/r)) * (180/pi);
end
    
%%
a   = (h^2)/(mu*(1-e^2));
T   = 2*pi*a^(3/2)/(mu^(1/2));
 
t_f = t_0+(72*3600);
 
%% The orbits no.

n       = round(t_f/T);
t_n     = ((t_f/T)-(round(t_f/T)-1))*T; % (s)
M_n     = 2*pi*t_n/T;                   % (rad)

E_n     = double(vpasolve(E_n-0.1*sin(E_n)==M_n,E_n));

theta_n = double(360+(vpasolve(tand(theta_n/2)==sqrt((1+e)/(1-e))*tan(E_n/2),theta_n)));
r_n     = (h^2)/(mu*(1+e*cosd(theta_n)));
r_x     = r_n*[cosd(theta_n); sind(theta_n); 0];

v_n     = (mu/h).*[-sind(theta_n); e+cosd(theta_n); 0];

omega_dot   = -(3*sqrt(mu)*(r_e^2)*J_2*cosd(i))/(2*((1-(e^2))^2)*(a^(7/2))); % rad/s
omega_dot_d = rad2deg(omega_dot);                                            % deg/s
omega_n     = omega+(omega_dot_d*Dt);
w_dot       = -(3*sqrt(mu)*(r_e^2)*J_2*(5/2*((sind(i))^2)-2)/(2*((1-(e^2))^2)*(a^(7/2))));
w_dot_d     = rad2deg(w_dot);
w_n         = w+(w_dot_d*Dt);

Q           =[cosd(omega_n)*cosd(w_n)-sind(omega_n)*sind(w_n)*cosd(i)    -cosd(omega_n)*sind(w_n)-sind(omega_n)*cosd(i)*cosd(w_n)     sind(omega_n)*sind(i)
              sind(omega_n)*cosd(w_n)+cosd(omega_n)*cosd(i)*sind(w_n)    -sind(omega_n)*sind(w_n)+cosd(omega_n)*cosd(i)*cosd(w_n)     -cosd(omega_n)*sind(i)
                        sind(i)*sind(w_n)                                                  sind(i)*cosd(w_n)                                         cosd(i)];
 
r_t         = Q*r_x;
v_t         = Q*v_n;

%% Result

disp('4.26 Calculate r and v precisely 72 hours later')
fprintf('\n Initial Condition \n       r = -2429.1I-4555.1J-4577.0K (km) \n       V = -4.7689I-5.6113J+3.0535K  (km/s)')
fprintf('\n\n')

fprintf(' 72 hrs Later')
fprintf('\n')
fprintf('               R_t = [ %g\n                       %g\n                      %g ] ', r_t(1,1) , r_t(2,1), r_t(3,1))
fprintf('km')
fprintf('\n\n')
fprintf('               V_t = [ %g\n                        %g\n                        %g ]', v_t(1,1) , v_t(2,1), v_t(3,1))
fprintf('  km\\s')
fprintf('\n')
