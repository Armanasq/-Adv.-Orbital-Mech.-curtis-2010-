clc
clear all
close all
disp('A.Asgharpoor          email: A.Asgharpoor@ut.ac.ir')
disp('FNST')
disp('===================================================================================')
disp('Adv. Orbital Mech.')
disp('HW.5 P5.20')
fprintf('\n')

%%
% phi - latitude of site (deg)
% t - vector of observation times t1, t2, t3 (s)
% theta - vector of local sidereal times for t1, t2, t3 (deg)
% R - matrix of site position vectors at t1, t2, t3 (km)
% rho - matrix of direction cosine vectors at t1, t2, t3
% fac1, fac2 - common factors
% r_old, v_old - the state vector without iterative improvement (km,km/s)
% r, v - the state vector with iterative improvement (km,km/s)

%%
global mu
deg     = pi/180;
mu      = 398600;
Re      = 6378;
f       = 1/298.26;

%% Data declaration for Problem 5.19:
alt     = 0;
phi     = 29*deg;
t       = [ 0       60          120];
ra      = [15.0394  25.7539     48.6055]    *deg;
dec     = [20.7487  30.1410     43.8910]    *deg;
theta   = [ 90      90.2507     90.5014]    *deg;

%% Equs 5.56 & 5.57
fac1    = Re/sqrt(1-(2*f - f*f)*sin(phi)^2);
fac2    = (Re*(1-f)^2/sqrt(1-(2*f - f*f)*sin(phi)^2) + alt)*sin(phi);
for i = 1:3
    R(i,1)      = (fac1 + alt)*cos(phi)*cos(theta(i));
    R(i,2)      = (fac1 + alt)*cos(phi)*sin(theta(i));
    R(i,3)      = fac2;
    rho(i,1)    = cos(dec(i))*cos(ra(i));
    rho(i,2)    = cos(dec(i))*sin(ra(i));
    rho(i,3)    = sin(dec(i));
end

%% Algorithms 5.5 & 5.6:
[r, v, r_old, v_old] = gauss(rho(1,:), rho(2,:), rho(3,:),R(1,:), R(2,:), R(3,:), t(1), t(2), t(3));


for i = 1:3
 fprintf('\n %9.4g %11.4f %19.4f %20.4f',t(i), ra(i)/deg, dec(i)/deg, theta(i)/deg)
end
        fprintf('\n\n                               Results')
        fprintf('\n             =================================================')
fprintf('\n')
fprintf('\n Without iterative improvement (P5.19)')
fprintf('\n     r       = [%g        %g           %g] (km)', r_old(1), r_old(2), r_old(3))
fprintf('\n     v       = [%g      %g         %g](km/s)', v_old(1), v_old(2), v_old(3))

fprintf('\n\n With iterative improvement (P5.20)')
fprintf('\n     r       = [%g        %g       %g] (km)', r(1), r(2), r(3))
fprintf('\n     v       = [%g       %g      %g] (km/s)', v(1), v(2), v(3))
fprintf('\n\n')
