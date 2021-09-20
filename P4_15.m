clc
clear all
close all
disp('A.Asgharpoor         email: A.Asgharpoor@ut.ac.ir')
disp('FNST')
disp('===================================================================================')
disp('Adv. Orbital Mech.')
disp('P4.15')
fprintf('\n')

%Problem 4.15 Adv. Orbital Mech. [curtis, 2010]

%Calculate r and v at perigee
e       = 1.5;
u       = 398600;
z_p     = 300;
r_e     = 6378;
i       = 35;
omega   = 130;
o       = omega;
w       = 115;
p_hat   = [1 ;	   0 ;      0];
r_x_hat = [1 ;     0 ;      0];
q_hat   = [0 ;     1;       0];
r_p     = z_p+r_e;
v_x_hat = [0 ;     1;       0];

%% a)Relative to the perifocal referenceframe.
r_x     = r_p*r_x_hat; % in direction of p vector
h       = sqrt(r_p*u*(1+e));
v_p     = h/r_p;
V       = v_p*q_hat;

%% b)Relative to the geocentric equatorial frame.

Q       = [-sind(o)*cosd(i)*sind(w)+cosd(o)*cosd(w)     -sind(o)*cosd(i)*cosd(w)-cosd(o)*sind(w)    sind(o)*sind(i)
            cosd(o)*cosd(i)*sind(w)+sind(o)*cosd(w)      cosd(o)*cosd(i)*cosd(w)-sind(o)*sind(w)   -cosd(o)*sind(i)
                       sind(i)*sind(w)                                  sind(i)*cosd(w)                    cosd(i)];
r_f     = Q*r_x;
v_x     = v_p*v_x_hat;
v       = Q*v_x;

%% Result
disp('4.15 Calculate r and v at perigee ')
fprintf('                 i = 35° \n                 \x03A9 = 130° \n                 \x03C9 = 115°')
fprintf('\n\n')
fprintf(' a) Relative to the perifocal reference frame')
fprintf('\n')
          
fprintf('               Rf = [ %g\n                         %g\n                         %g ] ', r_x(1,1) , r_x(2,1), r_x(3,1))
fprintf('km')
fprintf('\n\n')
fprintf('               V = [ %g\n                       %g\n                       %g       ]', V(1,1) , V(2,1), V(3,1))
fprintf('  km\\s')
fprintf('\n')
fprintf('\n\n')
fprintf(' b) Relative to the geocentric equatorial frame')
fprintf('\n\n')
        
fprintf('               r_x = [ %g\n                         %g\n                         %g ] ', r_f(1,1) , r_f(2,1), r_f(3,1))
fprintf('km')
fprintf('\n\n')
fprintf('               V_p = [ %g\n                         %g\n                         %g ]', v(1,1) , v(2,1), v(3,1))
fprintf('km\\s')
fprintf('\n')

