%% Problem 4.15 Adv. Orbital Mech. [curtis, 2010]
clc
clear all
close all
disp('A.Asgharpoor        email: A.Asgharpoor@ut.ac.ir')
disp('FNST')
disp('===================================================================================')
disp('Adv. Orbital Mech.')
disp('P4.27')
fprintf('\n')


% Dist. between successive ground tracks at the equator

%% intial Condition
    alt     = 180;      % (km)
    i       = 30;       %( degree)
    mu      = 398600;
    J_2     = 1.08263e-3;
    r_e     = 6378;
    r       = r_e+alt;
    T       = 2*pi*(r^(3/2))/sqrt(mu);
    omega_dot   = -(3*sqrt(mu)*J_2*(r_e^2)*cosd(i))/(2*r^(7/2)); % rad/s
    w_e     =(2*pi+(2*pi/365.26))/(24*3600); % rad/s
    d_lambda    =(w_e-omega_dot)*T;
    s       = r_e*d_lambda; % Km

%% Result

disp('4.27 Calculate Dist. between successive ground tracks at the equator')
fprintf('\n                           Initial Condition')
fprintf('\n             =================================================')
fprintf('\n\n                           alt  = 180 km ')
fprintf('\n                             i  = 30\x00B0 ')
fprintf('\n\n')
fprintf('\n\n                               Results')
fprintf('\n             =================================================')
formatSpec = '\n\n\n                         Dist. is %4.2f km';
        fprintf(formatSpec,s)
        fprintf('\n')

