%% Problem 4.7 Adv. Orbital Mech. [curtis, 2010]

clc
clear all
close all
disp('A.Asgharpoor       email: A.Asgharpoor@ut.ac.ir')
disp('FNST')
disp('===================================================================================')
disp('Adv. Orbital Mech.')
disp('P4.7')
fprintf('\n')

%%
r       = [-6634.2;      -1261.8;        -5530.9];
e       = [-0.40907;     -0.48751;       -0.63640];
u       = 398600;
k       = [0;                0;                 1];

%% Calculate magnitude of r and e vector 
mag_r   = sqrt(dot(r,r));
mag_e   = sqrt(dot(e,e));

%% Equ 4.13a
theta   = 360-acosd((dot(r,e)/(mag_r*mag_e)));

%% Calculate angular momentum
h       = cross(r,e)/norm(cross(r,e));
mag_h   = sqrt(mag_r*u*(1+mag_e*cos(theta)));
vec_h   = mag_h.*h;

%% Calculate inclination 
i       =acosd(dot(vec_h,k)/mag_h);

%% Result
disp('4.7 Calculate the inclination of the orbit')
fprintf('\n       r = -6634.2I -1261.8J -5230.9K (km) \n       e = -0.40907I -0.48751J -0.63640K ')
fprintf('\n\n')
        formatSpec = '                     i is %4.2f degeree';
        fprintf(formatSpec,i)
        fprintf('\n')
