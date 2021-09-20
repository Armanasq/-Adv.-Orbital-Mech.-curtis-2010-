%% Problem 4.13 Adv. Orbital Mech. [curtis, 2010]

clc
clear all
close all
disp('A.Asgharpoor         email: A.Asgharpoor@ut.ac.ir')
disp('FNST')
disp('===================================================================================')
disp('Adv. Orbital Mech.')
disp('P4.13')
fprintf('\n')

%% Calculate DCM
alpha   =300;
beta    =80;
gama    =30;

%% Rotation Matrix
R1  = [1        0          0
       0   cosd(gama)  sind(gama)
       0  -sind(gama)  cosd(gama)];

R2  = [cosd(beta)  0    -sind(beta)
           0       1       0     
       sind(beta)  0     cosd(beta)];

R3  = [cosd(alpha)   sind(alpha)  0
       -sind(alpha)  cosd(alpha)  0
            0            0        1];
       
Q   = R1*(R2*R3);

%% Calculate alpha2
alpha2  = 360+atand(Q(3,1)/(-Q(3,2)));

%% Calculate beta2
beta2   = acosd(Q(3,3));

%% Calculate gama2
gama2   = 360+atand(Q(1,3)/Q(2,3));

%% Result
disp('4.13 Calculate DCM as the yaw-pitch-roll sequence ')
fprintf('\n')
        fprintf('           \x03B1'' = %g\t',alpha2)
        fprintf('\n')
        fprintf('           \x03B2'' = %g\t',beta2)
        fprintf('\n')
        fprintf('           \x03B3'' = %g\t',gama2)
        fprintf('\n')
        
