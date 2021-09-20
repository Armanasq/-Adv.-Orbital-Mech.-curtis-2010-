clc
clear all
close all
disp('A.Asgharpoor          email: A.Asgharpoor@ut.ac.ir')
disp('FNST')
disp('===================================================================================')
disp('Adv. Orbital Mech.')
disp('HW.5 P5.29')
fprintf('\n')


% t - vector of observation times t1, t2, t3 (s)
% theta - vector of local sidereal times for t1, t2, t3 (deg)
% R - matrix of site position vectors at t1, t2, t3 (km)
% rho - matrix of direction cosine vectors at t1, t2, t3
% r_old, v_old - the state vector without iterative improvement (km,km/s)
% r, v - the state vector with iterative improvement (km,km/s)

%%
    global mu
    deg     = pi/180;
    mu      = 398600;
    
%% Data 
    t   = [ 0 300 600];
    R   = [ 5582.84     0           3073.90
            5581.50     122.122     3073.90
            5577.50     244.186     3073.90];
    rho = [ 0.846428    0           0.532504
            0.749290    0.463023    0.473470
            0.529447    0.777163    0.340152];

%% Algorithms 5.5 & 5.6
    [r, v, r_old, v_old]    = gauss(rho(1,:), rho(2,:), rho(3,:), R(1,:), R(2,:), R(3,:), t(1), t(2), t(3));

    for i = 1:3
        fprintf('\n         t       =   %g s:',     t(i))
        fprintf('\n         R       =   [%9.4f         %9.4f      %9.4f]',   R(i,1), R(i,2), R(i,3))
        fprintf('\n         rho     =   [%9.4f         %9.4f      %9.4f]',   rho(i,1), rho(i,2), rho(i,3))
        disp(' ')
    end
    
        fprintf('\n\n                               Results')
        fprintf('\n             =================================================')
        fprintf('\n Before Iterative Improvement')
        fprintf('\n         r       = [%g           %g        %g] (km) ', r_old(1), r_old(2), r_old(3))
        fprintf('\n         v       = [%.4f       %9.4f      %9.4f] (km/s)', v_old(1), v_old(2), v_old(3))
        fprintf('\n');
        
        fprintf('\n\n After Iterative Improvement')
        fprintf('\n         r       = [%g           %g        %g] (km) ', r(1), r(2), r(3))
        fprintf('\n         v       = [%.4f        %9.4f      %9.4f] (km/s)', v(1), v(2), v(3))
        fprintf('\n');
        