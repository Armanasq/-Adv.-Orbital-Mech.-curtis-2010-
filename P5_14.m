clc
clear all
close all
disp('A.Asgharpoor      ID:830398023    email: A.Asgharpoor@ut.ac.ir')
disp('FNST')
disp('===================================================================================')
disp('Adv. Orbital Mech.')
disp('P5.14')
fprintf('\n')

%Problem 4.15 Adv. Orbital Mech. [curtis, 2010]
%Using the Gibbs method to calculate the state vector of the satellite at the central observation time.

r_e     = 6378;
w_e     = 7.292115e-5;
f       = 0.0033528;
mu      = 398600;

phi     = -20;
alt     = 0.5;
t       = [0           2               4];
theta   = [60        60.5014     61.0027];
Az      = [165.932   145.970     2.40973];
El      = [8.81952   44.2734     20.7594];
Rn      = [1212.48   410.596     726.464];

R       = zeros(3,3);
V       = zeros(3,3);

for i=1:3    
[r,v] = rv_from_obs(Rn(i), Az(i), El(i) , 0, 0 , 0, alt , theta(i), phi);
R(i,:)= [r(1) r(2) r(3)] ;
V(i,:)= [v(1) v(2) v(3)] ;
end

r1=R(1,:);
r2=R(2,:);
r3=R(3,:);

V2=gibbs(r1, r2, r3);

%% Result
disp('P5.14 Using the Gibbs method to calculate the state vector of the satellite at the central observation time.') 

fprintf('\n  Time    Local Sidereal Time          Azimuth        Elevation angle     Range   ')
for     i = 1:3
    fprintf('\n %5.0f     %15.4f       %15.5f  %15.5f %15.3f ',t(i),theta(i), Az(i), El(i), Rn(i))
    
end



fprintf('\n')
fprintf('\n\n                               Results')
fprintf('\n             =================================================')
fprintf('\n\n  Time       Geocentric Position Vector. (Km) ')
fprintf('\n')

for i=1:3
fprintf('\n %5.1f       [%g     %g      %g] ', t(i) , R(i,1), R(i,2), R(i,3))
end

fprintf('\n')
    R2  = norm(r2);
    v2  = norm(V2);
    
fprintf('\n ')
formatspec = '	R2 is %4.2f km';
fprintf(formatspec, R2)
fprintf('\n')

formatspec = '	V2 is %4.4f km/s';
fprintf(formatspec, v2)
fprintf('\n')
