function [r,v] = rv_from_obs(rn, A, a, d_rho, d_A, d_a, alt, theta, phi)
    %  rho   range from observation [km]
    %  A     azimuth angle of the target [deg]
    %  a     angular elevation of the object [deg]
    %  drho  time rate of change of rho [km]
    %  dA    time rate of change of A [deg]
    %  da    time rate of change of a [deg]
    %  alt   height of the observation location [km]
    %  th    local sidereal time [deg]
    %  phi   latitude of the observer [deg]
    %	r     geocentric position vector
    %   v     geocentric velocity vector

    
%% Constants    
    Re      = 6378;
    wE      = 7.292115e-5;
    f       = 0.0033528;
    omega   = [0 0 wE];
    
%% Convert angles from degrees to radians
    deg     = pi/180;

    
%% Convert angles from degrees to radians
    A       = A *deg;
    dA    = d_A *deg;
    a       = a *deg;
    da    = d_a *deg;
    theta   = theta*deg;
    phi     = phi *deg;
    
%% Equ 5.56
    R       = [(Re/sqrt(1-(2*f - f*f)*sin(phi)^2) + alt)*cos(phi)*cos(theta)     (Re/sqrt(1-(2*f - f*f)*sin(phi)^2) + alt)*cos(phi)*sin(theta)     (Re*(1 - f)^2/sqrt(1-(2*f - f*f)*sin(phi)^2) + alt)*sin(phi)];

%% Equ 5.66
    dR      = cross(omega, R);

%% 5 Equ 5.83a
    dec     = asin(cos(phi)*cos(A)*cos(a) + sin(phi)*sin(a));

%% Equ 5.83b

    if (A > 0) && (A < pi)
        h   = 2*pi - acos((cos(phi)*sin(a) - sin(phi)*cos(A)*cos(a))/cos(dec));
    else
        h   =  acos((cos(phi)*sin(a) - sin(phi)*cos(A)*cos(a))/cos(dec));   
    end

%% Equ 5.83c
    RA      = theta - h;

%% Equ 5.57
    Rho     = [cos(RA)*cos(dec) sin(RA)*cos(dec) sin(dec)];

%% Equ 5.63
    r       = R + rn*Rho;

%% Equ 5.84
    d_Dec   = (-dA*cos(phi)*sin(A)*cos(a) + da*(sin(phi)*cos(a) - cos(phi)*cos(A)*sin(a)))/cos(dec);

%% Equ 5.85
    dRA     = wE + (dA*cos(A)*cos(a) - da*sin(A)*sin(a) + d_Dec*sin(A)*cos(a)*tan(dec)) /(cos(phi)*sin(a) - sin(phi)*cos(A)*cos(a));

%% Equs 5.69 & 5.72
    dRho    = [-dRA*sin(RA)*cos(dec) - d_Dec*cos(RA)*sin(dec) dRA*cos(RA)*cos(dec) - d_Dec*sin(RA)*sin(dec) d_Dec*cos(dec)];

%% Equ 5.64
    v       = dR + d_rho*Rho + rn*dRho;

end 