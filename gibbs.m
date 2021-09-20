function [V2] = gibbs(R1,R2,R3)
%% Input
    %  R1, R2, R3 Vectors
    
    % Output:
    %   V2    velocity
    
%% Constants    
    mu      = 398600;
    
%% Vrctors Mag.
    r1      = norm(R1);
    r2      = norm(R2);
    r3      = norm(R3);

%% Cross Products
    c12     = cross(R1,R2);
    c23     = cross(R2,R3);
    c31     = cross(R3,R1);


%% Equ 5.13
    N       = r1*c23 + r2*c31 + r3*c12;

%% Equ 5.14
    D       = c12 + c23 + c31;

%% Equ 5.21
    S       = R1*(r2 - r3) + R2*(r3 - r1) + R3*(r1 - r2);

%% Equ 5.22
    V2      = sqrt(mu/norm(N)/norm(D))*(cross(D,R2)/r2 + S);


end 