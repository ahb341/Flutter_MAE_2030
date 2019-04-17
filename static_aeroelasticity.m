l = 1; c = 1; % wing span (m); wing chord (m)
e = 0.5; % distance from aerodynamic center to elastic axis (m)

rho = 1.225; % air density
U = 1; % velocity (m/s)
q = (1/2)*rho*U^2;
S = 1; % wing area (m^2)


Cmac = 1; Sc = 1;
Mac = Cmac*q*Sc;

alpha = 0; % angle of attack (deg)

% Flat Plate
CLdot = 2*pi; % dCL/dalpha for a flat plate
Cmac0 = 0; CL0 = 0; % flat plate

CL = CL0 + CLdot*alpha; % coefficient of lift
L = CL*q*S;