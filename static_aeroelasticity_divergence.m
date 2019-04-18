%% Inputs/Constants
l = 64.44; c = 14.63; S = c*l; % wing span (m); wing chord (m); % wing area (m^2)
e = .13*c; % distance from aerodynamic center to elastic axis (m)
t = c/10; % airfoil thickness
J = (1/12)*c*t*(c^2+t^2); % assuming polar moment of inertia of rectangular cross section
G = 24*10^9; % Aluminum (Pa)
K_alpha = (pi/2)^2*G*J/l; % torsional spring stiffness (assuming = JG)

rho = 1.225; % air density
U = 224; % velocity (m/s)
q = (1/2)*rho*U^2;

alpha0 = pi/18; alphae = 0; % initial AoA (rad); twist due to spring (rad) 
alpha = alpha0 + alphae; % angle of attack (rad)

%% Flat Plate
CLa = 2*pi; % dCL/dalpha for a flat plate
C_MAC0 = 0; CL0 = 0; % flat plate

% CL = CL0 + CLa*alpha; % coefficient of lift
% L = CL*q*S; % lift (N)
% C_MAC = C_MAC0; % pitching moment coefficient about aerodynamic center
% M_AC = C_MAC*q*S*c; % moment about aerodynamic center
% My = M_AC + L*e - K_alpha*alphae; % moment about elastic axis or center
[UD, qD, alphae] = divergence(K_alpha,S,e,q,CLa,rho,alpha0);
fprintf('Divergence Velocity: %f m/s\n', UD);
fprintf('Divergence Dynamic Pressure: %f Pa\n', qD);
fprintf('Divergence Angle of Twist: %f deg\n', alphae*180/pi);

%% Functions
function [UD, qD, alphae] = divergence(K_alpha,S,e,q,CLa,rho,alpha0)
% solving for alphae when My = 0 (assuming C_MAC0 = 0 for simplicity)
% alphae = (q*S/K_alpha)*(e*CLdot*alpha0)/(1-q*(S*e/K_alpha)*CLdot);

% divergence condition: alphae -> inf as denominator -> 0
% 1-q*(S*e/K_alpha)*CLdot=0

% divergence dynamic pressure (solving divergence condition for q)
% note that only for e > 0 will divergence occur
qD = K_alpha/(S*e*CLa);
UD = sqrt(2*qD/rho); % divergence velocity (m/s)

% condensed equation for alphae
alphae = (q/qD)*alpha0/(1-q/qD);
end

function [a, b] = control_reversal()
end