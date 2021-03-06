% Sruti Vutukury, Aaron Brown
% MAE 2030, Spring 2019, Extra Credit Project
% Cornell University
%
% Dynamic Aeroelasticity
clear; clc;
clf
close all

%% Inputs
% Geometry
p.b = 10; p.c = 1; p.S = p.b*p.c; p.e_ac = 0.1; p.e_cg = 0.1;

% Properties
p.m = 0.5; p.g = 9.81; p.Kh = 100; p.Ka = 1000; p.Ch = 0; p.Ca = 0;
p.Ia = 1; p.Sa = p.m*p.e_cg;

% Aerodynamics
p.CLa = 2*pi; p.rho = 1.225; p.v = 1;

% n = 0.001;
% A = p.m*p.Ia-p.Sa^2;
% B = p.m*(p.Ka-p.q*p.S*p.e*p.CLa)+p.Kh*p.Ia-p.Sa*p.q*p.S*p.CLa;
% C = p.Kh*(p.Ka-p.q*p.S*p.e*p.CLa);
% while (B^2-4*A*C > 10^-1)
%     B^2-4*A*C
%     p.q = p.q + n;
%     A = p.m*p.Ia-p.Sa^2;
%     B = p.m*(p.Ka-p.q*p.S*p.e*p.CLa)+p.Kh*p.Ia-p.Sa*p.q*p.S*p.CLa;
%     C = p.Kh*(p.Ka-p.q*p.S*p.e*p.CLa);
% end
% p.q


%% Solve
tstart = 0; tend = 2; npointspers = 100;
tstart = 0; tend = 4; npointspers = 200;
ntimes = tend*npointspers+1; % total number of time points
t = linspace(tstart,tend,ntimes);

h0 = 0.1; hd0 = 0.1; alc0 = 5*pi/180; alcd0 = -1;
z0 = getZ0(h0,hd0,alc0,alcd0,p);

% ODE45
small = 1e-7;
options = odeset('RelTol', small, 'AbsTol', small);
f = @(t,z) nonLinearFlutterRHS(t,z,p);
[t,z] = ode45(f, t, z0, options);

h = z(:,1); hd = z(:,2); al = z(:,3); ald = z(:,4);
alc = z(:,1); alcd = z(:,2);

minh = min(h); maxh = max(h);
minal = min(al); maxal = max(al);
minalc = min(alc); maxalc = max(alc);

%% Plot
foil_str = 'naca0012.xlsx';
graph2(foil_str,t,p,h,al,alc);
%animate2(foil_str,t,p,h,al,alc);
%% Functions
function z0 = getZ0(h0,hd0,alc0,alcd0,p)
v = p.v; S = p.S; CLa = p.CLa; m = p.m; Ia = p.Ia; Sa = p.Sa;
e_ac = p.e_ac; g = p.g; e_cg = p.e_cg; Kh = p.Kh; Ch = p.Ch;
Ka = p.Ka; Ca = p.Ca; rho = p.rho;

al0 = alc0 - atan(-hd0/v); %al0 = 5*pi/180; 

% Calculate hdd0
q = (1/2)*rho*sqrt(v^2+hd0^2); L = q*S*CLa*al0;
A = (m*Ia/Sa) - Sa; B = -L*((Ia/Sa)+e_ac*cos(alc0));
C = m*g*((Ia/Sa)-e_cg*cos(alc0)); D = -(Ia/Sa)*Kh*h0;
E = -(Ia/Sa)*Ch*hd0; F = Ka*al0; G = Ca*alcd0; % should be ald0...but it is not defined yet
hdd0 = (1/A)*(B+C+D+E+F+G);

ald0 = alcd0 + hdd0/(v*(1+(hd0/v)^2)); %ald0 = -1;

z0 = [h0;hd0;al0;ald0;alc0;alcd0];
end

function zdot = nonLinearFlutterRHS(t,z,p)
h = z(1); hd = z(2);
al = z(3); ald = z(4);
alc = z(5); alcd = z(6);

m = p.m; g = p.g; Kh = p.Kh; Ch = p.Ch; Ka = p.Ka; Ca = p.Ca; Ia = p.Ia;
CLa = p.CLa; S = p.S; e_ac = p.e_ac; Sa = p.Sa; rho = p.rho; v = p.v;
e_cg = p.e_cg;

q = (1/2)*rho*(v^2+hd^2); L = q*S*CLa*al; Mz = L*e_ac*cos(alc);

A = (m*Ia/Sa) - Sa; B = -L*((Ia/Sa)+e_ac*cos(alc));
C = m*g*((Ia/Sa)-e_cg*cos(alc)); D = -(Ia/Sa)*Kh*h;
E = -(Ia/Sa)*Ch*hd; F = Ka*al; G = Ca*ald;

hdd = (1/A)*(B+C+D+E+F+G);
aldd = (Mz-Ka*al-Ca*ald+m*g*e_cg*cos(alc)-Sa*hdd)/Ia;
alcdd = aldd + (2*hd/(v^2*(((hd/v)^2)+1)^2))*(hdd/v)^2;

zdot = [hd;hdd;ald;aldd;alcd;alcdd];
end

%% Graph
function graph2(foil_str,t,p,h,al,alc)
minh = min(h); maxh = max(h); minal = min(al); maxal = max(al);
minalc = min(alc); maxalc = max(alc);

h0 = h(1);
%Initial Airfoil
% afX0 = [-p.c p.c]/2;
% afY0 = [-h0 -h0];
% af0 = [afX0;afY0];
airfoil = xlsread(foil_str);
x_a0 = p.c*airfoil(:,1); y_a0 = p.c*airfoil(:,2)+h0; % x and y coords of airfoil
af0 = [x_a0'; y_a0'];

figure(1);
plot(af0(1,:),af0(2,:)); 
grid on; axis equal;

figure(2);
plot(t,h,'r')
title('h(t)'); xlabel('t'); ylabel('h');
grid on; axis([t(1) t(end) minh maxh]);

figure(3);
plot(t,al,'g')
title('alpha(t)'); xlabel('t'); ylabel('alpha');
grid on; axis([t(1) t(end) minal maxal]);

figure(4);
plot(t,alc,'y')
title('alpha_chord(t)'); xlabel('t'); ylabel('alpha_chord');
grid on; axis([t(1) t(end) minalc maxalc]);

figure(5);
plot(al,h,'b')
title('h vs alpha'); xlabel('alpha'); ylabel('h');
grid on; axis([minal maxal minh maxh]);
end

%% Animate
function animate2(foil_str,t,p,h,al,alc)
minh = min(h); maxh = max(h); minal = min(al); maxal = max(al);
minalc = min(alc); maxalc = max(alc);

h0 = h(1);
%Initial Airfoil
% afX0 = [-p.c p.c]/2;
% afY0 = [-h0 -h0];
% af0 = [afX0;afY0];
airfoil = xlsread(foil_str);
x_a0 = p.c*airfoil(:,1); y_a0 = p.c*airfoil(:,2)+h0; % x and y coords of airfoil
af0 = [x_a0'; y_a0'];

fig = figure(1);

% plot/animate
for i = 1:length(t)
    % Create Rotation Matrix
    R = [cos(-alc(i)), -sin(-alc(i)); sin(-alc(i)), cos(-alc(i))];
    %Determine Airfiol End-Point Locations
    af_new = R*af0;
    
    %Plot Trajetory of Airfoil and Trajetory of Vertices
    subplot(4,1,1)
    plot(af_new(1,:), af_new(2,:),'k') %Plot Airfoil
    %plot(0,h(i),'ro','LineWidth',1) %Plot G
    title('Trajectory of Airfoil'); xlabel('x'); ylabel('y');
    grid on; axis equal;
    
    subplot(4,1,2)
    hold on
    plot(t(i),h(i),'r.')
    title('h(t)'); xlabel('t'); ylabel('h');
    grid on; axis([t(1) t(end) minh maxh]);
    hold off
    
    subplot(4,1,3)
    hold on
    plot(t(i),al(i),'g.')
    title('alpha(t)'); xlabel('t'); ylabel('alpha');
    grid on; axis([t(1) t(end) minal maxal]);
    hold off
    
    subplot(4,1,4)
    hold on
    plot(t(i),alc(i),'b.')
    title('alpha_c vs t'); xlabel('t'); ylabel('alpha_c');
    grid on; axis([t(1) t(end) minalc maxalc]);
    hold off
    
    %pause(.01) %uncomment to animate
    
    if ~ishghandle(fig)
        break
    end
end
end

