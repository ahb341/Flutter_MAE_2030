% Sruti Vutukury, Aaron Brown
% MAE 2030, Spring 2019, Extra Credit Project
% Cornell University
%
% Dynamic Aeroelasticity
clear; clc;
%% Inputs
% Geometry
p.b = 10; p.c = 1; p.S = p.b*p.c; p.e = 0.1;

% Properties
p.m = 1; p.g = 9.81; p.Kh = 100; p.Ka = 1000; p.Ch = 0; p.Ca = 0;
p.My = 1; p.Ia = 1; p.Sa = 0.1;

% Aerodynamics
p.L = 0.5; p.CLa = 2*pi;
%p.q = 0.1
p.q = 0.1; n = 0.001;
A = p.m*p.Ia-p.Sa^2;
B = p.m*(p.Ka-p.q*p.S*p.e*p.CLa)+p.Kh*p.Ia-p.Sa*p.q*p.S*p.CLa;
C = p.Kh*(p.Ka-p.q*p.S*p.e*p.CLa);
% while (B^2-4*A*C > 10^-1)
%     B^2-4*A*C
%     p.q = p.q + n;
%     A = p.m*p.Ia-p.Sa^2;
%     B = p.m*(p.Ka-p.q*p.S*p.e*p.CLa)+p.Kh*p.Ia-p.Sa*p.q*p.S*p.CLa;
%     C = p.Kh*(p.Ka-p.q*p.S*p.e*p.CLa);
% end
% p.q


%% Solve
tstart = 0; tend = 3; npointspers = 100;
ntimes = tend*npointspers+1; % total number of time points
t = linspace(tstart,tend,ntimes);

h0 = 0; hd0 = 0; al0 = 5*pi/180; ald0 = -1;
z0 = [h0;hd0;al0;ald0];

% ODE45
small = 1e-7;
options = odeset('RelTol', small, 'AbsTol', small);
f = @(t,z) detailedFlutterRHS(t,z,p);
[t,z] = ode45(f, t, z0, options);

h = z(:,1); hd = z(:,2); al = z(:,3); ald = z(:,4);
minh = min(h); maxh = max(h);
minal = min(al); maxal = max(al);

%% Plot
foil_str = 'naca0012.xlsx';
%graph(foil_str,t,al,h,p);
animate(foil_str,t,al,h,p);

%% More Detailed Flutter RHS Function
function zdot = detailedFlutterRHS(t,z,p)
h = z(1); hd = z(2);
al = z(3); ald = z(4);

m = p.m; Kh = p.Kh; Ch = p.Ch;
My = p.My; Ka = p.Ka; Ca = p.Ca; Ia = p.Ia;
%detailed:
q = p.q; CLa = p.CLa; S = p.S; e = p.e; Sa = p.Sa;

L = q*S*CLa*al; My = L*e;

aldd = (My-Ka*al+(Sa/m)*(Kh*h+L*sin(al)))/(Ia-Sa^2/m);
hdd = (-1/m)*(Sa*aldd+Kh*h+L*sin(al));
% hdd = (-1/m)*(Kh*h+Ch*hd+q*S*CLa*al);
% aldd = (1/Ia)*(q*S*e*CLa*al-Ka*al-Ca*ald);

zdot = [hd;hdd;ald;aldd];
end

%% Graph
function graph(foil_str,t,al,h,p)
minh = min(h); maxh = max(h); minal = min(al); maxal = max(al);
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
plot(al,h,'b')
title('h vs alpha'); xlabel('alpha'); ylabel('h');
grid on; axis([minal maxal minh maxh]);
end

%% Animate
function animate(foil_str,t,al,h,p)
minh = min(h); maxh = max(h); minal = min(al); maxal = max(al);
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
    R = [cos(-al(i)), -sin(-al(i)); sin(-al(i)), cos(-al(i))];
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
    plot(al(i),h(i),'b.')
    title('h vs alpha'); xlabel('alpha'); ylabel('h');
    grid on; axis([minal maxal minh maxh]);
    hold off
    
    pause(.01) %uncomment to animate
    
    if ~ishghandle(fig)
        break
    end
end
end

