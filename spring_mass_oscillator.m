threeMassSystem();
%threeMassSystem2();
%%  10.3.8
function threeMassSystem()
    p.m  = 1;  p.k = 1;                  
    tend = 2*pi;                         
    npointspers = 100;                    
    ntimes = tend*npointspers+1;
    tarray = linspace(0, tend, ntimes);
    x0 = [1,0,-1]; v0 = [0,0,0];     % Inititial conditions (ICs)
    z0 = [x0;v0]; 

    [tarray zarray]  = myMidpointSolver(@myrhs, tarray, z0, p);
    x = zeros(3,length(tarray)); v = zeros(3,length(tarray));
    for i = 1:length(tarray)
        x(1:3,i) = zarray{i,1}; v(1:3,i) = zarray{i,2};
    end
    plot (x(1,:),x(2,:)); axis 'equal'; shg;
end
function [tarray, zarray] = myMidpointSolver(eqn, tarray, z0, p)
    n = length(tarray);
    zarray = cell(n, 2);
    zarray{1,1} = z0(1,:); zarray{1,2} = z0(2,:);
    for i = 2:n
        h = tarray(i)-tarray(i-1); z = [zarray{i-1,1}; zarray{i-1,2}]; t = tarray(i-1);
        zdot_temp = eqn(t,z,p);
        z_temp = z + zdot_temp*h/2;
        zdot_new = eqn(t+h/2,z_temp,p);
        z_new = z + zdot_new*h;
        zarray{i,1} = z_new(1,:); zarray{i,2} = z_new(2,:);
    end
end
function zdot = myrhs(t,z,p)
    m = p.m; k = p.k;
    x = z(1,:); v = z(2,:);

    xdot = zeros(1,length(x));
    for i = 1:length(x)
      xdot(i) = v(i);
    end
    vdot(1) = (k/m)*((x(2)-x(1))-x(1));
    vdot(2) = (k/m)*(x(3)-2*x(2)+x(1));
    vdot(3) = (k/m)*(-2*x(3)+x(2));
 
    zdot = [xdot;vdot];
end

%% 10.3.7
function threeMassSystem2()
    m = 1; k = 1;
    M = [m 0 0;0 m 0;0 0 m]; K = k*[3 -1 -1;-1 3 -1;-1 -1 3];
    A = M\K;
    [V,D] = eig(A)
    N1 = V(:,1); N2 = V(:,2); N3 = V(:,3);
    lambda1 = D(1,1); lambda2 = D(2,2); lambda3 = D(3,3);
    omegan1 = sqrt(lambda1); omegan2 = sqrt(lambda2); omegan3 = sqrt(lambda3);
    
    t = linspace(0, 2*pi, 100);
    x1 = N1*sin(omegan1*t);
    x2 = N2*sin(omegan2*t);
    x3 = N3*sin(omegan3*t);
    figure
    plot(t,x1);
    figure
    plot(t,x2);
    figure
    plot(t,x3);
    
end