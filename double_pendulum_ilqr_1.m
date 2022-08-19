% Taylor Howell
% Research-Manchester
% June 19, 2018
% iLQR for control of a pendulum (simple) 
    % Midpoint discretization
    % "Iterative Linear Quadratic Regulator Design for Nonlinear Biological
    % Movement Systems" derivation

clear; clc; close all
%global tht1 tht2 omg1 omg2 T1 T2
 
% system 
n = 4; % dimensions of system
%m = [1;2]; % mass
m1=2;m2=3;
L1=2;L2=3;
%l = [1;2]; % length
g = 9.8; % gravity
mu = 0.01; % friction coefficient

xs=[0;0;0;0];

fc = @(x,u,dt)[x(3); x(4);
    ((-1).*L2.^2.*m2+(-1).*L1.*L2.*m2.*cos(x(2))).*(L1.^2.*L2.^2.*m1.* ...
  m2+L1.^2.*L2.^2.*m2.^2+(-1).*L1.^2.*L2.^2.*m2.^2.*cos(x(2)).^2).^( ...
  -1).*(u(2)+(-1).*L1.*L2.*m2.*x(3).^2.*sin(x(2))+(-1).*g.*L2.*m2.* ...
  sin(x(1)+x(2)))+L2.^2.*m2.*(L1.^2.*L2.^2.*m1.*m2+L1.^2.*L2.^2.* ...
  m2.^2+(-1).*L1.^2.*L2.^2.*m2.^2.*cos(x(2)).^2).^(-1).*(u(1)+(-1).* ...
  g.*L1.*(m1+m2).*sin(x(1))+L1.*L2.*m2.*x(4).*(2.*x(3)+x(4)).*sin( ...
  x(2))+(-1).*g.*L2.*m2.*sin(x(1)+x(2)));
  
  (L2.^2.*m2+L1.^2.*(m1+m2)+2.*L1.*L2.*m2.*cos(x(2))).*(L1.^2.* ...
  L2.^2.*m1.*m2+L1.^2.*L2.^2.*m2.^2+(-1).*L1.^2.*L2.^2.*m2.^2.*cos( ...
  x(2)).^2).^(-1).*(u(2)+(-1).*L1.*L2.*m2.*x(3).^2.*sin(x(2))+(-1).* ...
  g.*L2.*m2.*sin(x(1)+x(2)))+((-1).*L2.^2.*m2+(-1).*L1.*L2.*m2.*cos( ...
  x(2))).*(L1.^2.*L2.^2.*m1.*m2+L1.^2.*L2.^2.*m2.^2+(-1).*L1.^2.* ...
  L2.^2.*m2.^2.*cos(x(2)).^2).^(-1).*(u(1)+(-1).*g.*L1.*(m1+m2).*sin( ...
  x(1))+L1.*L2.*m2.*x(4).*(2.*x(3)+x(4)).*sin(x(2))+(-1).*g.*L2.* ...
  m2.*sin(x(1)+x(2)))];
        
dynamics_midpoint = @(x,u,dt) x + fc(x + fc(x,u)*dt/2,u)*dt;
%{
A_midpoint = @(x,dt) [(1 + g*cos(x(1))*(dt^2)/(2*l)) (dt - mu*(dt^2)/(2*m*l^2));
                      (g*cos(x(1) + x(2)*dt/2)*dt/l - mu*g*cos(x(1))*(dt^2)/(2*m*l^3)) (1 + g*cos(x(1) + x(2)*dt/2)*(dt^2)/(2*l) - mu*dt/(m*l^2) + (mu^2)*(dt^2)/(2*(m^2)*l^4))];

B_midpoint = @(x,dt) [(dt^2)/(2*m*l^2); 
                      (-mu*(dt^2)/(2*(m^2)*l^4) + dt/(m*l^2))];

%}

Amat=@(x,dt) [0,0,1,0;0,0,0,1;(-1).*g.*L2.*m2.*(L2.^2.*m2+L1.*L2.*m2.*cos(x(2)) ...
  ).*((-1).*L1.^2.*L2.^2.*m1.*m2+(-1).*L1.^2.*L2.^2.*m2.^2+L1.^2.* ...
  L2.^2.*m2.^2.*cos(x(2)).^2).^(-1).*cos(x(1)+x(2))+(-1).*L2.^2.* ...
  m2.*((-1).*L1.^2.*L2.^2.*m1.*m2+(-1).*L1.^2.*L2.^2.*m2.^2+L1.^2.* ...
  L2.^2.*m2.^2.*cos(x(2)).^2).^(-1).*((-1).*g.*L1.*(m1+m2).*cos( ...
  x(1))+(-1).*g.*L2.*m2.*cos(x(1)+x(2))),(-1).*L2.^3.*m2.^2.*((-1).* ...
  L1.^2.*L2.^2.*m1.*m2+(-1).*L1.^2.*L2.^2.*m2.^2+L1.^2.*L2.^2.* ...
  m2.^2.*cos(x(2)).^2).^(-1).*(L1.*x(4).*(2.*x(3)+x(4)).*cos(x(2))+( ...
  -1).*g.*cos(x(1)+x(2)))+(-1).*L2.*m2.*(L2.^2.*m2+L1.*L2.*m2.*cos( ...
  x(2))).*((-1).*L1.^2.*L2.^2.*m1.*m2+(-1).*L1.^2.*L2.^2.*m2.^2+ ...
  L1.^2.*L2.^2.*m2.^2.*cos(x(2)).^2).^(-1).*(L1.*x(3).^2.*cos(x(2))+ ...
  g.*cos(x(1)+x(2))),(-2).*L1.*L2.^3.*m2.^2.*x(4).*((-1).*L1.^2.* ...
  L2.^2.*m1.*m2+(-1).*L1.^2.*L2.^2.*m2.^2+L1.^2.*L2.^2.*m2.^2.*cos( ...
  x(2)).^2).^(-1).*sin(x(2))+(-2).*L1.*L2.*m2.*x(3).*(L2.^2.*m2+L1.* ...
  L2.*m2.*cos(x(2))).*((-1).*L1.^2.*L2.^2.*m1.*m2+(-1).*L1.^2.* ...
  L2.^2.*m2.^2+L1.^2.*L2.^2.*m2.^2.*cos(x(2)).^2).^(-1).*sin(x(2)),( ...
  -2).*L1.*L2.^3.*m2.^2.*(x(3)+x(4)).*((-1).*L1.^2.*L2.^2.*m1.*m2+( ...
  -1).*L1.^2.*L2.^2.*m2.^2+L1.^2.*L2.^2.*m2.^2.*cos(x(2)).^2).^(-1) ...
  .*sin(x(2));(-1).*g.*L2.*m2.*((-1).*L1.^2.*m1+(-1).*L1.^2.*m2+(-1) ...
  .*L2.^2.*m2+(-2).*L1.*L2.*m2.*cos(x(2))).*((-1).*L1.^2.*L2.^2.* ...
  m1.*m2+(-1).*L1.^2.*L2.^2.*m2.^2+L1.^2.*L2.^2.*m2.^2.*cos(x(2)) ...
  .^2).^(-1).*cos(x(1)+x(2))+(L2.^2.*m2+L1.*L2.*m2.*cos(x(2))).*(( ...
  -1).*L1.^2.*L2.^2.*m1.*m2+(-1).*L1.^2.*L2.^2.*m2.^2+L1.^2.*L2.^2.* ...
  m2.^2.*cos(x(2)).^2).^(-1).*((-1).*g.*L1.*(m1+m2).*cos(x(1))+(-1) ...
  .*g.*L2.*m2.*cos(x(1)+x(2))),L2.*m2.*(L2.^2.*m2+L1.*L2.*m2.*cos( ...
  x(2))).*((-1).*L1.^2.*L2.^2.*m1.*m2+(-1).*L1.^2.*L2.^2.*m2.^2+ ...
  L1.^2.*L2.^2.*m2.^2.*cos(x(2)).^2).^(-1).*(L1.*x(4).*(2.*x(3)+ ...
  x(4)).*cos(x(2))+(-1).*g.*cos(x(1)+x(2)))+(-1).*L2.*m2.*((-1).* ...
  L1.^2.*m1+(-1).*L1.^2.*m2+(-1).*L2.^2.*m2+(-2).*L1.*L2.*m2.*cos( ...
  x(2))).*((-1).*L1.^2.*L2.^2.*m1.*m2+(-1).*L1.^2.*L2.^2.*m2.^2+ ...
  L1.^2.*L2.^2.*m2.^2.*cos(x(2)).^2).^(-1).*(L1.*x(3).^2.*cos(x(2))+ ...
  g.*cos(x(1)+x(2))),(-2).*L1.*L2.*m2.*x(3).*((-1).*L1.^2.*m1+(-1).* ...
  L1.^2.*m2+(-1).*L2.^2.*m2+(-2).*L1.*L2.*m2.*cos(x(2))).*((-1).* ...
  L1.^2.*L2.^2.*m1.*m2+(-1).*L1.^2.*L2.^2.*m2.^2+L1.^2.*L2.^2.* ...
  m2.^2.*cos(x(2)).^2).^(-1).*sin(x(2))+2.*L1.*L2.*m2.*x(4).*( ...
  L2.^2.*m2+L1.*L2.*m2.*cos(x(2))).*((-1).*L1.^2.*L2.^2.*m1.*m2+(-1) ...
  .*L1.^2.*L2.^2.*m2.^2+L1.^2.*L2.^2.*m2.^2.*cos(x(2)).^2).^(-1).* ...
  sin(x(2)),2.*L1.*L2.*m2.*(x(3)+x(4)).*(L2.^2.*m2+L1.*L2.*m2.*cos( ...
  x(2))).*((-1).*L1.^2.*L2.^2.*m1.*m2+(-1).*L1.^2.*L2.^2.*m2.^2+ ...
  L1.^2.*L2.^2.*m2.^2.*cos(x(2)).^2).^(-1).*sin(x(2))];



Bmat=@(x,dt) [0,0;0,0;(-1).*L2.^2.*m2.*((-1).*L1.^2.*L2.^2.*m1.*m2+(-1).* ...
  L1.^2.*L2.^2.*m2.^2+L1.^2.*L2.^2.*m2.^2.*cos(x(2)).^2).^(-1),( ...
  L2.^2.*m2+L1.*L2.*m2.*cos(x(2))).*((-1).*L1.^2.*L2.^2.*m1.*m2+(-1) ...
  .*L1.^2.*L2.^2.*m2.^2+L1.^2.*L2.^2.*m2.^2.*cos(x(2)).^2).^(-1);( ...
  L2.^2.*m2+L1.*L2.*m2.*cos(x(2))).*((-1).*L1.^2.*L2.^2.*m1.*m2+(-1) ...
  .*L1.^2.*L2.^2.*m2.^2+L1.^2.*L2.^2.*m2.^2.*cos(x(2)).^2).^(-1),(( ...
  -1).*L1.^2.*m1+(-1).*L1.^2.*m2+(-1).*L2.^2.*m2+(-2).*L1.*L2.*m2.* ...
  cos(x(2))).*((-1).*L1.^2.*L2.^2.*m1.*m2+(-1).*L1.^2.*L2.^2.*m2.^2+ ...
  L1.^2.*L2.^2.*m2.^2.*cos(x(2)).^2).^(-1)];


% initial conditions
tht1=deg2rad(20);tht2=deg2rad(35);omg1=0.5;omg2=0.5;
x0 = [tht1;tht2;omg1;omg2];

% goal
xf = [pi; pi/2;3;3]; 

% costs
%Q = 1e-5*eye(4);
Q =  1e-1*eye(4);
Qf = 15*eye(4);
R = 1e-4*eye(2);

e_dJ = 1e-12;

% simulation
dt = 0.1;
tf = 1;
N = floor(tf/dt);
t = linspace(0,tf,N);
iterations = 100;

% initialization
u = zeros(2,N-1);
x = zeros(n,N);
x_prev = zeros(n,N);
x(:,1) = x0;


% first roll-out
for k = 2:N-1
        x(:,k) = dynamics_midpoint(x(:,k-1),u(:,k-1),dt);
        x(:,k-1);
        fc(x(:,k-1),u(:,k-1),dt);
       % fc() 
        
end

% original cost
J = 0;
for k = 1:N-1
    J = J + (x(:,k) - xf)'*Q*(x(:,k) - xf) + u(:,k)'*R*u(:,k);
end
disp('Original cost:')
J = 0.5*(J + (x(:,N) - xf)'*Qf*(x(:,N) - xf))





%% iterations of iLQR using Todorov derivation
for i = 1:iterations
%% Backward pass
    S = zeros(n,n,N);
    v = zeros(n,1,N);
    %vdim= zeros(n,2);
    S(:,:,N) = Qf;
    v(:,N) = Qf*(x(:,N) - xf);
    K = zeros(2,n,N-1);
    Kv = zeros(2,n,N-1);
   % pause;
    Ku = zeros(2,2,N-1);
    %disp('st1')
    for k = N-1:-1:1
        K(:,:,k) = (Bmat(x(:,k+1),dt)'*S(:,:,k+1)*Bmat(x(:,k+1),dt) + R)\Bmat(x(:,k+1),dt)'*S(:,:,k+1)*Amat(x(:,k+1),dt);
        Kv(:,:,k) = (Bmat(x(:,+1),dt)'*S(:,:,k+1)*Bmat(x(:,k+1),dt) + R)\Bmat(x(:,k+1),dt)';
        Ku(:,:,k) = (Bmat(x(:,k+1),dt)'*S(:,:,k+1)*Bmat(x(:,k+1),dt) + R)\R;
        S(:,:,k) = Amat(x(:,k+1),dt)'*S(:,:,k+1)*(Amat(x(:,k+1),dt) - Bmat(x(:,k+1),dt)*K(:,:,k)) + Q;
        v(:,k) = ((Amat(x(:,k+1),dt) - Bmat(x(:,k+1),dt)*K(:,:,k))'*v(:,k+1) - K(:,:,k)'*R*u(:,k) + Q*x(:,k));
       %v(:,k) = vdim(:,1); 
    end
    %disp('st2')
    % update control, roll out new policy, calculate new cost
    x_prev = x;
    J_prev = J;
    J = Inf;
    alpha = 1;
    iter = 0;

    
    while J > J_prev 
    %while iter <4
        x = zeros(n,N);
        x(:,1) = x0;
            %disp('st3')
        for k = 2:N
            u_(:,k-1) = u(:,k-1) -K(:,:,k-1)*(x(:,k-1) - x_prev(:,k-1)) + alpha*(-Kv(:,:,k-1)*v(:,k) - Ku(:,:,k-1)*u(:,k-1));
            x(:,k) = dynamics_midpoint(x(:,k-1),u_(:,k-1),dt);
        end
            %disp('st4')
        J = 0;
        for k = 1:N-1
            J = J + (x(:,k) - xf)'*Q*(x(:,k) - xf) + u_(:,k)'*R*u_(:,k);
        end
        J = 0.5*(J + (x(:,N) - xf)'*Qf*(x(:,N) - xf))
        alpha = alpha/2
        iter = iter + 1;
    end
    disp('New cost:')
    J
    u = u_;
    
    if abs(J - J_prev) < e_dJ
        disp(strcat('eps criteria met at iteration: ',num2str(i)))
        break
    end

    
end


%% Results
%{
% Animation
r = 1;
figure
U = [u(1) u];
for i = 1:N
    p1 = subplot(1,2,1);
    X = r*cos(x(1,i) - pi/2);
    Y = r*sin(x(1,i) - pi/2);
    plot([0 X],[0 Y],'k-')
    hold on
    plot(X,Y,'ko','MarkerFaceColor', 'k')
    xlabel('pendulum (simple)')
    axis([-1.5*r 1.5*r -1.5*r 1.5*r])
    axis square
    
    p2 = subplot(1,2,2);
    stairs(t(1:i),U(1:i))
    xlabel('t')54
    ylabel('u(t)')
    axis([0 tf min(u) max(u)])
    axis square

    drawnow
    %pause(dt)

    if i ~= N
        cla(p1);
        cla(p2);
    end
     
end

figure
hold on
plot(linspace(0,tf,N),x(1,:))
plot(linspace(0,tf,N),x(2,:))
legend('\theta','\theta_{d}')

%}




