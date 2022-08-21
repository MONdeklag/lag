clear; 
clc; 
close all

%n = 4; % dimensions of system
global m1 m2 L1 L2 g 
m1=2;
m2=3;
L1=2;L2=3;
g = 9.8; % gravity

Nx = 4;
Nu  = 2;
Tf = 1;
dt = 0.01;
Nt = round(Tf/dt)+1;

fc = @(x,u,dt)[x(3); x(4);

  L2.^2.*m2.*(L1.^2.*L2.^2.*m1.*m2+L1.^2.*L2.^2.*m2.^2+(-1).*L1.^2.* ...
  L2.^2.*m2.^2.*cos(x(1)+(-1).*x(2)).^2).^(-1).*(u(1)+(-1).*u(2)+(-1).* ...
  g.*L1.*(m1+m2).*sin(x(1))+(-1).*L1.*L2.*m2.*x(3).*x(4).*sin(x(1)+( ...
  -1).*x(2))+L1.*L2.*m2.*(x(3)+(-1).*x(4)).*x(4).*sin(x(1)+(-1).* ...
  x(2)))+(-1).*L1.*L2.*m2.*cos(x(1)+(-1).*x(2)).*(L1.^2.*L2.^2.*m1.* ...
  m2+L1.^2.*L2.^2.*m2.^2+(-1).*L1.^2.*L2.^2.*m2.^2.*cos(x(1)+(-1).* ...
  x(2)).^2).^(-1).*(u(2)+L1.*L2.*m2.*x(3).*(x(3)+(-1).*x(4)).*sin( ...
  x(1)+(-1).*x(2))+(-1).*L2.*m2.*((-1).*L1.*x(3).*x(4).*sin(x(1)+( ...
  -1).*x(2))+g.*sin(x(2))));


  (-1).*L1.*L2.*m2.*cos(x(1)+(-1).*x(2)).*(L1.^2.*L2.^2.*m1.*m2+ ...
  L1.^2.*L2.^2.*m2.^2+(-1).*L1.^2.*L2.^2.*m2.^2.*cos(x(1)+(-1).* ...
  x(2)).^2).^(-1).*(u(1)+(-1).*u(2)+(-1).*g.*L1.*(m1+m2).*sin(x(1))+(-1) ...
  .*L1.*L2.*m2.*x(3).*x(4).*sin(x(1)+(-1).*x(2))+L1.*L2.*m2.*(x(3)+( ...
  -1).*x(4)).*x(4).*sin(x(1)+(-1).*x(2)))+L1.^2.*(m1+m2).*(L1.^2.* ...
  L2.^2.*m1.*m2+L1.^2.*L2.^2.*m2.^2+(-1).*L1.^2.*L2.^2.*m2.^2.*cos( ...
  x(1)+(-1).*x(2)).^2).^(-1).*(u(2)+L1.*L2.*m2.*x(3).*(x(3)+(-1).* ...
  x(4)).*sin(x(1)+(-1).*x(2))+(-1).*L2.*m2.*((-1).*L1.*x(3).*x(4).* ...
  sin(x(1)+(-1).*x(2))+g.*sin(x(2))))];
        
dynamics_midpoint = @(x,u,dt) x + fc(x + fc(x,u)*dt/2,u)*dt;



% initial conditions
tht1=deg2rad(20);tht2=deg2rad(35);omg1=0.5;omg2=0.5;
x0 = [tht1;tht2;omg1;omg2];

% goal
xf = [pi;pi/3;0;0]; 

% costs
%Q = 1e-5*eye(4);
Q =  1e-5*eye(4);
Qf = 15*eye(4);
R = 5*1e-7*eye(2);
I = eye(4);

e_dJ = 1e-12;

% simulation
%dt = 0.1;
%tf = 1;
%N = floor(tf/dt);
%t = linspace(0,tf,N);
%iterations = 100;

% initialization
u = zeros(Nu,Nt-1);
x = zeros(Nx,Nt);
x_prev = zeros(Nx,Nt);
x(:,1) = x0;



% first roll-out
for k = 2:Nt-1
        x(:,k) = dynamics_midpoint(x(:,k-1),u(:,k-1),dt);
        %x(:,k-1)
        %fc(x(:,k-1) + fc(x(:,k-1),u(:,k-1))*dt/2,u(:,k-1))
        % fc() 
end

% original cost
J = 0;
for k = 1:Nt-1
    J = J + (x(:,k)-xf)'*Q*( x(:,k)-xf) + u(:,k)'*R*u(:,k);
end
disp('Original cost:')
J = 0.5*(J + (x(:,Nt)-xf)'*Qf*(x(:,Nt)-xf))



%%%%%%%%%%%%%%%% ILQR Algorithm  %%%%%%%%%%%%%%%%%%%%
p = ones(Nx,Nt);
P = zeros(Nx,Nx,Nt);
%d = ones(Nu,Nu,Nt-1);
d = ones(Nu,Nt-1);
K = zeros(Nu,Nx,Nt-1);
%pdim = ones(Nx,Nu,Nt);
dJ = 0.0;  % change in cost

xn = zeros(Nx,Nt);
un = zeros(Nu,Nt-1);
% func g(dx,du) is perturbation of val func
% grad- g/ hessian-G of change in value fun
gx = zeros(Nx);
gu = zeros(Nu);
Gxx = zeros(Nx,Nx);
Guu = zeros(Nu,Nu);
Gxu = zeros(Nx,Nu);
Gux = zeros(Nu,Nx);

iter = 0;
while max(abs(d(:))) >  1e-3
    
    iter = iter +  1 

 %%%%% Backward Pass %%%%%
    dJ = 0.0;
    p(:,Nt) = Qf*(x(:,Nt)-xf);     %%% P is vx
    P(:,:,Nt) = Qf;                %%% P is vxx
    mu_reg = 0;
    for k = (Nt-1):-1:1
          %Calculate derivatives of stage cost
           q = Q*( x(:,k)-xf);     % lx
           r = R*u(:,k);           % lu
            
            A = Astate(x(:,k),u(:,k),dt);
            B = Bstate(x(:,k),u(:,k),dt);

           %gradient of change in val fn
            gx = q + A'*p(:,k+1);% gx = dg/dx  
            gu = r + B'*p(:,k+1);% gu = dg/du
    
          %iLQR (Gauss-Newton) version
          %Hessian
             Gxx = Q + A'*(P(:,:,k+1))*A;
             Guu = R + B'*(P(:,:,k+1)+ mu_reg*eye(Nx))*B;
             Gxu = A'*P(:,:,k+1)*B;
             Gux = B'*(P(:,:,k+1) + mu_reg*eye(Nx))*A;     
             
             %beta = 0.1;
             log = issymmetric([Guu]);
             eigv = eig([Guu]);

          if any(eig(Guu)<0)
            mu_reg = mu_reg + 1;
            k = Nt-1;
            disp('regularized')
          end
        %{
              while (log==0) || all(eigv < 0) 
                    Gxx = Gxx + A'*beta*I*A
                    Guu = Guu + B'*beta*I*B
                    Gxu = Gxu + A'*beta*I*B
                    Gux = Gux + B'*beta*I*A
                    beta = 2*beta
                    %display("regularizing G")
                    display(beta)
                    log = issymmetric([Gxx Gxu; Gux Guu]);
                    eigv = eig([Gxx Gxu; Gux Guu]);
              end
         %}
            d(:,k) = Guu\gu;  % feedforward term
            K(:,:,k) = Guu\Gux; % feedback gain term
    
             p(:,k) = gx - K(:,:,k)'*gu + K(:,:,k)'*Guu*d(:,k) - Gxu*d(:,k);
             P(:,:,k) = Gxx + K(:,:,k)'*Guu*K(:,:,k) - Gxu*K(:,:,k) - K(:,:,k)'*Gux;
             dJ = dJ +  gu'*d(:,k);

          
    end
    disp('EOBP')
    
  
%%%% End of Backward Pass %%%%%
      
    %Forward rollout with line search
    xn(:,1) = x(:,1);
    alpha = 1.0;
   for k = 1:(Nt-1)
        un(:,k) = u(:,k) - alpha*d(:,k) - (K(:,:,k)*(xn(:,k)-x(:,k)));
        xn(:,k+1) = dynamics_midpoint(xn(:,k),un(:,k),dt);
    end
    disp('EOFP')
    Jn = 0;
    for k = 1:Nt-1
        Jn = Jn + (xn(:,k)-xf)'*Q*(xn(:,k)-xf) + un(:,k)'*R*un(:,k);
    end
   Jn = 0.5*(Jn + (xn(:,Nt)-xf)'*Qf*(xn(:,Nt)-xf))
    
    
 
    while isnan(Jn) || Jn > (J - 1e-2*alpha*dJ)
        alpha = 0.5*alpha
        for k = 1:(Nt-1)
            un(:,k) = u(:,k) - alpha*d(:,k) - (K(:,:,k)*(xn(:,k)-x(:,k)));
            xn(:,k+1) = dynamics_midpoint(xn(:,k),un(:,k),dt);
            %xn(:,k+1) = dynamics_rk4(xn(:,k),un(k)
        end
        %Jn = cost(xn,un);
        Jn = 0;
        for k = 1:Nt-1
            Jn = Jn + (xn(:,k) - xf)'*Q*(xn(:,k) - xf) + un(:,k)'*R*un(:,k);
        end
     Jn = 0.5*(Jn + (xn(:,Nt) - xf)'*Qf*(xn(:,Nt) - xf))
    end
    
 
    J = Jn;
    x = xn;
    u = un;
   %if iter > 5
    %   break
    %end
  end

   
    



















