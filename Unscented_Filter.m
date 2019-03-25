%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is written to filter out the readings 
%
% - Prabhjeet Singh Arora
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Generating Synthetic Readings %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = .1;                   % Timestep for discrete sensor reading/output discrete time 1/60 hr ie, 1min
tf = 240*60;                 % Total time 4hr
t = 0:dt:tf;                 % Discrete time scale
m = length(t);               % Length of the time scale
runs = 1;                   % monte carlo runs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% System of Manuevering Target %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = [0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1;0 0 0 0 0 0;0 0 0 0 0 0];      % Dynamic system model
B = [0 0;0 0;0 0;0 0;1 0;0 1];                                                      % Control model
G = [0 0;0 0;1 0;0 1;0 0;0 0];                                                      % Error model
[phi,gamw] = c2d(A,G,dt);                                                           % Error - dynamic model discretization
[phi,gam] = c2d(A,B,dt);                                                            % Control - dynamic model discretization

%%%%%%%%%%%%% Initial State %%%%%%%%%%%%%%%%%
x0 =  [10 15 11.11 8.33 0 0]';
%%%%%%%%%%%%% Control Input %%%%%%%%%%%%%%%%%
u = zeros(m,2);
u((m-1)/4,1) =  10^-2 *(1/6)*(5/18)*(1/6);                  % 10 kmph acceleration at 1hr in x-direction
u((m-1)/2,1) = -10^-2*(1/6)*(5/18)*(1/6);                  % 10 kmph deceleration at 2hr in x-direction
u((m-1)/4,2) =   10^-2*(1/6)*(5/18)*(1/6);                  % 10 kmph acceleration at 1hr in y-direction
u((m-1)/2,2) = - 10^-2*(1/6)*(5/18)*(1/6);                  % 10 kmph deceleration at 2hr in y-direction

%%%%%%%%%%%%%%%%%% Generating truth model and error model and Estimation %%%%%%%%%%%%%%%%%
%%%%% Creating store variables to store generated synthetic states %%%%%
x_true = zeros(6,m);    % To store true states
x_true(:,1) = x0;
x_store = zeros(6,m);   % To store states with process noise
x_store(:,1) = x0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Error Covariances %%%%%
q = 10^(-10);     % Process noise error covariance
r = 10^(1);     % Measurement noise error covariance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Qt = sqrt(q)*randn(2,m);
%%%%% Generating True States and Synthetic "Real" States with Process noise %%%%%
for i =2:m
    c = ode4(@fun,[t(i-1) t(i)],[x_store(:,i-1);u(i-1,:)';0;0],A,B,G);
    x_true(:,i) = c(end,1:6)';
    b = ode4(@fun,[t(i-1) t(i)],[x_store(:,i-1);u(i-1,:)';Qt(:,i-1)],A,B,G);
    x_store(:,i) = b(end,1:6)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Sensor Position at all time t %%%%%
%%% Sensor 1 %%%
X1_x = -10*ones(1,m)*10000;
X1_y =  5*ones(1,m)*10000;
%%% Sensor 2 %%%
X2_x =  5*ones(1,m)*10000;
X2_y = -5*ones(1,m)*10000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Sensor True Reading for True States %%%%%
%%%%% Sensor 1 %%%%%
x1d_x = (X1_x - x_true(1,:)).^2;
x1d_y = (X1_y - x_true(2,:)).^2;
%%%%% Sensor 2 %%%%%
x2d_x = (X2_x - x_true(1,:)).^2;
x2d_y = (X2_y - x_true(2,:)).^2;
%%%%% Two sensors for triangulation %%%%%
%%%%% Sensor output %%%%%
y_true = [(x1d_x + x1d_y).^0.5;(x2d_x + x2d_y).^0.5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Sensor True Reading for Synthetic "Real" States %%%%%
%%%%% Sensor 1 %%%%%
x1d_m_x = (X1_x - x_store(1,:)).^2;
x1d_m_y = (X1_y - x_store(2,:)).^2;
%%%%% Sensor 2 %%%%%
x2d_m_x = (X2_x - x_store(1,:)).^2;
x2d_m_y = (X2_y - x_store(2,:)).^2;
%%%%% Two sensors for triangulation %%%%%
%%%%% Sensor output %%%%%
y_m = [(x1d_m_x+x1d_m_y).^0.5;(x2d_m_x+x2d_m_y).^0.5];
%%%%% Measurement noise implementation %%%%%
y_m = y_m + sqrt(r)*randn(2,m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Figure showing the synthetic position wrt true position %%%%%
figure(1)
measure = plot(x_store(1,:),x_store(2,:),'r',x_true(1,:),x_true(2,:),'b',X1_x,X1_y,'.k',X2_x,X2_y,'.k');
%%%%% Plot properties %%%%%
hold on
grid on
%axis([0 250 0 250])
xlabel('x-position (m)','FontSize',18)
ylabel('y-position (m)','FontSize',18)
title('Model position','FontSize',20);
legend([measure],{'Synthetic "Real" Position','True Position','1st Sensor Position','2nd Sensor Position'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%          Unscented Filter          %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unscented Filter Parameters
alp  =  .1;
beta =  2;
L    =  10;
kap  =  3-L;
lam  =  alp^2*(L+kap)-L;
w0m  =  lam/(L+lam);
w0c  =  lam/(L+lam)+(1-alp^2+beta);
wim  =  1/(2*(L+lam));
yez  =  zeros(1,2*L);

Qu = q*eye(2);
Ru = r*eye(2);

pcovu = diag([1 1 1 1 1 1 diag(Qu)' diag(Ru)']);    
pcovu_store = zeros(6,m);
pcovu_store(:,1) = diag(pcovu(1:6,1:6));

xu0 = [10 15 11.11 8.33 0 0 zeros(1,2) zeros(1,2)]';
xu_store = zeros(6,m);
xu_store(:,1) = xu0(1:6);
xu = xu0;


%Hu = [(xe(1)-X1(1))/sqrt((X1(1)-xe(1))^2+(X2(1)-xe(2))^2) (xe(2)-X2(1))/sqrt((X1(1)-xe(1))^2+(X2(1)-xe(2))^2) 0 0 0 0];


for iu = 2:m
    psquare = chol(pcovu)';
    sigv = [sqrt(L+lam)*psquare -sqrt(L+lam)*psquare];
    xx0 = xu;                           % Mean Value
    xx = sigv + kron(xu,ones(1,2*L));  % Sigma Points
    
    xx_prop = phi*xx(1:6,:) + gamw*xx(7:8,:) + kron(gam*u(iu-1,:)',ones(1,2*L));  %propagate
    xx0_prop = phi*xx0(1:6) + gamw*xx0(7:8) + gam*u(iu-1,:)';  %propagate
    
    x_mean = w0m*xx0_prop + wim*sum(xx_prop,2);
    
    % Covariance after the propagation
    pp0  = w0c*(xx0_prop-x_mean)*(xx0_prop-x_mean)';
    pmat = xx_prop-kron(x_mean,ones(1,2*L));
    pcovu = pp0 + wim*pmat*pmat';
    
    
    yu0 = [((X1_x(iu)-xx0_prop(1))^2+(X1_y(iu)-xx0_prop(2))^2)^0.5;((X2_x(iu)-xx0_prop(1))^2+(X2_y(iu)-xx0_prop(2))^2)^0.5];
    for j = 1:2*L
        yuz(:,j) = [((X1_x(iu)-xx_prop(1,j))^2+(X1_y(iu)-xx_prop(2,j))^2)^0.5;((X2_x(iu)-xx_prop(1,j))^2+(X2_y(iu)-xx_prop(2,j))^2)^0.5];
    end
    
    yu     = w0m*yu0 + wim*sum(yuz,2);
    
    % Calculate pyy
    pyy0   = w0c*(yu-yu0)*(yu-yu0)';
    pyymat = yuz-yu;
    pyy    = pyy0 + wim*pyymat*pyymat';
    
    
    % Calculate pxy
    pxy0 = w0c*(xx0_prop-x_mean)*(yu0-yu)';
    pxy  = pxy0 + wim*pmat*pyymat';

    % Innovations Covarinace
    pvv = pyy + Ru;
    
    
    % Gain and Update
    gain = pxy*inv(pvv);
    pcovu = pcovu-gain*pvv*gain';
    pcovu_store(:,iu) = diag(pcovu);

    xu = x_mean + (gain*(y_m(:,iu)-yu));
    
    xu_store(:,iu) = xu;
    pcovu = [pcovu zeros(6,4);zeros(2,6) Qu zeros(2,2);zeros(2,8) Ru];
    xu = [xu;zeros(4,1)];
    %xu = phi*xu + gam*u(iu,:)';  %propagate
end
sig3_u = pcovu_store.^(0.5)*3;
    
%%%%% Figure showing the synthetic position wrt true position %%%%%
figure(2)
measure2 = plot(t,xu_store(1,:)-x_store(1,:),'.-b',t,sig3_u(1,:),'r',t,-sig3_u(1,:),'r');
%%%%% Plot properties %%%%%
hold on
grid on
axis([0 tf -.5 .5])
xlabel('x-position (km)','FontSize',18);
ylabel('y-position (km)','FontSize',18);
title(title1,'FontSize',20);
legend([measure2],{'Synthetic "Real" Position','True Position'});


function f2 = fun(t,x,A,B,G)
    xm = x(1:6);
    um = x(7:8);
    qm = x(9:10);
    xm_dot = A*xm+B*um+G*qm;
    f2 = [xm_dot;um;qm];
end