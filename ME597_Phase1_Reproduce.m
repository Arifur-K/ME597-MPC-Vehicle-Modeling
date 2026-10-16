clear; close all; clc;

%% A. Physical parameters and static axle loads
p.m = 1500;                   
p.Iz = 2500;                  
p.lf = 1.2; p.lr = 1.6;      
p.g = 9.81;                   
p.B = 12; p.C = 1.15;        % simplified tire-shape coefficients
p.Fzf = p.m*p.g*p.lr/(p.lf+p.lr);
p.Fzr = p.m*p.g*p.lf/(p.lf+p.lr);
Ts = 0.02;                   % s: candidate MPC update time
Tref = 0.001;                % s: high-resolution RK4 integration
fprintf('Static axle loads: Fzf = %.3f N; Fzr = %.3f N\n',p.Fzf,p.Fzr);
fprintf('Euler step Ts = %.4f s (%.0f Hz); RK4 step = %.4f s\n',...
        Ts,1/Ts,Tref);

%% B. Prediction tests about a straight-line nominal trajectory
% Physical inputs to nonlinear model: [rear-axle longitudinal force; steering]
% For the linear model the state is a perturbation about
%   xbar(t) = [20,0,0,20*t,0,0]'.
cas(1) = struct('label','Mild','mu',0.80,'deltaDeg',1,'Fbrake',-1200);
cas(2) = struct('label','Low friction / near saturation',...
                'mu',0.35,'deltaDeg',3,'Fbrake',-1000);
for j=1:numel(cas)
    q=cas(j);
    [t,xNL,xEuler,xLin,Ad,Bd] = validationCase(p,q,Ts,Tref);
    labels = {'v_x','v_y','r','e_y'};
    idx = [1 2 3 5];
    Enl = zeros(1,4); Eeuler = zeros(1,4);
    for k=1:4
        ii=idx(k);
        span = max(xNL(:,ii))-min(xNL(:,ii));
        Enl(k) = 100*max(abs(xLin(:,ii)-xNL(:,ii)))/span;
        Eeuler(k) = 100*max(abs(xEuler(:,ii)-xNL(:,ii)))/span;
    end
    fprintf('\n%s (mu=%.2f, delta=%.1f deg, Fx=%g N):\n',...
            q.label,q.mu,q.deltaDeg,q.Fbrake);
    for k=1:4
        fprintf('  %-3s: fixed-linear E=%.2f%%; nonlinear-Euler E=%.2f%%\n',...
                labels{k},Enl(k),Eeuler(k));
    end
    if j==1
        filename='fig1_mild';
    else
        filename='fig2_low_friction';
    end
    plotValidation(t,xNL,xLin,xEuler,q,filename);
    if j==1
        fprintf('  Discrete state-transition matrix Ad (6x6):\n'); disp(Ad);
        fprintf('  Discrete input matrix Bd (6x2):\n'); disp(Bd);
    end
end

%% C. Full nonlinear open-loop curved-path simulation
% Front steering is FIXED throughout this test (no steering feedback).
kap = 0.006;                 % m^-1 (constant curvature)
deltaFix = atan((p.lf+p.lr)*kap); % rad, low-slip kinematic approximation
Tf = 6.5;
x0 = [20;0;0;0;0;0];
[tOpen,xOpen,etaF,etaR,muOpen] = simulateOpenLoop(p,x0,Tf,Ts,kap,deltaFix);
fprintf('\nOpen-loop test: curvature = %.4f 1/m, fixed steering = %.4f deg\n',...
        kap,rad2deg(deltaFix));
fprintf('  Final v_x=%.4f m/s; final e_y=%.4f m\n',...
        xOpen(end,1),xOpen(end,5));
fprintf('  Maximum eta_f=%.4f; maximum eta_r=%.4f\n',...
        max(etaF),max(etaR));

% Compare Euler step sizes against independent 1-ms fourth-order RK4.
[tRK,xRK] = integrateRK4(@(tt,xx) openLoopRHS(tt,xx,p,kap,deltaFix),...
                         x0,Tf,Tref);
for step = [0.01 0.02 0.05]
    [te,xe] = integrateEuler(@(tt,xx) openLoopRHS(tt,xx,p,kap,deltaFix),...
                            x0,Tf,step);
    eyRef = interp1(tRK,xRK(:,5),te,'linear');
    fprintf('  dt = %.3f s: max |e_y(Euler)-e_y(RK4)| = %.5f m\n',...
            step,max(abs(xe(:,5)-eyRef)));
end
plotOpenLoop(tOpen,xOpen,etaF,etaR,muOpen);

%% D. Straightforward curvature-based reference-speed guide
for muTest=[0.35,0.8]
    vLimit = sqrt(muTest*p.g/abs(kap));
    fprintf('Curvature-only speed bound at mu=%.2f: %.3f m/s\n', muTest,vLimit);
end

function [t,xNL,xEuler,xLin,Ad,Bd] = validationCase(p,q,dt,dtFine)
    V=20; Tf=2;
    x0=[V;0;0;0;0;0];
    % rk4 is evaluated at a finer step and sampled at the candidate grid.
    [tFine,xFine] = integrateRK4(@(tt,xx) mildRHS(tt,xx,p,q),x0,Tf,dtFine);
    [t,xEuler] = integrateEuler(@(tt,xx) mildRHS(tt,xx,p,q),x0,Tf,dt);
    xNL = interp1(tFine,xFine,t,'linear');

    % Straight-path six-state Jacobian in coordinates
    % x=[vx,vy,r,s,ey,epsi]' about [V,0,0,Vt,0,0]'.
    Cf=q.mu*p.Fzf*p.B*p.C; Cr=q.mu*p.Fzr*p.B*p.C;
    Az=[ -(Cf+Cr)/(p.m*V), (-p.lf*Cf+p.lr*Cr)/(p.m*V)-V;
         (-p.lf*Cf+p.lr*Cr)/(p.Iz*V),...
          -(p.lf^2*Cf+p.lr^2*Cr)/(p.Iz*V) ];
    bz=[Cf/p.m;p.lf*Cf/p.Iz];
    A=zeros(6,6); B=zeros(6,2);
    A(2:3,2:3)=Az;
    A(4,1)=1;
    A(5,2)=1; A(5,6)=V; A(6,3)=1;
    B(1,1)=1/p.m; B(2:3,2)=bz;
    % Exact ZOH for the linearized continuous-time model via block expm.
    block=expm([A B; zeros(2,8)]*dt);
    Ad=block(1:6,1:6); Bd=block(1:6,7:8);
    xLin=zeros(numel(t),6); d=zeros(6,1);
    for k=1:numel(t)-1
        u=validationInput(t(k),q);
        d=Ad*d+Bd*u;
        nominal=[V;0;0;V*t(k+1);0;0];
        xLin(k+1,:)=(nominal+d).';
    end
    xLin(1,:)=x0.';
end
function dx=mildRHS(t,x,p,q)
    u=validationInput(t,q);
    dx=vehicleRHS(x,u,p,q.mu,0);
end
function u=validationInput(t,q)
    % Defined on half-open intervals, as in the written report:
    % steering steps at 0.5 s; braking is applied 0.8 <= t < 1.8 s.
    delta=0;
    if t>=0.5, delta=deg2rad(q.deltaDeg); end
    Fx=0;
    if (t>=0.8)&&(t<1.8), Fx=q.Fbrake; end
    u=[Fx;delta];
end
function dx=openLoopRHS(t,x,p,kap,delta)
    s=x(4);
    mu=frictionMap(s);
    Fx=0;
    if t>=2 && t<4, Fx=-600; end
    dx=vehicleRHS(x,[Fx;delta],p,mu,kap);
end
function mu=frictionMap(s)
    % A8: this effective friction value applies to BOTH axles at once.
    if s>=45 && s<=90
        mu=0.35;
    else
        mu=0.80;
    end
end
function dx=vehicleRHS(x,u,p,mu,kap)
    vx=x(1); vy=x(2); r=x(3); ey=x(5); ep=x(6);
    Fx=u(1); delta=u(2);
    if vx<=0 || 1-kap*ey<=0
        error('Vehicle model outside assumed domain (vx>0 and 1-kappa*ey>0).');
    end
    if abs(Fx)>mu*p.Fzr+1e-8
        error('Rear longitudinal force exceeds friction-circle limit.');
    end
    alphaF=delta-atan((vy+p.lf*r)/vx);
    alphaR=-atan((vy-p.lr*r)/vx);
    h=@(alpha) sin(p.C*atan(p.B*alpha));
    FyF=mu*p.Fzf*h(alphaF);
    FyR=sqrt(max(0,(mu*p.Fzr)^2-Fx^2))*h(alphaR);
    dvx=r*vy+(Fx-FyF*sin(delta))/p.m;
    dvy=-r*vx+(FyF*cos(delta)+FyR)/p.m;
    dr=(p.lf*FyF*cos(delta)-p.lr*FyR)/p.Iz;
    ds=(vx*cos(ep)-vy*sin(ep))/(1-kap*ey);
    dey=vx*sin(ep)+vy*cos(ep);
    dep=r-kap*ds;
    dx=[dvx;dvy;dr;ds;dey;dep];
end
function [t,x] = integrateEuler(f,x0,Tf,dt)
    N=round(Tf/dt);
    if abs(N*dt-Tf)>1e-10, error('Tf must be an integer multiple of dt.'); end
    t=(0:N)'*dt; x=zeros(N+1,numel(x0)); x(1,:)=x0.';
    for k=1:N
        x(k+1,:)=x(k,:)+dt*f(t(k),x(k,:).').';
    end
end
function [t,x] = integrateRK4(f,x0,Tf,dt)
    N=round(Tf/dt);
    if abs(N*dt-Tf)>1e-10, error('Tf must be an integer multiple of dt.'); end
    t=(0:N)'*dt; x=zeros(N+1,numel(x0)); x(1,:)=x0.';
    for k=1:N
        xx=x(k,:).'; tk=t(k);
        k1=f(tk,xx);
        k2=f(tk+dt/2,xx+dt*k1/2);
        k3=f(tk+dt/2,xx+dt*k2/2);
        k4=f(tk+dt,xx+dt*k3);
        x(k+1,:)=(xx+dt*(k1+2*k2+2*k3+k4)/6).';
    end
end
function [t,x,etaF,etaR,muSeries] = simulateOpenLoop(p,x0,Tf,dt,kap,delta)
    [t,x]=integrateEuler(@(tt,xx) openLoopRHS(tt,xx,p,kap,delta),x0,Tf,dt);
    N=numel(t); etaF=zeros(N,1); etaR=zeros(N,1); muSeries=zeros(N,1);
    for k=1:N
        mu=frictionMap(x(k,4)); muSeries(k)=mu;
        Fx=0; if t(k)>=2 && t(k)<4, Fx=-600; end
        vx=x(k,1); vy=x(k,2); r=x(k,3);
        alphaF=delta-atan((vy+p.lf*r)/vx);
        alphaR=-atan((vy-p.lr*r)/vx);
        h=@(alpha) sin(p.C*atan(p.B*alpha));
        FyF=mu*p.Fzf*h(alphaF);
        FyR=sqrt(max(0,(mu*p.Fzr)^2-Fx^2))*h(alphaR);
        etaF(k)=abs(FyF)/(mu*p.Fzf);
        etaR(k)=hypot(Fx,FyR)/(mu*p.Fzr);
    end
end
function plotValidation(t,nl,lin,euler,q,name)
    f=figure('Color','w','Position',[100 100 710 690]);
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    nexttile; hold on;
    plot(t,nl(:,1),'LineWidth',1.6);
    plot(t,lin(:,1),'--','LineWidth',1.6);
    plot(t,euler(:,1),':','LineWidth',1.6);
    ylabel('Forward speed v_x (m/s)'); grid on; xlim([0 2]);
    title(sprintf('%s: mu = %.2f, steering = %.1f deg', q.label, q.mu, q.deltaDeg), 'Interpreter', 'none');
    nexttile; hold on;
    plot(t,nl(:,5),'LineWidth',1.6);
    plot(t,lin(:,5),'--','LineWidth',1.6);
    plot(t,euler(:,5),':','LineWidth',1.6);
    ylabel('Cross-track error e_y (m)'); xlabel('Time (s)');
    grid on; xlim([0 2]);
    legend({'Nonlinear RK4 reference','Linear ZOH','Nonlinear Euler'},...
           'Location','southoutside','Orientation','horizontal');
    exportgraphics(f,[name '.png'],'Resolution',300);
    exportgraphics(f,[name '.pdf'],'ContentType','vector');
end
function plotOpenLoop(t,x,etaF,etaR,muSeries)
    f=figure('Color','w','Position',[100 100 740 920]);
    tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
    nexttile; plot(t,x(:,1),'LineWidth',1.6); grid on;
    ylabel('v_x (m/s)'); title('Longitudinal response'); xlim([0 t(end)]);
    nexttile; plot(t,x(:,5),'LineWidth',1.6); grid on;
    ylabel('e_y (m)'); title('Path-tracking response'); xlim([0 t(end)]);
    nexttile; hold on;
    plot(t,etaF,'LineWidth',1.6);plot(t,etaR,'LineWidth',1.6);
    yline(1,'k--','HandleVisibility','off'); grid on; ylim([0 1.1]);
    ylabel('Utilization \eta (-)'); title('Tire-force utilization');
    legend({'Front','Rear'},'Location','northeast'); xlim([0 t(end)]);
    nexttile; plot(t,muSeries,'LineWidth',1.6); grid on;
    ylabel('Friction \mu (-)');xlabel('Time (s)');
    title('Spatial terrain condition'); xlim([0 t(end)]);
    exportgraphics(f,'fig3_open_loop.png','Resolution',300);
    exportgraphics(f,'fig3_open_loop.pdf','ContentType','vector');
end