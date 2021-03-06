%% newmann.m -- an executable m-file for solving a partial differential equation
% Newman model
% 
%% Setup preferences for solving the problem.
%opts = pdeset('Eps', 1e-6, 'PDEflag', pdeflag, 'Ylim', [0.5,1]);
%opts = pdeset('Eps', 1e-5, 'Plot','on');%, 'Ylim', [-1,2]);

options = cheboppref();

% Print information to the command window while solving:
options.display = 'iter';

% Option for tolerance.
options.bvpTol = 1e-5;

% Option for damping.
options.damping = false;

% Specify the discretization to use. Possible options are:
%  'values' (default)
%  'coeffs'
%  A function handle (see 'help cheboppref' for details).
options.discretization = 'values';

% Option for determining how long each Newton step is shown.
options.plotting = 1;

chebfunpref.setDefaults('eps',1e-5)
chebfunpref.setDefaults('splitting',0)
chebfunpref.setDefaults('maxLength',999)
chebfunpref.setDefaults('extrapolate',1)
%chebfunpref.setDefaults('happinessCheck','Plateau')

%% coefficients
De = 1e-3;
Ds = 1e-12;
sigma = 1;
k  = 1;
kD = 1;
I  = 1;   % applied current
A  = 1;   % electrode surface/cell width*depth
t0 = .5;  % constant
F  = 1;   % faraday constant
n  = .3;  % porosity
alpha  = .5; % constant
K  = 1;   % reaction kinetic constant
T  = 300; % temperature
L  = 1e-3; % cell length
lambda = .4; % relative length of catode/anode
cmax = 1;  % concentration max

dt   = 1e0;
nt   = 100;%floor(ttot/dt);
ttot = dt*nt;%1e2;


%% Problem set-up
% Create an interval of the space domain...
dom = [0 L];
% Construct a chebfun of the space variable on the domain
x = chebfun(@(x) x, dom);
eps = L*1e-2;  
Om1 = (tanh((lambda*L-x)/eps)/2 +1/2);
Om2 = (tanh((x-lambda*L)/eps)/2 - tanh((x-(1-lambda)*L)/eps)/2);
Om3 = (tanh((x-(1-lambda)*L)/eps)/2+1/2);
%...and specify a sampling of the time domain:
t = 0:dt:ttot;

% overpotential
U0  = @(x,cs) (1-cs-Om2(x)).*Om1(x) + (cs).*Om3(x);
eta = @(x,cs,pe,ps) ps-pe-U0(x,cs);

% Butler-Volmer
% sq=chebfun(@(x) abs(x).^alpha, [-inf,inf]);
% sinj=chebfun(@(x) sinh(x), [-100,100],'splitting','on','maxLength',2000);
% %i0 = @(x,ce,cs) feval(chebfun(@(x) K*((cs.*(cmax*(1-Om2(x))-cs).*ce)).^alpha,dom),x);
% i0 = @(x,ce,cs) K*sq((cs.*(cmax*(1-Om2(x))-cs).*ce));
% j  = @(x,ce,cs,pe,ps,ps2) A*(Om1(x)+Om3(x)).*i0(x,ce,cs)...
%     .*sinj(alpha*F*eta(x,cs,pe,ps.*Om1(x)+ps2.*Om3(x)));
% %jx = @(x,ce,cs,pe,ps) feval(j(pe,ps,ce,cs),x);
% %jx = @(x,ce,cs,pe,ps) 1; % EXPLICIT
% jx = @(x,ce,cs,pe,ps,ps2) A*(Om1(x)+Om3(x))*...
%     K.*((cs.*(cmax*(1-Om2(x))-cs).*ce)).* ...
%     sinh(alpha*F*(ps.*Om1(x)+ps2.*Om3(x)-pe-((1-cs-Om2(x)).*Om1(x) + (cs).*Om3(x)))); %IMPLICIT
j=chebfun(@(t) sin(t),[0,20]);

% Make the right-hand side of the PDE.

N = chebop(@(r,C) [De*diff(ce,2) + (1-t0)*jx(x,ce,cs,pe,ps,ps2)/F/n - ce/dt;...
    -jx(x,ce,cs,pe,ps,ps2)/F/n + Ds*diff(cs,2) - cs/dt; ...
    -k*diff(pe,2) - kD*diff(diff(ce)./ce) - jx(x,ce,cs,pe,ps,ps2); ...
    -sigma*diff(ps,2) + jx(x,ce,cs,pe,ps,ps2);...
    -sigma*diff(ps2,2) + jx(x,ce,cs,pe,ps,ps2)], dom);
% pdeflag = [1  1  0  0]; % Zero when a variable is indep of time.

% Assign boundary conditions.
%bc.left = @(t,ce,cs,pe,ps) [diff(ce); diff(cs); diff(pe); ps-1];
%bc.right = @(t,ce,cs,pe,ps) [diff(ce); diff(cs); diff(pe); ps];
N.bc = @(x,ce,cs,pe,ps,ps2) [feval(diff(ce),0); feval(diff(cs),0); ...
    feval(diff(pe),0); ps(0)-1; feval(diff(ps2),(1-lambda)*L); ...
    feval(diff(ce),L); feval(diff(cs),L); ...
    feval(diff(pe),L); feval(diff(ps),lambda*L); ps2(L)];

% initial conditions.
ce0 = chebfun(@(x) .1, dom);
cs0 = .9*Om1 + .1*Om3;
pe0 = chebfun(@(x) 0, dom);
%ps0 = (L-x)/L;
ps0 = chebfun(@(x) 1, dom);
ps20 = chebfun(@(x) 0, dom);

ce1= 1;%sum(ce0+eps);
cs1= 1;%sum(cs0+eps);
ps1= 1;%sum(ps0+eps);
pe1= 1;%sum(pe0+eps);
j1 = 1;%sum(j(pe0,ps0,ce0,cs0,x)+eps);

ce=ce0;
cs=cs0;
pe=pe0;
ps=ps0;
ps2=ps20;
ttot=0;
res_ce=[ce];
res_cs=[cs];
res_pe=[pe];
res_ps=[ps];

for i=1:nt
    %cmax=max(abs(N.op(x,ce,cs,pe,ps)./([ce;cs;pe;ps]+eps)));
    %dt=min(1./cmax(1:2))/10;
    %ttot=ttot+dt;
    disp([i,dt,ttot]);
    N.init = [ce;cs;pe;ps;ps2];
    %% Solve!
    % Call solvebvp to solve the problem.
    % (With the default options, this is equivalent to u = N\rhs.)
    %jj=0; % IMPLICIT
    jj=(j(x,ce,cs,pe,ps,ps2)-jx(x,ce,cs,pe,ps,ps2))*dt; % EXPLICIT
%     figure(1)
%     hold on
%     plot(ce/ce1, 'r', 'LineWidth', 2)
%     axis auto
%     figure(2)
%     hold on
%     plot(cs/cs1, 'g', 'LineWidth', 2)
%     axis auto
%     figure(3)
%     hold on
%     plot(pe/pe1, 'b', 'LineWidth', 2)
%     axis auto
%     figure(4)
%     hold on
%     plot(ps/ps1, 'k', 'LineWidth', 2)
%     axis auto
%     figure(5)
%     hold on
%     %plot(jj/dt/j1, 'y', 'LineWidth', 2)
%     plot(jx(x,ce,cs,pe,ps)/j1, 'y', 'LineWidth', 2)
%     axis auto
%     drawnow
    [ce,cs,pe,ps,ps2] = solvebvp(N, ...
        [-ce-(1-t0)*jj/F/n;...
        jj/F/n-cs;...
        +jj;...
        -jj;...
        -jj]...
        /dt, options);
    res_ce(:,end+1)=ce;
    res_cs(:,end+1)=cs;
    res_pe(:,end+1)=pe;
    res_ps(:,end+1)=ps.*Om1+ps2.*Om3;
end

% figure(1)
% hold on
% plot(ce/ce1, 'r', 'LineWidth', 2)
% axis auto
% figure(2)
% hold on
% plot(cs/cs1, 'g', 'LineWidth', 2)
% axis auto
% figure(3)
% hold on
% plot(pe/pe1, 'b', 'LineWidth', 2)
% axis auto
% figure(4)
% hold on
% plot(ps/ps1, 'k', 'LineWidth', 2)
% axis auto
% figure(5)
% hold on
% plot(jj/dt/j1, 'y', 'LineWidth', 2)
% axis auto
% drawnow
% 

%% Call pde15s to solve the problem.
%[t, ce, cs, pe, ps] = pde15s(pdefun, t, sol0, bc, opts);

% %% Plot the solution components.
% figure
% waterfall(ce, t)
% xlabel('x'), ylabel('t'), title('ce')
% figure
% waterfall(pe, t)
% xlabel('x'), ylabel('t'), title('pe')
% figure
% waterfall(ps, t)
% xlabel('x'), ylabel('t'), title('ps')