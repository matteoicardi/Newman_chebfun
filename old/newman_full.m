%% newmann.m -- an executable m-file for solving a partial differential equation
% Newman model
% 

clf
clear all
clc

%% Setup preferences for solving the problem.
%opts = pdeset('Eps', 1e-6, 'PDEflag', pdeflag, 'Ylim', [0.5,1]);
%opts = pdeset('Eps', 1e-5, 'Plot','on');%, 'Ylim', [-1,2]);

%options = cheboppref();

% Print information to the command window while solving:
%options.display = 'iter';

% Option for tolerance.
%options.bvpTol = 1e-5;

% Option for damping.
%options.damping = false;

% Specify the discretization to use. Possible options are:
%  'values' (default)
%  'coeffs'
%  A function handle (see 'help cheboppref' for details).
%options.discretization = 'values';

% Option for determining how long each Newton step is shown.
%options.plotting = 1;

chebfuneps 1e-8
chebfun2eps 1e-8
%chebfunpref.setDefaults('eps',1e-8)
chebfunpref.setDefaults('splitting', 1)
%chebfunpref.setDefaults('maxLength',999)
%chebfunpref.setDefaults('extrapolate',1)
%chebfunpref.setDefaults('happinessCheck','Plateau')

%% load parameters
newman_param;

%% Time stepping
dt   = 3e-2*hour;
ttot = 3*hour;
nt   = floor(ttot/dt);


%% Problem set-up
% Create an interval of the space domain...
Om = [0 L];
Om1 = [0 L1];
Om2 = [L2 L];
Oms1 = [0 R1];
Oms2 = [0 R2];

% Construct a chebfun of the space variable on the domain
x = chebfun(@(x) x, Om);
x1 = chebfun(@(x) x, Om1);
x2 = chebfun(@(x) x, Om2);
r1 = chebfun(@(x) x, Oms1);
r2 = chebfun(@(x) x, Oms2);
xr1 = chebfun2(@(x,r) x, [Om1,Oms1]);
rx1 = chebfun2(@(x,r) r, [Om1,Oms1]);
xr2 = chebfun2(@(x,r) x, [Om2,Oms2]);
rx2 = chebfun2(@(x,r) r, [Om2,Oms2]);
eps = L*1e-2;  
% IOm1 = (tanh((lambda*L-x)/eps)/2 +1/2);
% IOm2 = (tanh((x-lambda*L)/eps)/2 - tanh((x-(1-lambda)*L)/eps)/2);
% IOm3 = (tanh((x-(1-lambda)*L)/eps)/2+1/2);
IOm1 = heaviside(x).*heaviside(-x+L1);
IOm2 = heaviside(x-L1).*heaviside(-x+L2);
IOm3 = heaviside(x-L2).*heaviside(-x+L);

k = k1*IOm1+ke*IOm2+k2*IOm3;
kD = kD1*IOm1+kDe*IOm2+kD2*IOm3;
De = De1*IOm1+Dee*IOm2+De2*IOm3;
nn = n1*IOm1+ne*IOm2+n2*IOm3;


% % mapping
% y=chebfun(@(x) atan(x)/pi+1/2,[-inf,inf]);
% yinv=chebfun(@(x) tan(pi*x-pi/2), [0,1]);
% dy = @(y) 1/(pi*(y.^2+1));
% ddy = @(y) -2*y./(pi*(y.^2+1).^2);
% dlny = @(y) 2/((2*atan(y)+pi).*(y.^2+1));
% ddlny = @(y) -4*(pi*y+2*y.*atan(y)+1)./((2*atan(y)+pi).*(y.^2+1)).^2;

%...and specify a sampling of the time domain:
tt = dt:dt:ttot;
VV = zeros(1,nt);
II = zeros(1,nt);
SoC = zeros(1,nt);
SoC2 = zeros(1,nt);

% % overpotential (Doyle Foller Newman 1993)
% cs = chebfun(@(cs) cs, [0,1]);
% U01  = (-4.656 + 88.669*cs.^2 - 401.119*cs.^4 + 342.909*cs.^6 - 462.471*cs.^8 + 433.434*cs.^10) ...
%        ./(-1 +18.933*cs.^2 - 79.532*cs.^4 + 37.311*cs.^6 - 73.083*cs.^8 + 95.96*cs.^10);
% U02  = @(cs) cs;

% overpotential (Ramos)
% U01  = chebfun(@(cs) RT*log(K1*(1-(limit01(cs)))./(limit01(cs)))/F,[-1,2]);
% U02  = chebfun(@(cs) RT*log(K2*(1-(limit01(cs)))./(limit01(cs)))/F,[-1,2]);
 U01  = @(y) -(RT*log(K1*(y)./(1-y))/F);
 U02  = @(y) (RT*log(K2*(1-y)./y)/F);

% Butler-Volmer
% cs_sq = chebfun(@(cs) sqrt((limit01(cs)).*(1-(limit01(cs))))*cmax,[-1 2]);
% ce_sq = chebfun(@(ce) sqrt((limit01(ce))*cmax),[-1 2]);
cs_sq = @(y) (sqrt(y.*(1-y))*cmax);
ce_sq = @(y) (sqrt(y*cmax));
i01 = chebfun2(@(ce,cs) K1*(sqrt((cs).*(1-(cs)).*(ce)))*cmax^1.5,[0 1 0 1]);
i02 = chebfun2(@(ce,cs) K2*(sqrt((cs).*(1-(cs)).*(ce)))*cmax^1.5,[0 1 0 1]);
j1  = @(ce,cs,pe,ps) 3*n1*K1.*cs_sq(cs).*restrict(ce_sq(ce),Om1)...
    .*sinh(alpha*F*(ps-restrict(pe,Om1)-U01(cs))/RT)/R1;
j2  = @(ce,cs,pe,ps) 3*n2*K2.*cs_sq(cs).*restrict(ce_sq(ce),Om2)...
    .*sinh(alpha*F*(ps-restrict(pe,Om2)-U02(cs))/RT)/R2;
j  = @(ce,cs,pe,ps) ...
    3*ce_sq(ce).*cs_sq(cs).*(n1*K1*sinh(alpha*F*(ps-U01(cs)-pe)/RT).*IOm1/R1 ...
    + n2*K2*sinh(alpha*F*(ps-U02(cs)-pe)/RT).*IOm3/R2);

% initial conditions.
ce = chebfun(@(x) .05, Om);
Cs1 = chebfun2(@(x,r) .5, [Om1,Oms1]);
Cs2 = chebfun2(@(x,r) .5, [Om2,Oms2]);
pe = chebfun(@(x) 0, Om);
cs1 = chebfun(@(t) feval(Cs1,t,R1)',Om1);
cs2 = chebfun(@(t) feval(Cs1,t,R2)',Om2);
ps1 = chebfun(@(x) U01(mean(cs1)), Om1);
ps2 = chebfun(@(x) U02(mean(cs2)), Om2);
% joined functions
cs = chebfun(@(x) cs1(x).*IOm1(x)+cs2(x).*IOm3(x)+.5*IOm2(x),Om); 
ps = chebfun(@(x) ps1(x).*IOm1(x)+ps2(x).*IOm3(x),Om); 

res_ce=[ce];
res_cs1=[cs1];
res_cs2=[cs2];
res_pe=[pe];
res_ps1=[ps1];
res_ps2=[ps2];



% Define operators
pde_ps1 = chebop(Om1);
pde_ps1.rbc = 'neumann';
pde_ps2 = chebop(Om2);
pde_ps2.lbc = 'neumann';

pde_e = chebop(Om);
pde_e.bc='neumann';
pde_pe = chebop(Om);
pde_pe.bc='neumann';
%pde_pe.bc=@(x,pe) pe(0);
pde_ce = chebop(Om);
pde_ce.bc='neumann';

% pde_cs1 = chebop2([Om1,Oms1]);
% pde_cs1.lbc =  @(y,Cs) diff(Cs);
% pde_cs1.rbc =  @(y,Cs) diff(Cs);
% pde_cs1.dbc =  0;%@(x,Cs) diff(Cs);
% pde_cs2 = chebop2([Om2,Oms2]);
% pde_cs2.lbc = @(y,Cs) diff(Cs);
% pde_cs2.rbc = @(y,Cs) diff(Cs); 
% pde_cs2.dbc = @(x,Cs) diff(Cs);

pde_cs1 = chebop(Om1);
pde_cs1.lbc =  @(Cs) diff(Cs);
pde_cs1.rbc =  @(Cs) diff(Cs);
pde_cs2 = chebop(Om2);
pde_cs2.lbc = @(Cs) diff(Cs);
pde_cs2.rbc = @(Cs) diff(Cs); 

init=0;
% Time loop and update operators accordingly
for i=1:nt
    disp([num2str(i), ' of ',num2str(nt) , '\n',...
        num2str(tt(i)), ' of ', num2str(ttot)])
%     figure(1)
%     plot(ps1)
%     hold on
    if mod(i,5)==0
        figure(1)
        plot(ps)
        hold on
        plot(pe)
        figure(2)
        plot(ce)
        hold on
        plot(cs)
        figure(3)
        plot(j(ce,cs,pe,ps)./F./cmax)
        hold on
        figure(4)
        hold off
        plot(SoC2(1:i-1),VV(1:i-1))
        drawnow
        %I=0;
    end
    if max(cs1)>.95 || max(cs2)>.95 || max(1-cs1)>.95 || max(1-cs2)>.95
        I = -I;
        init = 0;
        disp('Charge/Discharge finished, invert current')
    end
    pde_ps1.op = @(ps) -sigma1*diff(ps,2)+j1(ce,cs1,pe,ps);
    pde_ps1.lbc = @(ps) sigma1*diff(ps)+I/A1 ;
    if init
        pde_ps1.init = ps1;
    end
    ps1 = pde_ps1 \ 0;
    
%     plot(ps2)
    pde_ps2.op = @(ps) -sigma2*diff(ps,2)+j2(ce,cs2,pe,ps);
    pde_ps2.rbc = @(ps) sigma2*diff(ps)+I/A2;
    if init
        pde_ps2.init = ps2;
    end
    ps2 = pde_ps2 \ 0;
    ps = chebfun(@(x) ps1(x).*IOm1(x)+ps2(x).*IOm3(x),Om);
    
%     pde_e.op= @(ce,pe) [De*diff(ce,2) + (1-t0)*j(ce,cs,pe,ps)/F/n/cmax - ce/dt;...
%         -k*diff(pe,2)-diff(kD*diff(ce)./ce) - j(ce,cs,pe,ps)];
%     ce = pde_e \ [(ce/dt); 0];
%     plot(pe)
    pde_pe.op= @(pe) -k.*diff(pe,2)-kD.*diff(diff(ce)./ce) - j(ce,cs,pe,ps);
    if init
        pde_pe.init = pe;
    end
    pe = pde_pe \0;
    
%     figure(2)
%     plot(ce)
%     hold on
    pde_ce.op= @(ce) De.*diff(ce,2) + ...
        (1-t0).*j(ce,cs,pe,ps)/F/cmax./(1-nn) - ce/dt;
       % (- 15.99*ce.^6 + 53.59*ce.^5 - 70.62*ce.^4 + 46.64*ce.^3 - 16.62*ce.^2 + 3.945*ce + 0.05083)*...
    pde_ce.init = ce;
    ce = pde_ce \(-ce/dt);

    pde_cs1.op= @(cs) Ds1*diff(cs,2) - ...
        (1-t0)*j1(ce,cs1,pe,ps1)/F/cmax/n1 - cs/dt;
       % (- 310*cs.^8 + 1.24e3*cs.^7 - 2.06e3*cs.^6 + 1.83e3*cs.^5 - 943*cs.^4 + 286*cs.^3 - 51*cs.^2 + 5.8*cs + 0.0266)*...
    pde_cs1.init = cs1;
    cs1 = pde_cs1 \(-cs1/dt);

    pde_cs2.op= @(cs) Ds2*diff(cs,2) - ...
        (1-t0)*j2(ce,cs2,pe,ps2)/F/cmax/n2 - cs/dt;
       % (- 310*cs.^8 + 1.24e3*cs.^7 - 2.06e3*cs.^6 + 1.83e3*cs.^5 - 943*cs.^4 + 286*cs.^3 - 51*cs.^2 + 5.8*cs + 0.0266)*...
    pde_cs2.init = cs2;
    cs2 = pde_cs2 \(-cs2/dt);
    
    cs = chebfun(@(x) cs1(x).*IOm1(x)+cs2(x).*IOm3(x)+.5*IOm2(x),Om);

    % QoI
    VV(i) = ps1(0)-ps2(L);
    II(i) = I;
    SoC(i) = mean(cs1);
    SoC2(i) = mean(cs2);
    
    disp([VV(i), II(i), SoC(i), SoC2(i)])

    init = 1;
    
    continue
    
    pde_cs1.op = @(x,y,Cs) diff(Cs,2,1)  + diff(Cs,2,2) ... + 2*diff(Cs,1,2)./y...
        - Cs/dt/Ds1 + Cs1/dt/Ds1;
    pde_cs1.ubc = @(x,Cs) diff(Cs)+ R1*j1(ce,cs1,pe,ps1)/(3*n1*F)/cmax;
    Cs1 = pde_cs1 \ 0;
    
    pde_cs2.op = @(x,y,Cs) diff(Cs,2,1)  + diff(Cs,2,2) ... + 2*diff(Cs,1,2)./y ...
        - Cs/dt/Ds2 + Cs2/dt/Ds2;
    pde_cs2.ubc = @(x,Cs) diff(Cs)+R2*j2(ce,cs2,pe,ps2)/(3*n2*F)/cmax;
    Cs2 = pde_cs2 \ 0;
    cs1 = chebfun(@(t) feval(Cs1,t,R1)',Om1);
    cs2 = chebfun(@(t) feval(Cs1,t,R2)',Om2);
    cs = chebfun(@(x) cs1(x).*IOm1(x)+cs2(x).*IOm3(x),Om);
    

end


    figure(1)
    plot(ps)
    hold on
    plot(pe)
    figure(2)
    plot(ce)
    hold on
    plot(cs)
    figure(3)
    plot(j(ce,cs,pe,ps)./F./cmax)
    hold on
