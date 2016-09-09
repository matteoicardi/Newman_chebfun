
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
% xr1 = chebfun2(@(x,r) x, [Om1,Oms1]);
% rx1 = chebfun2(@(x,r) r, [Om1,Oms1]);
% xr2 = chebfun2(@(x,r) x, [Om2,Oms2]);
% rx2 = chebfun2(@(x,r) r, [Om2,Oms2]);
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
tt = 0:dt:ttot;
VV = zeros(1,nt);
II = zeros(1,nt);
SoC = zeros(1,nt);
SoC2 = zeros(1,nt);


% Butler-Volmer
% cs_sq = chebfun(@(cs) sqrt((limit01(cs)).*(1-(limit01(cs))))*cmax,[-1 2]);
% ce_sq = chebfun(@(ce) sqrt((limit01(ce))*cmax),[-1 2]);
cs_sq = @(y) (sqrt(y.*(1-y))*cmax);
ce_sq = @(y) (sqrt(y*cmax));
%cs_sq = @(y) (y.*(1-y)*cmax^2);
%ce_sq = @(y) (y*cmax);
%i01 = chebfun2(@(ce,cs) K1*(sqrt((cs).*(1-(cs)).*(ce)))*cmax^1.5,[0 1 0 1]);
%i02 = chebfun2(@(ce,cs) K2*(sqrt((cs).*(1-(cs)).*(ce)))*cmax^1.5,[0 1 0 1]);
j1  = @(ce,cs,pe,ps) 3*n1*K1.*cs_sq(cs).*restrict(ce_sq(ce),Om1)...
    .*sinh(alpha*F*(ps-restrict(pe,Om1)-U01(cs))/RT)/R1;
j2  = @(ce,cs,pe,ps) 3*n2*K2.*cs_sq(cs).*restrict(ce_sq(ce),Om2)...
    .*sinh(alpha*F*(ps-restrict(pe,Om2)-U02(cs))/RT)/R2;
j  = @(ce,cs,pe,ps) ...
    3*ce_sq(ce).*cs_sq(cs).*(n1*K1*sinh(alpha*F*(ps-U01(cs)-pe)/RT).*IOm1/R1 ...
    + n2*K2*sinh(alpha*F*(ps-U02(cs)-pe)/RT).*IOm3/R2);

%% initial conditions.
ce0=.1;
cs10=.5;
cs20=.5;
ce = chebfun(@(x) ce0, Om);
Cs1 = chebfun2(@(x,r) cs10, [Om1,Oms1]);
Cs2 = chebfun2(@(x,r) cs20, [Om2,Oms2]);
pe = chebfun(@(x) 0, Om);
cs1 = chebfun(@(t) feval(Cs1,t,R1)',Om1);
cs2 = chebfun(@(t) feval(Cs1,t,R2)',Om2);
ps1 = chebfun(@(x) U01(mean(cs1)), Om1);
ps2 = chebfun(@(x) U02(mean(cs2)), Om2);
% joined functions
cs = chebfun(@(x) cs1(x).*IOm1(x)+cs2(x).*IOm3(x)+.5*IOm2(x),Om); 
ps = chebfun(@(x) ps1(x).*IOm1(x)+ps2(x).*IOm3(x),Om); 

res_ce=chebmatrix(ce);
res_cs1=chebmatrix(cs1);
res_cs2=chebmatrix(cs2);
res_pe=chebmatrix(pe);
res_ps1=chebmatrix(ps1);
res_ps2=chebmatrix(ps2);



%% Define operators
pde_ps1 = chebop(Om1);
%pde_ps1.rbc = 'neumann';
pde_ps2 = chebop(Om2);
%pde_ps2.lbc = 'neumann';

pde_e = chebop(Om);
pde_e.bc='neumann';
pde_pe = chebop(Om);
pde_pe.bc='neumann';
pde_pe.bc=@(x,u) [u(L/2)];
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

%% additional dimension
v1 = 1/(2*R1)*r1.^2;
v2 = 1/(2*R2)*r2.^2;
r=roots(diff(chebfun(@(x) sphbes(0,x),[0 100])));
r=r(r>1e-5);
ef1=chebmatrix(x.^0,Oms1);
ef2=chebmatrix(x.^0,Oms2);
w1 = r1.^2;
w2 = r2.^2;
for i=1:M-1
    ef1(i,:)=chebfun(@(x) sphbes(0,r(i)*x/R1),Oms1);
    ef1(i,:) = ef1(i,:)/sqrt(sum(w1.*ef1(i,:).^2));
    ef2(i,:)=chebfun(@(x) sphbes(0,r(i)*x/R2),Oms2);
    ef2(i,:) = ef2(i,:)/sqrt(sum(w2.*ef2(i,:).^2));
end
ef1(M,:) = chebfun(@(x) 1, Oms1)/sqrt(sum(w1));
ef2(M,:) = chebfun(@(x) 1, Oms2)/sqrt(sum(w2));

% compute mass and stiffness matrix
def1 = diff(ef1);
def2 = diff(ef2);
for ind=1:M
    vcoef1(ind) = sum(w1.*v1.*ef1(ind));
    vcoef2(ind) = sum(w2.*v2.*ef2(ind));
    for kind=1:M
        dcoef1(ind,kind) = sum(w1.*def1(ind).*def1(kind));
        dcoef2(ind,kind) = sum(w2.*def2(ind).*def2(kind));
    end
    ncoef1(ind) = 3*sum(w1.*ef1(ind))/R1;
    ncoef2(ind) = 3*sum(w2.*ef2(ind))/R2;
end
dcoef2(dcoef2/norm(dcoef2)<1e-8)=0;
dcoef1(dcoef1/norm(dcoef1)<1e-8)=0;

dtwarn = max(dt*norm(dcoef1*Ds1,'inf'),dt*norm(dcoef2*Ds2,'inf'));
if (dtwarn>.99)
    disp('WARNING, time step too big, reduce it by ')
    dtwarn
end
m1=chebmatrix(r1,Oms1);
m2=chebmatrix(r1,Oms1);
for i=1:M
    coef1 = sum(w1.*ef1(i)*cs10);
    m1(i,:) = chebfun(@(x) coef1, Om1);
    coef2 = sum(w2.*ef2(i)*cs20);
    m2(i,:) = chebfun(@(x) coef2, Om2);
end

dj1=0;
dj2=0;
