%% newmann_solid_diffusion2.m -- an executable m-file for solving a partial differential equation
% Automatically created in CHEBGUI by user icardi.
% Created on February 17, 2016 at 19:47.

%% Problem description.
% Solving
%   u_t = (0.1*x.^2*u')'/x,
% for x in [0.1,10] and t in [0,100], subject to
%   dirichlet at x = 0.1
% and
%   u'=sin(2*t)*10 at x = 10

%% Problem set-up
% Create an interval of the space domain...
R = 1;
dom = [0 R];
%...and specify a sampling of the time domain:
t = 0:.1:20;

Ds=5e-2;
j = chebfun(@(tt) sin(tt).*.10,[t(1) t(end)]);
M=10;


%% Equations
x = chebfun('x',dom);
v = 1/(2*R)*x.^2;

% Make the right-hand side of the PDE.
pdefun = @(t,x,u) Ds*(diff(u,2)) + Ds*2*diff(u)./(x+1e-5);
%pdefun = @(t,x,u) Ds*diff(x.^2.*diff(u))./(x+1e-5).^2;
L=chebop(dom);
%L.op = @(x,u) Ds*(diff(u,2)) + Ds* 2*diff(u)./(x+1e-5);
L.op = @(x,u) Ds*diff(x.^2.*diff(u))./(x+1e-5).^2;

% Assign boundary conditions.
L.lbc = 'neumann';
L.rbc = 'neumann';
bc.left = 'neumann';
bc.right = @(t,u) diff(u)-j(t);

u0 = v*j(0);

% cartesian coordinates
%[ef,ev] = eigs(L,M);
w=x.^0;
% spherical coordinates
r=roots(diff(chebfun(@(x) sphbes(0,x),[0 100])));
r=r(r>1e-5);
ef=chebmatrix(x.^0,dom);
w = x.^2;
for i=1:M
    ef(i,:)=chebfun(@(x) sphbes(0,r(i)*x/R),dom);
    ef(i,:) = ef(i,:)/sqrt(sum(w.*ef(i,:).^2));
end
ef(M,:)=x.^0;
ef(M,:) = ef(M,:)/sqrt(sum(w.*ef(M,:).^2));

% compute mass and stiffness matrix
def = diff(ef);
for i=1:M
    vcoef(i) = sum(w.*v.*ef(i));
    for k=1:M
        dcoef(i,k) = sum(w.*def(i).*def(k));
    end
    ncoef(i) = sum(w.*ef(i))/R;
end
% spherical only
ncoef = 3*ncoef;
dcoef(dcoef/norm(dcoef)<1e-8)=0;

%% Call pde23t to solve the problem.
opts = pdeset('Eps', 1e-5, 'Ylim', [-3,3]);
[t, u] = pde23t(pdefun, t, u0, bc, opts);

%% Plot the solution.
%waterfall(u, t)
%xlabel('x'), ylabel('t')

%% compare with spectral method
tspan = [t(1),t(end)];
m=ones(M,length(t));
mref=ones(M,length(t));

for i=1:length(t)
    for k=1:M
        mref(k,i)=sum(w.*(u(:,i)-v*j(t(i))).*ef(k));
    end
end
m(:,1)=mref(:,1);

myfun = @(t,y) eig_sph_diff(t,y,{ncoef,vcoef,dcoef,j(t),feval(diff(j),t),Ds});

%% solve
sol = ode23t(myfun, tspan, m(:,1)');
m = deval(sol,t);

um=u*0;
for i=1:length(t)
    uu = 0;
    for k=1:M
        uu = uu+m(k,i)*ef{k};
    end
    um(:,i)=uu+v*j(t(i));
end

%% plot
figure(1)
hold off
err=(mref'-m');
plot(t,err,'-')

figure(2)
hold off
waterfall(u,t)
hold on
waterfall(um,t)