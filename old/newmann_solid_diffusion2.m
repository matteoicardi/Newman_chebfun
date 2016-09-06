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


Ds=5e-2;
j = @(tt) sin(tt).*.10;
R = 1;

%% Problem set-up
% Create an interval of the space domain...
dom = [0 R];
%...and specify a sampling of the time domain:
t = 0:.1:20;


% Make the right-hand side of the PDE.
pdefun = @(t,x,u) Ds*(diff(u,2));% + 2*diff(u)./(x+1e-15));

% Assign boundary conditions.
bc.left = 'neumann';
bc.right = @(t,u) diff(u)-j(t);

% Construct a chebfun of the space variable on the domain,
x = chebfun(@(x) x, dom);
% and of the initial condition.
u0 = chebfun(1,dom);

%% Setup preferences for solving the problem.
opts = pdeset('Eps', 1e-5, 'Ylim', [-3,3]);

%% Call pde23t to solve the problem.
[t, u] = pde23t(pdefun, t, u0, bc, opts);

%% Plot the solution.
%waterfall(u, t)
%xlabel('x'), ylabel('t')

%% compare with qmom

M = 16;

% compute legendre polynomials
yl=ones(12,1);
for i=1:M
    yl(i)=R^(i)/(i+1);
end
global abl abml
[abl,abml]=chebyshev(floor(M/2),yl);

tspan = [t(1),t(end)];
m=ones(M,length(t));
mref=ones(M,length(t));

for i=1:length(t)
    for k=0:M-1
        mref(k+1,i)=sum((x.^k).*u(:,i));
    end
end

m(:,1)=mref(:,1);

myfun = @(t,y) mom_sph_diff(t,y,[M,Ds,R,j(t)]);

sol = ode23t(myfun, tspan, m(:,1));
solt = deval(sol,t)';

figure(1)
hold off
err=(mref'-solt)./mref';
plot(t,err(:,1:floor(M/2)),'-')
hold on
%plot(t,deval(sol,t)','*')
cst=feval(u,R);

for i=1:size(t)
    [ab,abm]=chebyshev(floor(M/2),solt(i,:));
    ab=real(ab);
    abm=abs(abm);
    ml=floor(M/2);
    while ml>1
        %q=real(radau(ml-1,ab,R));
        %q=real(lobatto(ml-2,ab,0,R));
        q=real(gauss(ml,ab));
        if min(q(:,1))<0 || max(q(:,1))<R
            ml=ml-1;
        else
            break
        end
    end

    coeff = (polyeval(abl,abml,q(:,1))*q(:,2))/R;
    Cs = sum(polyeval(abl,abml,R).*coeff);

    if (min(q(:,1))<0)
        disp(q)
    end
     figure(2)
     hold on
     plot(t(i),(Cs-cst(i))/cst(i),'*')
end