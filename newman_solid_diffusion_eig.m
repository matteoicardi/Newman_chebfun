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
clf
clear all
close all

%% Problem set-up
% Create an interval of the space domain...
R = 1;
dom = [0 R];
%...and specify a sampling of the time domain:
t = linspace(0,10,1000)*2*pi;

%load jj
Ds = 1e-2;
omega = 1;
j = chebfun(@(tt) cos(omega*tt),[t(1) t(end)]);
%j = chebfun(jj',dom,'equi');
M=5;


%% Equations
x = chebfun('x',dom);

% Make the right-hand side of the PDE.
pdefun = @(t,x,u) Ds*(diff(u,2)) + Ds*2*diff(u)./(x+1e-8);
%pdefun = @(t,x,u) Ds*diff(x.^2.*diff(u))./(x+1e-5).^2;
L=chebop(dom);
L.op = @(x,u) Ds*(diff(u,2)) + Ds* 2*diff(u)./(x+1e-8);
%L.op = @(x,u) Ds*diff(x.^2.*diff(u))./(x+1e-8).^2;

%% Assign boundary conditions.
%L.lbc = 'dirichlet';
%L.rbc = @(u) [Ds*diff(u)-j(0),u];
L.bc = @(x,u) [Ds*feval(diff(u),R)-j(0);feval(u,R)];

bc.left = 'neumann';
bc.right = @(t,u) Ds*diff(u)-j(t);

% function to transform into homogeneous Neumann
v = 1/(2*R)*x.^2 - 3*R/10;

u0 = L\0; %0.5*x.^0 + v*j(0);
%u0 = 0 *x;

%% Eigenfunctions
% cartesian coordinates
%[ef,ev] = eigs(L,M);
w=x.^0;
% spherical coordinates
% compute spherical bessel functions
r=roots(diff(chebfun(@(x) sphbes(0,x),[0 200])));
r=r(r>1e-5);
ef=chebmatrix(x.^0,dom);
w = x.^2;  % weight for integrating in spherical coord
for i=1:M
    ef(i,:)=chebfun(@(x) sphbes(0,r(i)*x/R),dom);
    ef(i,:) = ef(i,:)/sqrt(sum(w.*ef(i,:).^2));
end
ef(M,:)=x.^0;
ef(M,:) = ef(M,:)/sqrt(sum(w.*ef(M,:).^2));

% % Modify v to be orthogonal to ef
% for i=1:M
% v=v-sum(w.*v.*ef{i})*ef{i};
% end

%% compute spectral mass and stiffness matrix
def = diff(ef);
for i=1:M
    vcoef(i) = sum(w.*v.*ef(i));
    %for k=1:M
    %    dcoef(i,k) = sum(w.*def(i).*def(k));
    %end
    dcoef(i,i) = sum(w.*def(i).*def(i));
    ncoef(i) = sum(w.*ef(i))/R;
end
% spherical only
ncoef = 3*ncoef;
%dcoef(dcoef/norm(dcoef)<1e-8)=0;

%% Call pde23t to solve the problem.
opts = pdeset('Eps', 1e-10, 'Ylim', [-3,3]);
[t, u] = pde23t(pdefun, t, u0, bc, opts);


% %% compute reference spectral coefficients
% tspan = [t(1),t(end)];
% m=ones(M,length(t));
% mref=ones(M,length(t));
% 
% uneu = u;
% for i=1:length(t)
%     uneu(:,i) = u(:,i)-v*j(t(i))/Ds;
%     for k=1:M
%         mref(k,i)=sum(w.*(uneu(:,i)).*ef(k));
%     end
% end
% m(:,1)=mref(:,1);
% 
% %% solve for spectral coefficients
% % Modified Laplacian
% myfun = @(t,y) eig_sph_diff(t,y,{ncoef,vcoef,dcoef,j(t),feval(diff(j),t),Ds});
% 
% sol = ode23t(myfun, tspan, m(:,1)');
% m = deval(sol,t);
% 
% umneu=u*0;
% um=u*0;
% for i=1:length(t)
%     uu = m(end,i)*ef{end};
%     for k=1:M-1
%         uu = uu+m(k,i)*ef{k};
%         err(k,i)=max(abs((uu+v*j(t(i))/Ds)-u(:,i)));
%         err2(k,i)=norm(((uu+v*j(t(i))/Ds)-u(:,i)));
%     end
%     umneu(:,i)=uu;
%     um(:,i)=uu+v*j(t(i))/Ds;
% end
% 
% mcheck=m;
% for i=1:length(t)
%     for k=1:M
%         mcheck(k,i)=sum(w.*(umneu(:,i)).*ef(k));
%     end
% end
% 
% 
% % %% plot
% % figure(1)
% % hold off
% % errm=(mref'-m');
% % plot(t,errm,'-')
% % 
% % figure(2)
% % hold off
% % waterfall(u,t)
% % hold on
% % waterfall(um,t)