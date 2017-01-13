%% newmann.m -- an executable m-file for solving a partial differential equation
% Automatically created in CHEBGUI by user icardi.
% Created on February 10, 2016 at 20:24.

%% Problem description.
% Solving
%   ce_t = ce" - ce,
%   pe" + (log(ce))' + ce = 0,
%   ps" - ce = 0,
% for x in [-1,1] and t in [0,1], subject to
%   ce = 1, pe = 1, ps = 1 at x = -1
% and
%   ce = 1, pe = 1, ps = 1 at x = 1

%% Problem set-up
% Create an interval of the space domain...
dom = [-1 1];
%...and specify a sampling of the time domain:
t = 0:.1:1;

% Make the right-hand side of the PDE.
pdefun = @(t,x,ce,pe,ps) [diff(ce,2)-ce; -diff(pe,2)-diff(log(ce))-ce; -diff(ps,2)+ce];
pdeflag = [1  0  0]; % Zero when a variable is indep of time.

% Assign boundary conditions.
bc.left = @(t,ce,pe,ps) [ce-1; pe-1; ps-1];
bc.right = @(t,ce,pe,ps) [ce-1; pe-1; ps-1];

% Construct a chebfun of the space variable on the domain,
x = chebfun(@(x) x, dom);
% and of the initial conditions.
ce0 = 1;
pe0 = 1 - .5.*cos(.5.*pi.*x);
ps0 = 1;
sol0 = [ce0, pe0, ps0];

%% Setup preferences for solving the problem.
opts = pdeset('Eps', 1e-6, 'PDEflag', pdeflag, 'Ylim', [0.5,1]);

%% Call pde15s to solve the problem.
[t, ce, pe, ps] = pde15s(pdefun, t, sol0, bc, opts);

%% Plot the solution components.
figure
waterfall(ce, t)
xlabel('x'), ylabel('t'), title('ce')
figure
waterfall(pe, t)
xlabel('x'), ylabel('t'), title('pe')
figure
waterfall(ps, t)
xlabel('x'), ylabel('t'), title('ps')