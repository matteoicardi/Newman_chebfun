%% Setup preferences for solving the problem.

eps = 1e-3;   % dimensional tolerance for solvers and for smoothed domains
releps = 1e-6;  % relative tolerance for solvers

chebfunpref.setDefaults('factory')
chebfunpref.setDefaults('chebfuneps', releps);
chebfunpref.setDefaults({'cheb2Prefs','chebfun2eps'}, releps^2);
chebfunpref.setDefaults('splitting', 0);
%chebfunpref.setDefaults('maxLength',999)
%chebfunpref.setDefaults('extrapolate',1)
%chebfunpref.setDefaults('happinessCheck','Plateau')

cheboppref.setDefaults('factory')
cheboppref.setDefaults('bvpTol', releps);
cheboppref.setDefaults('maxIter', 500);
cheboppref.setDefaults('damping', 0);
%cheboppref.setDefaults('discretization', 'values');
%cheboppref.setDefaults('discretization', 'coeffs');

%% Load Physical parameters
newman_param_sym  % fully symmetric paramters
%newman_param      % asymmetrical cell

%% Time stepping and discretization
dt   = hour/10000;
ttot = hour/1;
nt   = 120;%floor(ttot/dt);
% modes in the additional dimension
M = 10;

%% Problem set-up
% applied current  % I_ref is the 1C rate
cycles = 4;
I = chebfun(@(t) 10*I_ref*sin(2*pi*cycles*t/ttot),[0,ttot]);
%I = chebfun(@(t)    I_ref    ,[0,ttot]);

% initialise all equations with last time step
% for some reason it is better not to initialise the chebfun operators
init  = 0;

% plot while solving
draw  = 1;
drawstep = 1;  % interval for plotting

% store all functions at all time steps
store = 1;

% update after each equation or at the end of each time step
update_after = 0;