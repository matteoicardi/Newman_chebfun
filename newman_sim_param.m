%% Setup preferences for solving the problem.
chebfuneps 1e-15
chebfun2eps 1e-15
%chebfunpref.setDefaults('eps',1e-8)
chebfunpref.setDefaults('splitting', 1)
%chebfunpref.setDefaults('maxLength',999)
%chebfunpref.setDefaults('extrapolate',0)
%chebfunpref.setDefaults('happinessCheck','Plateau')
cheboppref.setDefaults('bvpTol', 1e-15);
cheboppref.setDefaults('maxIter', 1000);
%cheboppref.setDefaults('damping', true);
%cheboppref.setDefaults('discretization', 'values');

%% Time stepping and discretization
dt   = hour/1000;
ttot = hour;
nt   = floor(ttot/dt);
% modes in the additional dimension
M = 1;

%% Problem set-up
% applied current
I_1C =  chebfun(@(t) I_ref,[0,ttot]);
cycles = 3;
%I = chebfun(@(t) I_ref*sin(2*pi*cycles*t/ttot),[0,ttot]);
I = chebfun(@(t) I_ref/100,[0,ttot]);

% initialise all equations with last time step
% TODO for some reason it is better not to initialise the chebfun operators
init  = 0;

% plot while solving
draw  = 1;

% store all functions at all time steps
store = 1;

% update after each equation or at the end of each time step
update_after = 0;