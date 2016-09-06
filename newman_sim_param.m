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
%I = chebfun(@(t) I_ref*cos(0*2*pi*cycles*t/ttot),[0,ttot]);
I = chebfun(@(t) I_ref,[0,ttot]);

init=0;
draw=0;
