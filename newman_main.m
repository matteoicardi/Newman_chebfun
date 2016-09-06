%% newmann.m -- an executable m-file for solving a partial differential equation
% Newman model
% 

close all
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

chebfuneps 1e-6
chebfun2eps 1e-6
%chebfunpref.setDefaults('eps',1e-8)
chebfunpref.setDefaults('splitting', 1)
%chebfunpref.setDefaults('maxLength',999)
%chebfunpref.setDefaults('extrapolate',0)
%chebfunpref.setDefaults('happinessCheck','Plateau')
cheboppref.setDefaults('bvpTol', 1e-6);

%% load parameters
newman_param_sym;

%% load simulation macroscale parameters (timestep, applied current, etc.)
newman_sim_param;

%% setup all functions, operators, equations, etc.. and initialise
newman_setup;


%% Time loop and update operators accordingly
i=1;
while i<=nt+1
    newman_timeloop
end

%% Plotting examples
newman_plot