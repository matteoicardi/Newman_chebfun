%% newmann.m -- an executable m-file for solving a partial differential equation
% Newman model
% 

close all
clear all
clc
clf

%% load parameters
newman_param_sym;

%% load simulation macroscale parameters (timestep, applied current, etc.)
newman_sim_param;

%% setup all functions, operators, equations, etc.. and initialise
newman_setup;


%% Time loop and update operators accordingly
i=1;
while i<=nt
    %newman_coupled_timeloop
    newman_reduced_timeloop
    %newman_timeloop
    %diffusion_timeloop
end

%% Plotting examples
%newman_plot
newman_reduced_plot