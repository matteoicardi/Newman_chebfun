%% newmann.m -- an executable m-file for solving a partial differential equation
% Newman model
% 

close all
clear all
clc
clf

%% load simulation parameters
newman_sim_param;

%% setup all functions, operators, equations, etc.. and initialise
newman_setup;


%% Time loop and update operators accordingly
i=1;
while i<=nt
    %newman_coupled_timeloop    % coupled newman
    %newman_reduced_timeloop    % only solid equations in half cell
    %newman_timeloop            % full newman
    diffusion_timeloop         % only solid diffusion
end

%% Plotting examples
%newman_plot
newman_reduced_plot