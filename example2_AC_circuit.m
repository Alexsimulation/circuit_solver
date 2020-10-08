% ------------------------------------------------------------------------
% CIRCUIT_SOLVER
% ------------------------------------------------------------------------
% EXAMPLE 2: AC circuit
% ------------------------------------------------------------------------
% By Alexis Angers
% https://github.com/Alexsimulation
% ------------------------------------------------------------------------
clear all; close all; clc;

% --- USER INPUT ---

% Define circuit edges: e1 -> starting nodes, e2 -> end nodes
% If you aren't sure how to define a circuit, see example 1
% Each column is an edge, n0 are the starting nodes, n1 the end nodes
n0 = [1, 2, 3, 4, 2, 3, 4, 1, 5];
n1 = [2, 3, 4, 1, 5, 6, 6, 5, 6];

% Impedance vector z, where z_n is the resistance in edge n
z = [0, 2+1i, 5-3i, 0, 1, 1+1i, 0, 4, 1];

% Voltage source vector Vs, where Vs_n is the tension source in edge n
Vs = [120*exp(0.2i), 0, 0, 0, 0, 0, 0, 0, 0];

% --- END OF USER INPUT ---

% Create circuit
C = circuit(n0, n1, z, Vs);

% Get the currents, potentials and dissipated powers
% - You can change 'IVP' to any combination you want
% - Use the .numeric tab to just get the numeric data, without it you'll 
%   get both the display-ready values and the numerical data
results = C.get_results('IVP').numeric;

% Print results in console, note that angles are in degrees
C.print_results();

% Display circuit with current value at each edge
% Note: The circuit plot isn't that great at the moment, try, you'll see
C.plot_circuit('I','3D');