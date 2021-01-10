% ------------------------------------------------------------------------
% CIRCUIT_SOLVER
% ------------------------------------------------------------------------
% EXAMPLE 3: DC circuit with diodes
% ------------------------------------------------------------------------
% By Alexis Angers
% https://github.com/Alexsimulation
% ------------------------------------------------------------------------
clear all; close all; clc;

% --- USER INPUT ---

% Define circuit edges: e1 -> starting nodes, e2 -> end nodes
% Each column is an edge, n0 are the starting nodes, n1 the end nodes
n0 = [1, 2, 4, 3, 4];
n1 = [2, 3, 3, 4, 1];

% Resistance vector r, where r_n is the resistance in edge n
r = [0, 2, 5, 10, 0];

% Voltage source vector Vs, where Vs_n is the tension source in edge n
Vs = [20, 0, 0, 0, 0];

% Diodes vector, where edges with diodes are specified, always forward
% direction
D = 3;

% --- END OF USER INPUT ---

% Create circuit
C = circuit(n0, n1, r, Vs, D);

% Get the currents, potentials and dissipated powers
% - You can change 'IVP' to any combination you want
% - Use the .numeric tab to just get the numeric data, without it you'll 
%   get both the display-ready values and the numerical data
results = C.get_results('IVP').numeric;

% Print results in console
C.print_results();

% Display circuit with current value at each edge
% Note: The circuit plot isn't that great at the moment, try, you'll see
C.plot_circuit('I');