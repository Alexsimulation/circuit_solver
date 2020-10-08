classdef circuit
    % CIRCUIT Electric circuit
    %   C = circuit(e1, e2, Impedance,Source) constructs a circuit with
    %   geometry described by the edges. edges must be row vectors,
    %   and both impedance and source must be real or complex valued
    %   vectors.
    %
    %   circuit properties:
    %       Graph            - Graph containing the circuit geometry
    %       Impedance        - Vector of edges impedance
    %       Source           - Vector of edges voltage sources
    %       Currents         - Vector of edges currents
    %       PotentialDiffs   - Vector of edges potential differences 
    %       Powers           - Vector of edges dissipated powers
    %
    %   circuit methods:
    %       plot_circuit     - Plots the circuit using dioGraph methods
    %       get_results      - Returns the currents, potentials and powers
    %       print_results    - Prints information about the circuit
    %
    %   By Alexis Angers
    %   https://github.com/Alexsimulation
    %
    %   See also DIOGRAPH, SOLVE_CIRCUIT
    
    properties
        Graph
        Impedance
        Source
        Currents
        PotentialDiffs
        Powers
    end
    methods
        function C = circuit(e1,e2,Impedance,Source)
            if size(Impedance,1) ~= length(Impedance)
                Impedance = Impedance.';
            end
            if size(Source,1) ~= length(Source)
                Source = Source.';
            end
            C.Graph = dioGraph(e1,e2);
            C.Impedance = Impedance;
            C.Source = Source;
            [I,V,P] = circuit.solve_circuit(C.Graph,Impedance,Source);
            C.Currents = I;
            C.PotentialDiffs = V;
            C.Powers = P;
        end
        
        function [] = plot_circuit(C,varargin)
            if isempty(varargin)
                C.Graph.plot_graph();
            else
                switch varargin{1}
                    case 'I'
                        u = abs(C.Currents);
                    case 'V'
                        u = abs(C.PotentialDiffs);
                    case 'P'
                        u = abs(C.Powers);
                end
                is3D = false;
                if length(varargin) > 1
                    if strcmp(varargin{2},'3D')
                        is3D = true;
                    end
                end
                C.Graph.plot_graph(u,is3D)
                title('Circuit plot');
            end
        end
        
        function results = get_results(C,res)
            I = circuit.polar_form(C.Currents,'deg');
            V = circuit.polar_form(C.PotentialDiffs,'deg');
            P = C.Powers;
            Z = C.Impedance;
            Edges = [1:length(I)]';
            results.display = [table(Edges),table(Z)];
            for n = 1:length(res)
                switch res(n)
                    case 'I'
                        results.display = [results.display, table(I)];
                    case 'V'
                        results.display = [results.display, table(V)];
                    case 'P'
                        results.display = [results.display, table(P)];
                end
            end
            I = C.Currents;
            V = C.PotentialDiffs;
            P = C.Powers;
            Z = C.Impedance;
            Edges = [1:length(I)]';
            results.numeric = [table(Edges),table(Z)];
            for n = 1:length(res)
                switch res(n)
                    case 'I'
                        results.numeric = [results.numeric, table(I)];
                    case 'V'
                        results.numeric = [results.numeric, table(V)];
                    case 'P'
                        results.numeric = [results.numeric, table(P)];
                end
            end
        end
        
        function [] = print_results(C)
            results = C.get_results('IVP');
            fprintf('\n');
            disp(results.display);
        end
        
    end
    
    methods(Static)
        
        function [varargout] = solve_circuit(Graph,Impedance,Source)
            %SOLVE_CIRCUIT solves an electrical circuit problem.
            %   I = SOLVE_CIRCUIT(Graph,Impedance,Source) finds the currents in a
            %   circuit. The circuit is described by a directed graph, an impedance
            %   vector and a potential source vector. 
            %   
            %   [I,V] = SOLVE_CIRCUIT(Graph,Impedance,Source) also returns the
            %   potential differences across each edges.
            %
            %   [I,V,P] = SOLVE_CIRCUIT(Graph,Impedance,Source) also returns the
            %   power dissipated in each edges.
            %   
            %   Graph := A dioGraph object.
            %   I,V,P := Collumn vectors of complex values.
            %
            %   See also dioGRAPH

            % Incidence matrix A
            A = Graph.incidence_matrix();

            % Impedance matrix Z
            Z = diag(Impedance);

            % Find null space matrix N for kirchoff's loop rule N'*Z*y = s -> 
            % first system, the null space matrix of the incidence matrix contains 
            % all loops of graph G
            N = null(A','r');

            % Close system by using kirchhoff's nodal rule A'*y = 0 -> second system
            K = rref(A');
            K = K(1:end-1,:);

            % Assemble full system matrix D and the full solution vector d
            D = [N'*Z; K];
            d = [N'*Source; zeros(size(D,1)-length(N'*Source),1)];

            % Solve for currents I,  potential differences V and power P
            I = D^-1*d;
            V = Z*I+Source;
            P = Impedance.*(abs(I).^2);

            % Return results
            varargout{1} = I;
            varargout{2} = V;
            varargout{3} = P;
        end
        
        function out = polar_form(A,varargin)
            if all(imag(A) == 0)
                out = A;
            else
                absA = abs(A);
                if ~isempty(varargin) & (varargin{1} == 'deg')
                    phaseA = angle(A)*180/pi;
                else
                    phaseA = angle(A);
                end
                out = arrayfun(@(x, y) sprintf('%0.3f < %0.3f', x, y), absA, phaseA, 'uni', 0);
            end
        end
        
    end
end