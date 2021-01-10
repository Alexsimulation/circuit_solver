classdef dioGraph
    %DIOGRAPH Directed graph with edge order
    %   G = dioGraph(e1,e2) constructs a directed graph with edges specified
    %   by the node pairs (S,T). S and T must both be numeric, string vectors,
    %   or cell arrays of character vectors. S and T must have the same number
    %   of elements or be scalars.
    %   
    %   dioGraph properties:
    %       Edges            - Table containing edge information.
    %       Nodes            - Table containing node information.
    %
    %   dioGraph methods:
    %       incidence_matrix - Returns the incidence matrix of the graph
    %       plot_graph       - Plots the graph.
    %
    %   By Alexis Angers
    %   https://github.com/Alexsimulation
    %
    
    properties
        Edges
        Nodes
    end
    
    methods
        
        function G = dioGraph(e1,e2)
            if size(e1,1) ~= length(e1)
                e1 = e1.';
            end
            if size(e2,1) ~= length(e2)
                e2 = e2.';
            end
            EndNodes = [e1, e2];
            G.Edges = table(EndNodes);
        end
        
        function A = incidence_matrix(G)
            A = zeros(size(G.Edges.EndNodes,1),max(max(G.Edges.EndNodes)));
            for n = 1:size(A,1)
                for v = 1:size(A,2)
                    if G.Edges.EndNodes(n,1) == v
                        A(n,v) = -1;
                    elseif G.Edges.EndNodes(n,2) == v
                        A(n,v) = 1;
                    end
                end
            end
        end
        
        function [] = plot_graph(G,varargin)
            is3D = false;
            e1 = G.Edges.EndNodes(:,1);
            e2 = G.Edges.EndNodes(:,2);
            figure('units','normalized','outerposition',[0 0 1 1]);
            if isempty(varargin)
                H = digraph(e1,e2,1:length(e1));
            else
                if isempty(varargin{1})
                    H = digraph(e1,e2,1:length(e1));
                else
                    H = digraph(e1,e2,varargin{1});
                end
                if length(varargin) > 1
                    is3D = varargin{2};
                end
            end
            if is3D
                p = plot(H,'Layout','force3','EdgeLabel',H.Edges.Weight);
            else
                p = plot(H,'Layout','force','EdgeLabel',H.Edges.Weight);
            end
            axis equal;
            p.EdgeColor = [0, 0, 0];
            p.NodeColor = [0, 0, 0];
            p.MarkerSize = 7;
            p.LineWidth = 3;
        end
        
    end
end