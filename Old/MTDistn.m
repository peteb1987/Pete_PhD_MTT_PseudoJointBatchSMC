classdef MTDistn < handle
    %MTDISTN A container object for a MF IPPF MTT posterior distribution.
    
    properties
        clusters            % Cell array of TrackDistn objects
        weight             % Particle weights
    end
    
    methods
        
        % Constructor
        function obj = MTDistn(clusters, weights)
            obj.clusters = clusters;
            
            if nargin == 1
                N = length(clusters{1}.particles);
                obj.weight = ones(N, 1)/N;
            elseif nargin == 2
                obj.weight = weights;
            else
                assert(false);
            end
            
        end
        
        
        % Copy
        function new = Copy(obj)
            c = cell(size(obj.clusters));
            for k = 1:length(c)
                c{k} = obj.clusters{k}.Copy;
            end
            new = MTDistn(c, obj.weight);
        end
        
        
        
        
    end
    
end

