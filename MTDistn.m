classdef MTDistn < handle
    %MTDISTN A container object for a MF IPPF MTT posterior distribution.
    
    properties
        clusters            % Cell array of TrackDistn objects
        N                   % Number of clusters
    end
    
    methods
        
        % Constructor
        function obj = MTDistn(clusters)
            obj.clusters = clusters;
            obj.N = length(clusters);
        end
        
        
        % Copy
        function new = Copy(obj)
            c = cell(size(obj.clusters));
            for k = 1:length(c)
                c{k} = obj.clusters{k}.Copy;
            end
            new = MTDistn(c);
        end
        
        
        
        
    end
    
end

