classdef TrackGroupDistn < handle
    %TRACKGROUPDISTN A track cluster particle distribution.
    
    properties
        particles           % Cell array of TrackSet objects
        weights             % Particle weights
        members             % Array indicating which targets are represented by this cluster
        N                   % Number of targets in the group
    end
    
    methods
        
        % Constructor
        function obj = TrackGroupDistn(members, particles, weights)
            obj.particles = particles;
            obj.members = members;
            obj.N = length(obj.members);
            
            if nargin == 2
                N = length(particles);
                obj.weights = ones(N, 1)/N;
            elseif nargin == 3
                obj.weights = weights;
            else
                assert(false);
            end
        end
        
        
        % Copy
        function new = Copy(obj)
            p = cell(size(obj.particles));
            for k = 1:length(p)
                p{k} = obj.particles{k}.Copy;
            end
            new = TrackGroupDistn(obj.members, p, obj.weights);
        end
        
        
        
        
    end
    
end

