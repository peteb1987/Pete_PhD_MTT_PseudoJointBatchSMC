classdef TrackGroupDistn < handle
    %TRACKGROUPDISTN A track cluster particle distribution.
    
    properties
        particles           % Cell array of TrackSet objects
        weight             % Particle weights
        members             % Array indicating which targets are represented by this cluster
    end
    
    methods
        
        % Constructor
        function obj = TrackGroupDistn(members, particles, weights)
            obj.particles = particles;
            obj.members = members;
            
            if nargin == 2
                N = length(particles);
                obj.weight = ones(N, 1)/N;
            elseif nargin == 3
                obj.weight = weights;
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
            new = TrackGroupDistn(obj.members, p, obj.weight);
        end
        
        
        
        
    end
    
end

