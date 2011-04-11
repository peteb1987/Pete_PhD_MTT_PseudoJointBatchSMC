classdef TrackGroupDistn < handle
    %TRACKGROUPDISTN A track cluster particle distribution.
    
    properties
        particles           % Cell array of TrackSet objects
        weights             % Particle weights
        members             % Array indicating which targets are represented by this cluster
        N                   % Number of targets in the group
        locked              % A counter giving the number of time steps before this cluster may break
    end
    
    methods
        
        % Constructor
        function obj = TrackGroupDistn(members, particles, weights, locked)
            obj.particles = particles;
            obj.members = members;
            obj.N = length(obj.members);
            
            if nargin == 2
                N = length(particles);
                obj.weights = log(ones(N, 1)/N);
            elseif nargin > 2
                obj.weights = weights;
            else
                assert(false);
            end
            
            if nargin == 4
                obj.locked = locked;
            else
                obj.locked = 0;
            end
            
        end
        
        
        % Copy
        function new = Copy(obj)
            p = cell(size(obj.particles));
            for k = 1:length(p)
                p{k} = obj.particles{k}.Copy;
            end
            new = TrackGroupDistn(obj.members, p, obj.weights, obj.locked);
        end
        
        
        
        
    end
    
end

