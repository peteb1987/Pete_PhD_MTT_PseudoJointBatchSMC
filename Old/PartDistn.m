classdef PartDistn < handle
    %PARTDISTN A particle approximation of a multi-target state
    
    properties
        particles       % Cell array of TrackSet objects
        weight          % Vector of particle weights, of equal length
        prev_ppsl       % Vector of proposal probabilities for the previous step (required in the weight evaluation of the subsequent step)
        prev_post       % Vector of posterior probabilities for the previous step (required in the weight evaluation of the subsequent step)
    end
    
    methods
        
        % Constructor
        function obj = PartDistn(particles, weight, prev_ppsl, prev_post)
            if nargin == 4
                obj.particles = particles;
                obj.weight = weight;
                obj.prev_ppsl = prev_ppsl;
                obj.prev_post = prev_post;
                assert(length(particles)==length(weight), 'Weight array is the wrong size');
                assert(length(particles)==length(prev_ppsl), 'Previous proposal probability array is the wrong size');
                assert(length(particles)==length(prev_post), 'Previous posterior probability array is the wrong size');
            elseif nargin == 2
                obj.particles = particles;
                obj.weight = weight;
                assert(length(particles)==length(weight), 'Particle array and weight array have different sizes');
            else
                obj.particles = particles;
                obj.weight = ones(length(particles),1)/length(particles);
                obj.prev_ppsl = zeros(length(particles),1);
                obj.prev_post = zeros(length(particles),1);
            end
        end
        
        
        
        % Copy
        function new = Copy(obj)
            p = cell(size(obj.particles));
            for k = 1:length(p)
                p{k} = obj.particles{k}.Copy;
            end
            new = PartDistn(p, obj.weight, obj.prev_ppsl, obj.prev_post);
        end
        
        
        
        
    end
    
end

