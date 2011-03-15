function [ TrueState, TargSpecOut ] = GenerateTargetMotion(TargSpec)
%GENERATETARGETMOTION Use a target specification to generate multi-target
%state sets for each frame

global Par;

% Input:    TargSpec: structure array specifying target existence and dynamics
%                       .birth - birth time
%                       .death - death time
%                       .state - initial birth state
%                       .acc - array of deterministic accelerations

% Output:   TrueState: cell array of track objects

% Initialise State cell array
TrueState = cell(Par.NumTgts, 1);

% Generate tracks
for j = 1:Par.NumTgts
    
    birth = TargSpec(j).birth;
    death = TargSpec(j).death;
    num = death - birth;
    
    % Create track arrays
    state = cell(num, 1);
    assoc = zeros(num, 1);
    
    % First frame
    state{1} = TargSpec(j).state;
    
    % Loop through frames
    for k = 2:num
        
%         [bng, range] = cart2pol(state{k-1}(1), state{k-1}(2));
%         if range > 0.75*Par.Xmax
%             magn = -0.01*(range - 0.75*Par.Xmax);
%             [acc(1), acc(2)] = pol2cart(bng, magn);
%             TargSpec(j).acc(k, :) = acc;
%         end

        % Calculate expected state
        exp_state = Par.A * state{k-1} + Par.B * TargSpec(j).acc(k, :)';
        
        % Sample state from Gaussian
        state{k} = mvnrnd(exp_state', Par.Q)';

        % Kill if outside scene
        if any(abs(state{k}(1:2))>Par.Xmax)
            state(k:end) = [];
            assoc(k:end) = [];
            death = birth + k - 1;
            num = death - birth;
            TargSpec(j).death = death;
            break;
        end
        
        % Limit velocity
        state{k}(3:4) = min( max( state{k}(3:4), -Par.Vmax), Par.Vmax);
        
    end
    
    TrueState{j} = Track(birth, death, state, assoc);
    
end

TargSpecOut = TargSpec;

end