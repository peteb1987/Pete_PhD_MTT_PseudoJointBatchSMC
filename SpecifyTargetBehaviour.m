function [ TargSpec ] = SpecifyTargetBehaviour
%SPECIFYTARGETBEHAVIOUR Generates a set of parameters specifying targets
%present in a scene

global Par;

% Output:   %TargSpec: structure array specifying target existence and dynamics            
                        % .birth - birth time
                        % .death - death time
                        % .state - initial birth state
                        % .acc - array of deterministic accelerations to simulate manoeuvring targets.

% Number of targets
N = Par.NumTgts;

% Number of time steps
T = Par.T;

% Initialise cell array for targets
TargSpec = repmat(struct('birth', 0, 'death', 0, 'state', zeros(4, 1), 'acc', []), N, 1);

% Set default parameters
for j = 1:N
    TargSpec(j).birth = 1;
    if Par.FLAG_DyingTargs
        TargSpec(j).death = unidrnd(T);
    else
        TargSpec(j).death = T + 1;
        [~]=unidrnd(T);
    end
    num = TargSpec(j).death - TargSpec(j).birth;
    TargSpec(j).state = zeros(4, 1);
%     if Par.FLAG_ObsMod == 0
%         TargSpec(j).state(1) = unifrnd(-0.5*Par.Xmax, 0.5*Par.Xmax);
%         TargSpec(j).state(2) = unifrnd(-0.5*Par.Xmax, 0.5*Par.Xmax);
%     elseif Par.FLAG_ObsMod == 1
      rng = unifrnd(0.15*Par.Xmax, 0.25*Par.Xmax);
%       rng = unifrnd(0.3*Par.Xmax, 0.5*Par.Xmax);
%         rng = unifrnd(0.10*Par.Xmax, 0.50*Par.Xmax);
        bng = unifrnd(-pi, pi);
        TargSpec(j).state(1) = rng*cos(bng);
        TargSpec(j).state(2) = rng*sin(bng);
%     end
    TargSpec(j).state(3) = unifrnd(-Par.Vmax, Par.Vmax);
    TargSpec(j).state(4) = unifrnd(-Par.Vmax, Par.Vmax);
    TargSpec(j).acc = zeros(num, 2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Manually overwrite individual target values if desired              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TargSpec(1).death = 11;
% TargSpec(1).birth = 5;
% TargSpec(2).birth = 20;
% TargSpec(3).birth = 35;

% TargSpec(1).state = [-100, 100, 3, 0]';
% TargSpec(2).state = [-100, 90, 3, 0]';
% TargSpec(3).state = [-100, 110, 3, 0]';
% TargSpec(4).state = [-100, 80, 3, 0]';
% TargSpec(5).state = [-100, 120, 3, 0]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% End of manual overwrites                                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end %SpecifyTargetBehaviour