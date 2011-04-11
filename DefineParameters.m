% Define constant parameters for target tracker execution

global Par;

% Set random seed
Par.rand_seed = 0;

% 0, with 5 targets, T=50, is standard test
% 5, with 2 targets, T=20, set next to each other is track-coallescing test
% 2, with 5 targets, T=50, is challenging (tgt 1 turns a sharp corner with no detections

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Flags                                                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Par.FLAG_ObsMod = 0;            % 0 = cartesian, 1 = polar,
Par.FLAG_PseudoJoint = true;   % Use joint tracking for colliding targets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Scene parameters                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Par.T = 20;                             % Number of frames
Par.P = 1; P = Par.P;                   % Sampling period
Par.Xmax = 500;                         % Scene limit (half side or radius depending on observation model)
Par.Vmax = 10;                          % Maximum velocity


Par.UnifVelDens = 1/(2*Par.Vmax)^2;             % Uniform density on velocity

if Par.FLAG_ObsMod == 0
    Par.UnifPosDens = 1/(2*Par.Xmax)^2;         % Uniform density on position
    Par.ClutDens = Par.UnifPosDens;             % Clutter density in observation space
elseif Par.FLAG_ObsMod == 1
    Par.UnifPosDens = 1/(pi*Par.Xmax^2);        % Uniform density on position
    Par.ClutDens = (1/Par.Xmax)*(1/(2*pi));     % Clutter density in observation space
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Scenario parameters                                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Par.NumTgts = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Target dynamic model parameters                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Par.ProcNoiseVar = 1;                                                      % Gaussian process noise variance (random accelerations)
Par.A = [1 0 P 0; 0 1 0 P; 0 0 1 0; 0 0 0 1];                              % 2D transition matrix using near CVM model
Par.B = [P^2/2*eye(2); P*eye(2)];                                          % 2D input transition matrix (used in track generation when we impose a deterministic acceleration)
Par.Q = Par.ProcNoiseVar * ...
    [P^3/3 0 P^2/2 0; 0 P^3/3 0 P^2/2; P^2/2 0 P 0; 0 P^2/2 0 P];          % Gaussian motion covariance matrix (discretised continous random model)
%     [P^4/4 0 P^3/2 0; 0 P^4/4 0 P^3/2; P^3/2 0 P^2 0; 0 P^3/2 0 P^2];      % Gaussian motion covariance matrix (piecewise constant acceleration discrete random model)
Par.ExpBirth = 0.5;                                                        % Expected number of new targets in a frame (poisson deistributed)
Par.PDeath = 0.01;                                                          % Probability of a (given) target death in a frame

Par.Qchol = chol(Par.Q);                                                   % Cholesky decompostion of Par.Q

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Observation model parameters                                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Par.ExpClutObs = 100;                      % Number of clutter objects expected in scene - 1000 is dense for Xmax=500, 160 for Xmax=200
Par.PDetect = 0.75;                         % Probability of detecting a target in a given frame

if Par.FLAG_ObsMod == 0
    Par.ObsNoiseVar = 1;                % Observation noise variance
    Par.R = Par.ObsNoiseVar * eye(2);   % Observation covariance matrix
    Par.C = [1 0 0 0; 0 1 0 0];         % 2D Observation matrix
elseif Par.FLAG_ObsMod == 1
    Par.BearingNoiseVar = 1E-4;                                 % Bearing noise variance
    Par.RangeNoiseVar = 1;                                      % Range noise variance
    Par.R = [Par.BearingNoiseVar 0; 0 Par.RangeNoiseVar];       % Observation covariance matrix
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Algorithm parameters                                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Par.L = 5;                              % Length of rolling window
Par.NumPart = 500;                      % Number of particles

Par.PRemove = 0.05;                     % Probability of removing target in a given particle
Par.PReplace = 0.9;                     % Probability of putting a target back

Par.ResamThresh = 0.1;                  % Resampling threshold as a proportion of maximum

Par.Vlimit = 1.5*Par.Vmax;                % Limit above which we do not accept velocity (lh=0)
Par.KFInitVar = 1E-20;                  % Variance with which to initialise Kalman Filters (scaled identity matrix)

if Par.FLAG_ObsMod == 0
    Par.BirthExclusionRadius = 0;
else
    Par.BirthExclusionRadius = 100;
end