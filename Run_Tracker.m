% Base script for multi-frame multi-target tracker using batch SMC with
% independence assumptions and collision detection

% Clear the workspace (maintaining breakpoints)
clup
dbstop if error

% Define all the necessary parameters in a global structure.
DefineParameters;

% Set a standard random stream (for repeatability)
s = RandStream('mt19937ar', 'seed', Par.rand_seed);
RandStream.setDefaultStream(s);

% Specify target behaviour
TargSpec = SpecifyTargetBehaviour;

for i=1:5-Par.NumTgts
    [~]=unidrnd(50);
    [~] = unifrnd(0.15*Par.Xmax, 0.25*Par.Xmax);
    [~] = unifrnd(-pi, pi);
    [~] = unifrnd(-pi, pi);
    [~] = unifrnd(-pi, pi);
end

% Generate target motion
[TrueState, TargSpec] = GenerateTargetMotion(TargSpec);

% Generate observations from target states
[Observs, detections] = GenerateObs(TrueState);

% Plot states and observations
fig = PlotTrueState(TrueState);
PlotObs(Observs, detections);

% Run tracker
[ Distns, ObsTargIndexes, ESS_post, ESS_pre, num_resamples ] = MultiTargetTrack(detections, Observs, {TargSpec(:).state} );

% Plot final estimates
PlotTracks(Distns{Par.T}, fig);

% Analyse associations
[ass, count, present] = AnalyseAss( detections, Distns{Par.T}, Par.T);
% 
% % Plot ESS
% % figure, plot(ESS_post), ylim([0 Par.NumPart])
% % figure, plot(ESS_pre), ylim([0 Par.NumPart])
% disp(['Particles resampled ' num2str(num_resamples) ' times']);