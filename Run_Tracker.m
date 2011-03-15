% Base script for multi-frame multi-target tracker using batch SMC with
% independence assumptions and collision detection

% Clear the workspace (maintaining breakpoints)
clup
dbstop if error

% Set a standard random stream (for repeatability)
s = RandStream('mt19937ar', 'seed', 0);
RandStream.setDefaultStream(s);

% Define all the necessary parameters in a global structure.
DefineParameters;

% Specify target behaviour
TargSpec = SpecifyTargetBehaviour;

% Generate target motion
[TrueState, TargSpec] = GenerateTargetMotion(TargSpec);

% Generate observations from target states
[Observs, detections] = GenerateObs(TrueState);

% Plot states and observations
fig = PlotTrueState(TrueState);
PlotObs(Observs, detections);

% Run tracker
[ Distns, ESS_post, ESS_pre, num_resamples ] = MultiTargetTrack( Observs, {TargSpec(:).state} );

% % Plot final estimates
% PlotTracks(Distns{Par.T}, fig);
% 
% % Analyse associations
% [ass, count] = AnalyseAss( detections, Distns);
% 
% % Plot ESS
% % figure, plot(ESS_post), ylim([0 Par.NumPart])
% % figure, plot(ESS_pre), ylim([0 Par.NumPart])
% disp(['Particles resampled ' num2str(num_resamples) ' times']);