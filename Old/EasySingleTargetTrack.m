function [ Distns, ESS, ESS_pre_resam, num_resamples ] = EasySingleTargetTrack( Observs )
%EASYSINGLETARGETTRACK Runs a batch SMC tracking algorithm for a single
% target with no missed observations or clutter

global Par;

% Initialise particle array (an array of particle arrays)
Distns = cell(Par.T, 1);
ESS = zeros(Par.T, 1);
ESS_pre_resam = zeros(Par.T, 1);
num_resamples = 0;

% Initialise particle set
init_track = Track(0, 1, {Par.TargInitState{1}-[Par.TargInitState{1}(3:4)' 0 0]'}, 0);
init_track_set = TrackSet({init_track});
part_set = repmat({init_track_set}, Par.NumPart, 1);
InitEst = PartDistn(part_set);

% Loop through time
for t = 1:Par.T
    
    tic;
    
    disp('**************************************************************');
    disp(['*** Now processing frame ' num2str(t)]);
    
    if t==1
        [Distns{t}, ESS(t), ESS_pre_resam(t), resample] = BatchSMC(t, t, InitEst, Observs);
    elseif t<Par.L
        [Distns{t}, ESS(t), ESS_pre_resam(t), resample] = BatchSMC(t, t, Distns{t-1}, Observs);
    else
        [Distns{t}, ESS(t), ESS_pre_resam(t), resample] = BatchSMC(t, Par.L, Distns{t-1}, Observs);
    end

    if resample
        num_resamples = num_resamples + 1;
    end
    
    if resample
        disp('*** Resampled in this frame');
    end
    disp(['*** Frame ' num2str(t) ' processed in ' num2str(toc) ' seconds']);
    disp('**************************************************************');

end




end



function [Distn, ESS, ESS_pre, resample] = BatchSMC(t, L, Previous, Observs)
% Execute a step of the SMC batch sampler

global Par;

% Runs an SMC on a batch of frames. We can alter states in t-L+1:t. Frame
% t-L is also needed for filtering, transition density calculations, etc.

% Initialise particle array and weight array
Distn = Previous.Copy;

% Extend the first track with a blank state (it will be overwritten by the
% KF output) and known association
dummy = Distn.particles{1}.tracks{1}.Copy;
dummy.Extend(t, [0 0 0 0]', 1);

% Generate a list of associated observations
obs = ListAssocObservs(t, L, dummy, Observs);

% Create a temporary array for un-normalised weights
weight = zeros(Par.NumPart, 1);


post_arr = zeros(Par.NumPart, 1);
ppsl_arr = zeros(Par.NumPart, 1);

% Loop through particles
for ii = 1:Par.NumPart
    
    Part = Distn.particles{ii};
        
    % Filter the observations with a Kalman filter
    [KFMean, KFVar] = KalmanFilter(obs, Part.tracks{1}.GetState(t-L), 1E-2*eye(4));
%     [KFMean, KFVar] = KalmanFilter(obs, Part.tracks{1}.GetState(t-L), eye(4));
%     [KFMean, KFVar] = KalmanFilter(obs, Part.tracks{1}.GetState(t-L), Par.Q);
    
    % Propose a new track from the KF distribution
    [NewTrack, ppsl_prob] = SampleKalman(KFMean, KFVar);
    
    % Caluclate the probability of the artificial distribution term
    [~, prev_ppsl_prob] = SampleKalman(KFMean(1:L-1, 1), KFVar(1:L-1, 1), Part.tracks{1});
    
    % Calculate the outgoing posterior probability term
    prev_post_prob = Posterior(t-1, L-1, Part, Observs);
    
    % Update the track
    Part.tracks{1}.Update(t, NewTrack, ones(L,1));
    
    % Update the weights
    post_prob = Posterior(t, L, Part, Observs);
    weight(ii) = Distn.weight(ii) ...
               + (post_prob - prev_post_prob) ...
               + (prev_ppsl_prob - ppsl_prob);

    post_arr(ii) = (post_prob - prev_post_prob);
    ppsl_arr(ii) = (prev_ppsl_prob - ppsl_prob);
           
	% Store probabilities for the next time step
%     Distn.prev_ppsl(ii) = ppsl_prob;
%     Distn.prev_post(ii) = post_prob;
    
    if isnan(weight(ii))
        weight(ii) = -inf;
    end
    
end

assert(~all(isinf(weight)), 'All weights are zero');

disp(['Mean/Variance of proposal term: ' num2str([mean(ppsl_arr), var(ppsl_arr)])]);
disp(['Mean/Variance of posterior term: ' num2str([mean(post_arr), var(post_arr)])]);

% Normalise weights
weight = exp(weight);
weight = weight/sum(weight);
weight = log(weight);

% Attach weights to particles
for ii = 1:Par.NumPart
    Distn.weight(ii) = weight(ii);
end

% Test effective sample size
ESS_pre = CalcESS(weight);
assert(~isnan(ESS_pre), 'Effective Sample Size is non defined (probably all weights negligible)');

% PlotTracks(Distn)
% uiwait

if (ESS_pre < 0.5*Par.NumPart)
    % Resample
    Distn = SystematicResample(Distn, weight);
    ESS = Par.NumPart;
    resample = true;
   
else
    resample = false;
    ESS = ESS_pre;
    
end

% PlotTracks(Distn)
% uiwait

end