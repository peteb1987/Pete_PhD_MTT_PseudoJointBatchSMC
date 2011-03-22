function [ Distns, ObsTargIndexes, ESS_post, ESS_pre, num_resamples ] = MultiTargetTrack( Observs, InitState )
%MULTITARGETTRACK Runs a batch SMC tracking algorithm for multiple
% targets with missed observations and clutter

global Par;

% Initialise particle array and diagnostics (an array of particle arrays)
Distns = cell(Par.T, 1);
ObsTargIndexes = cell(Par.T, 1);
Area = cell(Par.T, 1);
ESS_post = cell(Par.T, Par.NumTgts);
ESS_pre = cell(Par.T, Par.NumTgts);
num_resamples = zeros(Par.NumTgts, 1);

% Initialise particle set - separate
init_track = cell(Par.NumTgts, 1);
init_track_distn = cell(Par.NumTgts, 1);
for j = 1:Par.NumTgts
    init_track{j} = Track(0, 1, {InitState{j}-[InitState{j}(3:4)' 0 0]'}, 0);
    init_track_distn{j} = TrackGroupDistn(j, repmat({TrackSet(j, init_track(j))}, Par.NumPart, 1));
end
InitEst = MTDistn(init_track_distn);

% % Initialise particle set - together
% init_track = cell(Par.NumTgts, 1);
% for j = 1:Par.NumTgts
%     init_track{j} = Track(0, 1, {InitState{j}-[InitState{j}(3:4)' 0 0]'}, 0);
% end
% init_track_distn = TrackGroupDistn(1:Par.NumTgts, repmat({TrackSet(1:Par.NumTgts, init_track)}, Par.NumPart, 1));
% InitEst = MTDistn({init_track_distn});


% Loop through time
for t = 1:Par.T
    
    tic;
    
    disp('**************************************************************');
    disp(['*** Now processing frame ' num2str(t)]);
    
    if t==1
        [Distns{t}, ESS_post{t}, ESS_pre{t}, resamples, ObsTargIndexes, Area] = BatchSMC(t, t, InitEst, Observs, ObsTargIndexes, Area);
    elseif t<Par.L
        [Distns{t}, ESS_post{t}, ESS_pre{t}, resamples, ObsTargIndexes, Area] = BatchSMC(t, t, Distns{t-1}, Observs, ObsTargIndexes, Area);
    else
        [Distns{t}, ESS_post{t}, ESS_pre{t}, resamples, ObsTargIndexes, Area] = BatchSMC(t, Par.L, Distns{t-1}, Observs, ObsTargIndexes, Area);
    end

    num_resamples = num_resamples + resamples;
    
    assoc = [];
    for c = 1:Distns{t}.N
        for j = 1:Distns{t}.clusters{c}.N
            assoc = [assoc, Distns{t}.clusters{c}.particles{round(Par.NumPart/2)}.tracks{j}.GetAssoc(t-min(t,Par.L)+1)];
        end
    end
    disp(['*** Particle ' num2str(round(Par.NumPart/2)) ' associations at t-L+1: ' num2str(assoc)]);

    disp(['*** Frame ' num2str(t) ' processed in ' num2str(toc) ' seconds']);
    disp('**************************************************************');

    if round(t/10)==t/10
        PlotTracks(Distns{t});
        pause(1);
    end
    
end




end



function [Distn, ESS_post, ESS_pre, resamples, ObsTargIndexes, Area] = BatchSMC(t, L, Previous, Observs, ObsTargIndexes, Area)
% Execute a step of the SMC batch sampler

% Runs an SMC on a batch of frames. We can alter states in t-L+1:t. Frame
% t-L is also needed for filtering, transition density calculations, etc.

global Par;

ESS_post = zeros(Previous.N,1);
ESS_pre = zeros(Previous.N,1);
resamples = zeros(Previous.N,1);

% Initialise particle array and weight array
Distn = Previous.Copy;

% Loop through clusters
for c = 1:Distn.N
    % Loop through particles
    for ii = 1:Par.NumPart
        % Project all tracks with an ML prediction
        Distn.clusters{c}.particles{ii}.ProjectTracks(t);
    end
end

% Detect relevant observations
[ObsTargIndexes{t}] = DetectNearbyObservations(t, Observs, Distn);%, Area{t}
disp(['*** ' num2str(cellfun(@length, ObsTargIndexes{t})') ' observations in respective target vicinities']);

% post_prob_array = zeros(Par.NumPart, 1);
% state_ppsl_array = zeros(Par.NumPart, 1);
% jah_ppsl_array = zeros(Par.NumPart, 1);
curr_arr = zeros(Par.NumPart, 1);
% prev_arr = zeros(Par.NumPart, 1);

% Loop through clusters
for c = 1:length(Distn.clusters)
    
    weights = zeros(Par.NumPart, 1);
    
    % Loop through particles
    for ii = 1:Par.NumPart
        
        Cluster = Distn.clusters{c}.particles{ii};
        
%         prev_post_prob = Posterior(t-1, L-1, Cluster, Observs);
%         prev_jah_ppsl = Cluster.SampleAssociations(t-1, L-1, Observs, true);
%         prev_state_ppsl = Cluster.SampleStates(t-1, L-1, Observs, true);

        Cluster_OTI = cellfun(@(x) x(c), ObsTargIndexes(1:t));

        % Sample and update associations
        jah_ppsl = Cluster.SampleAssociations(t, L, Observs, Cluster_OTI, false);
    
        % Sample and update states
        state_ppsl = Cluster.SampleStates(t, L, Observs, false);
        
        % Calculate new posterior
        post_prob = Posterior(t, L, Cluster, Observs, Cluster_OTI);%, Area
    
        % Update the weight
        weights(ii) = Distn.clusters{c}.weights(ii) ...
                   + (post_prob - sum(state_ppsl) - jah_ppsl);%...
%                    - (prev_post_prob - sum(prev_state_ppsl) - prev_jah_ppsl);
               
        curr_arr(ii) = (post_prob - sum(state_ppsl) - jah_ppsl);
%         prev_arr(ii) = (prev_post_prob - sum(prev_state_ppsl) - prev_jah_ppsl);
               
%         post_prob_array(ii) = post_prob;
%         state_ppsl_array(ii) = sum(state_ppsl(:));
%         jah_ppsl_array(ii) = jah_ppsl;
        
        if isnan(weights(ii))
            weights(ii) = -inf;
        end
        
        if isinf(weights(ii))
            disp(['Uh-oh: Zero weight in cluster ' num2str(c) ', particle ' num2str(ii)]);
        end
        
    end
    
    assert(~all(isinf(weights)), 'All weights are zero');
    
    % Normalise weights
    max_weight = max(weights); max_weight = max_weight(1); weights = weights - max_weight;
    weights = exp(weights); weights = weights/sum(weights);  weights = log(weights);
    
    % Attach weights to particles
    Distn.clusters{c}.weights = weights;
    
    % Calculate effective sample size for diagnostics
    ESS_pre(c) = CalcESS(weights);
    assert(~isnan(ESS_pre(c)), 'Effective Sample Size is non defined (probably all weights negligible)');
    
    if (ESS_pre(c) < Par.ResamThresh*Par.NumPart)
        [Distn.clusters{c}, weights] = ConservativeResample(Distn.clusters{c}, weights);
        resamples(c) = resamples(c) + 1;
        ESS_post(c) = CalcESS(weights);
        disp(['*** Target Cluster' num2str(c) ': Effective Sample Size = ' num2str(ESS_pre(c)) '. RESAMPLED (Conservative). ESS = ' num2str(ESS_post(c))]);
    else
        disp(['*** Target Cluster' num2str(c) ': Effective Sample Size = ' num2str(ESS_pre(c))]);
        ESS_post(c) = ESS_pre(c);
    end
    
end

end

