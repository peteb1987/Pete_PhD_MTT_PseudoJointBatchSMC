function [ Distns, ESS, ESS_pre_resam, num_resamples ] = MultiTargetTrack( Observs, InitState )
%MULTITARGETTRACK Runs a batch SMC tracking algorithm for multiple
% targets with missed observations and clutter

global Par;

% Initialise particle array and diagnostics (an array of particle arrays)
Distns = cell(Par.T, 1);
ESS = cell(Par.T, Par.NumTgts);
ESS_pre_resam = cell(Par.T, Par.NumTgts);
num_resamples = zeros(Par.NumTgts, 1);

% Initialise particle set
init_track = cell(Par.NumTgts, 1);
for j = 1:Par.NumTgts
    init_track{j} = Track(0, 1, {InitState{j}-[InitState{j}(3:4); 0; 0]}, 0);
end

%%% REWRITING HERE %%%


% Loop through time
for t = 1:Par.T
    
    tic;
    
    disp('**************************************************************');
    disp(['*** Now processing frame ' num2str(t)]);
    
    if t==1
        [Distns{t}, ESS{t}, ESS_pre_resam{t}, resample, Observs] = BatchSMC(t, t, InitEst, Observs, Par.NumPart);
    elseif t<Par.L
        [Distns{t}, ESS{t}, ESS_pre_resam{t}, resample, Observs] = BatchSMC(t, t, Distns{t-1}, Observs, ESS{t-1});
    else
        [Distns{t}, ESS{t}, ESS_pre_resam{t}, resample, Observs] = BatchSMC(t, Par.L, Distns{t-1}, Observs, ESS{t-1});
    end

    if resample
        num_resamples = num_resamples + 1;
    end
    
    assoc = [];
    for j=1:Par.NumTgts
        assoc = [assoc, Distns{t}.clusters{j}.particles{round(Par.NumPart/2)}.tracks{1}.GetAssoc(t-min(t,Par.L)+1)];
    end
    disp(['*** Particle ' num2str(round(Par.NumPart/2)) ' associations at t-L+1: ' num2str(assoc)]);
    
%     if resample
%         disp(['*** ESS = ' num2str(ESS_pre_resam(t)) ' : Resampled in this frame']);
%     else
%         disp(['*** ESS = ' num2str(ESS_pre_resam(t))]);
%     end
%     disp(['*** Tracks detected in first particle: ' num2str(Distns{t}.particles{1}.N)]);
    disp(['*** Frame ' num2str(t) ' processed in ' num2str(toc) ' seconds']);
    disp('**************************************************************');

    if round(t/10)==t/10
        PlotTracks(Distns{t});
        pause(5);
    end
    
end




end



function [Distn, ESS, ESS_pre, resample, PrunedObservs] = BatchSMC(t, L, Previous, Observs, ESS_prev)
% Execute a step of the SMC batch sampler

global Par;

ESS = zeros(Par.NumTgts,1);
ESS_pre = zeros(Par.NumTgts,1);
resample = false;

% Runs an SMC on a batch of frames. We can alter states in t-L+1:t. Frame
% t-L is also needed for filtering, transition density calculations, etc.

% Initialise particle array and weight array
Distn = Previous.Copy;

% Create a temporary array for un-normalised weights
weight = zeros(Par.NumPart, 1);

% % Identify potential birth sites
% if (~Par.FLAG_TargInit)&&(t>2)
%     BirthSites = FindBirthSites(t, Observs);
% else
%     BirthSites = cell(0, 1);
% end
% disp(['*** ' num2str(size(BirthSites, 1)) ' birth sites in this frame']);

% Loop through clusters
for c = 1:length(Distn.clusters)
    
    % Loop through particles
    for ii = 1:Par.NumPart
        
        Cluster = Distn.clusters{c}.particles{ii};
        
        % Project all tracks with an ML prediction
        Cluster.ProjectTracks(t);
        
    end
    
end

% Prune observations
Observs(t) = PruneObservs(t, Observs(t), Distn);
PrunedObservs = Observs;
% for tt = t-L+1:t
%     Observs(tt) = PruneObservs(tt, Observs(tt), Distn);
% end
disp(['*** ' num2str(Observs(t).N) ' observations survive gating in frame ' num2str(t)]);


% post_prob_array = zeros(Par.NumPart, 1);
% state_ppsl_array = zeros(Par.NumPart, 1);
% jah_ppsl_array = zeros(Par.NumPart, 1);
curr_arr = zeros(Par.NumPart, 1);
% prev_arr = zeros(Par.NumPart, 1);


% Loop through clusters
for c = 1:length(Distn.clusters)
    
    % Loop through particles
    for ii = 1:Par.NumPart
        
        Cluster = Distn.clusters{c}.particles{ii};
        
%         % Project all tracks with an ML prediction
%         Cluster.ProjectTracks(t);
        
%         prev_post_prob = Posterior(t-1, L-1, Cluster, Observs);
%         prev_jah_ppsl = Cluster.SampleAssociations(t-1, L-1, Observs, true);
%         prev_state_ppsl = Cluster.SampleStates(t-1, L-1, Observs, true);

        % Sample and update associations
        jah_ppsl = Cluster.SampleAssociations(t, L, Observs, false);
    
        % Sample and update states
        state_ppsl = Cluster.SampleStates(t, L, Observs, false);
        
        % Calculate new posterior
        post_prob = Posterior(t, L, Cluster, Observs);
    
        % Update the weight
        weight(ii) = Distn.clusters{c}.weight(ii) ...
                   + (post_prob - sum(state_ppsl) - jah_ppsl);%...
%                    - (prev_post_prob - sum(prev_state_ppsl) - prev_jah_ppsl);
               
        curr_arr(ii) = (post_prob - sum(state_ppsl) - jah_ppsl);
%         prev_arr(ii) = (prev_post_prob - sum(prev_state_ppsl) - prev_jah_ppsl);
               
%         post_prob_array(ii) = post_prob;
%         state_ppsl_array(ii) = sum(state_ppsl(:));
%         jah_ppsl_array(ii) = jah_ppsl;
        
        if isnan(weight(ii))
            weight(ii) = -inf;
        end
        
        if isinf(weight(ii))
            disp(['Uh-oh: Zero weight in cluster ' num2str(c) ', particle ' num2str(ii)]);
        end
        
    end
    
    assert(~all(isinf(weight)), 'All weights are zero');
    
    % Normalise weights
    max_weight = max(weight); max_weight = max_weight(1); weight = weight - max_weight;
    weight = exp(weight); weight = weight/sum(weight);  weight = log(weight);
    
    % Attach weights to particles
    Distn.clusters{c}.weight = weight;
    
    % Calculate effective sample size for diagnostics
    ESS_pre(c) = CalcESS(weight);
    assert(~isnan(ESS_pre(c)), 'Effective Sample Size is non defined (probably all weights negligible)');
    
    if (ESS_pre(c) < 0.0*Par.NumPart)%&&(ESS_prev(c) < 0.05*Par.NumPart)
        Distn.clusters{c} = SystematicResample(Distn.clusters{c}, weight);
        resample = true;
        ESS(c) = Par.NumPart;
        disp(['*** Target Cluster' num2str(c) ': Effective Sample Size = ' num2str(ESS_pre(c)) '. RESAMPLED (Systematic). ESS = ' num2str(ESS(c))]);
    elseif (ESS_pre(c) < 0.1*Par.NumPart)%&&(ESS_prev(c) < 0.05*Par.NumPart)
        [Distn.clusters{c}, weight] = ConservativeResample(Distn.clusters{c}, weight);
        resample = true;
        ESS(c) = CalcESS(weight);
        disp(['*** Target Cluster' num2str(c) ': Effective Sample Size = ' num2str(ESS_pre(c)) '. RESAMPLED (Conservative). ESS = ' num2str(ESS(c))]);
    else
        resample = false;
        disp(['*** Target Cluster' num2str(c) ': Effective Sample Size = ' num2str(ESS_pre(c))]);
        ESS(c) = ESS_pre(c);
    end
    
end

% % Find total particle weight
% for ii = 1:Par.NumPart
%     num_assign = 0;
%     for c = 1:length(Distn.clusters)
%         for j = 1:Distn.clusters{c}.particles{ii}.N
%             if Distn.clusters{c}.particles{ii}.tracks{j}.GetAssoc(t) ~= 0
%                 num_assign = num_assign + 1;
%             end
%         end
%     end
%     
%     num_clut = Observs(t).N - num_assign;
%     assoc_prob = log( poisspdf(num_clut, Par.ExpClutObs) );
%     
%     weight(ii) = Distn.weight(ii) + assoc_prob;
%     
% end
% 
% % Normalise weights
% max_weight = max(weight); max_weight = max_weight(1); weight = weight - max_weight;
% weight = exp(weight); weight = weight/sum(weight);  weight = log(weight);
% 
% % Attach weights to particles
% for ii = 1:Par.NumPart
%     Distn.weight(ii) = weight(ii);
% end
% 
% % Calculate effective sample size for diagnostics
% ESS_pre = CalcESS(weight);
% assert(~isnan(ESS_pre), 'Effective Sample Size is non defined (probably all weights negligible)');
% disp(['*** Combined MTT: Effective Sample Size = ' num2str(ESS_pre)]);

end






% % Loop through particles
% for ii = 1:Par.NumPart
%     
%     Part = Distn.particles{ii};
%     
%     % % Calculate outgoing posterior probability term
%     % prev_post_prob = Posterior(t-1, L-1, Part, Observs);
%     
%     % Extend all tracks with an ML prediction
%     Part.ProjectTracks(t);
%     
%     % Sample associations
%     jah_ppsl = Part.SampleAssociations(t, L, Observs, BirthSites);
%     
%     state_ppsl = zeros(Par.NumTgts, 1);
%     NewTracks = cell(Par.NumTgts, 1);
%     
%     % Loop through targets
%     for j = 1:Part.N
%         
%         % Only need examine those which are present after t-L
%         if Part.tracks{j}.death > t-L+1
%             
%             % How long should the KF run for?
%             last = min(t, Part.tracks{j}.death - 1);
%             first = max(t-L+1, Part.tracks{j}.birth+1);
%             num = last - first + 1;
%             
%             % Draw up a list of associated hypotheses
%             obs = ListAssocObservs(last, num, Part.tracks{j}, Observs);
%             
%             % Run a Kalman filter the target
%             [KFMean, KFVar] = KalmanFilter(obs, Part.tracks{j}.GetState(first-1), Par.KFInitVar*eye(4));
%             
%             % Sample Kalman filter
%             [NewTracks{j}, state_ppsl(j)] = SampleKalman(KFMean, KFVar);
%             
%             % Update distribution
%             Part.tracks{j}.Update(last, NewTracks{j}, []);
%             
%         end
%         
%     end
%     
%     % Calculate new posterior
%     post_prob = Posterior(t, L, Part, Observs);
% %     post_prob = Posterior(t-1, L-1, Part, Observs);
%     
%     % Update the weight
%     weight(ii) = Distn.weight(ii) ...
%                + (post_prob)... - prev_post_prob) ...
%                - (sum(state_ppsl) + jah_ppsl);
%            
% 	  Diagnostics.post_arr(ii) = post_prob;
% %     Diagnostics.prev_post_arr(ii) = prev_post_prob;
%     Diagnostics.state_ppsl_arr(ii) = sum(state_ppsl);
%     Diagnostics.jah_ppsl_arr(ii) = jah_ppsl;
%     Diagnostics.weight_arr(ii) = weight(ii);
%            
%     if isnan(weight(ii))
%         weight(ii) = -inf;
%     end
%     
%     if isinf(weight(ii))
%         disp(['zero weight in particle ' num2str(ii)]);
%     end
%     
% end
% 
% assert(~all(isinf(weight)), 'All weights are zero');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% DIAGNOSTICS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Diagnostics.post_arr = sort(Diagnostics.post_arr, 'descend');
% % % Diagnostics.prev_post_arr = sort(Diagnostics.prev_post_arr, 'descend');
% % Diagnostics.state_ppsl_arr = sort(Diagnostics.state_ppsl_arr, 'descend');
% % Diagnostics.jah_ppsl_arr = sort(Diagnostics.jah_ppsl_arr, 'descend');
% % Diagnostics.weight_arr = sort(Diagnostics.weight_arr, 'descend');
% % 
% % Diagnostics.post_arr(1) - Diagnostics.post_arr(2)
% % Diagnostics.prev_post_arr(1) - Diagnostics.prev_post_arr(2)
% % Diagnostics.state_ppsl_arr(1) - Diagnostics.state_ppsl_arr(2)
% % Diagnostics.jah_ppsl_arr(1) - Diagnostics.jah_ppsl_arr(2)
% % Diagnostics.weight_arr(2) - Diagnostics.weight_arr(2)
% 
% % std(Diagnostics.post_arr)
% % std(Diagnostics.prev_post_arr)
% % std(Diagnostics.state_ppsl_arr)
% % std(Diagnostics.jah_ppsl_arr)
% % std(Diagnostics.weight_arr)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Normalise weights
% max_weight = max(weight);
% max_weight = max_weight(1);
% weight = weight - max_weight;
% weight = exp(weight);
% weight = weight/sum(weight);
% weight = log(weight);
% 
% % Attach weights to particles
% for ii = 1:Par.NumPart
%     Distn.weight(ii) = weight(ii);
% end
% 
% % Test effective sample size
% ESS_pre = CalcESS(weight);
% assert(~isnan(ESS_pre), 'Effective Sample Size is non defined (probably all weights negligible)');
% 
% % PlotTracks(Distn)
% % uiwait
% 
% if (ESS_pre < 0.5*Par.NumPart)
%     % Resample
%     ResamDistn = SystematicResample(Distn, weight);
%     if Par.FLAG_ResamMove
%         Distn = MoveMCMC(t, L, ResamDistn, Distn, Observs);
%     else
%         Distn = ResamDistn;
%     end
%     ESS = Par.NumPart;
%     resample = true;
%    
% %     Distn = SecondarySampling(t, L, Distn, Observs);
%     
% else
%     resample = false;
%     ESS = ESS_pre;
%     
% end
% 
% 
% 
% % PlotTracks(Distn)
% % uiwait
% 
% end