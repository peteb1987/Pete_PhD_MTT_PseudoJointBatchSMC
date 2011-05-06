function [ Distns, ObsTargIndexes, ESS_post, ESS_pre, num_resamples ] = MultiTargetTrack( detections, Observs, InitState )
%MULTITARGETTRACK Runs a batch SMC tracking algorithm for multiple
% targets with missed observations and clutter

global Par;

% Initialise particle array and diagnostics (an array of particle arrays)
Distns = cell(Par.T, 1);
SearchDistns = cell(Par.T, 1);
ObsTargIndexes = cell(Par.T, 1);
Area = cell(Par.T, 1);
ESS_post = cell(Par.T, Par.NumTgts);
ESS_pre = cell(Par.T, Par.NumTgts);
num_resamples = zeros(Par.NumTgts, 1);

if ~Par.FLAG_InitTargs
    
    % Initialise particle set - no targets
    init_track_distn = cell(0, 1);
    InitEst = MTDistn(init_track_distn);
    
else
    
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
    
end

% Initialise search track
search_track = Track(0, 0, [], []);
search_track_set = TrackSet(0, {search_track});
SearchTrack = TrackGroupDistn(0, repmat({search_track_set}, Par.NumPart, 1)); 

% Loop through time
for t = 1:Par.T
    
    tic;
    
    disp('**************************************************************');
    disp(['*** Now processing frame ' num2str(t)]);
    
    if t==1
        [Distns{t}, ESS_post{t}, ESS_pre{t}, resamples, ObsTargIndexes, SearchDistns{t}] = BatchSMC(t, t, InitEst, Observs, ObsTargIndexes, SearchTrack);
    elseif t<Par.L
        [Distns{t}, ESS_post{t}, ESS_pre{t}, resamples, ObsTargIndexes, SearchDistns{t}] = BatchSMC(t, t, Distns{t-1}, Observs, ObsTargIndexes, SearchDistns{t-1});
    else
        [Distns{t}, ESS_post{t}, ESS_pre{t}, resamples, ObsTargIndexes, SearchDistns{t}] = BatchSMC(t, Par.L, Distns{t-1}, Observs, ObsTargIndexes, SearchDistns{t-1});
    end
    
%     num_resamples = num_resamples + resamples;
    
    assoc = [];
    for c = 1:Distns{t}.N
        for j = 1:Distns{t}.clusters{c}.N
            get_ass = cellfun(@(x) x.tracks{j}.GetAssoc(t-min(t,Par.L)+1), Distns{t}.clusters{c}.particles);
            mode_ass = mode(get_ass);
            assoc = [assoc, mode_ass];
%             assoc = [assoc, Distns{t}.clusters{c}.particles{round(Par.NumPart/2)}.tracks{j}.GetAssoc(t-min(t,Par.L)+1)];
        end
    end
    disp(['*** Correct associations at frame ' num2str(t-min(t,Par.L)+1) ': ' num2str(detections(t-min(t,Par.L)+1,:))]);
    disp(['*** Modal associations at frame ' num2str(t-min(t,Par.L)+1) ': ' num2str(assoc)]);
    
    disp(['*** Frame ' num2str(t) ' processed in ' num2str(toc) ' seconds']);
    disp('**************************************************************');
    
    if mod(t, 100)==0
        PlotTracks(Distns{t});
%         plot(Observs(t).r(:, 2).*cos(Observs(t).r(:, 1)), Observs(t).r(:, 2).*sin(Observs(t).r(:, 1)), 'x', 'color', [1,0.75,0.75]);
%         saveas(gcf, ['Tracks' num2str(t) '.eps'], 'epsc2');
%         close(gcf)
        pause(1);
    end
    
end




end



function [Distn, ESS_post, ESS_pre, resamples, ObsTargIndexes, SearchTrack] = BatchSMC(t, L, Previous, Observs, ObsTargIndexes, SearchTrackIn)
% Execute a step of the SMC batch sampler

% Runs an SMC on a batch of frames. We can alter states in t-L+1:t. Frame
% t-L is also needed for filtering, transition density calculations, etc.

global Par;

ESS_post = zeros(Previous.N,1);
ESS_pre = zeros(Previous.N,1);
resamples = zeros(Previous.N,1);

% Initialise particle array and weight array
Distn = Previous.Copy;
SearchTrack = SearchTrackIn.Copy;

% Loop through clusters
for c = 1:Distn.N
    % Loop through particles
    for ii = 1:Par.NumPart
        % Project all tracks with an ML prediction
        Distn.clusters{c}.particles{ii}.ProjectTracks(t);
    end
end

ProjectedPrevious = Distn.Copy;

% Detect relevant observations - NOT USED. FUDGED REMOVAL
[ObsTargIndexes{t}] = DetectNearbyObservations(t, Observs, Distn);%, Area{t}
% disp(['*** ' num2str(cellfun(@length, ObsTargIndexes{t})') ' observations in respective target vicinities']);
Cluster_OTI = cellfun(@(x) x(1), ObsTargIndexes(1:t));

post_prob_array = zeros(Par.NumPart, 1);
state_ppsl_array = zeros(Par.NumPart, 1);
jah_ppsl_array = zeros(Par.NumPart, 1);
curr_arr = zeros(Par.NumPart, 1);
prev_arr = zeros(Par.NumPart, 1);
% pred_like_arr = zeros(Par.NumPart, 1);
d_arr = zeros(Par.NumPart, 1);

% Collision detection loop
clusters_done = false(Distn.N, 1);
ass_used = cell(1,L);
while ~all(clusters_done)
    
    % Loop through clusters
    for c = find(~clusters_done)'
        
        weights = zeros(Par.NumPart, 1);
        
        % Loop through particles
        for ii = 1:Par.NumPart
            
            Cluster = Distn.clusters{c}.particles{ii};
            
%             d = unidrnd(L);
%             d = L;
            
            prev_post_prob = Posterior(t-1, L-1, Cluster, Observs);
%             prev_jah_ppsl = Cluster.SampleAssociations(t-1, d-1, Observs, Cluster_OTI, true);
%             prev_state_ppsl = Cluster.SampleStates(t-1, d-1, Observs, true);
            
            % Calculate artfical density using PDAF approximation
%             PDAF_prob = 0;
%             for d = 1:L-1
            PDAF_prob = PDAF_prob + PDAFEstimatePosterior(t-1, d, Cluster, Observs);
%             end

%             Cluster_OTI = cellfun(@(x) x(c), ObsTargIndexes(1:t));
            
            d = unidrnd(L);
            d_arr(ii) = d;
            
            % Sample and update associations
            jah_ppsl = Cluster.SampleAssociations(t, d, Observs, Cluster_OTI, false);
            
            % Sample and update states
            state_ppsl = Cluster.SampleStates(t, d, Observs, false);
            
            % Calculate new posterior
            post_prob = Posterior(t, L, Cluster, Observs, Cluster_OTI);%, Area

            % Update the weight
            weights(ii) = Distn.clusters{c}.weights(ii) ...
                + (post_prob - sum(state_ppsl) - jah_ppsl) ...
                - (prev_post_prob - PDAF_prob);
%                 - (prev_post_prob - sum(prev_state_ppsl) - prev_jah_ppsl);
%                 - PL_est;

            
            curr_arr(ii) = (post_prob - sum(state_ppsl) - jah_ppsl);
%             prev_arr(ii) = (prev_post_prob - sum(prev_state_ppsl) - prev_jah_ppsl);
%             prev_arr(ii) = (prev_post_prob - PDAF_prob);

            post_prob_array(ii) = post_prob;
            state_ppsl_array(ii) = sum(state_ppsl(:));
            jah_ppsl_array(ii) = jah_ppsl;
            
            if isnan(weights(ii))
                weights(ii) = -inf;
            end
            
            if isinf(weights(ii))
%                 disp(['Uh-oh: Zero weight in cluster ' num2str(c) ', particle ' num2str(ii)]);
            end
            
        end
        
        assert(~all(isinf(weights)), 'All weights are zero');
        
        save_weights = weights;
        
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
%             [Distn.clusters{c}] = SystematicResample(Distn.clusters{c}, weights);
            resamples(c) = resamples(c) + 1;
            ESS_post(c) = CalcESS(weights);
            disp(['*** Target Cluster' num2str(c) ': Effective Sample Size = ' num2str(ESS_pre(c)) '. RESAMPLED (Conservative). ESS = ' num2str(ESS_post(c))]);
        else
            [Distn.clusters{c}, weights] = LowWeightRemoval(Distn.clusters{c}, weights);
            ESS_post(c) = CalcESS(weights);
            disp(['*** Target Cluster' num2str(c) ': Effective Sample Size = ' num2str(ESS_pre(c))]);
        end
        
    end
    
    % Detect Collisions
    [clusters_done, new_groups, new_groups_ind, ass_used] = DetectCollisions(t, L, Distn);
    
    if Par.FLAG_PseudoJoint
        if ~all(clusters_done)
            [Distn, clusters_done] = JoinClusters(Distn, ProjectedPrevious, new_groups, new_groups_ind);
        end
    else
        clusters_done = true;
    end
    
end

% Search track processing
if Par.FLAG_UseSearchTrack && t >= Par.BirthWindow
    
    % Find birth sites
    BirthSites = FindBirthSites( t, L, Observs, ass_used );
    
    if size(BirthSites, 1) > 0
        
        weights = zeros(Par.NumPart, 1);
        
        post_prob_array = zeros(Par.NumPart, 1);
        state_ppsl_array = zeros(Par.NumPart, 1);
        assoc_ppsl_array = zeros(Par.NumPart, 1);
        
        % Loop through search track particles
        for ii = 1:Par.NumPart
            
            SearchTrack.particles{ii}.ProjectTracks(t);
            
            % Propose associations
            assoc_ppsl = SearchTrack.particles{ii}.SampleSearchAssociations(t, L, Observs, Cluster_OTI, BirthSites );
            
            % Propse state
            state_ppsl = SearchTrack.particles{ii}.SampleStates(t, L, Observs, false);
            
            % Calculate new posterior
            post_prob = Posterior(t, L, SearchTrack.particles{ii}, Observs, Cluster_OTI);
            
            % Birth term
            if SearchTrack.particles{ii}.tracks{1}.birth > (t-L)
                birth_prob = log(poisspdf(1, Par.ExpBirth));
            else
                birth_prob = log(poisspdf(0, Par.ExpBirth));
            end
            
            % Update the weight
            weights(ii) = SearchTrack.weights(ii) + ...
                (post_prob + birth_prob - sum(state_ppsl) - assoc_ppsl);
            
            post_prob_array(ii) = post_prob;
            state_ppsl_array(ii) = sum(state_ppsl(:));
            assoc_ppsl_array(ii) = assoc_ppsl;
            
            if isnan(weights(ii))
                weights(ii) = -inf;
            end
            
            if isinf(weights(ii))
                disp(['Uh-oh: Zero weight in search track, particle ' num2str(ii)]);
            end
            
        end
        
        % Normalise weights
        max_weight = max(weights); max_weight = max_weight(1); weights = weights - max_weight;
        weights = exp(weights); weights = weights/sum(weights);  weights = log(weights);
        
        %     % Attach weights to particles
        %     SearchTrack.weights = weights;
        
        % Resample
%         PlotTracks( MTDistn({SearchTrack}) )
%         SearchTrack = SystematicResample(SearchTrack, weights);
        [SearchTrack, weights] = ConservativeResample(SearchTrack, weights);
%         PlotTracks( MTDistn({SearchTrack}) )
        
        origins=zeros(Par.NumPart,1); for ii=1:Par.NumPart, origins(ii)=SearchTrack.particles{ii}.members; end
        origins(origins==0)=[];
        most_common = mode(origins);
        proportion = sum(most_common==origins)/Par.NumPart;
        disp(['*** Most promising birth site has ' num2str(100*proportion) '% of the particles: ' num2str(most_common)]);
        
        % Promote search track if it has enough of the particles
        if (most_common ~= 0) && (proportion > Par.SearchPromoteThresh)
            
            % Count tracks already present
            id = 1;
            for c = 1:Distn.N
                for j = 1:Distn.clusters{c}.N
                    id = id + 1;
                end
            end
            
            % Copy the search track for promotion
            Promoted = SearchTrack.Copy;
            
            % Set the track id in each particle and remove those with the wrong id (set weight to 0)
            Promoted.members = id;
            for ii = 1:Par.NumPart
                if Promoted.particles{ii}.members ~= most_common
                    weights(ii) = -inf;
                end
                Promoted.particles{ii}.members = id;
            end
            
            % Normalise weights
            max_weight = max(weights); max_weight = max_weight(1); weights = weights - max_weight;
            weights = exp(weights); weights = weights/sum(weights);  weights = log(weights);
            [Promoted, weights] = LowWeightRemoval(Promoted, weights);
            
            Distn.clusters = [Distn.clusters; {Promoted}];
            Distn.N = Distn.N + 1;
            
%             % Remove the used particles from the search track
%             for ii = 1:Par.NumPart
%                 if SearchTrack.particles{ii}.members == most_common
%                     weights(ii) = -inf;
%                 end
%             end
%             [SearchTrack, weights] = LowWeightRemoval(SearchTrack, weights);
            
            % Reinitialise search track
            search_track = Track(0, 0, [], []);
            search_track_set = TrackSet(0, {search_track});
            SearchTrack = TrackGroupDistn(0, repmat({search_track_set}, Par.NumPart, 1));
            
        end
        
    end
    
end

if Par.FLAG_PseudoJoint
    % Detect Separations
    Distn = SeparateClusters(t, L, Distn);
end

end

