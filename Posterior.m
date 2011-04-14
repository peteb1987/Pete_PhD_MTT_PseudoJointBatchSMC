function [ post ] = Posterior(t, L, Set, Observs, Cluster_OTI)%, Area 
%POSTERIOR Calculate posterior probability of a TrackSet

global Par;

% Initialise probabilities
like = zeros(Set.N, 1);
trans = zeros(Set.N, 1);
clut = zeros(L, 1);
assoc = zeros(L, 1);
birth = zeros(L, 1);

% Loop through targets
for j = 1:Set.N
    
    end_time = min(t, Set.tracks{j}.death-1);
    start_time = max(t-L+1, Set.tracks{j}.birth);
    
    % Loop through window
    for tt = start_time:end_time
        
        % Get states
        state = Set.tracks{j}.GetState(tt);
        
        % Calculate likelihood
        if any(abs(state(3:4))>Par.Vlimit)||any(abs(state(1:2))>2*Par.Xmax)
            like(j) = -inf;
        else
            ass = Set.tracks{j}.GetAssoc(tt);
            if ass~=0
                if Par.FLAG_ObsMod == 0
                    like(j) = like(j) + log( mvnpdfFastSymm(Observs(tt).r(ass, :), state(1:2)', Par.ObsNoiseVar) );
                elseif Par.FLAG_ObsMod == 1
                    [bng, rng] = cart2pol(state(1), state(2));
                    if (Observs(tt).r(ass, 1) - bng) > pi
                        bng = bng + 2*pi;
                    elseif (Observs(tt).r(ass, 1) - bng) < -pi
                        bng = bng - 2*pi;
                    end
                    like(j) = like(j) + log( mvnpdf(Observs(tt).r(ass, :), [bng rng], diag(Par.R)') );
                end
            end
        end
        
        % Calculate transition density
        if tt==Set.tracks{j}.birth
            trans(j) = trans(j) + log(Par.UnifPosDens*Par.UnifVelDens);
        else
            prev_state = Set.tracks{j}.GetState(tt-1);
            trans(j) = trans(j) + log( (1-Par.PDeath) * mvnpdfQ(state', (Par.A * prev_state)') );
            % trans(j) = trans(j) + log( (1-Par.PDeath) * mvnpdf(state', (Par.A * prev_state)', Par.Q) );
        end
        
    end
    
    if (Set.tracks{j}.death <= t) && (Set.tracks{j}.death > t-L)
        trans(j) = trans(j) + log(Par.PDeath);
    end
    
end

% Clutter and association terms
for tt = t-L+1:t        %start_time:end_time
    k = tt - (t-L);
    
    % Association prior
    num_unassigned = Observs(tt).N;%length(ObsTargIndexes{t}{c});%sum(cellfun(@length, ObsTargIndexes{t}));%
    obs_assigned = [];
    for j = 1:Set.N
        if Set.tracks{j}.Present(tt)
            ass = Set.tracks{j}.GetAssoc(tt);
            if any(ass==obs_assigned)
                assoc(k) = -inf;
                break
            elseif ass>0
                assoc(k) = assoc(k) + log(Par.PDetect);
%                 assoc(k) = assoc(k) + log(Par.PDetect * num_unassigned);
%                 assoc(k) = assoc(k) + log(Par.PDetect/Par.NumTgts);
%                 assoc(k) = assoc(k) + log(Par.PDetect/num_unassigned);
                obs_assigned = [obs_assigned; ass];
                num_unassigned = num_unassigned - 1;
            elseif ass==0
%                  assoc(k) = assoc(k) + log((1-Par.PDetect) * Par.ExpClutObs*Par.ClutDens*sum(cell2mat(Area{t})));
                 assoc(k) = assoc(k) + log((1-Par.PDetect) * Par.ExpClutObs);
%                 assoc(k) = assoc(k) + log((1-Par.PDetect) * Par.ExpClutObs*Par.ClutDens*Area{t}{c});
%                 assoc(k) = assoc(k) + log(1-Par.PDetect);
            else
                error('Invalid association');
            end
        else
            assoc(k) = assoc(k) + log(Par.ExpClutObs);
        end
    end
    
    % Clutter term
%     num_targ = length(obs_assigned);
%     num_clut = Observs(tt).N - num_targ;
    num_missed = Set.N - length(obs_assigned);
    clut(k) = num_missed * log(Par.ClutDens);

    birth(k) = 0;
%     % Birth term
%     num_births = 0;
%     for j = 1:Set.N
%         if tt==Set.tracks{j}.birth
%             num_births = num_births + 1;
%         end
%     end
%     birth(k) = log(poisspdf(num_births, Par.ExpBirth));
    
end    

% if any(isinf(like))
%     disp('Zero likelihood');
% end
% if any(isinf(trans))
%     disp('Zero transition density');
% end
% if any(isinf(clut))
%     disp('Zero clutter probability');
% end
% if any(isinf(assoc))
%     disp('Zero association probability');
% end
% if any(isinf(birth))
%     disp('Zero birth probability');
% end

% Combine likelihood and transition density terms
post = sum(like) + sum(trans) + sum(clut) + sum(assoc) + sum(birth);
% post = sum(like) + sum(trans) + sum(assoc) + sum(birth);

end