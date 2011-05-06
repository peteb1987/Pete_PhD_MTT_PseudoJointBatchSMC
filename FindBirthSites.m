function BirthSites = FindBirthSites( t, L, Observs, ass_used )
%FINDBIRTHSITES Locate points where a target could have been born in this
% frame.

global Par;

B = Par.BirthWindow;

BirthSites = cell(0,1);

% Convert observations into cartesians
obs_x = cell(B, 1);
obs_y = cell(B, 1);
for k = 1:B
    tt = t-B+k;
    if Par.FLAG_ObsMod == 0
        obs_x{k} = Observs(tt).r(:,1);
        obs_y{k} = Observs(tt).r(:,2);
    elseif Par.FLAG_ObsMod == 1
        [obs_x{k}, obs_y{k}] = pol2cart(Observs(tt).r(:,1), Observs(tt).r(:,2));
        obs_x{k}(Observs(tt).r(:,2)<Par.BirthExclusionRadius) = inf;
        obs_y{k}(Observs(tt).r(:,2)<Par.BirthExclusionRadius) = inf;
    end
    
    % Remove all observations already used in a track
    k_L = tt - (t-L);
    if k_L > 0
        obs_x{k}(ass_used{k_L}) = inf;
        obs_y{k}(ass_used{k_L}) = inf;
    end
    
end

% Loop backwards through time looking for in-range obersvation pairs
ind = cell(B, 1);
for tt = t:-1:t-B+2
    k = tt - (t-B);
    
    if tt == t
        use = 1:length(obs_x{k});
    else
        use = unique(ind{k+1}(:,2));
    end
    
    % Calculate distance between points
    dist_ = sqrt( bsxfun(@minus, obs_x{k}, obs_x{k-1}').^2 + bsxfun(@minus, obs_y{k}, obs_y{k-1}').^2 );
    dist = inf*ones(size(dist_)); dist(use, :) = dist_(use, :);
    
    % Find points in range
    in_range = dist < (Par.P * Par.Vlimit);
    [later, earlier] = find(in_range);
    ind{k} = [later, earlier];
    
end

% Loop through tree constructing possibilities
pos_start_pts = unique(ind{2}(:,2));
for i = 1:length(pos_start_pts)
    BirthSites = TreeBranch(2, ind, pos_start_pts(i), BirthSites);
end

% Construct a score for each site
score = zeros(length(BirthSites), 1);
pos_track = Track(t-B+1, t+1, cell(B,1), zeros(B,1));
pos_track_set = TrackSet(1,{pos_track});

for ii = 1:length(BirthSites)
    
    x0 = zeros(4, 1);
    [x0(1), x0(2)] = pol2cart(Observs(t-B+1).r(BirthSites{ii}(1), 1), Observs(t-B+1).r(BirthSites{ii}(1), 2));
    [x1, y1] = pol2cart(Observs(t-B+2).r(BirthSites{ii}(2), 1), Observs(t-B+2).r(BirthSites{ii}(2), 2));
    x0(3) = x1 - x0(1);
    x0(4) = y1 - x0(2);
    
    pos_track.assoc = BirthSites{ii};
    
    % Draw up a list of associated hypotheses
    obs = ListAssocObservs(t, B, pos_track, Observs);
    
    % Run a Kalman filter the target
    [KFMean, KFVar] = KalmanFilter(obs, x0, Par.KFInitVar*eye(4));
    pos_track.state = KFMean;
    
    % Calculate the posterior probability of the track stub
    [ score(ii) ] = Posterior(t, B, pos_track_set, Observs, []);
    
%     if Par.FLAG_ObsMod == 1
%         [~, start_range] = cart2pol(KFMean{1}(1), KFMean{1}(2));
%         score(ii) = score(ii) + B*log(start_range);
%     end
    
end

% Sort sites and select most probable
sort_arr = [score, cumsum(ones(length(score), 1))];
sort_arr = sortrows(sort_arr, -1);
BirthSites = BirthSites(sort_arr(1:min(Par.NumBirthSites,length(BirthSites)), 2));

disp(['Found ' num2str(length(BirthSites)) ' birth sites']);

end

function BirthSites = TreeBranch(k, ind, hist, BirthSites)

recent = hist(end);
list = find(ind{k}(:,2)==recent);

for i = list'
    if k == length(ind)
        BirthSites = [BirthSites; {[hist, ind{k}(i,1)]}];
    else
        BirthSites = TreeBranch(k+1, ind, [hist, ind{k}(i,1)], BirthSites);
    end
end

end

