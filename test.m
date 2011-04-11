B = 4;
t = 5;
score = zeros(length(BirthSites), 1);

for ii = 1:length(BirthSites)
    
    x0 = zeros(4, 1);
    [x0(1), x0(2)] = pol2cart(Observs(t-B+1).r(BirthSites{ii}(1), 1), Observs(t-B+1).r(BirthSites{ii}(1), 2));
    [x1, y1] = pol2cart(Observs(t-B+2).r(BirthSites{ii}(2), 1), Observs(t-B+2).r(BirthSites{ii}(2), 2));
    x0(3) = x1 - x0(1);
    x0(4) = y1 - x0(2);
    
    pos_track = Track(t-B+1, t+1, cell(B,1), BirthSites{ii});
    
    % Draw up a list of associated hypotheses
    obs = ListAssocObservs(t, B, pos_track, Observs);
    
    % Run a Kalman filter the target
    [KFMean, KFVar] = KalmanFilter(obs, x0, Par.KFInitVar*eye(4));
    pos_track.state = KFMean;
    
    [ score(ii) ] = Posterior(t, B, TrackSet(1,{pos_track}), Observs, []);
    
end