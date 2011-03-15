function [ Distn ] = SecondarySampling( t, L, DistnIn, Observs )
%SECONDARYSAMPLING Do a secondary step resampling states only

global Par;

Distn = DistnIn.Copy;

for ii = 1:Par.NumPart

    Part = Distn.particles{ii};
    
    % Loop through targets
    for j = 1:Part.N
        
        % Only need examine those which are present after t-L
        if Part.tracks{j}.death > t-L+1
            
            % How long should the KF run for?
            last = min(t, Part.tracks{j}.death - 1);
            first = max(t-L+1, Part.tracks{j}.birth+1);
            num = last - first + 1;
            
            % Draw up a list of associated hypotheses
            obs = ListAssocObservs(last, num, Part.tracks{j}, Observs);
            
            % Run a Kalman filter the target
            [KFMean, KFVar] = KalmanFilter(obs, Part.tracks{j}.GetState(first-1), Par.KFInitVar*eye(4));
            
            % Sample Kalman filter
            [NewTracks{j}, state_ppsl(j)] = SampleKalman(KFMean, KFVar);
            
            % Update distribution
            Part.tracks{j}.Update(last, NewTracks{j}, []);
            
        end
        
    end
    
    % Calculate new posterior
    post_prob = Posterior(t, L, Part, Observs);
    
    % Update the weight
    weight(ii) = Distn.weight(ii) ...
               + (post_prob)... - prev_post_prob) ...
               - (sum(state_ppsl));
    
end
    
% Normalise weights
max_weight = max(weight);
max_weight = max_weight(1);
weight = weight - max_weight;
weight = exp(weight);
weight = weight/sum(weight);
weight = log(weight);

% Attach weights to particles
for ii = 1:Par.NumPart
    Distn.weight(ii) = weight(ii);
end

% Test effective sample size
ESS_pre = CalcESS(weight);
if ESS_pre < 0.5*Par.NumPart
    Distn = SystematicResample(Distn, weight);
end
    
end

