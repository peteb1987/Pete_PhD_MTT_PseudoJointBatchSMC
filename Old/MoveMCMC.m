function [ Distn ] = MoveMCMC( t, L, DistnIn, PreResamDistn, Observs )
%MOVEMCMC Carry out an MCMC move on each particle

global Par;

% Copy distribution
Distn = DistnIn.Copy;

% Loop through particles
for ii = 1:Par.NumPart
    
    Part = Distn.particles{ii};

    for count = 1:5
    
    % Loop through targets
    for j = 1:Part.N
        
        ass_old = Distn.particles{ii}.tracks{j}.GetAssoc(t);
        
        % Propose a change to the association
        k = unidrnd(Par.NumPart);
        ass_new = PreResamDistn.particles{k}.tracks{j}.GetAssoc(t);
        
        if (ass_new ~= ass_old)%&&(ass_new ~= 0)
        
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
                [~, reverse_state_ppsl] = SampleKalman(KFMean, KFVar, Part.tracks{j});
                
                % Build the new track
                NewTrack = Part.tracks{j}.Copy;
                NewTrack.SetAssoc(t, ass_new);
                
                % Draw up a list of associated hypotheses
                obs = ListAssocObservs(last, num, NewTrack, Observs);
                
                % Run a Kalman filter the target
                [KFMean, KFVar] = KalmanFilter(obs, NewTrack.GetState(first-1), Par.KFInitVar*eye(4));
                
                % Sample Kalman filter
                [NewTrack, state_ppsl] = SampleKalman(KFMean, KFVar);
                
                % Posteriors
                reverse_post_prob = Posterior(t, L, Part, Observs);
                PartPpsl = Part.Copy;
                PartPpsl.tracks{j}.Update(last, NewTrack, []);
                post_prob = Posterior(t, L, PartPpsl, Observs);
                
                % Test acceptance
                acc_prob = (post_prob - reverse_post_prob) ...
                         + (reverse_state_ppsl - state_ppsl);
                     
                if log(rand) < acc_prob
                    Distn.particles{j} = PartPpsl.Copy;
                end
                
            end
            
        end
        
    end
    
    end
    
end

