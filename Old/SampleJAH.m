function [ jah_ppsl ] = SampleJAH( t, L, Set, Observs, FLAG_ns )
%SAMPLEJAH Probabilistically generates a joint association hypothesis for a
% single frame.

global Par;

jah_ppsl = zeros(Set.N, L);

R = Par.R;
Q = Par.Q;
invR = inv(Par.R);
invQ = inv(Par.Q);
A = Par.A;

% Generate a random permutation order
order = randperm(Set.N);

% List of used associations
used_ass = cell(L,1);

% Loop through targets
for j = order

    d = inf;
    
    % Re-propose point
%     z = unidrnd(L);
    z = L;
    
%     % Get state
%     x = Set.tracks{j}.GetState(t-L);
    

%%%%%%%%%%%%%%%%%%%%%%%

    % How long should the KF run for?
    last = t-z;
    first = max(t-(z+L)+1, Set.tracks{j}.birth+1);
    num = last - first + 1;
    
%     if (num==0)%||(rand<0.5)
        x = Set.tracks{j}.GetState(t-z);
%     else
%         
%         % Draw up a list of associated hypotheses
%         obs = ListAssocObservs(last, num, Set.tracks{j}, Observs);
%         
%         % Run a Kalman filter the target
%         [KFMean, KFVar] = KalmanFilter(obs, Set.tracks{j}.GetState(first-1), Par.KFInitVar*eye(4));
%         
%         x = KFMean{end};
%         
%     end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Loop through time
    for tt = t:-1:t-z+1
        
        k = tt - (t-z);
        N = Observs(tt).N;
        
        % Calculate deterministic prediction and jacobian
        p_x = (A^k) * x;
        [p_bng, p_rng] = Cart2Pol(p_x(1:2));
        p_rngsq = p_rng^2;
        J = [-p_x(2)/p_rngsq, p_x(1)/p_rngsq, 0, 0; p_x(1)/p_rng, p_x(2)/p_rng, 0, 0];
        
        % Calculate the mean and variance
        if d>L
            mu = [p_bng; p_rng];
            S = R + k*J*(k*Q)*J';
        else
            nV = (R+next_J*(d*Q)*next_J');
            invSig = J'*(R\J) + ((A^d)'*next_J'/nV)*next_J*(A^d) + invQ/k;
            invS = invR - ((R\J)/invSig)*J'/R;
            S = inv(invS);
            mu = [p_bng; p_rng] - J*p_x + invS \ ( (((R\J)*(invSig\(A^d)')*next_J'/nV)*(next_y-next_p_pol+next_J*next_p_x)) + (((R\(J/invSig))/(k*Q))*(A^k)*x) );
        end
        
        S = (S+S')/2;
        
%         if (d>z)&&(Set.tracks{j}.GetAssoc(t-z)==0)&&(t-z>0)
%             S = 10*S;
%         end
        
        thresh1 = 5*sqrt(S(1,1));
        thresh2 = 5*sqrt(S(2,2));
        
        % Precalculate some things to speed up the next loop, which is over
        % 100's of observations within a loop over 100's of particles (i.e.
        % the bottleneck of the whole program)
        innov = bsxfun(@minus, Observs(tt).r, mu')';
        wrap_around = innov(1,:)>pi;
        innov(1, wrap_around) = 2*pi - innov(1, wrap_around);
        test1 = abs(innov(1, :)) < thresh1;
        test2 = abs(innov(2, :)) < thresh2;
        indexes = find(test1&test2);
        test = zeros(1, Observs(tt).N);
        for i = indexes
            test(i) = ((innov(:,i)'/S)*innov(:,i) < 16);
        end
        validated = find(test);
%         disp(length(validated));
        
        % Calculate weights
        ppsl_weights = zeros(N+1, 1);
        for i = validated
            ppsl_weights(i) = (Par.PDetect) * mvnpdf(Observs(tt).r(i, :), mu', S);
        end
        
%         % Calculate weights
%         ppsl_weights = zeros(N+1, 1);
%         for i = 1:N
% %             obs_cart  = Pol2Cart(Observs(tt).r(i, 1), Observs(tt).r(i, 2));
% %             if Dist(obs_cart, p_x)<k*Par.Vmax
%             innov = Observs(tt).r(i, :)' - mu;
%             if innov(1)>pi
%                 innov(1) = 2*pi - innov(1);
%             end
%             if (innov(1)<thresh1)&&(innov(2)<thresh2)&&((innov'/S)*innov < 16)
%                 ppsl_weights(i) = (Par.PDetect) * mvnpdf(Observs(tt).r(i, :), mu', S);
%             else
%                 ppsl_weights(i) = 0;
%             end
%         end
        
        % Clutter
        ppsl_weights(N+1) = Par.ClutDens * (1-Par.PDetect);
        
        % Remove used ones
        ppsl_weights(used_ass{k}) = 0;
        
        % Normalise
        ppsl_weights = ppsl_weights/sum(ppsl_weights);
        
        % Set minimum value for clutter as 1-PDetect
        if (d>L)&&(ppsl_weights(N+1)<(1-Par.PDetect))
            ppsl_weights(1:N) = Par.PDetect*ppsl_weights(1:N)/sum(ppsl_weights(1:N));
            ppsl_weights(N+1) = 1-Par.PDetect;
        end
        
        if ~FLAG_ns
            % Sample
            ass = randsample(N+1, 1, true, ppsl_weights);
        else
           ass = Set.tracks{j}.GetAssoc(tt);
           if ass == 0
               ass = N+1;
           end
        end
        
        % Probability
        jah_ppsl(j, k) = log(ppsl_weights(ass));
        
        % Assign it
        if ass==N+1
            ass = 0;
        else
            used_ass{k} = [used_ass{k}; ass];
        end
        
        if ~FLAG_ns    
            Set.tracks{j}.SetAssoc(tt, ass);
        end
        
        % Store numbers for next step
        if ass==0
            d=d+1;
        else
            next_y = Observs(tt).r(ass, :)';
            next_p_pol = [p_bng; p_rng];
            next_p_x = p_x;
            next_J = J;
            d=1;
        end
                
    end
    
end

end
        
%     range_squ = range^2;
%     range = rng;
%     jac = [-x(2)/range_squ, x(1)/range_squ, 0, 0; x(1)/range, x(2)/range, 0, 0];
%     mean_obs = [bng; rng];
%     var_obs = Par.R + jac * Par.Q * jac';
%     
%     % Find the marginal likelihood of each observation, given the previous state
%     for i = 1:N
%         
%         obs_cart  = Pol2Cart(Observs(t).r(i, 1), Observs(t).r(i, 2));
%         mean_cart = Pol2Cart(mean_obs(1), mean_obs(2));
%         if Dist(obs_cart, mean_cart)<Par.Vmax
%             ppsl_weights(i) = (Par.PDetect/Observs(t).N) * mvnpdf(Observs(t).r(i, :), mean_obs', var_obs);
%             %                 ppsl_weights(i) = mvnpdfFastSymm(obs_cart, mean_cart, Par.AuctionVar);
%         else
%             ppsl_weights(i) = 0;
%         end
%         
%     end
%     
%     % Clutter
%     %         ppsl_weights(N+1) = Par.UnifPosDens;
%     ppsl_weights(N+1) = Par.ClutDens * (1-Par.PDetect);
%     
%     % Remove used ones
%     ppsl_weights(used_ass) = 0;
%     
%     % Normalise
%     ppsl_weights = ppsl_weights/sum(ppsl_weights);
%     
%     % Sample
%     ass = randsample(N+1, 1, true, ppsl_weights);
%     
%     % Probability
%     jah_ppsl(j) = log(ppsl_weights(ass));
%     
%     % Assign it
%     if ass==N+1
%         ass = 0;
%     else
%         used_ass = [used_ass; ass];
%     end
%     Set.tracks{j}.SetAssoc(t, ass);
%     
% end
% 
% end

% function d = Dist(x1, x2)
% d = sqrt((x1(1)-x2(1))^2+(x1(2)-x2(2))^2);
% end
% 
% function diff = BngDist(b1, b2)
% diff=abs(b1-b2);
% if diff>pi
%     diff = 2*pi - diff;
% end
% end