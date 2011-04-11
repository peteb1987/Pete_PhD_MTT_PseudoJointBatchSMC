function [ jah_ppsl ] = SampleJAH( t, L, Set, Observs, Cluster_OTI, FLAG_ns )
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
    
    % Look for missed-detection associations
    first_miss = 0;
    for tt = t-L+1:t-1
        if Set.tracks{j}.GetAssoc(tt) == 0
            first_miss = tt;
            break
        end
    end
    
    % With some probability, end the track
    if (rand < Par.PRemove) && (Set.tracks{j}.death == t+1) && (first_miss > 0)
        Set.tracks{j}.EndTrack(first_miss)
        jah_ppsl(j, L) = log(Par.PRemove);
    end
        
    % If the track ended before the window, do nothing, otherwise propose associations
    if Set.tracks{j}.death <= t-L+1
        jah_ppsl(j, L) = 0;
        
    else
        
        d = inf;
        
        z = L;
        
        x = Set.tracks{j}.GetState(t-z);
        
        last = min(t, Set.tracks{j}.death - 1);
        
        % Loop through time
        for tt = last:-1:t-z+1
            
            k = tt - (t-z);
            N = Observs(tt).N;
            
            % Calculate deterministic prediction and jacobian
            p_x = (A^k) * x;
            if Par.FLAG_ObsMod == 0
                C = [1 0 0 0; 0 1 0 0];
            elseif Par.FLAG_ObsMod == 1
                [p_bng, p_rng] = cart2pol(p_x(1), p_x(2));
                p_rngsq = p_rng^2;
                J = [-p_x(2)/p_rngsq, p_x(1)/p_rngsq, 0, 0; p_x(1)/p_rng, p_x(2)/p_rng, 0, 0];
            end
            
            % Calculate the mean and variance
            if Par.FLAG_ObsMod == 0
                if d>L
                    mu = C*p_x;
                    S = R + C*(k*Q)*C';
                else
                    nV = R+C*(d*Q)*C';
                    invSig = C'*(R\C) + ((A^d)'*C'/nV)*C*(A^d) + invQ/k;
                    invS = invR - ((R\C)/invSig)*C'/R;
                    S = inv(invS);
                    mu = invS \ ( (((R\C)*(invSig\(A^d)')*C'/nV)*next_y) + (((R\(C/invSig))/(k*Q))*(A^k)*x) );
                end
            elseif Par.FLAG_ObsMod == 1
                if d>L
                    mu = [p_bng; p_rng];
                    S = R + J*(k*Q)*J';
                else
                    nV = (R+next_J*(d*Q)*next_J');
                    invSig = J'*(R\J) + ((A^d)'*next_J'/nV)*next_J*(A^d) + invQ/k;
                    invS = invR - ((R\J)/invSig)*J'/R;
                    S = inv(invS);
                    mu = [p_bng; p_rng] - J*p_x + invS \ ( (((R\J)*(invSig\(A^d)')*next_J'/nV)*(next_y-next_p_pol+next_J*next_p_x)) + (((R\(J/invSig))/(k*Q))*(A^k)*x) );
                end
            end
            
            S = (S+S')/2;
            
            thresh1 = 5*sqrt(S(1,1));
            thresh2 = 5*sqrt(S(2,2));
            
            % Precalculate some things to speed up the next loop, which is over
            % 100's of observations within a loop over 100's of particles (i.e.
            % the bottleneck of the whole program)
            ind = Cluster_OTI{tt};
            innov = bsxfun(@minus, Observs(tt).r(ind,:), mu')';
            if Par.FLAG_ObsMod == 1
                wrap_around1 = innov(1,:)>pi; innov(1, wrap_around1) = innov(1, wrap_around1) - 2*pi;
                wrap_around2 = innov(1,:)<-pi; innov(1, wrap_around2) = innov(1, wrap_around2) + 2*pi;
            end
            test1 = abs(innov(1, :)) < thresh1;
            test2 = abs(innov(2, :)) < thresh2;
            indexes = find(test1&test2);
            test = false(1, length(ind));
            for i = indexes
                test(i) = ((innov(:,i)'/S)*innov(:,i) < 16);
            end
            validated = ind(test);
            
            % Calculate weights
            ppsl_weights = zeros(N+1, 1);
            for i = validated
                ppsl_weights(i) = (Par.PDetect) * mvnpdf(Observs(tt).r(i, :), mu', S);
            end
            
            % Clutter
            ppsl_weights(N+1) = Par.ClutDens * (1-Par.PDetect);
            
            % Remove used ones
            ppsl_weights(used_ass{k}) = 0;
            
            % Normalise
            ppsl_weights = ppsl_weights/sum(ppsl_weights);
            
            % Set minimum value for clutter as 1-PDetect
            %         if (d>L)&&(ppsl_weights(N+1)<(1-Par.PDetect))
            if (ppsl_weights(N+1)<(1-Par.PDetect))
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
            jah_ppsl(j, k) = jah_ppsl(j, k) + log(ppsl_weights(ass));
            
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
                d=1;
                if Par.FLAG_ObsMod == 0
                    
                elseif Par.FLAG_ObsMod == 1
                    next_p_pol = [p_bng; p_rng];
                    next_p_x = p_x;
                    next_J = J;
                end
            end
            
        end
        
    end
    
end

end

