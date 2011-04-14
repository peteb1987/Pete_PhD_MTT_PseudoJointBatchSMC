function [ ppsl_prob ] = SampleSearchJAH( t, L, Set, Observs, Cluster_OTI, BirthSites )
% Sample associations for the search track, i.e. some are drawn from a
% shortlist of "birth sites".

global Par;

R = Par.R;
Q = Par.Q;
invR = inv(Par.R);
invQ = inv(Par.Q);
A = Par.A;

B = Par.BirthWindow;

ProbNewSite = 9/20;
ProbContinue = 9/20;
ProbNoTarg = 2/20;


u = rand;
if u < ProbNewSite
    
    % Propose a new birth site
    i = unidrnd(length(BirthSites));
    
    Set.tracks{1}.birth = t-B+1;
    Set.tracks{1}.death = t+1;
    Set.tracks{1}.num = B;
    Set.tracks{1}.state = cell(B, 1);
    Set.tracks{1}.assoc = zeros(B, 1);
    Set.members = (t-1)*Par.NumBirthSites+i;
    
    for tt = t-B+1:t
        k = tt - (t-B);
        
        ass = BirthSites{i}(k);
        Set.tracks{1}.SetAssoc(tt, ass);
        Set.tracks{1}.SetState(tt, inf*ones(4, 1));
        
    end
    
    init_state = zeros(4, 1);
    if Par.FLAG_ObsMod == 0
        init_state(1:2) = Observs(t-B+1).r(BirthSites{i}(1), :)';
        next_pos = Observs(t-B+2).r(BirthSites{i}(2), :)';
    elseif Par.FLAG_ObsMod == 1
        init_state(1:2) = pol2cart( Observs(t-B+1).r(BirthSites{i}(1), 1), Observs(t-B+1).r(BirthSites{i}(1), 2) );
        next_pos = pol2cart( Observs(t-B+2).r(BirthSites{i}(2), 1), Observs(t-B+2).r(BirthSites{i}(2), 2) );
    end
    init_state(3:4) = (next_pos - init_state(1:2))/Par.P;
    Set.tracks{1}.SetState(t-B+1, init_state);
    
    ppsl_prob = log(1/length(BirthSites)) + log(ProbNewSite);
    
elseif u < ProbNewSite + ProbNoTarg
    
    % Propose that no target is born
    Set.tracks{1}.birth = 0;
    Set.tracks{1}.death = 0;
    Set.tracks{1}.num = 0;
    Set.tracks{1}.state = [];
    Set.tracks{1}.assoc = [];
    Set.members = 0;
    
    ppsl_prob = 0 + log(ProbNoTarg);
    
else
    
    % Propose a new association for time t (not the whole window), in the normal way
    
    if Set.members ~= 0
        
        x = Set.tracks{1}.GetState(t-1);
        
        N = Observs(t).N;
        
        % Calculate deterministic prediction and jacobian
        p_x = A*x;
        if Par.FLAG_ObsMod == 0
            C = [1 0 0 0; 0 1 0 0];
        elseif Par.FLAG_ObsMod == 1
            [p_bng, p_rng] = cart2pol(p_x(1), p_x(2));
            p_rngsq = p_rng^2;
            J = [-p_x(2)/p_rngsq, p_x(1)/p_rngsq, 0, 0; p_x(1)/p_rng, p_x(2)/p_rng, 0, 0];
        end
        
        % Calculate the mean and variance
        if Par.FLAG_ObsMod == 0
            mu = C*p_x;
            S = R + C*Q*C';
        elseif Par.FLAG_ObsMod == 1
            mu = [p_bng; p_rng];
            S = R + J*Q*J';
        end
        
        S = (S+S')/2;
        
        thresh1 = 5*sqrt(S(1,1));
        thresh2 = 5*sqrt(S(2,2));
        
        % Precalculate some things to speed up the next loop, which is over
        % 100's of observations within a loop over 100's of particles (i.e.
        % the bottleneck of the whole program)
        ind = Cluster_OTI{t};
        innov = bsxfun(@minus, Observs(t).r(ind,:), mu')';
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
            ppsl_weights(i) = (Par.PDetect) * mvnpdf(Observs(t).r(i, :), mu', S);
        end
        
        % Clutter
        ppsl_weights(N+1) = Par.ClutDens * (1-Par.PDetect);
        
        % Normalise
        ppsl_weights = ppsl_weights/sum(ppsl_weights);
        
        % Set minimum value for clutter as 1-PDetect
        if (ppsl_weights(N+1)<(1-Par.PDetect))
            ppsl_weights(1:N) = Par.PDetect*ppsl_weights(1:N)/sum(ppsl_weights(1:N));
            ppsl_weights(N+1) = 1-Par.PDetect;
        end
        
        % Sample
        ass = randsample(N+1, 1, true, ppsl_weights);
        
        % Probability
        ppsl_prob = log(ppsl_weights(ass));
        
        % Assign it
        if ass==N+1
            ass = 0;
        end
        Set.tracks{1}.SetAssoc(t, ass);
        
    else
        
        % No target, no change
        ppsl_prob = 0;
        
    end
    
    ppsl_prob = ppsl_prob + log(ProbContinue);
    
end

end
