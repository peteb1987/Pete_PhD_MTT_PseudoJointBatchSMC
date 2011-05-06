function [ PL_est ] = PDAFPredLike(t, L, Cluster, Observs)
% Run a PDAF through the window and estimate the prediction likelihood,
% i.e. p(y_{t-L+1:t}|x_{t-L})

global Par;

PL_est = 0;

% Loop over targets in cluster
for j = 1:Cluster.N
    
    init_state = Cluster.tracks{j}.GetState(t-L);
    init_covar = Par.KFInitVar*eye(4);
    
    state_arr = cell(L+1, 1);
    covar_arr = cell(L+1, 1);
    
    state_arr{1} = init_state;
    covar_arr{1} = init_covar;
    
    % Loop through time
    for tt = t-L+1:t
        k = tt - (t-L-1);
        
        [state_arr{k}, covar_arr{k}, pred_like] = PDAFFrame(tt, state_arr{k-1}, covar_arr{k-1}, Observs(tt));
        PL_est = PL_est + pred_like;
        
    end
    
end



end

function [est_state, est_covar, pred_like] = PDAFFrame(t, prev_state, prev_covar, Observs)

global Par;

% Predict forwards a step
state = Par.A*prev_state;
covar = Par.A*prev_covar*Par.A' + Par.Q;

% Initialise stuff
GatedObs = Observs;
x1 = state(1); x2 = state(2);

if Par.FLAG_ObsMod == 0
    
    C = Par.C;
    gate_mean = state(1:2);
    
elseif Par.FLAG_ObsMod == 1
    
    C = zeros(2, 4);
    C(1,1) = -x2/(x1^2+x2^2);
    C(1,2) = x1/(x1^2+x2^2);
    C(2,1) = x1/sqrt(x1^2+x2^2);
    C(2,2) = x2/sqrt(x1^2+x2^2);
    
    [bng, rng] = cart2pol(x1, x2);
    gate_mean = [bng; rng];
    
end
    
gate_covar = C * covar * C' + Par.R;
gate_covar = 0.5*(gate_covar+gate_covar');

% Gate observations

thresh1 = 10*sqrt(gate_covar(1,1));
thresh2 = 10*sqrt(gate_covar(2,2));

innov = bsxfun(@minus, GatedObs.r, gate_mean')';
if Par.FLAG_ObsMod == 1
    wrap_around1 = innov(1,:)>pi; innov(1, wrap_around1) = innov(1, wrap_around1) - 2*pi;
    wrap_around2 = innov(1,:)<-pi; innov(1, wrap_around2) = innov(1, wrap_around2) + 2*pi;
end
test1 = abs(innov(1, :)) < thresh1;
test2 = abs(innov(2, :)) < thresh2;
indexes = find(test1&test2);
test = false(1, GatedObs.N);
for i = indexes
    test(i) = ((innov(:,i)'/gate_covar)*innov(:,i) < 100);
end
GatedObs.r = GatedObs.r(test, :);
GatedObs.N = size(GatedObs.r, 1);

% for i = GatedObs.N:-1:1
%     innov = GatedObs.r(i, :)' - gate_mean;
%     dist = sqrt( (innov' / gate_covar) * innov );
%     if dist > 4 % 4 gives 0.999
%         GatedObs.r(i, :) = [];
%         GatedObs.N = GatedObs.N - 1;
%     end
% end

% Run a KF on each gated observation
mean_arr = zeros(4, GatedObs.N);
covar_arr = zeros(4, 4, GatedObs.N);
assoc_arr = zeros(GatedObs.N, 1);
% Loop through valid associations
for i = 1:GatedObs.N
    
    % Calculate association probabilities
    assoc_arr(i) = mvnpdf(GatedObs.r(i, :), gate_mean', gate_covar);
    
    % Kalman filter
    [m, p] = KalmanFilter( {GatedObs.r(i, :)'}, prev_state, prev_covar );
    mean_arr(:, i) = m{1};
    covar_arr(:, :, i) = p{1};
    
end

assoc_clut = Par.ClutDens * (1-Par.PDetect) / Par.PDetect;
mean_clut = state;
covar_clut = covar;

% Normalise
norm_const = sum(assoc_arr) + assoc_clut;
assoc_arr = assoc_arr / norm_const;
assoc_clut = assoc_clut / norm_const;

% Combine
innov = zeros(2, GatedObs.N);
tot_innov = zeros(2,1);
for i = 1:GatedObs.N
    innov(:, i) = GatedObs.r(i, :)' - gate_mean;
    tot_innov = tot_innov + assoc_arr(i) * innov(:, i);
end
weight = covar * C' / gate_covar;
est_state = state + weight * tot_innov;

middle_bit = - tot_innov*tot_innov';
for i = 1:GatedObs.N
    middle_bit = middle_bit + assoc_arr(i)*innov(:,i)*innov(:,i)';
end

if GatedObs.N > 0
    est_covar = assoc_clut*covar + (1-assoc_clut)*covar_arr(:,:,1)...
        + weight * middle_bit * weight';
else
    est_covar = covar;
end

est_covar = 0.5*(est_covar+est_covar');

pred_like = Par.ClutDens / (GatedObs.N+1);
for i = 1:GatedObs.N
    pred_like = pred_like + mvnpdf(GatedObs.r(i,:), (C*est_state)', Par.R+C*est_covar*C') / (GatedObs.N+1);
end
pred_like = log(pred_like);

end