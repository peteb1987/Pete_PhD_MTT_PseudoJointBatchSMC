function [ NewTrack, ppsl_prob ] = SampleKalman( KFMean, KFVar, track )
%SAMPLEKALMAN Sample from a set of Kalman estimates to generate a new track

% If a new track is supplied, just calculate its probability

global Par;

% Sample backwards from the end, as in Doucet, Briers, Senecal

L = size(KFMean,1); 

if L == 0
    NewTrack = [];
    ppsl_prob = 0;
    return
end

for k = 1:L
    KFVar{k} = (KFVar{k}+KFVar{k}')/2;
end

if nargin == 3
    NewTrack = track.state(end-L+1:end);
    
elseif nargin == 2
    NewTrack = cell(L, 1);
    
    % Last point
    NewTrack{end} = mvnrnd(KFMean{end}', KFVar{end})';
    
%     % % % FUDGE TO STOP OUT-OF-RANGE PROPOSALS
%     if any(abs(NewTrack{k}(3:4))>Par.Vlimit)
%         NewTrack{k}(3:4) = min(max(NewTrack{k}(3:4), -Par.Vlimit), Par.Vlimit);
%         disp('Proposal limited');
%     end
    
end

prob = zeros(size(KFMean));
prob(end) = mvnpdf(NewTrack{end}', KFMean{end}', KFVar{end});

% Loop though time
for k = L-1:-1:1
    
    norm_var = inv(Par.A' * (Par.Q \ Par.A) + inv(KFVar{k}));
    norm_mean = norm_var * (Par.A' * (Par.Q \ NewTrack{k+1}) + (KFVar{k} \ KFMean{k})); %#ok<MINV>
    
    norm_var = (norm_var+norm_var')/2;
    
    if nargin == 2
        NewTrack{k} = mvnrnd(norm_mean', norm_var)';
%         NewTrack{k} = mvnrnd(KFMean{k}', KFVar{k})';
        
%         % % % FUDGE TO STOP OUT-OF-RANGE PROPOSALS
%         if any(abs(NewTrack{k}(3:4))>Par.Vlimit)
%             NewTrack{k}(3:4) = min(max(NewTrack{k}(3:4), -Par.Vlimit), Par.Vlimit);
%             disp('Proposal limited');
%         end

    end
    
    prob(k) = mvnpdf(NewTrack{k}', norm_mean', norm_var);
%     prob(k) = mvnpdf(NewTrack{k}', KFMean{k}', KFVar{k});

end

ppsl_prob = sum(log(prob));

end