function [ EstState, EstVar ] = KalmanFilter( obs, init_state, init_var )
%KALMANFILTER Kalman filter a set of observations to give a track estimate

%[ BackState, BackVar ]

global Par;

if size(obs, 2) == 0
    EstState = cell(0,1);
    EstVar = cell(0,1);
    return
end

L = length(obs(:,1));

PredState = cell(L, 1);
PredVar = cell(L, 1);

EstState = cell(L, 1);
EstVar = cell(L, 1);

% Set C depending on model
if Par.FLAG_ObsMod == 0
    C = Par.C;

elseif Par.FLAG_ObsMod == 1
    % Use EKF approximation
    C = zeros(2, 4);
    
end

% Loop through time
for k = 1:L
    
    % Prediction step
    
    if k==1
        PredState{1} = Par.A * init_state;
        PredVar{1} = Par.A * init_var * Par.A' + Par.Q;
    else
        PredState{k} = Par.A * EstState{k-1};
        PredVar{k} = Par.A * EstVar{k-1} * Par.A' + Par.Q;
    end
    
    % Update step
    
    if length(obs{k})==2
        
        % Observation associated with target
        
        if (Par.FLAG_ObsMod == 1)
            
            % Linearisation
            x1 = PredState{k}(1);
            x2 = PredState{k}(2);
            C(1,1) = -x2/(x1^2+x2^2);
            C(1,2) = x1/(x1^2+x2^2);
            C(2,1) = x1/sqrt(x1^2+x2^2);
            C(2,2) = x2/sqrt(x1^2+x2^2);
            
        end
        
        % Innovation
        if Par.FLAG_ObsMod == 0
            y = obs{k} - C * PredState{k};
        elseif Par.FLAG_ObsMod == 1
            [bng, rng] = cart2pol(PredState{k}(1), PredState{k}(2));
            y = obs{k} - [bng; rng];
            if y(1) > pi
                y(1) = y(1) - 2*pi;
            elseif y(1) < -pi
                y(1) = y(1) + 2*pi;
            end
        end
        s = C * PredVar{k} * C' + Par.R;
        gain = PredVar{k} * C' / s;
        EstState{k} = PredState{k} + gain * y;
        EstVar{k} = (eye(4)-gain*C) * PredVar{k};
        
    elseif isempty(obs{k})
        
        % No observation
        
        EstState{k} = PredState{k};
        EstVar{k} = PredVar{k};
        
    else
        
        error('Observation model allows for 0 or 1 observations only');
        
    end

end

% BackState = cell(L, 1);
% BackVar = cell(L, 1);
% BackState{k} = EstState{k}; BackVar{k} = EstVar{k};
% % Smoothing
% for k = L-1:-1:1
%     
%     G = EstVar{k}*Par.A'/PredVar{k+1};
%     BackState{k} = EstState{k} + G * (BackState{k+1}-PredState{k+1});
%     BackVar{k} = EstVar{k} + G * (BackVar{k+1}-PredVar{k+1}) * G';
%     
% end

end