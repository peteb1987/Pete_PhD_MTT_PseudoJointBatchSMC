classdef TrackSet < handle
    %TRACKSET A set of track objects, comprising a multi-target state over
    %the duration of a scene
    
    properties
        tracks          % Cell array of track objects
        N               % Number of tracks
        
    end
    
    methods
        
        % Constructor
        function obj = TrackSet(tracks)
            obj.tracks = tracks;
            obj.N = length(tracks);
        end %Constructor
        
        % Copy
        function new = Copy(obj)
            t = cell(size(obj.tracks));
            for k = 1:length(t)
                t{k} = obj.tracks{k}.Copy;
            end
            new = TrackSet(t);
        end
        
        
        
        % Remove Track
        function RemoveTrack(obj, j)
            obj.tracks(j) = [];
            obj.N = obj.N - 1;
        end
        
        
        
        % Add Track
        function AddTrack(obj, NewTrack)
            obj.tracks = [obj.tracks; {NewTrack}];
            obj.N = obj.N + 1;
        end
        
        
        
        %ProjectTracks - Projects each track in a TrackSet forward by one frame.
        function ProjectTracks(obj, t)
            
            global Par;
            
            % Loop through targets
            for j = 1:obj.N
                
                if obj.tracks{j}.death == t
                    % Only extend it if it dies in the current frame. If
                    % not, its expired completely.
                    
                    % % Check that the track ends at t-1
                    % assert(obj.tracks{j}.death==t, 'Track to be projected does not end at t-1');
                    
                    % Get previous state
                    prev_state = obj.tracks{j}.GetState(t-1);
                    
                    % Project it forward
                    state = Par.A * prev_state;
                    
                    % Extend the track
                    obj.tracks{j}.Extend(t, state, 0);
                    
                end
                
            end
            
        end
        
        
        
        % SampleStates - Propose changes to state and calculate
        % proposal probability
        function ppsl_prob = SampleStates(obj, t, L, Observs, FLAG_ns)
            
            global Par;
            
            ppsl_prob = zeros(obj.N, 1);
            NewTracks = cell(obj.N, 1);
            
            % Loop through targets
            for j = 1:obj.N
                
                % Only need examine those which are present after t-L
                if obj.tracks{j}.death > t-L+1
                    
                    % How long should the KF run for?
                    last = min(t, obj.tracks{j}.death - 1);
                    first = max(t-L+1, obj.tracks{j}.birth+1);
                    num = last - first + 1;
                    
                    % Draw up a list of associated hypotheses
                    obs = ListAssocObservs(last, num, obj.tracks{j}, Observs);
                    
                    % Run a Kalman filter the target
                    [KFMean, KFVar] = KalmanFilter(obs, obj.tracks{j}.GetState(first-1), Par.KFInitVar*eye(4));
                    
                    % Sample Kalman filter
                    if ~FLAG_ns
                        [NewTracks{j}, ppsl_prob(j)] = SampleKalman(KFMean, KFVar);
                    else
                        track = obj.tracks{j}.Copy;
                        track.state(end) = [];
                        [NewTracks{j}, ppsl_prob(j)] = SampleKalman(KFMean, KFVar, track);
                    end
                    
                    if ~FLAG_ns
                        % Update distribution
                        obj.tracks{j}.Update(last, NewTracks{j}, []);
                    end
                    
                end
                
            end
            
            
            
        end
        
        
        
        %SampleAssociations - Propose changes to associations and calculate
        % proposal probability
        function ppsl_prob = SampleAssociations(obj, t, L, Observs, FLAG_ns)%, BirthSites )
            
            global Par;
            
            ppsl_prob = 0;
            
            jah_ppsl = SampleJAH(t, L, obj, Observs, FLAG_ns);
            ppsl_prob = ppsl_prob + sum(jah_ppsl(:));

            
%             % Propose target birth
%             if ~Par.FLAG_TargInit
%                 
%                 for j = 1:obj.N
%                     for tt = t-2:t
%                         for k = length(BirthSites):-1:1
%                             if any(BirthSites{k}==obj.tracks{j}.GetAssoc(tt))
%                                 BirthSites(k)=[];
%                             end
%                         end
%                     end
%                 end
%                 if ~isempty(BirthSites)&&(rand<Par.PAdd)
%                     
%                     ppsl_prob = ppsl_prob + log(Par.PAdd);
%                     
%                     % Select a random birth site
%                     k = unidrnd(length(BirthSites));
%                     
%                     % Construct a new track start point
%                     NewStates = cell(3,1);
%                     NewStates{1} = zeros(4,1);
%                     NewStates{3} = zeros(4,1);
%                     if Par.FLAG_ObsMod == 0
%                         NewStates{3}(1:2) = Observs(t).r( BirthSites{k}(3), : )';
%                         vel = (Observs(t).r( BirthSites{k}(3), : ) - Observs(t-1).r( BirthSites{k}(2), : )) / Par.P;
%                         NewStates{3}(3:4) = vel';
%                     elseif Par.FLAG_ObsMod == 2
%                         
%                         % Set final state (for projecting forward)
%                         xt = Pol2Cart(Observs(t).r( BirthSites{k}(3), 1 ), Observs(t).r( BirthSites{k}(3), 2 ));
%                         xt_1 = Pol2Cart(Observs(t-1).r( BirthSites{k}(2), 1 ), Observs(t-1).r( BirthSites{k}(2), 2 ));
%                         NewStates{3}(1:2) = xt;
%                         NewStates{3}(3:4) = (xt-xt_1)/Par.P;
%                         
%                         % Set first state (for initialising KF)
%                         xt = Pol2Cart(Observs(t-2).r( BirthSites{k}(1), 1 ), Observs(t-2).r( BirthSites{k}(1), 2 ));
%                         xt_1 = Pol2Cart(Observs(t-1).r( BirthSites{k}(2), 1 ), Observs(t-1).r( BirthSites{k}(2), 2 ));
%                         NewStates{1}(1:2) = xt;
%                         NewStates{1}(3:4) = (xt_1-xt)/Par.P;
%                         
%                     end
%                     
%                     % Create and add the new track
%                     NewTrack = Track(t-2, t+1, NewStates, BirthSites{k}');
%                     obj.AddTrack(NewTrack);
%                     
%                 end
%                 
%             end
%             
%             % Propose target death
%             if ~Par.FLAG_NoDeath
%                 
%                 if (t>2)
%                     for j = 1:obj.N
%                         if (obj.tracks{j}.Present(t)) ...
%                                 && (obj.tracks{j}.GetAssoc(t)==0) ...
%                                 && (obj.tracks{j}.GetAssoc(t-1)==0) ...
%                                 && (obj.tracks{j}.GetAssoc(t-2)==0)
%                             if rand < Par.PRemove
%                                 
%                                 ppsl_prob = ppsl_prob + log(Par.PRemove);
%                                 
%                                 tt = t-2;
%                                 while (obj.tracks{j}.GetAssoc(tt)==0)
%                                     if (tt==0)||(tt==t-L), break; end
%                                     tt = tt-1;
%                                 end
%                                 tt = tt+1;
%                                 obj.tracks{j}.EndTrack(tt);
%                             end
%                             
%                         end
%                     end
%                 end
%                 
%             end
            
%             jah_ppsl = SampleJAH(t, obj, Observs);
%             ppsl_prob = ppsl_prob + sum(jah_ppsl);
            
%             % Dig out target states
%             states = cell(obj.N, 1);
%             for j = obj.N:-1:1
%                 if obj.tracks{j}.Present(t)
%                     states{j} = obj.tracks{j}.GetState(t)';
% %                     states{j} = mvnrnd((Par.A * obj.tracks{j}.GetState(t-1))', Par.Q);
% %                     ppsl_prob = ppsl_prob + log( mvnpdf(states{j}, (Par.A * obj.tracks{j}.GetState(t-1))', Par.Q) );
%                 else
%                     states(j) = [];
%                 end
%             end
%             
%             % Fudge to make sure its the right way round
%             if isempty(states)
%                 states = cell(0, 1);
%             end
%             
%             % Generate a list of observations
%             if Par.FLAG_ObsMod == 0
%                 obs = Observs(t).r;
%             elseif Par.FLAG_ObsMod == 2
%                 obs = Observs(t).r;
% %                 obs = zeros(size(Observs(t).r));
% %                 obs(:,1) = Observs(t).r(:,2).*cos(Observs(t).r(:,1));
% %                 obs(:,2) = Observs(t).r(:,2).*sin(Observs(t).r(:,1));
%             end
%             
%             % Auction algorithm for ML associations
%             assoc = AuctionAssoc( states, obs );
%             
%             % Set associations
%             i = 0;
%             for j = 1:obj.N
%                 if (obj.tracks{j}.Present(t))
%                     i = i + 1;
%                     obj.tracks{j}.SetAssoc(t, assoc(i));
%                 end
%             end
            
        end
        
    end
    
end

