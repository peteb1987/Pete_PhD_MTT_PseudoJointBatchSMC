function [clusters_done, new_groups, new_groups_ind, all_ass_used] = DetectCollisions( t, L, Distn )
%DETECTCOLLISIONS Detect target "collisions", i.e. where two targets which
% are currently being tracked independently (IPPF-style) are associated
% with the same observation.

all_ass_used = cell(1, L);

clusters_done = true(Distn.N, 1);

% Work out current group structure
groups = cellfun(@(x){x.members}, Distn.clusters);
new_groups = cell(0,1);
new_groups_ind  = cell(0,1);

% Loop through time
for tt = t-L+1:t
    
    k = tt - (t-L);
    
    ass_used = cell(1, Distn.N);
    
    % Loop through clusters
    for c = 1:Distn.N
        ass = [];
        
        % Loop through targets and fetch association for each particle
        for j = 1:Distn.clusters{c}.N
            ass =  [ass, cellfun(@(x) x.tracks{j}.GetAssoc(tt), Distn.clusters{c}.particles)'];
        end
        
        % List unique, non-zero associations
        ass_used{c} = unique(ass);
        ass_used{c}(ass_used{c}==0)=[];
        
        % Loop through previous targets and detect collisions
        for cc = 1:c-1
            for i = 1:length(ass_used{c})
                if any(ass_used{c}(i)==ass_used{cc})
                    disp(['Collision between ' num2str(cc) ' and ' num2str(c) ' at observation ' num2str(ass_used{c}(i)) ' in frame ' num2str(tt)]);
                    clusters_done(c) = false;
                    clusters_done(cc) = false;
                    
                    % Create new cluster
                    new_cluster = [groups{cc}; groups{c}];
                    new_cluster_ind = [cc; c];
                    
                    % Check to see if there is already a new cluster with one of these clusters in
                    placed = false;
                    for g = 1:length(new_groups)
                        if length(new_groups_ind{g})==length(new_cluster_ind)
                            if all(new_groups_ind{g}==new_cluster_ind)
                                % Group already exists
                                placed = true;
                            end
                        end
                        
                        if ~isempty(intersect(new_groups_ind{g}, new_cluster_ind))&&~placed
                            % Groups overlaps with another
                            new_groups{g} = unique([new_groups{g}; new_cluster]);
                            new_groups_ind{g} = unique([new_groups_ind{g}; new_cluster_ind]);
                            placed = true;
                        end
                    end
                    if ~placed
                        new_groups = [new_groups; new_cluster];
                        new_groups_ind = [new_groups_ind; new_cluster_ind];
                    end
                    
                end
            end
        end

    end
    
    all_ass_used{k} = unique(cell2mat(ass_used));
    
end

end

