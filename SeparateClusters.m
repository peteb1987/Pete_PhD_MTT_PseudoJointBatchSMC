function Distn = SeparateClusters( t, L, Distn )
%SEPARATECLUSTERS Detect when targets in a cluster have moved apart and
%separate them. (Reintroduce independence assumption)

global Par;

indep = cell(Distn.N, 1);

% Loop through clusters
for c = 1:Distn.N
    
    done = false;
    
    % Don't consider breaking a cluster until after L time steps have passed
    if Distn.clusters{c}.locked > 0
        Distn.clusters{c}.locked = Distn.clusters{c}.locked - 1;
        break
    end
    
    indep{c} = true(Distn.clusters{c}.N, 1);
    
    % Loop through window
    for tt = t-L+1:t
        
        ass_used = cell(Distn.clusters{c}.N, 1);
        
        % Loop through targets
        for j = 1:Distn.clusters{c}.N
            
            % Fetch associations
            ass = cellfun(@(x) x.tracks{j}.GetAssoc(tt), Distn.clusters{c}.particles);
            ass_used{j} = unique(ass);
            ass_used{j}(ass_used{j}==0)=[];
            
            % Loop through targets so far looking for collisions
            for jj = 1:j-1
                if ~isempty(intersect(ass_used{j}, ass_used{jj}))
                    indep{c}(j) = false; indep{c}(jj) = false;
                    if all(indep{c}==false)
                        done = true;
                        break
                    end
                end
            end
            
            if done, break, end
            
        end
        
        if done, break, end
        
    end
    
    % Now break out each independent target
    for j = find(indep{c})'
        
        % Test that this isn't the only track left
        if Distn.clusters{c}.N > 1
            
            % Create a new cluster
            NewCluster = SeparateParticles(j, Distn.clusters{c});
            
            % Add it to the distn
            Distn.clusters = [Distn.clusters; {NewCluster}];
            Distn.N = Distn.N + 1;
            
            % Remove the track from the old cluster
            Cluster = Distn.clusters{c};
            for ii = 1:Par.NumPart
                Cluster.particles{ii}.tracks(j) = [];
                Cluster.particles{ii}.N = Cluster.particles{ii}.N - 1;
                Cluster.particles{ii}.members(j) = [];
            end
            Cluster.members(j) = [];
            Cluster.N = Cluster.N - 1;
            
        end
        
    end
    
end

end

function NewCluster = SeparateParticles(j, OldCluster)

global Par;

NewCluster = OldCluster.Copy;

NewCluster.members = NewCluster.members(j);
NewCluster.N = 1;

for ii = 1:Par.NumPart
    NewCluster.particles{ii}.tracks([1:j-1 j+1:end]) = [];
    NewCluster.particles{ii}.N = 1;
    NewCluster.particles{ii}.members = OldCluster.particles{ii}.members(j);
end

end
