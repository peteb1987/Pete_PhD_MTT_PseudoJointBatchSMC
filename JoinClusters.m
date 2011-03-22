function [Distn, clusters_done] = JoinClusters(IPPFDistn, UnprocessedDistn, new_groups, new_groups_ind)
%JOINCLUSTERS Lump two clusters together to form a joint distribution

global Par;

% Copy current distribution
Distn = IPPFDistn.Copy;

clusters_done = true(Distn.N, 1);

% Loop through clusters and delete those to be merged
for c = Distn.N:-1:1
    if any(c==cell2mat(new_groups_ind))
        Distn.clusters(c) = [];
        Distn.N = Distn.N - 1;
        clusters_done(c) = [];
    end
end

% Loop through new groups and create new clusters from old distn
for c_new = 1:length(new_groups)
    
    members = new_groups{c_new};
    particles = [];
    
    for c_old = new_groups_ind{c_new}'
        
        % Resample the old distn so we can combine particles
        OldCluster = SystematicResample(UnprocessedDistn.clusters{c_old}, UnprocessedDistn.clusters{c_old}.weights);
        
        if isempty(particles)
            particles = OldCluster.particles;
        else
            particles = JoinParticles(particles, OldCluster);
        end
        
    end
    
    NewCluster = TrackGroupDistn(members, particles);
    Distn.clusters = [Distn.clusters; {NewCluster}];
    Distn.N = Distn.N + 1;
    clusters_done = [clusters_done; false];
    
end

end

function particles = JoinParticles(particles, OldCluster)

global Par;

for ii = 1:Par.NumPart
    particles{ii}.tracks = [particles{ii}.tracks; OldCluster.particles{ii}.tracks];
    particles{ii}.N = particles{ii}.N + OldCluster.particles{ii}.N;
    particles{ii}.members = [particles{ii}.members; OldCluster.particles{ii}.members];
end

end
