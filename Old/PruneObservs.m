function [ Pruned ] = PruneObservs(t, Observs, Distn )
%PRUNEOBSERVS Remove all observations which are not near a target to save
%checking them every time

global Par;

Ns = length(Distn.clusters);
States = zeros(Ns, 2);
Spread = zeros(Ns, 2);

% Create a grid for crude observation validation
Nbox = ceil(2*Par.Xmax / Par.Vlimit);
grid = zeros(Nbox);

% Loop through targets and fetch state
for j = 1:Ns
    Parts = zeros(Par.NumPart, 4);
    for ii = 1:Par.NumPart
        Parts(ii, :) = Distn.clusters{j}.particles{ii}.tracks{1}.GetState(t)';
%         x = ceil(Parts(ii, 1) / Par.Vlimit) + Nbox/2;
%         y = ceil(Parts(ii, 2) / Par.Vlimit) + Nbox/2;
%         grid(x-1:x+2, y-1:y+1) = 1;
    end
    Parts(:, 3:4) = [];
    
    
    
    States(j, :) = mean(Parts);
    Spread(j, :) = 0.5*range(Parts);
end

% Loop through observations and delete obvious clutter
for i = Observs.N:-1:1
    obs_cart = Pol2Cart(Observs.r(i, 1), Observs.r(i, 2));
%     x = ceil(obs_cart(1) / Par.Vlimit) + Nbox/2;
%     y = ceil(obs_cart(2) / Par.Vlimit) + Nbox/2;
%     if grid(x, y) == 0
%         Observs.r(i, :) = [];
%         Observs.N = Observs.N - 1;
%     end
    
    
    keep = false;
    for j = 1:Ns
%         innov = obs_cart' - States(j, 1:2)';
%         rng = sqrt(Spread(j,1)^2+Spread(j,2)^2);
%         P = (rng + Par.Vlimit)^2 * eye(2);
        if (abs(obs_cart(1)-States(j,1))<(Spread(j,1)+Par.Vlimit))&&(abs(obs_cart(2)-States(j,2))<(Spread(j,2)+Par.Vlimit))
%         if innov' / P * innov < 1
            keep = true;
            break;
        end
    end
    if ~keep
        Observs.r(i, :) = [];
        Observs.N = Observs.N - 1;
    end
    
    
end

Pruned = Observs;

end