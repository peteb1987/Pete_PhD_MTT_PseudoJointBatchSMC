function [ResamDistn, new_weight] = LowWeightRemoval( Distn, weight )
%Remove particles with very low weights and replace by systematic
%resampling

global Par;

Np = Par.NumPart;

% Change the weights to linear
weight = exp(weight);

% Copy distribution
ResamDistn = Distn.Copy;

% One child per input particle
parent = cumsum(ones(Np, 1));

% Get rid of low weights
thresh = max(weight) / exp(Par.ResampleLowWeightThresh);
parent(weight<thresh)=[];
weight(weight<thresh)=[];
Nct = length(parent);

% Systematically downsample
u = ceil( Nct * cumsum( (1/Np)*ones(Np,1) ) - unifrnd(0, 1/Np) );
samp = parent(u);
new_weight = weight(u);

new_weight = new_weight/sum(new_weight);
new_weight = log(new_weight);

% Construct new, resampled array
ResamDistn.particles = cell(Par.NumPart,1);
for ii = 1:Par.NumPart
    ResamDistn.particles{ii} = Distn.particles{samp(ii)};
    ResamDistn.weights(ii) = new_weight(ii);
end

% "Hard copy" the new distribution (we just copied the pointers in the loop, for speed)
ResamDistn = ResamDistn.Copy;

end