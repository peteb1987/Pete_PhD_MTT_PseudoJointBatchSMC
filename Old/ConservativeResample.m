function [ResamDistn, new_weight] = ConservativeResample( Distn, weight )
%SYSTEMATICRESAMPLE Resample Distn systematically from weight

global Par;

Np = Par.NumPart;

% Change the weights to linear
weight = exp(weight);

% Copy distribution
ResamDistn = Distn.Copy;

% Create empty vector for offspring count
Nchild = zeros(Np, 1);

% Enumerate offspring
for ii = 1:Np
    Nchild(ii) = max(1,floor(  1 *  weight(ii)*Np));
end
Nct = sum(Nchild);

% Generate an array of selection indices and reweight
parent = zeros(Nct,1);
interm_weight = zeros(Nct,1);
jj=1;
for ii = 1:Np
    w = weight(ii)/Nchild(ii);
    parent(jj:jj+Nchild(ii)-1) = ii;
    interm_weight(jj:jj+Nchild(ii)-1) = w;
    jj = jj+Nchild(ii);
end

% Get rid of zero weights
thresh = max(interm_weight) / exp(30);
parent(interm_weight<thresh)=[];
interm_weight(interm_weight<thresh)=[];
Nct = length(parent);

% % Get rid of the least likely
% order = [cumsum(ones(Nct,1)) interm_weight];
% order = sortrows(order, -2);
% to_keep = order(1:Np, 1);
% samp = parent(to_keep);
% new_weight = interm_weight(to_keep);

% Systematically downsample
u = ceil( Nct * cumsum( (1/Np)*ones(Np,1) ) - unifrnd(0, 1/Np) );
samp = parent(u);
new_weight = interm_weight(u);

new_weight = new_weight/sum(new_weight);
new_weight = log(new_weight);

% Construct new, resampled array
ResamDistn.particles = cell(Par.NumPart,1);
for ii = 1:Par.NumPart
    ResamDistn.particles{ii} = Distn.particles{samp(ii)};
    ResamDistn.weight(ii) = new_weight(ii);
end

% "Hard copy" the new distribution (we just copied the pointers in the loop, for speed)
ResamDistn = ResamDistn.Copy;

% ResamDistn.weight = log(ones(Par.NumPart, 1)/Par.NumPart);

end