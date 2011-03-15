function ResamDistn = SystematicResample( Distn, weight )
%SYSTEMATICRESAMPLE Resample Distn systematically from weight

global Par;

% Change the weights to linear
weight = exp(weight);

% Copy distribution
ResamDistn = Distn.Copy;

% Create empty vector for offspring count
N = zeros(Par.NumPart, 1);

% Generate random index array
u = zeros(Par.NumPart, 1);
u(1) = rand/Par.NumPart;
for ii = 2:Par.NumPart
    u(ii) = u(1) + (ii-1)/Par.NumPart;
end    

% Generate cumulative weight array
w_sum = cumsum(weight);

% Enumerate offspring
for ii = 1:Par.NumPart
    if ii>1
        N(ii) = sum((u < w_sum(ii))&(u > w_sum(ii-1)));
    else
        N(ii) = sum(u < w_sum(ii));
    end
end

assert(sum(N)==Par.NumPart, 'Wrong number of resampled children');

% Construct new, resampled array
ResamDistn.particles = cell(Par.NumPart,1);
jj = 1;
for ii = 1:Par.NumPart
    for k = 1:N(ii)
        ResamDistn.particles{jj} = Distn.particles{ii};
        jj = jj + 1;
    end
end

% "Hard copy" the new distribution (we just copied the pointers in the loop, for speed)
ResamDistn = ResamDistn.Copy;

ResamDistn.weight = log(ones(Par.NumPart, 1)/Par.NumPart);

end