function ESS = CalcESS( weight )
%CALCESS Calculates effective sample size from array of (logarithmic) weights
weight = exp(weight);
ESS = 1/( sum( weight.^2 ) );
end

