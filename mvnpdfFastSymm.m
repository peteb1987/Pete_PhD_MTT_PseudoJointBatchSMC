function [ pdf ] = mvnpdfFastSymm( X, Mu, Var )
%MVNPDFFASTSYMM This is an alternative to mvnpdf when the covariance matrix
%is a multiple of the identity matrix, i.e. iso-probability surfaces are
%hyperspheres. Var is entered as a scalar.

A = -(1/(2*Var))*sum((X-Mu).^2)-(length(X)/2)*log(2*pi*Var);

pdf = exp(A);

end