function y = mvnpdfQ( X, Mu )
%MVNPDFQ A slim version of mvnpdf for the process covariance, for which we
% pre-calculate the cholesky decomposition and omit options/validity tests.

global Par;

R = Par.Qchol;
d = length(R);

X0 = X - Mu;
xRinv = X0 / R;
logSqrtDetSigma = sum(log(diag(R)));
quadform = sum(xRinv.^2, 2);
y = exp(-0.5*quadform - logSqrtDetSigma - d*log(2*pi)/2);

end