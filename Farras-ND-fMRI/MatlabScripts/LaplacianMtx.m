function H = LaplacianMtx(MapC,Tcrit)

%
%
%
%
%

degC = diag(sum(MapC,2)); % Connectivity strength matrix
beta = 1;
% Unnormalized Laplacian, aka "generator matrix":
A = degC - MapC;
% Normalized Laplacian matrix:
iDelta = diag(sqrt(1./(diag(degC + eps))));
% H = iDelta.^2 * H;
A = iDelta * A * iDelta;

H = expm(-beta*A*Tcrit);
%save LaplaceMtx H corCoef pval;