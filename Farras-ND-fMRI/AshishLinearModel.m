function [H,corCoef,pval] = AshishLinearModel(MapC,corr_graphs,inds)

% Estimate the (normalized) Laplacian matrix
% MapC: Structural connectivity matrix
% corr_graphs: Data-generated functional connectivity matrix
% H, Err, and pval are saved in the current directory

nt = 100;

beta = 1;
%corr_graphs1 = corr_graphs(:);
tvec = linspace(0,10, nt);
degC = diag(sum(MapC,2)); % Connectivity strength matrix

% Unnormalized Laplacian, aka "generator matrix":
A = degC - MapC;
% Normalized Laplacian matrix:
iDelta = diag(sqrt(1./(diag(degC + eps))));
% H = iDelta.^2 * H;
A = iDelta * A * iDelta;

for ii=1:nt
    Hhat = expm(-beta*A*tvec(ii));
    [Err,pval{ii}] = corrcoef(Hhat(inds),corr_graphs(inds));
    corCoef(ii) = Err(1,2);
end

erR = find(corCoef == max(corCoef));
H = expm(-beta*A*tvec(erR));
%save LaplaceMtx H corCoef pval;
