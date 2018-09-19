function sC = structFromFunc(meanfC,meanMapC,bT,thrsh,fig)

%
%
%
%
%
%

Cd = diag(sum(meanMapC,1)) - meanMapC; % "\Delta - C", Laplacian
%Cd = eye(size(meanMapC)) - diag(1./sum(meanMapC,1))*meanMapC; % "I - Delta^-1 C", Laplacian
D1 = diag(sqrt(1./sum(meanMapC,1)));
%Cd = eye(size(meanMapC)) - D1*meanMapC*D1; % "I - D^-1/2 C D^-1/2", Laplacian

%% Normalize
%Cd = (Cd - min(Cd(:))) / (max(Cd(:)) - min(Cd(:))); 
%% Threshold from "prior" knowledge
%mask = (Cd1 > 0.05);

fC2 = (meanfC + meanfC')/2; % Exact symmetry
%idx = (abs(fC2) > 0.10); % Some thresholding of functional connectivity
%fC2 = fC2 .* idx;

% Regularize fC2 via SVD
[U,S,V] = svd(fC2);
indx = (diag(S) > thrsh*max(diag(S)));
% Compute estimated structural Laplacian 
sL = U(:,indx) * diag(-log(diag(S(indx,indx)))/bT) * V(:,indx)';
sC = -sL;
for ii=1:size(sC,1)
sC(ii,ii) = 0;
end

% Normalize
%sC = (sC - min(sC(:))) / (max(sC(:)) - min(sC(:)));
ssort = sort(sC(:), 'Ascend');
maxsc = ssort(round(0.995*length(ssort)));
minsc = 0.1*maxsc;
sC = min(sC,maxsc);
sC = max(sC,minsc);
sC = (sC - minsc) / (maxsc - minsc);

if fig == 1
figure; 
subplot(1,3,1); 
imagesc(meanfC); 
axis tight off square
title(['thrsh = ' num2str(thrsh)]);
subplot(1,3,2); 
imagesc(meanMapC); 
axis tight off square
subplot(1,3,3); 
imagesc(sC); 
axis tight off square
end
 
