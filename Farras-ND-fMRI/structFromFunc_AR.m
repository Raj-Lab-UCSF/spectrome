function sC1 = structFromFunc_AR(fCall, sCall, thrsh)

meanfC = mean(fCall, 3);
meanMapC = mean(sCall, 3);

Cd = diag(sum(meanMapC,1)) - meanMapC; % "\Delta - C", Laplacian
%Cd = eye(size(meanMapC)) - diag(1./sum(meanMapC,1))*meanMapC; % "I - Delta^-1 C", Laplacian
D1 = diag(sqrt(1./sum(meanMapC,1)));
%Cd = eye(size(meanMapC)) - D1*meanMapC*D1; % "I - D^-1/2 C D^-1/2", Laplacian
nsubj = size(sCall,3);

for k = 1:nsubj
fCsubj = fCall(:,:,k);
sCsubj = sCall(:,:,k);

fC2 = (fCsubj + fCsubj')/2; % Exact symmetry
%idx = (abs(fC2) > 0.10); % Some thresholding of functional connectivity
%fC2 = fC2 .* idx;
 
% nOW USE svd OF sC
[U,S,V] = svd(Cd);
q = diag(U'*fC2*V);
minq = thrsh*max(abs(q));
q(q<minq) = 1;
qq = -log(q);
sL = U*diag(qq)*V';
sC = -sL;
for ii=1:size(sC,1)
sC(ii,ii) = 0;
end

% Normalize
ssort = sort(sC(:), 'Ascend');
maxsc = ssort(round(0.995*length(ssort)));
minsc = 0.1*maxsc;
sC = min(sC,maxsc);
sC = max(sC,minsc);
sC = (sC - minsc) / (maxsc - minsc);

%% Apply mask
%sC1 = sC1 .* mask;

sCall_est(:,:,k) = sC; 

figure; subplot(1,3,1); imagesc(fCsubj); title(['thrsh = ' num2str(thrsh)]);
subplot(1,3,2); imagesc(sCsubj); subplot(1,3,3); imagesc(sC); 

end
meansC_est = mean(sCall_est, 3);

figure; subplot(1,3,1); imagesc(meanfC); title(['thrsh = ' num2str(thrsh)]);
subplot(1,3,2); imagesc(meanMapC); subplot(1,3,3); imagesc(meansC_est); 
 
