function [corMax,maxIdx] = Honey09NonLinearSim(MapC,c,corr_graphs,inds,kk)

% Honey nonlinear simulation
% Cite paper:
%
% -MapC is a sparsified structural connectivity matrix. 
% -c is a free parameter 
% -kk is the kkth subject
%
% Outputs: 
%   C: estimated functional connectivity matrix obtained via the
%   Honey et al. nonlinear model.
%   parameter: the free global connectivity parameter c
%   timeElapsed: The time that took the program to estimate the functional
%   connectivity matrix C.
%
% Pearson correlation for all values of c are saved in corInfo.mat in addition to the index of the highest correlation.

%addpath('/home/farras/Documents/MATLAB/myToolboxes/Sporns');

disp(['Honey nonlinear model, subject ' num2str(kk)])

for ii=1:length(c)
    corMax = -1e12;
    disp(['Honey nonlinear model, c = ' num2str(c(ii)) ', subject ' num2str(kk)])
    ticID = tic;
    rn = ['Honey_C' num2str(c(ii)*100)];
    C = NMM_sim_function(rn,MapC,c(ii));
    savefile = ['HoneyMatrix_C' num2str(c(ii)*100) '.mat'];
    parameter = c(ii);
    timeElapsed = toc(ticID);
    
    [corCoef,Pval.Honey{ii}] = corrcoef(C(inds),corr_graphs(inds));
    corInfo.corr(ii) = corCoef(1,2);
    if corCoef(1,2) >= corMax,
        corMax = corCoef(1,2);
        maxIdx = ii;    
    end
    save(savefile,'C','parameter','timeElapsed');
end

corInfo.corMax = corMax;
corInfo.maxIdx = maxIdx;
save('corInfo.mat','corInfo');
    
    
    
