function main(subjects,Fldr)

% Simulate and generate figures in the mansucript titled 
% "Network diffusion accurately models the relationship between structural and 
% functional brain connectivity networks," 
% F. Abdelnour and A. Raj, submitted to PNAS. September 2012.
%
% subjects is a struct of folders, each containing one subject. Each folder
% contains subfolders DTI, fMRI, and T1scan.
% Fldr is the directory containing the subjects.
% 
% The program assumes that the structural connectivity matrices have been
% obtained.
% Atlas T1 file assumed to be at folder
% /"subject"/T1scan/Atlased/"atlasFileName"
% Structural connectivity matrix is assumed to be at
% ../"subjects"/FiberTracts116/Matrix_ACP116.mat
%
% Note: The files from which the functional connectivity matrices are to be
% estimated are assumed to have been preprocessed using DPARSF:
% "DPARSF: a MATLAB toolbox for "pipeline" data analysis of resting-state
% fMRI,"Yan Chao-Gan and Zang Yu-Feng, Frontiers in Systems Neuroscience.
%
% Figs: Folder where generated figures of the paper are saved. All figures
% are saved in *.fig and *.eps formats.
% matFiles: Folder where various *.mat files are saved:
% 
% -subjectInfo: Located in each subject's folder. Contains each subject's estimated matrices and other
% information.
% 
% -CorPval: correlation and pvals from all three models relative to empirical functional networks
% for all methods and subjets.
% -meanLinErrOverAllSubjects: Overall mean error due to each subject,
% linear model.
% -nonLinearCorr: Mean error over all subjects when nonlinear model is
% used.
% -normTble: Correlation table 8x4 for all subjects and methods (Table 1 in paper):
% -CorrCoefIntra: Struct containing correlation coefficients and
% corresponding p-values for the case of mean matrices with no cross-hemispohere
% connectivity.
% -CorrCoefFull: Struct containing correlation coefficients and
% corresponding p-values for the case of mean matrices with cross-hemisphere connectivity.
% -CorrCoefIntraHemi: Table of Pearson correlations of individual subjects
% for intra hemisphere case.
% -CorCoefMedTcrit: Table 3 in the manuscript
%
% Functions/packages called:
%
% SPM5
% CONN
% myprocess_fMRI
% CalcPaths_Bayesian_generic
% brain_network_plotting2
% plotScatter2
% GummiBrain
% lob_info,
% Honey09NonLinearSim
% Honey09LinearSimulation
% AshishLinearModel,
% LaplacianMtx
% NMM_sim_function

% Number of nodes in each network:
numROI = 90;

% Folder where figures are saved...
Figs = [Fldr filesep 'Figures'];
if isdir(Figs) == 0
   mkdir([Fldr filesep 'Figures']);
end

% Folder where various *.mat files are saved...
matFiles = [Fldr filesep 'matFiles'];
if isdir(matFiles) == 0
   mkdir([Fldr filesep 'matFiles']);
end

% For the linear simulation, we run iterations for each subject.    
% numIt is the number of random instances for each subject. The final 
% matrix is obtained by taking the average of all the numIt matrices.  
numIt = 5;

graphstr = {'MapC','fC','Lin','NonLin','H'};
atFileName = {'fs009a1001_Atlas'}; % The altas for GummiBrain figures. Needs to be supplied by user.

savefig_flg = 1;
savemovie_flg = 0;

% Parameter c values where the nonlinear model is evaluated:
c = [0.02:0.05:0.32];

simRange = 3; % Range over which to simulate the linear model
Err = zeros(10*2*simRange+1,10*2*simRange+1,numIt);
% graphMatrix = zeros(numROI,numROI,numIt);
% minErrPosL1 = zeros(1,2,numIt);
% minErr = zeros(1,numIt);

%% Find structural connectivity:

cd(Fldr);
for i = 1:length(subjects)
    PatDir = [Fldr filesep subjects{i}];
    % Need to supply T1 file name below...
    T1FileName = [Fldr filesep subjects{i} filesep 'T1scan' filesep "filename"];
    atlassize = 116;
    diff_dir = [PatDir filesep 'DTI'];
    diff_ext = '*.img';
    % Modify lines 9, 10, 28...
    [time_for_streamlines, time_for_findfibre, time_for_conn]=CalcPaths_Bayesian_generic(PatDir,T1FileName,atlassize,diff_dir,diff_ext);
end

%% Run SPM fMRI Connectivity Toolbox processing
cd(Fldr);
for i = 1:length(subjects)
    myprocess_fMRI(['..' filesep subjects{i}]);
    cd ..
end

% Calculate functional correlation using SPM connectivity toolbox
%func_graphs = cell(size(subjects1));
for i = 1:length(subjects)
    cd(subjects{i})
    F = load([Fldr filesep subjects{i} '/conn_process/results/firstlevel/ANALYSIS_01/resultsROI_Subject001_Condition001.mat']);
    fg = F.Z(1:116,1:116); 
    fg(isnan(fg)) = 0;
    func_graphs = fg;
    save func_graphs.mat func_graphs;
    cd ..
end

% Calculate functional correlation using simple correlation
%corr_graphs = cell(size(subjects));
for i = 1:length(subjects)
    cd(subjects{i});
    F = load([Fldr filesep subjects{i} '/conn_process/results/preprocessing/ROI_Subject001_Condition001.mat']);
    RT = zeros(116,length(F.data{4}));
    for j = 1:116 
        RT(j,:) = F.data{j+3};
    end
    corr_graphs = corrcoef(RT');
    save corr_graphs.mat corr_graphs;
    cd ..
end

%% The following loop generates the Pearson correlation information for all
% subjects and three models: linear, nonlinear, and Laplacian.

for jj = 1:length(subjects)
    disp(['Loading and sparsifying structural connectivity matrices of subject ', num2str(subjects{jj})])
    % Open structural connectivity matrix...
    cd([Fldr filesep subjects{jj}]);
    ACPmatrix = [Fldr filesep subjects{jj} filesep 'FiberTracts116/Matrix_ACP116.mat'];
    load(ACPmatrix); % Contains matrix MapACP
    
    % Obtain structural map,
    % Sparsify and nomalize Mapc
    nzs = nonzeros(MapACP);
    % Set elements less than $0.1*\sigma$ to zero
    L = (MapACP > 0.1*std(nzs)); 
    MapC = L .* MapACP;
    % Mean of nonzeros after thresholding
    MapC = MapC(1:numROI,1:numROI);
    NormMean = mean(nonzeros(MapC));
    MapC = MapC/NormMean; % Normalize ACP matrix
    clear NormMean;
    
    % Normalize each MapC row to 1
    Norm = sum(MapC,2);
    Norm = repmat(Norm,1,size(MapC,2));
    MapC = MapC ./ Norm;
    clear Norm;
    
    Mtx.MapC = MapC;
    % Obtain functional map
    load([Fldr filesep subjects{jj} filesep 'corr_graphs.mat']);
    corr_graphs = corr_graphs(1:numROI,1:numROI);
    inds = find(abs(corr_graphs(:)) > 0.05*max(abs(corr_graphs(:))));
    
    %% We next find the Pearson correlation for each model, each subject, in addition to
    % the Pearson corrleation of the structure and empirical functional
    % connectivity.
    
    %% Structure vs true function correlation
    [corCoef,Pval.Struct{jj}] = corrcoef(Mtx.MapC(inds),corr_graphs(inds)); % per subject
    Cor.Struct(jj) = corCoef(1,2);                                          % per subject
    
    %% Galan/Honey linear model    
 
    disp(['Honey linear model, subject ' num2str(subjects{jj})])
    meanGraph = zeros(numROI,numROI,numIt);
    % Run Honey/Galan linear model...
    for ii = 1:numIt
        [Err(:,:,ii),meanGraph(:,:,ii)] = Honey09LinearSimulation(Mtx.MapC,corr_graphs,inds,simRange);
    end
    % Mean estimation error
    Error.meanErr = mean(Err,3);
    Mtx.meanGraph = mean(meanGraph,3);
    Error.lin = norm(Mtx.meanGraph(inds) - corr_graphs(inds));
    [corCoef,Pval.LinHoney{jj}] = corrcoef(Mtx.meanGraph(inds),corr_graphs(inds));
    Cor.LinHoney(jj) = corCoef(1,2);
    disp(['End of linear model simulation... '])
    
    
    %% Laplacian linear model
    % Resulting matrix H is saved in the current folder as LaplaceMtx.mat
    
    disp(['Ashish linear model, subject ' num2str(subjects{jj})])
    % Run Ashish's normalized Laplacian matrix model
    [Mtx.H,ErrH{jj},pval] = AshishLinearModel(Mtx.MapC,corr_graphs,inds);
    [corCoef,Pval.Ashish{jj}] = corrcoef(Mtx.H(inds),corr_graphs(inds));
    Cor.Ashish(jj) = corCoef(1,2);
    
    disp(['End of Laplacian linear model simulation...'])
    errLap = ErrH{jj};
   
     %% Honey nonlinear model
     % The nonlinear model is run over a range of values of c for subject jj given 
     % the corresponding structural connectivity matrix MapC. Then the matrix giving
     % the highest Pearson correlation is chosen.
     disp(['Honey nonlinear model, subject ' num2str(subjects{jj})])
     [~,maxIdx] = Honey09NonLinearSim(MapC,c,corr_graphs,inds,jj);
     rn = ['HoneyMatrix_C' num2str(c(maxIdx)*100) '.mat'];
     load('rn','C');
     Mtx.Cnl = C;
     clear C;
     [corCoef,Pval.Honey{jj}] = corrcoef(Mtx.Cnl(inds),corr_graphs(inds));
     Cor.Honey(jj) = corCoef(1,2);
     
     %% Now save all info in subject's folder:
     save subjectInfo Mtx Error maxIdx errLap % <--- make it a struct!!
end
cd(matFiles);
save CorPval Cor Pval; 

% End of outer loop. Find mean network for struct, data, and the estimated functional connectivity networks:

mMatrix.MapC = zeros(numROI,numROI);
mMatrix.fC = zeros(numROI,numROI);
mMatrix.Lin = zeros(numROI,numROI);
mMatrix.NonLin = zeros(numROI,numROI);
mMatrix.H = zeros(numROI,numROI);

cd(Fldr);
for ii=1:length(subjects)
    load([Fldr filesep subjects{jj} filesep 'corr_graphs.mat']);
    corr_graphs = corr_graphs(1:numROI,1:numROI);
    load('subjectInfo','Mtx');
    mMatrix.fC = mMatrix.fC + corr_graphs;
    mMatrix.MapC = mMatrix.MapC + Mtx.MapC;
    mMatrix.Lin = mMatrix.Lin + Mtx.meanGraph;
    mMatrix.H = mMatrix.H + Mtx.H;
    mMatrix.NonLin = mMatrix.NonLin + Mtx.Cnl;
end

mMatrix.MapC = mMatrix.MapC / length(subjects);
mMatrix.fC = mMatrix.fC / length(subjects);
mMatrix.Lin = mMatrix.Lin / length(subjects);
mMatrix.NonLin = mMatrix.NonLin / length(subjects);
mMatrix.H = mMatrix.H / length(subjects);

%% Generate Figure 1, mean over all subjects' linear method recon error
% matrix
cd(Fldr);
temp = zeros(numROI,numROI);
for ii=1:length(subjects)
    cd([Fldr filesep subjects{ii}]);
    load('subjectInfo','Error');
    temp = temp + Error.lin;
end
cd(Fldr);
temp = temp/length(subjects);
meanLinErrOverAllSubjects = temp;
clear temp;
% Mean L2 error matrix over all subjects resulting from the linear model.
cd(matFiles);
save meanLinErrOverAllSubjects.mat meanLinErrOverAllSubjects;

% Generate Fig 1 in the manuscript, the error over all subjects for the linear model. 
clf;
cd(Figs);
figure(1); imagesc(meanLinErrOverAllSubjects); axis tight square;
saveas(gcf,'meanLinErrOverAllSubjects.fig');
orient portrait;
print('-depsc','meanLinErrOverAllSubjects');

%% Generate Figure 2A (Pearson correlation for all subjects over all values of c in nonlinear approach):

cd(Fldr);
nonLinearCorr = zeros(length(subjects),length(c));
for jj=1:length(subjects)
    for kk=1:length(c)
        cd([Fldr filesep subjects{jj}]);
        load('corInfo');
        nonLinearCorr(jj,:) = corInfo.corr;
    end
end
cd(matFiles);
% Save matrix containing Pearson correlation values for the nonlinear model, for all subjects and range of paramter c.  
save nonLinearCorr.mat nonLinearCorr;
clf;
cd(Figs);
figure(1); 
plot(c,nonLinearCorr);
title('Nonlinear model Pearson coefficients');
xlabel('c parameter','FontSize',10);
ylabel('Pearson correlation','FontSize',10);
% Generate figure:
figure(1); 
saveas(gcf,'nonLinearCorr.fig');
orient portrait;
print('-depsc','nonLinearCorr');

%% Generate Figure 2C, the Pearson coefficient as a function of the product \beta*t
clf;
cd(Figs);
tvec = linspace(0,10, 100);
figure(1);
for ii=1:length(subjects)
    plot(tvec,ErrH{ii});
    hold on;
end
hold off
saveas(gcf,'laplacianCorrelationCurves.fig');
orient portrait
print('-depsc','laplacianCorrelationCurves');


%% Generate network figures of Figure 3 as well as the network of Figure 5
clf;
cd(Figs);

L = (abs(mMatrix.MapC) < max(max(-diag(diag(mMatrix.MapC)) + mMatrix.MapC))*0.5 == 0);
meanMapC2 = mMatrix.MapC .* L;
NormMean = mean(nonzeros(meanMapC2));
meanMapC2 = meanMapC2/NormMean;

tt = (numROI - 1) / max(mMatrix.MapC(:));
CColorMat = mMatrix.MapC * tt;
CColorMat = round(CColorMat)+1; 

L = (abs(mMatrix.fC) < max(max(-diag(diag(mMatrix.fC)) + mMatrix.fC))*0.5 == 0);
meanfC2 = mMatrix.fC .* L;
NormMean = mean(nonzeros(meanfC2));
meanfC2 = meanfC2/NormMean; 

L = (abs(mMatrix.Lin) < max(max(-diag(diag(mMatrix.Lin)) + mMatrix.Lin))*0.5 == 0);
meanLin2 = mMatrix.Lin .* L;
NormMean = mean(nonzeros(meanLin2));
meanLin2 = meanLin2/NormMean; 

L = (abs(mMatrix.NonLin) < max(max(-diag(diag(mMatrix.NonLin)) + mMatrix.NonLin))*0.5 == 0);
meanNonLin2 = mMatrix.NonLin .* L;
NormMean = mean(nonzeros(meanNonLin2));
meanNonLin2 = meanNonLin2/NormMean; 

L = (abs(mMatrix.H) < max(max(-diag(diag(mMatrix.H)) + mMatrix.H))*0.5 == 0);
meanH2 = mMatrix.H .* L;
NormMean = mean(nonzeros(meanH2));
meanH2 = meanH2/NormMean; 

graph = {meanMapC2,meanfC2,meanLin2,meanNonLin2,meanH2};

atlasFileName = [Fldr filesep subjects{1} filesep 'T1scan' filesep 'Atlased' filesep atFileName{1} '.hdr']; % atlasFileName = '';
hh = spm_vol(atlasFileName);
[at,~] = spm_read_vols(hh);    

ticID = tic;
for ii=1:length(graph)
    graphSparse = graph{ii}; 
    figstr = ['meanNetwork2_' graphstr{ii}];
    brain_network_plotting2(at, graphSparse, ones(size(graphSparse,1),1), lobes, figstr,CColorMat,savefig_flg,savemovie_flg);
end
tElapsed = toc(ticID)

%% Figure 4: Scatter plots
clf;
cd(Figs);

plotScatter2(mMatrix,'scatterPlotCorr',mMatrix.fC);

%% Generate color maps given activity seed for four different regions.
% Figures 5, 6, 7, and 8
clf;
cd(Figs);
precuneus{1} = mMatrix.MapC(:,67);

FrontalSupR{1} = mMatrix.MapC(:,3); % view(30,30);
FrontalSupR{2} = mMatrix.fC(:,3);
FrontalSupR{3} = mMatrix.Lin(:,3);
FrontalSupR{4} = mMatrix.NonLin(:,3);
FrontalSupR{5} = mMatrix.H(:,3);

HippocampusR{1} = mMatrix.MapC(:,37); % view(180,0); lightangle(120,30);
HippocampusR{2} = mMatrix.fC(:,37);
HippocampusR{3} = mMatrix.Lin(:,37);
HippocampusR{4} = mMatrix.NonLin(:,37);
HippocampusR{5} = mMatrix.H(:,37);

CingulumPostR{1} = mMatrix.MapC(:,35); % view(0,30);
CingulumPostR{2} = mMatrix.fC(:,35);
CingulumPostR{3} = mMatrix.Lin(:,35);
CingulumPostR{4} = mMatrix.NonLin(:,35);
CingulumPostR{5} = mMatrix.H(:,35);

PlotHemi = ones(numROI,1);
PlotHemi(2:2:end) = 0;


% Precuneus
clf;
p = precuneus{1};
p = p - mean(p);
stdp = std(p);
p = max(p , -1.5*stdp);
p = min(p , 1.5*stdp);
p1 = (p - min(p))*(numROI - 1)/(max(p) - min(p)) + 1;
p1 = round(p1);

disp(['Currently running ' graphstr{1}]);
GummiBrain(p1,atlasFileName,colormap(jet(numROI)),PlotHemi);
%view(180,0);
%lightangle(120,30);
colorbar;
saveas(gcf,['colormapPrecuneusR_RH' graphstr{1} '.fig']);
orient portrait;
print('-depsc',['colormapPrecuneusR_RH' graphstr{1}]);


% Hippocampus
clf;
for ii=2:length(graphstr)
    p = HippocampusR{ii};
    p = p - mean(p);
    stdp = std(p);
    p = max(p , -1.5*stdp);
    p = min(p , 1.5*stdp);
   
    p1 = (p - min(p))*(numROI - 1)/(max(p) - min(p)) + 1;
    p1 = round(p1); 
    
    disp(['Currently running ' graphstr{ii}]);
    GummiBrain(p1,atlasFileName,colormap(jet(numROI)),PlotHemi);
    view(180,0);
    lightangle(120,30);
    colorbar;
    saveas(gcf,['colormapHippocampusR_RH' graphstr{ii} '.fig']);
    orient portrait;
    print('-depsc',['colormapHippocampusR_RH' graphstr{ii}]);  
end    

% Frontal Sup R
clf;
for ii=2:length(graphstr)
    p = FrontalSupR{ii};
    p = p - mean(p);
    stdp = std(p);
    p = max(p , -1.5*stdp);
    p = min(p , 1.5*stdp);
    
    p1 = (p - min(p))*(numROI - 1)/(max(p) - min(p)) + 1;
    p1 = round(p1); 
    
    disp(['Currently running ' graphstr{ii}]);
    GummiBrain(p1,atlasFileName,colormap(jet(numROI)),PlotHemi);
    lightangle(30,30);
    colorbar;
    saveas(gcf,['colormapFrontalSupR_RH' graphstr{ii} '.fig']);
    orient portrait;
    print('-depsc',['colormapFrontalSupR_RH' graphstr{ii}]);  
end        
    
% Cingulum Post R
clf;
for ii=2:length(graphstr)
    p = CingulumPostR{ii};
    p = p - mean(p);
    stdp = std(p);
    p = max(p , -1.5*stdp);
    p = min(p , 1.5*stdp);
    
    p1 = (p - min(p))*(numROI - 1)/(max(p) - min(p)) + 1;
    p1 = round(p1); 
    
    disp(['Currently running ' graphstr{ii}]);
    GummiBrain(p1,atlasFileName,colormap(jet(numROI)),PlotHemi);
    view(0,30);
    %lightangle(120,30);
    colorbar;
    saveas(gcf,['colormapCingulumPostR_RH' graphstr{ii} '.fig']);
    orient portrait;
    print('-depsc',['colormapCingulumPostR_RH' graphstr{ii}]);  
end

%% Tables 1, 2: Compute the Pearson coefficient for mean connectivity matrices with and
% without cross-hemisphere connectivities (near end of paper):

% Without cross-hemisphere connectivities:
tmeanfC = mMatrix.fC; tmeanfC(1:numROI/2,(numROI/2 + 1):end) = 0; tmeanfC((numROI/2 + 1):end,1:numROI/2) = 0;
tmeanLin = mMatrix.Lin; tmeanLin(1:numROI/2,(numROI/2 + 1):end) = 0; tmeanLin((numROI/2 + 1):end,1:numROI/2) = 0;
tmeanNonLin = mMatrix.NonLin; tmeanNonLin(1:numROI/2,(numROI/2 + 1):end) = 0; tmeanNonLin((numROI/2 + 1):end,1:numROI/2) = 0;
tmeanH = mMatrix.H; tmeanH(1:numROI/2,(numROI/2 + 1):end) = 0; tmeanH((numROI/2 + 1):end,1:numROI/2) = 0; 
% Sparsify
nzIndex = find(abs(tmeanfC(:)) > 0.05*max(abs(tmeanfC(:))));
% Compute correlation coefficients:
[corCoef,rho] = corrcoef(tmeanMapC(nzIndex),tmeanfC(nzIndex));
CorrCoefIntra.MapC = corCoef(1,2);
CorrCoefIntra.MapCrho = rho;
[corCoef,rho] = corrcoef(tmeanLin(nzIndex),tmeanfC(nzIndex));
CorrCoefIntra.Lin = corCoef(1,2);
CorrCoefIntra.Linrho = rho;
[corCoef,rho] = corrcoef(tmeanNonLin(nzIndex),tmeanfC(nzIndex));
CorrCoefIntra.NonLin = corCoef(1,2);
CorrCoefIntra.NonLinrho = rho;
[corCoef,rho] = corrcoef(tmeanH(nzIndex),tmeanfC(nzIndex));
CorrCoefIntra.H = corCoef(1,2);
CorrCoefIntra.Hrho = rho;

% Including cross-hemisphere connectivities:
% Sparsify the full matrix
nzIndex2 = find(abs(mMatrix.fC(:)) > 0.05*max(abs(mMatrix.fC(:))));
%
[corCoef,rho] = corrcoef(mMatrix.MapC(nzIndex2),mMatrix.fC(nzIndex2));
CorrCoefFull.MapC = corCoef(1,2);
CorrCoefFull.MapCrho = rho;
[corCoef,rho] = corrcoef(mMatrix.Lin(nzIndex2),mMatrix.fC(nzIndex2));
CorrCoefFull.Lin = corCoef(1,2);
CorrCoefFull.Linrho = rho;
[corCoef,rho] = corrcoef(mMatrix.NonLin(nzIndex2),mMatrix.fC(nzIndex2));
CorrCoefFull.NonLin = corCoef(1,2);
CorrCoefFull.NonLinrho = rho;
[corCoef,rho] = corrcoef(mMatrix.H(nzIndex2),mMatrix.fC(nzIndex2));
CorrCoefFull.H = corCoef(1,2);
CorrCoefFull.Hrho = rho;

meanCorrIntra = [CorrCoefIntra.MapC CorrCoefIntra.Lin CorrCoefIntra.NonLin CorrCoefIntra.H];
meanCorrFull = [CorrCoefFull.MapC CorrCoefFull.Lin CorrCoefFull.NonLin CorrCoefFull.H];

%% Correlation table 8x4 for all subjects and methods (Table 1 in paper):
normTble = [Cor.Struct' Cor.LinHoney' Cor.Honey' Cor.Ashish'];
normTble = [normTble; mean(normTble,1)];
normTble = [normTble; meanCorrFull];
normTble = [(1:length(subjects))' normTble];
cd(matFiles);
% Save Table 1
save normTble normTble;

save CorrCoefIntra CorrCoefIntra;
save CorrCoefFull CorrCoefFull;

%% Table 2, Compute the Pearson coefficient for **individual subjects** connectivity matrices with and
% without cross-hemisphere connectivities (near end of paper):
% 09.04.2012

cd(Fldr);
tbl = zeros(length(subjects),5);
tbl(:,1) = (1:length(subjects))';
for ii=1:length(subjects)
   cd([Fldr filesep subjects{ii}]);
   
   load corr_graphs;
   load subjectInfo;
   
   fC = corr_graphs(1:numROI,1:numROI);
   % Set inter-hemisphere elements equal to zero...
   fC(1:numROI/2,(numROI/2 + 1):end) = 0; fC((numROI/2 + 1):end,1:numROI/2) = 0;
   tMapC = Mtx.MapC; tMapC(1:numROI/2,(numROI/2 + 1):end) = 0; tMapC((numROI/2 + 1):end,1:numROI/2) = 0;
   tmeanGraph = Mtx.meanGraph; tmeanGraph(1:numROI/2,(numROI/2 + 1):end) = 0; tmeanGraph((numROI/2 + 1):end,1:numROI/2) = 0;
   tH = Mtx.H; tH(1:numROI/2,(numROI/2 + 1):end) = 0; tH((numROI/2 + 1):end,1:numROI/2) = 0;
   tCnl = Mtx.Cnl; tCnl(1:numROI/2,(numROI/2 + 1):end) = 0; tCnl((numROI/2 + 1):end,1:numROI/2) = 0;
   % Sparsify
   nzIndex = find(abs(fC(:)) > 0.05*max(abs(fC(:))));
   %
   [corCoef,rho] = corrcoef(tMapC(nzIndex),fC(nzIndex));
   CorrCoefIntraHemi(ii).MapC = corCoef(1,2); tbl(ii,2) = corCoef(1,2);
   CorrCoefIntraHemi(ii).MapCrho = rho;
   
   [corCoef,rho] = corrcoef(tmeanGraph(nzIndex),fC(nzIndex));
   CorrCoefIntraHemi(ii).meanGraph = corCoef(1,2); tbl(ii,3) = corCoef(1,2);
   CorrCoefIntraHemi(ii).meanGraphrho = rho;
   
   [corCoef,rho] = corrcoef(tCnl(nzIndex),fC(nzIndex));
   CorrCoefIntraHemi(ii).Cnl = corCoef(1,2); tbl(ii,4) = corCoef(1,2);
   CorrCoefIntraHemi(ii).Cnlrho = rho;
   
   [corCoef,rho] = corrcoef(tH(nzIndex),fC(nzIndex));
   CorrCoefIntraHemi(ii).H = corCoef(1,2); tbl(ii,5) = corCoef(1,2);
   CorrCoefIntraHemi(ii).Hrho = rho;
   
   clear corr_graphs nzIndex fC tMapC tmeanGraph tH tCnl
end

tbl = [tbl; mean(tbl,1); meanCorrIntra];
CorrCoefIntraHemiTable = tbl;
CorrCoefIntraHemiTable = [(1:length(subjects))' CorrCoefIntraHemiTable];
cd(matFiles);
save CorrCoefIntraHemi CorrCoefIntraHemi;
save CorrCoefIntraHemiTable CorrCoefIntraHemiTable

%% Table 3: Find Pearson coefficient for fixed \beta t_{crit} for all subjects:
medTcrit = 2.0;
 %beta = 1;
 for ii=1:length(subjects)
    cd([Fldr filesep 'DPARSFWorkingDir' filesep 'Subjects' filesep subjects{ii}])
    load('subjectInfo','Mtx');
    load corr_graphs corr_graphs;
    corr_graphs = corr_graphs(1:numROI,1:numRO);
    Hg = LaplacianMtx(Mtx.MapC,medTcrit); % medTcrit
    [ErrHg,pvalg{ii}] = corrcoef(Hg(inds),corr_graphs(inds));
    corCoefMedTcrit(ii) = ErrHg(1,2); 
 end
corCoefMedTcrit = corCoefMedTcrit';

corCoefMedTcrit = [(1:length(subjects))' Cor.Ashish' corCoefMedTcrit];
cd(matFiles);
save corCoefTcrit corCoefMedTcrit pvalg;


