function plotScatter2(Mtx,filename,corr_graphs)
% This function will create scatter plots and correlation tables

lobe_info;
lobes = lobes(1:90);
repLobes = repmat(lobes,90,1);

n = 2;
m = 2;

%npanels = flags.HoneyLinear + flags.HoneyNonLinear + flags.AshishLinear + 1;
figure;
subplot(n, m, 1); scatter(Mtx.MapC(:), corr_graphs(:),[],repLobes,'x'); lsline;
    title(sprintf('SC-TrueFC')); axis tight square;
    xlabel('True fMRI FC','FontSize',8)
    ylabel('SC','FontSize',8)
kk=1;
%if flags.HoneyLinear == 1
    kk=kk+1;
    subplot(n, m, kk); 
     %
    test = (Mtx.meanGraph < 1 & Mtx.meanGraph >= 0.0);
    aTest = Mtx.meanGraph .* test;
    %
    scatter(aTest(:), corr_graphs(:),[],repLobes,'x'); lsline; 
    axis tight;
    %scatter(Mtx.meanGraph(inds), corr_graphs(inds),[],repLobes(inds),'x'); lsline;
    title(sprintf('Linear model')); axis tight square;
    xlabel('True fMRI FC','FontSize',9)
    ylabel('FC linear model','FontSize',9)
%end
%if flags.HoneyNonLinear == 1
    kk=kk+1;
    subplot(n, m, kk); 
    %
    test = (Mtx.Cnl < 1 & Mtx.Cnl > -0.2);
    aTest = Mtx.Cnl .* test;
    %
    scatter(aTest(:), corr_graphs(:),[],repLobes,'x'); lsline;
    title(sprintf('Nonlinear model')); axis tight square;
    xlabel('True fMRI FC','FontSize',9)
    ylabel('FC nonlinear model','FontSize',9) 
%end
%if flags.AshishLinear == 1
    kk=kk+1;
    subplot(n, m, kk); 
    %
    test = (Mtx.H < 0.16);
    aTest = Mtx.H .* test;
    %
    scatter(aTest(:), corr_graphs(:),[],repLobes,'x'); lsline;
    %scatter(Mtx.H(inds), corr_graphs(inds),[],repLobes(inds),'x'); lsline;
    title(sprintf('Laplacian model')); axis tight square;
    xlabel('True fMRI FC','FontSize',9)
    ylabel('FC laplacian model','FontSize',9) 
%end

saveas(gcf,[filename '.fig']);
orient portrait
print('-depsc',filename);