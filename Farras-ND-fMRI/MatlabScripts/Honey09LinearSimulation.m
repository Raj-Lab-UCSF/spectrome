function [Err,meanGraph] = Honey09LinearSimulation(MapC,corr_mtx,inds,simRange,Timepts)

%
% MapC: structural connectivity matrix
% Use the sparsified MapC matrix. Maybe also normalize it row-wise...
% corr_mtx: data-obtained functional connectivity matrix
% corr_conn: Simulated functional connectivity matrix
% inds: index of nonzero entries 
% cite paper....

NumROI = size(MapC,1);

if nargin == 4
    Timepts = 170; % Original number of time points 180, then first 10 time points removed.
end
delta_t = 1;

Err = zeros(10*2*simRange+1,10*2*simRange+1);
minErr = 1e6;
dlambda = 1.0;

% Noise stats N(0,1)
m_noise = 0;
s_noise = 1;

for alpha = -10*simRange:10*simRange
    for beta = 0*simRange+1:20*simRange+1
        u = zeros(NumROI,Timepts);
        if ~(alpha == 10 && beta == 0)
            A = (1 - (alpha/10)*delta_t)*eye(size(MapC)) + (beta/10)*MapC*delta_t;
            A = A/norm(A) * dlambda;
            for ii=1:Timepts-1
                noise = s_noise.*randn(size(A,1),1)*delta_t + m_noise;
                u(:,ii+1) = A*u(:,ii) + noise;
            end
            corr_conn = corrcoef(u');
            corr_conn = corr_conn/max(corr_conn(:));
            corr_mtx = corr_mtx/max(corr_mtx(:)); % DOUBLE CHECK on corr_mtx
            Err(alpha+10*simRange+1,beta) = norm(corr_conn(inds) - corr_mtx(inds),1);
            
            %crossCor = corrcoef(corr_conn(:),corr_mtx(:));
            %Err(alpha+10*simRange+1,beta+10*simRange+1) = crossCor(1,2);
            if  Err(alpha+10*simRange+1,beta) < minErr % replace beta+10*simRange+1 with beta
                minErr = Err(alpha+10*simRange+1,beta);
                minErrPosition = [alpha+10*simRange+1 beta];
                prm.alpha = (minErrPosition(1) - 10*simRange-1)/10;
                prm.beta = (minErrPosition(2))/10;
                graphMatrix = corr_conn;
            end
        end
    end
end

%meanGraph = mean(graphMatrix,3);
meanGraph = graphMatrix;

save meanGraphDataNew.mat Err graphMatrix minErrPosition minErr prm;
