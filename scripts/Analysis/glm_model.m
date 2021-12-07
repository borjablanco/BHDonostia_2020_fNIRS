% Create TDR regressors (based on hmrDeconv_DriftSS function)

% TO DO: Clean the script
% TO DO: Create figures at different steps to understand the process

function [tdr_fw, tdr_bw, lstInc, Am] = glm_model(data)

% Initialize variables for deconvolution
s = data.s(:, [2,6]);
paramsBasis = [0 6 11]; % for gamma (tau sigma T)
trange = [-5 25];
t = data.time;
nT = length(t);
dt = t(2)-t(1);
nPre = round(trange(1)/dt);
nPost = round(trange(2)/dt);
nTpts = size(data.OD,1);
tHRF = (1*nPre*dt:dt:nPost*dt)';
ntHRF=length(tHRF);    

% According to Homer2-hmrDeconvHRF_DriftSS.m (0==noise and 1==clean)
% Have to modify for the analysis then
% Change to motion = 0 and clean = 1
tInc = data.tInc_auto;
tInc(tInc==0) = 2;
tInc(tInc==1) = 0;
tInc(tInc==2) = 1;
tIncAuto = tInc;

%figure; 
%plot(data.OD); hold on;
%plot(data.tInc_auto)

% Find stim onsets for each condition
lstCond = find(sum(s>0,1)>0);
nCond = length(lstCond); %size(s,2);
onset=zeros(nT,nCond);
nTrials = zeros(nCond,1);
for iCond = 1:nCond
    lstT=find(s(:,lstCond(iCond))==1);
    lstp=find((lstT+nPre)>=1 & (lstT+nPost)<=nTpts);
    lst=lstT(lstp);
    nTrials(iCond)=length(lst);
    onset(lst+nPre,iCond)=1;
end

% Create Modified Gamma function
tau = paramsBasis(1);
sigma = paramsBasis(2);
T = paramsBasis(3);
tbasis = (exp(1)*(tHRF-tau).^2/sigma^2) .* exp( -(tHRF-tau).^2/sigma^2 );

% Make zero baseline values
lstNeg = find(tHRF<0);
tbasis(lstNeg,1) = 0;

% I think this part is not doing anything
% is the same as previous two lines
if tHRF(1)<tau
    tbasis(1:round((tau-tHRF(1))/dt),1) = 0;
end

% Convolve basis function and task
% But it is not a box, just a line
foo = conv(tbasis,ones(round(T/dt),1)) / round(T/dt);
tbasis = foo(1:ntHRF,1);

% Create regressors for each condition
dA=zeros(nT,nCond);

for iCond=1:nCond
    
    clmn = conv(onset(:,iCond), tbasis);
    clmn=clmn(1:nT);
    dA(:,iCond)=clmn;
    
end

% Censoring for motion periods
% and baseline regressors
idxMA = find(diff(tIncAuto)==1);  % number of motion artifacts

nMA = length(idxMA);
nMC = nMA+1;
Amotion = zeros(nT,nMC);

% BUG in original script - if the first segment is marked as noise the first regressor
% becames all zeros making the matrix poor conditioned
% I had to add some lines of code for cases where the beginning of
% the experiment has been marked as motion

if tIncAuto(1) == 0
    Amotion(1:idxMA(1),1) = 1;
    for ii=2:nMA
        Amotion((idxMA(ii-1)+1):idxMA(ii),ii) = 1;
    end
    Amotion(:,2) = Amotion(:,1) + Amotion(:,2);
    Amotion(:,1) = [];
else
    Amotion(1:idxMA(1),1) = 1;
    for ii=2:nMA
        Amotion((idxMA(ii-1)+1):idxMA(ii),ii) = 1;
    end
end
Amotion((idxMA(nMA)+1):end,end) = 1;
lstInc = find(tIncAuto==1);

% Final design matrix
A = [dA(:,:) Amotion]; % Figure

% ----------------------------------
% GLM solution HOMER2

% Solve GLM
% ATA=A(lstInc,:)'*A(lstInc,:);
% pinvA=ATA\At(lstInc,:)'; % design matrix
% ytmp = y(lstInc,:); % data
% foo = pinvA*squeeze(ytmp); % betas

% % lstML = all channels
% % conc, can be removed because we are doing HbO and HbR separately


% for iCond=1:nCond
%     tb(:,lstML,conc,iCond)=foo([1:nB]+(iCond-1)*nB,lstML,conc);
%     %                yavg(:,lstML,conc,lstCond(iCond))=tbasis*tb(:,lstML,conc,lstCond(iCond));
%     if size(tbasis,3)==1
%         yavg(:,lstML,conc,iCond)=tbasis*tb(:,lstML,conc,iCond);
%     else
%         yavg(:,lstML,conc,iCond)=tbasis(:,:,conc)*tb(:,lstML,conc,iCond);
%     end
% end
% 
% % reconstruct y and yresid (y is obtained just from the HRF)
% % and R
% yresid(lstInc,conc,lstML) = ytmp - permute(At(lstInc,:)*foo(:,lstML,conc),[1 3 2]);
% ynew(lstInc,conc,lstML) = permute(dA(lstInc,:,conc)*foo(1:(nB*nCond),lstML,conc),[1 3 2]) + yresid(lstInc,conc,lstML);
% 
% yfit = permute(At(lstInc,:)*foo(:,lstML,conc),[1 3 2]);
% for iML=1:length(lstML)
%     yRtmp = corrcoef(ytmp(:,1,iML),yfit(:,1,iML));
%     yR(lstML(iML),conc) = yRtmp(1,2);
% end
% 
% %get error
% if glmSolveMethod==1 %  OLS  ~flagUseTed
%     pAinvAinvD = diag(pinvA*pinvA');
%     yest(:,lstML,conc) = At * foo(:,lstML,conc);
%     yvar(1,lstML,conc) = std(squeeze(y(:,conc,lstML))-yest(:,lstML,conc),[],1).^2; % check this against eq(53) in Ye2009
%     for iCh = 1:length(lstML)
%         bvar(:,lstML(iCh),conc) = yvar(1,lstML(iCh),conc) * pAinvAinvD;
%         for iCond=1:nCond
%             %                yavgstd(:,lstML(iCh),conc,iCond) = diag(tbasis*diag(bvar([1:nB]+(iCond-1)*nB,lstML(iCh),conc))*tbasis').^0.5;
%             if size(tbasis,3)==1
%                 yavgstd(:,lstML(iCh),conc,iCond) = diag(tbasis*diag(bvar([1:nB]+(iCond-1)*nB,lstML(iCh),conc))*tbasis').^0.5;
%             else
%                 yavgstd(:,lstML(iCh),conc,iCond) = diag(tbasis(:,:,conc)*diag(bvar([1:nB]+(iCond-1)*nB,lstML(iCh),conc))*tbasis(:,:,conc)').^0.5;
%             end
%             ysum2(:,lstML(iCh),conc,iCond) = yavgstd(:,lstML(iCh),conc,iCond).^2 + nTrials(iCond)*yavg(:,lstML(iCh),conc,iCond).^2;
%         end
%     end
%     
%     

% ------------------------------

% Design matrix after censoring
Am = A(lstInc,:);

% Plot design matrices
%figure
%subplot(121); imagesc(A); title('Design matrix before censoring')
%subplot(122); imagesc(Am); title('Design matrix after censoring')
%colormap gray

% Task regressors for each condition
tdr_fw = Am(:,1);
tdr_bw = Am(:,2);


end
