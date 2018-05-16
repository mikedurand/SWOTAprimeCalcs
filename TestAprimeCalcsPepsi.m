%% Computing cross-sectional area changes from height and width observations
% Mike Durand, durand.8@osu.edu
% 
% May 11, 2018

clear all; close all

addpath('./PepsiSrc');
%% 3. Setup: Read data and create synthetic observations
% So far, I've just looked at example cases on the Sacramento Upstream (cf=13) 
% and Garonne Downstream  (cf=4) from the Pepsi 1. 
%%
cf=4;
pathtoncfiles='./Pepsi1/';
Files=dir([pathtoncfiles '*.nc']);
numbfiles=size(Files,1);
ReadRiver;
FilterRiver;
R=Rivers(cf).Reaches; 
disp(['Working with: ' Rivers(cf).Name ' data.'])
tObs=ceil(Rivers(cf).orbit./86400)';
nObsCycle=length(tObs);
clear Rivers

[nR,nt]=size(R.W);
t=1:nt;

% SWOT observation
nCycle=floor(nt/21);
nReg=3; %number of regions
for i=2:nCycle
    tObs=[tObs tObs(1:nObsCycle)+(i-1)*21];
end
nObs=length(tObs);

%% 4. Visualize height-width relationship 
r=3;

disp(['Working with reach ' num2str(r)])

figure(1)
plot(R.H(r,tObs),R.W(r,tObs),'o','LineWidth',2); 
set(gca,'FontSize',14)
xlabel('Water surface elevation, m')
ylabel('River width, m')
title(['Fig. 1: "True" Height-width data for reach ' num2str(r)])
grid on


%% 5. Fitting funtions to the width & height data

ReDoFits=false;
if ReDoFits
    stdH=0.1;
    stdW=10;
    Hobs=R.H(:,tObs)+stdH.*randn(nR,nObs);
    Wobs=R.W(:,tObs)+stdW.*randn(nR,nObs);
    Hbar=median(Hobs,2);
    Wbar=median(Wobs,2);
    
    disp('Re-computing fits. This takes a while...')
    
    for i=1:nR  %fit H-W functions
        disp(['Computing fit for reach ' num2str(i) '/' num2str(nR)])
        [p{i},xbreak{i},iout{i},hval{i},wval{i}]=CalcWHFitsEIV(Hobs(i,:),Wobs(i,:),nReg,stdH,stdW);        
    end
    save('WHFits.mat','p','xbreak','iout','hval','wval','Hobs','Hbar','Wobs','Wbar','stdH','stdW')
else
    load WHFits.mat
end

m_zz=nan(2,2,nR);  %m_zz is the W-H observation covariance matrix for each reach
for i=1:nR    
    m_zz(:,:,i)=cov(Wobs(i,:),Hobs(i,:));
end

figure(2)
errorbar(Hobs(r,:),Wobs(r,:),stdW*ones(nObs,1),stdW*ones(nObs,1),stdH*ones(nObs,1),stdH*ones(nObs,1),'LineStyle','none','LineWidth',2,'Marker','o'); hold on
plot(hval{r},wval{r},'LineWidth',2)
a=axis;
for i=1:length(xbreak{r})
    plot(xbreak{r}(i)*ones(2,1),a(3:4),'g-')
end
plot(Hobs(r,iout{r}),Wobs(r,iout{r}),'kx','MarkerSize',12,'LineWidth',2)
hold off;
set(gca,'FontSize',14)
xlabel('Water surface elevation, m')
ylabel('River width, m')
title(['Fig. 2: Height-width fits & obs. for reach ' num2str(r)])
grid on

%% 6. Performing the$A \prime$calculations 
dAHbar=CalculatedAEIV(Hbar(r),Wbar(r),xbreak{r},p{r},nReg,0,stdW^2,stdH^2,m_zz(:,:,r),nObs); %sample calculation

Htest=linspace(xbreak{r}(1)+.001,xbreak{r}(4)-.001,100);

for i=1:length(Htest)
    [dAtest(i),Wtest(i)] = CalculatedAEIV(Htest(i),[],xbreak{r},p{r},nReg,dAHbar,stdW^2,stdH^2,m_zz(:,:,r),nObs); %sample calculation
end

[~,iSort]=sort(R.H(r,:));
[Hs,iUnique]=unique(R.H(r,:),'sorted');
As=R.A(r,iUnique);
AHbarhat_true=interp1(Hs,As,Hbar(r));

epsilonAbar=AHbarhat_true-median(R.A(r,:));

disp(['epsilon_Abar='  num2str(epsilonAbar) ' m^2'])
disp(['Abar='  num2str(median(R.A(r,:))) ' m^2'])

dAtrue=R.A(r,:)-median(R.A(r,:)); %true dA for comparison

figure(3)
plot(Reaches.H(r,tObs),dAtrue(tObs),'o',Htest,dAtest,Hbar(r),0,'ko','LineWidth',2)
set(gca,'FontSize',14)
xlabel('Height, m')
ylabel('A'', m^2')
legend('True','SWOT calibrated','Hbar,0','Location','Best')
grid on;
title(['Fig. 3: Illustration of A''(H) for reach ' num2str(r)])


%% 7. Timeseries estimation of $A\prime$ and preliminary uncertainty calculations

for i=1:nObs
    [dAhat(i),What(i),Hhat(i),dAunc(i)] = CalculatedAEIV(Hobs(r,i),Wobs(r,i),xbreak{r},p{r},nReg,dAHbar,stdW^2,stdH^2,m_zz(:,:,r),nObs);
end

figure(4)
plot(tObs,dAtrue(tObs),'s-','LineWidth',2); hold on;
errorbar(tObs,dAhat,dAunc,'s-','LineStyle','none','LineWidth',2,'MarkerSize',10); hold off;
set(gca,'FontSize',14)
xlabel('Obs Time')
ylabel('A'', m^2')
title(['Fig. 4: A'' timeseries for reach ' num2str(r)])
legend('True','SWOT')

Err.AprimeErrStd=std(dAhat-dAtrue(tObs));
Err.AprimeAvgUnc=mean(dAunc);
%% 8. Taking a quick look at the constrained widths and heights

figure(5)

plot([min(Wobs(r,~iout{r})) max(Wobs(r,~iout{r}))],[min(Wobs(r,~iout{r})) max(Wobs(r,~iout{r}))],'k-'); 
hold on;
plot(R.W(r,tObs(~iout{r})),Wobs(r,~iout{r}),'ro','LineWidth',2,'MarkerSize',10); 
plot(R.W(r,tObs(~iout{r})),What(~iout{r}),'bx','LineWidth',2,'MarkerSize',10); 
hold off
set(gca,'FontSize',14)
xlabel('True Width, m')
ylabel('SWOT Width, m')
legend('1:1','Observations','Constrained Obs','Location','Best')
title(['Fig. 5: Observed and constrained width for reach ' num2str(r)])
grid on;


Err.WidthObsStdDev=std(R.W(r,tObs(~iout{r}))-Wobs(r,~iout{r})); %close to stdW
Err.WidthConstrainedStdDev=std(R.W(r,tObs(~iout{r}))-What(~iout{r})); 
Err.HeightObsStdDev=std(R.H(r,tObs(~iout{r}))-Hobs(r,~iout{r})); %close to stdW
Err.HeightConstrainedStdDev=std(R.H(r,tObs(~iout{r}))-Hhat(~iout{r})); 


disp(Err)
