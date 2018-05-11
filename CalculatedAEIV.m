function [dAhat,What,Hhat,dAunc] = CalculatedAEIV(Hobs,Wobs,Hbp,p,nReg,dAHbar,sigma_ee,sigma_uu,m_ZZ,nObs)

%0) Integrate function: this will be removed as it's a trivial computation
for i=1:nReg
    pi{i}=polyint(p{i});
end

%1) Check whether out of boudns
OutOfBound=Hobs>max(Hbp) | Hobs<min(Hbp);

if OutOfBound    
    What=nan;
    dAunc=nan;
    
    if Hobs>max(Hbp)
        dAhat=polyval(pi{end},Hbp(end))-polyval(pi{end},Hbp(end-1))+dAHbar;
        What=Wobs;
        dAunc=sqrt(sigma_uu*Wobs^2+2*sigma_ee*(Hobs-Hbp(end))^2);
    else
        dAhat=-dAHbar-(Hbp(1)-Hobs)*(Wobs+p{1}(1)*Hbp(1)+p{1}(2))/2;
        What=Wobs;
        dAunc=sqrt(sigma_uu*Wobs^2+2*sigma_ee*(Hobs-Hbp(1))^2);
    end
    
    return
end

%2) Find region
if Hobs>=max(Hbp)
    iReg=nReg;
else
    iReg=find(Hobs>=Hbp, 1, 'last' ); %region of the curve that H falls into
end

%3) Compute height estimate Hhat using observed W
Hhat=EstimateH(Wobs,Hobs,p{iReg},sigma_ee,sigma_uu);

%4) Check region of updated Hhat. If it's changed, then update Hhat again
if Hhat>=max(Hbp)
    iRegHat=nReg;
else
    iRegHat=find(Hhat>=Hbp, 1, 'last' ); %region of the curve that H falls into
end

if iRegHat~=iReg
    iReg=iRegHat;
    Hhat=EstimateH(Wobs,Hobs,p{iReg},sigma_ee,sigma_uu);
end

%5) compute estimated W to go with estimated H. Note that this is 1.3.18 in Fuller
What=polyval(p{iReg},Hhat);

%6) estimate polynomial uncertainty. Shortcut: m_ZZ averaged over all of the sub-domains. 
% VarBeta=ComputeSlopeUncertainty(p{iReg}(1),sigma_uu,sigma_ee,m_ZZ,nObs);

%7) calculate dA: integrate the polynomial fits and evaluate at Hhat
for i=1:iReg
    dA(i)=polyval(pi{i},min(Hhat,Hbp(i+1)))-polyval(pi{i},Hbp(i));
end

dAhat=sum(dA)-dAHbar;

%8) estimate dA uncertainty: this needs to be improved in the future. Tried
%an EIV-based approach, but abandoned it. Implemented an OLS approach, but
%it is not as rigorous as it should be.

% EIV approach: this is unreliable, and needs a more thorough derivation.
% Leaving it in case hooks are helpful in future. MD 5/11/18
% dAunc=sqrt(VarBeta/4*(Hhat^2-Hbp(iReg)^2)^2 + p{iReg}(1)^2/4 * sigma_uu);

% 
dAunc=Wobs*sqrt(sigma_uu)*sqrt(2);

end

function Hhat=EstimateH(Wobs,Hobs,p,sigma_ee,sigma_uu)

if ~isempty(Wobs)
    %note this implements eqn. 1.3.17 in Fuller, assuming sigma_eu=0
    sigma_vv=sigma_ee+p(1)^2*sigma_uu;
    sigma_uv=-p(1)*sigma_uu;
    v=Wobs-p(2)-p(1)*Hobs;
    Hhat=Hobs-sigma_vv^-1*sigma_uv*v;
else
    Hhat=Hobs;
end

end

function lambda=RootsOfVarianceEqn(m_ZZ,sigma_ee,sigma_uu)

%comptue the roots lambda of equation 1.3.26 in Fuller

m_YY=m_ZZ(1,1); m_XX=m_ZZ(1,1); m_XY=m_ZZ(1,2);

a=sigma_ee*sigma_uu;
b=-(sigma_uu*m_YY+sigma_ee*m_XX);
c=m_YY*m_XX-m_XY^2;

lambda=roots([a b c]);

end

function VarBeta=ComputeSlopeUncertainty(betahat1,sigma_uu,sigma_ee,m_ZZ,nObs)

%This does not take into account that EIV objective function is minimized
%over regions... likely overestimates uncertainty.

lambda=RootsOfVarianceEqn(m_ZZ,sigma_ee,sigma_uu);

if max(lambda)>1
    %use eqn. 1.3.31
    sigma_uv=-betahat1*sigma_uu;
    sigma_vv=sigma_ee+betahat1^2*sigma_uu;
    H2prime=[0; 1;]-[1; -betahat1;]*sigma_vv^-1*sigma_uv;
    Sigma=[sigma_ee 0; 0 sigma_uu;];
    m_xx=H2prime'*(m_ZZ-Sigma)*H2prime;
    VarBeta=1/(nObs-1)*m_xx^-2*(m_xx*sigma_vv+sigma_uu*sigma_vv);
else
    %use OLS uncertainty. via Statistics by Johnson & Bhattacharyya §11.5
    Sxx=m_ZZ(2,2);
    VarBeta=sigma_ee/Sxx/(nObs-1);
end


end