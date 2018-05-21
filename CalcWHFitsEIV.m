function [p,xbreak,iout,hval,wval]=CalcWHFitsEIV(Hobs,Wobs,nReg,stdH,stdW)

%%
x=Hobs;
y=Wobs;
xbreak0=linspace(min(x),max(x)+.1,nReg+1);
p0hat=[];

% find outliers
iout=isoutlier(y) & isoutlier(x);
if any(iout) && all(max(x(iout))>=x) && all(max(y(iout))>=y)
    HighOut=true;
    xbreak0(nReg)=max(x(~iout));
    if xbreak0(2)>xbreak0(3)
       xbreak0(2)=mean(xbreak0([1 3]));
    end
else
    HighOut=false;
end

% calculate preliminary fits
for i=1:nReg
    iReg=x>=xbreak0(i) & x<xbreak0(i+1);
    if sum(iReg)>2
        p0{i}=polyfit(x(iReg),y(iReg),1);
    elseif i>1
        p0{i}=p0{1};
    else
        p0{i}=[100 min(y)];
    end
    p0hat=[p0hat p0{i}];
end
   
% fit analysis and unpacking
[xbreak,sse,phat]=WH_EIV(x,y,p0hat,stdH^2,stdW^2,xbreak0,nReg,HighOut);

for i=1:nReg
    p{i}=phat((i-1)*2+[1 2] );
end

%%
hval=linspace(min(x),max(x),100);
wval=[];
for i=1:nReg
    iVal=hval>=xbreak(i) & hval<xbreak(i+1);
    wval=[wval polyval(p{i},hval(iVal))];
end

end