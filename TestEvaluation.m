function [hval,wval]=TestEvaluation(x,xbreak,p,nReg)

hval=linspace(min(x),max(x),100);
wval=[];
for i=1:nReg
    iVal=hval>=xbreak(i) & hval<xbreak(i+1);
    wval=[wval polyval(p{i},hval(iVal))];
end