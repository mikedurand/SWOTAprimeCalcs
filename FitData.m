function [p,xbreak,iout,hval,wval,UseReg] = FitData(ReDoFits,nR,Wobs,stdW,Hobs,stdH,nReg)

if ReDoFits   
    disp('Re-computing fits. This takes a while...')    
    for i=1:nR  %fit H-W functions
        disp(['Computing fit for reach ' num2str(i) '/' num2str(nR)])
        SpecialCases{i} = CheckSpecialCases(Wobs(i,:),stdW,Hobs(i,:),stdH);
        if SpecialCases{i}.LowWidthVar
            RectangularFit;
        elseif SpecialCases{i}.LowHeightVar
            SimpleLinearFit;
        else
            [p{i},xbreak{i},iout{i},hval{i},wval{i}]=CalcWHFitsEIV(Hobs(i,:),Wobs(i,:),nReg,stdH,stdW);       
            UseReg{i}=[1 1 1];
        end
    end
    save('WHFits.mat','p','xbreak','iout','hval','wval','UseReg')
else
    load WHFits.mat
end