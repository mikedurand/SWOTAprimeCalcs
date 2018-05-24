p{i}{1}=[0 mean(Wobs(i,:))]; p{i}{2}=[]; p{i}{3}=[]; p{i}{4}=[];
xbreak{i}=[min(Hobs(i,:)) max(Hobs(i,:)+1E-4) nan nan];
iout{i}=false(size(Wobs(i,:)));
[hval{i},wval{i}]=TestEvaluation(Hobs(i,:),xbreak{i},p{i},nReg);
UseReg{i}=[1 0 0 ];
