function [xbreak,sse_end,phat_end]=WH_EIV(xall,yall,p0,sigma_uu,sigma_ee,xbreak0,nReg,HighOut)
    
    if ~HighOut
        Aeq=[1 0 0 0;
             0 0 0 1;];
        beq=[xbreak0(1);
             xbreak0(4);];
    else
        Aeq=[1 0 0 0;  %constrain solution for index 1,2 and 4. 
             0 0 1 0;  
             0 0 0 1;];          
         beq=[xbreak0(1);
              xbreak0(3);
              xbreak0(4);];             
    end
    
    Aie=[1 -1 0 0;  %constrain xbreak to be monotonic
       0 1 -1 0;
       0 0 1 -1;];
    bie=zeros(3,1);

    OptCon1=optimoptions('fmincon','Display','notify');       
    xbreak=fmincon(@RegSSE,xbreak0,Aie,bie,Aeq,beq,[],[],[],OptCon1);
    [sse_end,phat_end]=RegSSE(xbreak);
    
    function [sse,phat]=RegSSE(xbreak)

    A=zeros(nReg-1,nReg*2); b=zeros(nReg-1,1);    
    for j=1:nReg-1
        c1=(j-1)*2+1;
        A(j,c1:c1+3)=[xbreak(j+1) 1 -xbreak(j+1) -1];
    end   
    
    plb=-inf(size(p0));
    for j=1:nReg
        c1=(j-1)*2+1;
        plb(c1)=0;
    end
    
    OptCon=optimoptions('fmincon','Display','notify');    
    [phat,sse]=fmincon(@Reg_EIV_SSE,p0,[],[],A,b,plb,[],[],OptCon);  

        function sse = Reg_EIV_SSE(pall)

        sseReg=nan(nReg,1);        
        for i=1:nReg
            iReg=xall>=xbreak0(i) & xall<xbreak0(i+1);
            p=[pall( (i-1)*2 +1 ) pall( (i-1)*2+2 )];
            x=xall(iReg);
            y=yall(iReg);
            sseReg(i) = EIV_SSE(p);            
        end                                   
        sse=sum(sseReg);

            function sse = EIV_SSE(p)
            sse=(sigma_ee+p(1)^2*sigma_uu)^(-1) * sum( (p(1).*x+p(2)-y).^2 );
            end

        end

    end
end