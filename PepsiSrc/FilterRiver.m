%filters out cross-sections that are not in good reaches
Nodes=Rivers(cf).Nodes;
Reaches=Rivers(cf).Reaches;

Rivers(cf).Nodes=[];
Rivers(cf).Reaches=[];

numbnodes=size(Nodes.Z); 
numbnodes=numbnodes(1);
goodnodes= Nodes.xs_rch == Rivers(cf).gdrch(1);
for cr=2:length(Rivers(cf).gdrch)
    dummy=Nodes.xs_rch == Rivers(cf).gdrch(cr);
    goodnodes= goodnodes|dummy;
end
dummy=1:length(Nodes.xs_rch);
goodnodes=dummy(goodnodes);
Rivers(cf).Nodes.Z=Nodes.Z(goodnodes,:);
Rivers(cf).Nodes.xs_rch=Nodes.xs_rch(goodnodes,1);
Rivers(cf).Nodes.X=Nodes.X(goodnodes,1);
Rivers(cf).Nodes.W=Nodes.W(goodnodes,:);
Rivers(cf).Nodes.Q=Nodes.Q(goodnodes,:);
Rivers(cf).Nodes.H=Nodes.H(goodnodes,:);
Rivers(cf).Nodes.A=Nodes.A(goodnodes,:);
Rivers(cf).Nodes.P=Nodes.P(goodnodes,:);
Rivers(cf).Nodes.n=Nodes.n(goodnodes,:);
%filters out not good reachesRivers(cf).Nodes.n=ncread([pathtoncfiles Files(cf).name],'/XS_Timeseries/n');
Rivers(cf).Reaches.W=Reaches.W(Rivers(cf).gdrch,:);
Rivers(cf).Reaches.Q=Reaches.Q(Rivers(cf).gdrch,:);
Rivers(cf).Reaches.H=Reaches.H(Rivers(cf).gdrch,:);
Rivers(cf).Reaches.S=Reaches.S(Rivers(cf).gdrch,:);
Rivers(cf).Reaches.A=Reaches.A(Rivers(cf).gdrch,:);
Rivers(cf).Reaches.P=Reaches.P(Rivers(cf).gdrch,:);
for cr=1:length(Rivers(cf).gdrch)
    Rivers(cf).Reaches.x(cr)=mean(Nodes.X(Nodes.xs_rch == Rivers(cf).gdrch(cr)));
end

clear Nodes Reaches