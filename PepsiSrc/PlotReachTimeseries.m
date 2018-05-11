%temporary: define time indices
[~,nt]=size(Rivers(cf).Reaches.Q);
t=1:nt;

figure
subplot(2,2,1)
title([Rivers(cf).Name ' - reaches'])
plot(t,Rivers(cf).Reaches.Q')
xlabel('Time Index')
ylabel('Reach averaged discharge in m^3/s')
subplot(2,2,2)
plot(t,Rivers(cf).Reaches.H)
xlabel('Time Index')
ylabel('Reach averaged elevations (m)')
subplot(2,2,3)
plot(t,Rivers(cf).Reaches.W)
xlabel('Time Index')
ylabel('Reach averaged widths (m)')
subplot(2,2,4)
plot(t,Rivers(cf).Reaches.S*100000)
xlabel('Time Index')
ylabel('Reach averaged slopes (cm/km)')
subplot(2,2,1)
title([Rivers(cf).Name ' - reaches'])