% Random walk and mean square displacement
%
Nwalk=1000;
Nstep=1000;
% Generate Nwalk random walkers that move Nstep steps randomly to the left
% (x=-1) or to the right (x=+1)
x=2*randi(2,Nstep,Nwalk)-3;
% Sum up the steps to measure how they move (displacements/ trajectories)
y=cumsum(x,1);
% plot all the trajectories
%figure(1), plot(y)
% calculate the mean square displacement
%i.e. the mean of all the walkers' squared displacements
msd=mean(y.^2,2);
% plot the msd, slope is 2Dt
figure(2), plot(msd)
[counts_Ns,centers_Ns]=hist(y(Nstep,:),30);
[counts_Nsh,centers_Nsh]=hist(y(Nstep/2,:),30);
[counts_Nsq,centers_Nsq]=hist(y(Nstep/4,:),30);
figure(3), plot(centers_Ns,counts_Ns,'-or',centers_Nsh,counts_Nsh,'-db',centers_Nsq,counts_Nsq,'-sg')
figure(4), plot(centers_Ns,counts_Ns/Nwalk,'-or',centers_Nsh*sqrt(2),counts_Nsh/Nwalk,'-db',centers_Nsq*sqrt(4),counts_Nsq/Nwalk,'-sg')
% Calculate second moment of distribution = 2Dt
twoDt_Ns = sum(centers_Ns.^2.*counts_Ns)/Nwalk;
twoD_Ns=twoDt_Ns/Nstep
twoDt_Nsh = sum(centers_Nsh.^2.*counts_Nsh)/Nwalk; 
twoD_Nsh=twoDt_Nsh/(Nstep/2)
twoDt_Nsq = sum(centers_Nsq.^2.*counts_Nsq)/Nwalk; 
twoD_Nsq=twoDt_Nsq/(Nstep/4)
