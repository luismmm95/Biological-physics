%Intialize
clear all
%Sett st�rrelsen p� boksen
xsize=1;
ysize=1;
%Velg antall tidssteg, antall partikler og hvor langt partiklene skal
%flytte seg for hvert tidssteg
timesteps=389; % antall tidssteg
numpart=1000; % antall partikler
stepsize=xsize/500; %1/500 av boksst�rrelsen f�r det til � ligne p� mikroskopbildene
%Initialiser vektorene
xlist_all=zeros(timesteps,numpart);
ylist_all=zeros(timesteps,numpart);
poslist=zeros(timesteps*numpart,3);
tvec=ones(numpart,1);
%Sett tilfeldige startposisjoner som er fordelt over et areal 4x boksarealet. 
%Dette gj�res for at det alltid skal vises omtrent like mange partikler i
%boksen.
xy=2*xsize*(rand(numpart,2)-.5)+ones(numpart,2)*xsize/2;
xlist_all(1,:)=xy(:,1)';
ylist_all(1,:)=xy(:,2)';
%Finn de partiklene som er inni boksen
pnum=find(xy(:,1)<xsize&xy(:,1)>0&xy(:,2)<ysize&xy(:,2)>0);
xyin=xy(pnum,:);
plistlen=length(pnum);
%Lagre posisjonene til de som er inni boksen
poslist(1:plistlen,:)=[xyin 0*ones(size(pnum))];
%Vis posisjonene til de som er inni boksen
figure(1), plot(xy(:,1),xy(:,2),'ko','MarkerFaceColor','k'),axis([0 xsize 0 ysize])
for k=1:timesteps
    %Generer normalfordelte tilfeldige forflytninger med middelverdi 0 og
    %standardavvik 1
    xymove=stepsize*(randn(numpart,2));
    xy=xy+xymove;
    xlist_all(k,:)=xy(:,1)';
    ylist_all(k,:)=xy(:,2)';
    %Vis posisjonene til de som er inni boksen
    figure(1), plot(xy(:,1),xy(:,2),'ko','MarkerFaceColor','k'),axis([0 xsize 0 ysize]),
    %Lagre posisjonene til de som er inni boksen
    pnum=find(xy(:,1)<xsize&xy(:,1)>0&xy(:,2)<ysize&xy(:,2)>0);
    xyin=xy(pnum,:);
    padd=length(pnum);
    poslist(plistlen:plistlen+padd-1,:)=[xyin k*ones(size(pnum))];
    plistlen=plistlen+padd;
end
%Analyser partikkelbanene til alle partiklene (ikke bare dem som er
%inni boksen)
%Trekk fra startposisjonen til hver partikkel og beregn kvadratisk
%forflytning til alle partiklene
x0=xlist_all-ones(timesteps,1)*xlist_all(1,:);
xsq=x0.^2;
y0=ylist_all-ones(timesteps,1)*ylist_all(1,:);
ysq=y0.^2;
figure(2)
msd_x=mean(xsq,2)/stepsize^2;
msd_y=mean(ysq,2)/stepsize^2;
msd_xy=mean(xsq+ysq,2)/stepsize^2;
t=1:timesteps;
plot(t,msd_x,'b'),hold on,plot(t,msd_y,'g'),plot(t,msd_xy,'r'),hold off
xlabel('timesteps')
ylabel('Mean square displacement, stepsize^2')
px=t(:)\msd_x(:);
py=t(:)\msd_y(:);
pxy=t(:)\msd_xy(:);
legend(['<x^2>, d<x^2>/dt=' num2str(px)],...
    ['<y^2>, d<y^2>/dt=' num2str(py)],...
    ['<x^2+y^2>, d<x^2+y^2>/dt=' num2str(pxy)],'Location','NorthWest')

hold on, plot(t,px*t,'b'), plot(t,py*t,'g'), plot(t,pxy*t,'r'), hold off
title('Virrevandring med normalfordelte tilfeldige forflytninger')
axis tight


