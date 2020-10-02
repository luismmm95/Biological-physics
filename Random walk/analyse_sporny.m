% gjenbruk figure 1 fra finn_spor!
figure(1), hold on

tr = trackmatrix_new;
numpart=max(tr(:,4));
trny=[];
trnum = [];

%kontrollering av spor for � fjerne de partiklene som sitter fast
%loop gjennom alle partiklene for � kontrollere sporene.
figure(2), clf
hold on
m=1;
grense = 20; %sett inn et passende kriterium her, for � fjerne partikler som sitter fast

% 
figure(2)
for n=1:numpart 
    I=find(tr(:,4)==n);
%     %de tre linjene nedenfor gir muligheten til � sjekke spor manuelt
%     figure(1), plot(tr(I,1),tr(I,2))
%      godta=input('Godta/Forkast sporet? g/f: ','s');
%     if godta=='g'
%         trny=[trny;tr(I,:)];
%         trnum(m)=n;
%         m=m+1;
%         figure(2),plot(tr(I,1),tr(I,2),'g') %plot godkjente spor i gr�nt
%     else
%         figure(2)
%         plot(tr(I,1),tr(I,2),'r') % og ikke godkjente i r�dt
%     end
    % med de neste tre linjene tar man bare med partikler som har beveget
    % seg en viss lengde i x- og y-retning (alternativ til den manuelle
    % kontrollen)
    ddxx = max(tr(I,1))-min(tr(I,1));
    ddyy = max(tr(I,2))-min(tr(I,2));
    if ddxx > grense && ddyy > grense
        figure(2), plot(tr(I,1),tr(I,2))
        trny=[trny;tr(I,:)];
        trnum(m)=n;
        m=m+1;
        plot(tr(I,1),tr(I,2),'g') %plot godkjente spor i gr�nt
    else 
        plot(tr(I,1),tr(I,2),'r.') %plot godkjente spor i gr�nt
    end
   
   %plot(trny(I,1))
end

% Overskriv Figure 1 fra finn_spor med sporene
figure(1)
for n=1:numpart 
    I=find(tr(:,4)==n);
    ddxx = max(tr(I,1))-min(tr(I,1));
    ddyy = max(tr(I,2))-min(tr(I,2));
    if ddxx > grense && ddyy > grense
        plot(tr(I,1),tr(I,2),'g') %plot godkjente spor i grønt
    else 
        plot(tr(I,1),tr(I,2),'r.','MarkerSize',5) % og ikke godkjente i rødt
    end
end

figure(2),title 'Partikkelspor. Grønt spor: godkjent, rødt: forkastet'
hold off,axis equal
tr=trny;
numpart=length(trnum);

% kalkulere forflytningene
m=max(tr(:,3));
t=1:m;
sumdx2=zeros(m,1); sumdx4=zeros(m,1);
sumdy2=zeros(m,1); sumdy4=zeros(m,1); sumdxy4=zeros(m,1);
xstart = zeros(numpart,1); ystart = zeros(numpart,1);
xstop = zeros(numpart,1); ystop = zeros(numpart,1);
xstd = zeros(numpart,1); ystd = zeros(numpart,1);
for n=1:numpart % loop gjennom alle partiklene
    I=find(tr(:,4)==trnum(n));
    % finne hvor langt partikkelen har flyttet seg fra startposisjonen ved 
    % hvert tidssteg
    dx=tr(I,1)-tr(I(1),1);
    dy=tr(I,2)-tr(I(1),2);
    sumdx2=sumdx2+dx.^2;
    sumdy2=sumdy2+dy.^2;
    sumdx4=sumdx4+dx.^4;
    sumdy4=sumdy4+dy.^4;
    sumdxy4=sumdxy4+(dx.^2+dy.^2).^2;
    xstart(n) = tr(I(1),1);
    ystart(n) = tr(I(1),2);
    xstop(n) = tr(I(end),1);
    ystop(n) = tr(I(end),2);
    xstd(n) = std(diff(dx))^2;
    ystd(n) = std(diff(dy))^2;
    %fprintf('%d %f %f   %f %f\n',n,std(tr(I,1)),std(tr(I,2)),mean(tr(I,1)-tr(I(1),1)),mean(tr(I,2)-tr(I(1),2)))
end
msdx2=sumdx2/numpart;
msdy2=sumdy2/numpart;
msdxy2=(msdx2+msdy2);

msdx4 = sumdx4/numpart;
msdy4 = sumdy4/numpart;
msdxy4 = sumdxy4/numpart;

errdx = sqrt(msdx4-msdx2.^2)/sqrt(numpart);
errdy = sqrt(msdy4-msdy2.^2)/sqrt(numpart);
errdxy = sqrt(msdxy4-msdxy2.^2)/sqrt(numpart);

fprintf('%d tidssteg, %d partikler\n',m,numpart)

% Midlet forflytning alle sporet partikkler
fprintf('\nMidlet total lineær forflyttning: (%d partikkler)\n',numpart)
fprintf('      x-retning: %.1f +- %.1f piksler\n'  ,mean(xstop-xstart),std(xstop-xstart)/sqrt(length(xstop)))
fprintf('      y-retning: %.1f +- %.1f piksler\n\n',mean(ystop-ystart),std(ystop-ystart)/sqrt(length(ystop)))

% Midlet forflytning^2 alle sporet partikkler
fprintf('\nMidlet total kvadratisk forflyttning per tidssteg: (%d partikkler)\n',numpart)
fprintf('      x-retning: %.2f +- %.2f piksler^2/tidssteg\n'  ,mean((xstop-xstart).^2)/m,std((xstop-xstart).^2)/sqrt(length(xstop))/m)
fprintf('      y-retning: %.2f +- %.2f piksler^2/tidssteg\n\n',mean((ystop-ystart).^2)/m,std((ystop-ystart).^2)/sqrt(length(ystop))/m)

% Midlet standard avvik^2 alle sporet partikkler
fprintf('\n Midlet standard avvik^2 per tidssteg\n',numpart)
fprintf('      x-retning: %.3f +- %.3f piksler^2/tidssteg\n'  ,mean(xstd),std(xstd)/sqrt(numpart))
fprintf('      y-retning: %.3f +- %.3f piksler^2/tidssteg\n'  ,mean(ystd),std(ystd)/sqrt(numpart))


% Plotte start og stopp punkter
figure(1), plot(xstart,ystart,'ro')
plot(xstop,ystop,'ko')
title('Partikkelspor med start (r�d sirkel) og stopp (svart sirkel)')
%label('x [piksler]'), ylabel('y [piksler]')

% Plotte resultater 
figure(3), hold off
plotx=plot(msdx2,'b');
xband = plotx.XData';
hold on;
% Plot error in y as function of x as a transparent blue band
f=fill([xband;flip(xband)],[msdx2-errdx;flip(msdx2+errdx)],'b','EdgeColor','None');
f.FaceAlpha=0.3;

plot(msdy2,'g')
% Plot error in y as function of x as a transparent green band
f=fill([xband;flip(xband)],[msdy2-errdy;flip(msdy2+errdy)],'g','EdgeColor','None');
f.FaceAlpha=0.3;

plot(msdxy2,'r')
% Plot error in y as function of x as a transparent red band
f=fill([xband;flip(xband)],[msdxy2-errdxy;flip(msdxy2+errdxy)],'r','EdgeColor','None');
f.FaceAlpha=0.3;

% finne stigningstall
% tvinge løsningen til å gå gjennom origo.
px = t(:)\msdx2(:);
py = t(:)\msdy2(:);
pxy = t(:)\msdxy2(:);

xlabel('antall tidssteg')
ylabel('Midlere kvadratisk forflytning, piksler^2')
hold on, plot(t,px*t,'b'), plot(t,py*t,'g'), plot(t,pxy*t,'r'), hold off
title('Midlere forflytning av Brownske partikler')
axis tight

legend('<x^2>','\sigma_{<x^2>}','<y^2>','\sigma_{<y^2>}', '<x^2+y^2>','\sigma_{<x^2+y^2>}',...
    ['d<x^2>/dt=' num2str(px(1))],['d<y^2>/dt=' num2str(py(1))],['d<x^2+y^2>/dt=' num2str(pxy(1))],...
    'Location','NorthWest')

