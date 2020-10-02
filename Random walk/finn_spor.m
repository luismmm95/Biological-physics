% Lasting av bilder
filename = 'tmp/30s_40x_lav.avi'
%sjekk om filen ligger i Matlab search path:
if exist(filename,'file')
    movieobj = VideoReader(filename); % lager et <1x1 mmreader> object
    get(movieobj);
else %les inn manuelt
    disp(['Finner ingen videofil med navn: ' filename])
    disp 'Last inn manuelt:'
    [filename,filepath]=uigetfile('*.avi','Velg videofilen du vil analysere');
    filename=[filepath filename];
    movieobj=VideoReader(filename);
    get(movieobj);
end

% Info fra movieobj
nFrames = movieobj.NumberOfFrames; % f�r tak i antall frames 
width = movieobj.Width; % f�r tak i px bredden
height = movieobj.Height; % f�r tak i px h�yden
frames=read(movieobj,1);
a=double(frames(:,:,2,1));
figure(1)
imagesc(a), colorbar;
dim_part_c=input('Zoom inn p� noen av partiklene i bildet og angi en typisk diameter (udelelig med 4!!): ','s');
dim_part=str2num(dim_part_c);

% diameteren til partikkelen b�r v�re et partall som ikke er delelig med 4, 
% hvis ikke gir centrd() mange feilmeldinger. Sjekk om dim_part er odd:
if rem(dim_part,2)
    dim_part = dim_part+1; % n� er dim_part er partall
end
if rem(dim_part,4)==0
    error('partikkeldiameteren skal ikke v�re delelig med 4!!!')
    return;   
end
fprintf('dim_part = %f\n',dim_part)

b=bpass(a,1,dim_part);
imagesc(b), colorbar;
threshold_c=input('Zoom inn p� noen av partiklene i bildet og angi en intensitetsverdi "threshold" for � skille partikler fra bakgrunn: ','s');
threshold=str2num(threshold_c);

fprintf('\n')
fprintf('************************************************************\n')
fprintf('Slett ikke figuren - den gjenbrukes i analyse_spor!\n')
fprintf('************************************************************\n')

poslist=[]; %Initialiserer posisjonslisten 

m=1;
n=1;
% Samler opp statistikk
h = waitbar(0,'Finner partikler'); %dette tar litt tid, s� her er noe � se p� imens...
for k=1:nFrames
  
    frames=read(movieobj,k);
    a=double(frames(:,:,2,1));
    b=bpass(a,1,dim_part);
    pk = pkfnd(b,threshold,dim_part/2);
    cnt = cntrd(b,pk,dim_part/2);  % finner midtpunktet til hver partikkel
    % (ignorer feilmeldingen dere f�r her)
    
%     %Test
%     figure(k),imagesc(b),hold on,plot(cnt(:,1),cnt(:,2),'go'),hold off
    
    np=size(cnt,1); %number of particles
    
    poslist(m:m+np-1,1:3)=[cnt(:,1:2) n*ones(np,1)];
    
    m=m+np;
    n=n+1;
    waitbar(k/nFrames,h,['Finner partikler. ' num2str(floor(100*k/nFrames)) '% fullf�rt'])
end
close(h)

trackmatrix = track(poslist,10);

trackmatrix_new = [];
k=1;
m=1;
% Beholder bare partiklene som er med i alle bildene
for i=1:trackmatrix(end)
    T = find(trackmatrix(:,4) == i);
    tsize = size(T);
    if(tsize(1) == nFrames)
       trackmatrix_new(m:m+nFrames-1,1:4) = [trackmatrix(T,1:3) k*ones(nFrames,1)];
       m=m+nFrames;
       k=k+1;
    end
    
end

save('trackmatrix','trackmatrix')
save('trackmatrix_new','trackmatrix_new')