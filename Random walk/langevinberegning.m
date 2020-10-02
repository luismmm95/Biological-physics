% Eksempel på sammenligning mellom teori og eksperiment for 
% ligning 6 (første leddet) i Brownske bevegelser

%mu=8.9e-4 % dynamic viscosity of water at 25 C
mu=9.5e-4 % dynamic viscosity of water at 22 C
%mu=8.1e-4 % ved 29 C
%mu = 7.8e-4 % ved 31 C
boltz=1.38064852e-23
temp=273.15+22 % temperatur i kelvin
r=0.99e-6/2   % kuleradius

fps = 15 % frames per scond - graffen viser per tidssteg!
xyfaktor = 2 % både x og y, teori er bare x

teori = boltz*temp/(mu*r*3*pi)

% 432 piksler per 0.10 mm (100 um)
kalib=100/432*1e-6; 

stigningstallet = 1.8936

dx2dt = stigningstallet/xyfaktor*kalib^2*fps

fprintf('dx^2/dt målt %e\n',dx2dt)
fprintf('dx^2/dt teori %e\n',teori)
fprintf('målt/teori %e\n',dx2dt/teori)
