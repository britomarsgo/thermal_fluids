%%%%% Subroutine computes thermophysical properties of air and residual gases mixture 

function [spp,h,cp,cv,R,u,v,Y,Ru,F,Fs,cpp,cvp,vp]=propres1(p,T,phi,f)

% p=100e3; % Pressão em Pascal
% T=900; % TEmperatura em Kelvin
% phi=0.8 % mistura
% f=0.5; % sem gases residuais presentes;

[A]=propar(T);

SVal=4.184e3/8.31434;
SVec=SVal*[1e-3 1e-6 1e-9 1e-12 1e3 1 1];

% alpha=8;beta=18;gamma=0;delta=0; % ISOOCTANO
% Afuel=[-0.55313 181.62 -97.787 20.402 -0.03095 -60.751 27.2162/SVal].*SVec;

% alpha=8;beta=18;gamma=0;delta=0;
% Afuel=[6.678e-1 8.398e-2 -3.334e-5 0 0 -3.058e+4 2.351e+1];

% case'diesel '%Ferguson
alpha=14.4;
beta=24.9;
gamma=0;
delta=0;

Afuel=[7.9710 1.1954E-01 -3.6858E-05 0 0 -1.9385E+04 -1.7879];

Ru=8314.34;%J/kmol.K
M=[44.01 18.02 28.008 32.000 28.01 2.018]';%kg/kmol - ESTES VALORES ESTÂO EM FERGUSON- Matriz Coluna Massa Molecular [CO2 H20 N2 O2 CO H2]
MinMol=1e-25;
eps=0.210/(alpha +0.25*beta-0.5*gamma);

if phi<=1.0% Stoichiometric or lean mixture
% Partindo da Equação de Ferguson ( Ver revisão Bibliografica) CO=0 e H2=0
%nu=[b c d e 0 0]';%[CO2 H2O N2 O2 CO H2 ] numero de moles totais para cada produto
nu=[alpha*phi*eps beta*phi*eps/2 0.79+delta*phi*eps/2 0.21*(1-phi) 0 0]'; %[CO2 H2O N2 O2 CO H2 ] numero de moles totais para cada produto
dcdT=0;

else % Rich mixture
z=1000/T;
K=exp(2.743+z*(-1.761+z*(-1.611+z*0.2803)));
dKdT=-K*(-1.761+z*(-3.222+z*0.8409))/1000;
a=1-K;
b=0.42-phi*eps*(2-gamma)+K*(0.42*(phi-1)+alpha*phi*eps);
c=-0.42*alpha*phi*eps*(phi-1)*K;
nu5=(-b+sqrt(b^2-4*a*c))/2/a;
dcdT=dKdT*(nu5^2-nu5*(0.42*(phi-1)+alpha*phi*eps)+ 0.42*alpha*phi*eps*(phi-1))/(2*nu5*a+b);
nu=[alpha*phi*eps-nu5 0.42-phi*eps*(2*alpha-gamma)+nu5 0.79+delta*phi*eps/2 0 nu5 0.42*(phi-1)-nu5]';
end
    
 % Molar and mass fraction of the combustion products
tmoles=sum(nu); % Numero total de moles dos produtos / moles de ar
Y=nu/tmoles ;% Fração molar de cada produto divido pelo numero total de moles dos produtos
Mres=sum(Y.*M); % Peso molecular dos produtos Kg/kmol 
X=(M/Mres).*Y; % Fração em Massa de cada Especie

% Molar fraction of reacting species
fuel=eps*phi/(1+eps*phi); % FRAÇÃO MOLAR DE COMBUSTIVEL em uma mistura com ar fuel-air
o2=0.21/(1+eps*phi); % FRAÇÃO MOLAR DE OXIGÊNIO
n2=0.79/(1+eps*phi); % FRAÇÃO MOLAR DE NITROGÊNIO
Mf=fuel*(12.01*alpha+1.008*beta+16*gamma+14.01*delta);% Kg/kmol
Mar=32*o2+28.02*n2; % Kg/kmol
Mfa=Mf+Mar;% Kg/kmol

Fs=Mf/Mar; % Razão combustivel /ar estequimetrica
% Fs=eps*(12.01*alpha+1.008*beta+16*gamma+14.01*delta)
F=phi*Fs; % Razão Combustivel / ar real

% Properties of air as an ideal gas[ range from 273 - 1800 K] paper Liberec 
Tar0=[1 T T^2 T^3]';
Rar=Ru/Mar;% J/kg.K;
C=[28.11 0.1967e-2 0.4802e-5 -1.966e-9]';
cpar=sum(C.*Tar0);
cpar=(cpar/Mfa)*1000; % J/kg.K
har=cpar*T;
uar=har-Rar*T;
sar=-Rar*log(p/101.325e3);
var=Rar*T/p;

%Thermodynamic properties of products 
Tcp0=[1 T T^2 T^3 T^4]';% EQUAÇÂO 3.22 
Th0=[1 T/2 T^2/3 T^3/4 T^4/5 1/T]';% EQUAÇÂO 3.23
Ts0=[log(T) T T^2/2 T^3/3 T^4/4 1]';% EQUAÇÂO 3.24
cp0=A(1:6,1:5)*Tcp0;
h0=A(1:6,1:6)*Th0 ;
s0=A(1:6,[1:5 7])*Ts0 ;

%COmbustion products properties
Rp=Ru/sum(Y.*M);% CONSTANTE UNIVERSAL PARA PRODUTOS DA COMBUSTÃO
hp=sum(h0.*Y);
cpp=sum(cp0.*Y);
sp=sum((s0.*Y));
hp=Rp*T*hp;
up=hp-Rp*T;
vp=Rp*T/p;
sp=Rp*(-log(p/101.325e3)+sp);
cpp=Rp*cpp;
cvp=cpp-Rp;

% Fuel properties
Mfuel=12.01*alpha+1.008*beta+16.000*gamma+14.01*delta ;% MASSA DE COMBUSTIVEL
a0=Afuel(1);b0=Afuel(2);c0=Afuel(3);d0=Afuel(6);e0=Afuel(7); % INDICIES PARA DIESEL FERGUSON
cpfuel=Afuel(1:5)*[1 T T^2 T^3 1/T^2]'; % CALOR ESPECIFICO DO COMBUSTIVEL
hfuel=Afuel(1:6)*[1 T/2 T^2/3 T^3/4 -1/T^2 1/T]';% ENTALPIA DO COMBUSTIVEL
s0fuel=Afuel([1:5, 7])*[log(T) T  T^2/2 T^3/3 -1/T^2/2 1]';% ENTROPIA DO COMBUSTIVEL

Rfuel=Ru/Mfuel;
hfuel=Rfuel*T*hfuel;
ufuel=hfuel-Rfuel*T;
vfuel=Rfuel*T/p;
s0fuel=Rfuel*(-log(p/101.325e3)+s0fuel);
cpfuel=Rfuel*cpfuel;

%Properites of the mixture air + residual products of combustion
R=Rar*(1-f)+Rp*f;
h=har*(1-f)+hp*f;
u=uar*(1-f)+up*f;
v=var*(1-f)+vp*f;
cp=cpar*(1-f)+cpp*f;
spp=sp*(1-f)+sp*f;
cv=cp-R;

end
