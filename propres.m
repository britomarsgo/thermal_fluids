% SUBROTINA PARA CALCULO DAS PROPRIEDADES TERMODINAMICAS NA COMBUSTÃO


% 1º CALCULA PROPRIEDADES DO [(AR + AR RESIDUAL) + COMBUSTIVEL PARA FASE
% ANTES DA  COMBUSTÃO]


% 2º RETORNA OS VALORES DAS PROPRIEDADES TERMODINAMICAS PARA OS PRODUTOS DA
% COMBUSTÃO


function [MW Rmix hmix umix vmix smix cpmix cvmix Rp hp up vp sp cpp cvp ]=propres(p,T,phi,f) ;

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
table=[-1 1 0 0 1 -1]';
M=[44.01 18.02 28.008 32.000 28.01 2.018]';%kg/kmol - ESTES VALORES ESTÂO EM FERGUSON- Matriz Coluna Massa Molecular [CO2 H20 N2 O2 CO H2]

MinMol=1e-25;
dlvlT=1;% DERIVADAS PARCIAIS
dlvlp=-1; % DERIVADAS PARCIAIS
a=(alpha +0.25*beta-0.5*gamma);
eps=0.210/a;

if phi<=1.0% Mistura Estequimetrica ou pobre
% Partindo da Equação de Ferguson ( Ver revisão Bibliografica) CO=0 e H2=0
%nu=[b c d e 0 0]';%[CO2 H2O N2 O2 CO H2 ] numero de moles totais para cada produto
nu=[alpha*phi*eps beta*phi*eps/2 0.79+delta*phi*eps/2 0.21*(1-phi) 0 0]'; %[CO2 H2O N2 O2 CO H2 ] numero de moles totais para cada produto
dcdT=0;
else % MISTURA RICA
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
    
 % FRAÇÃO MOLAR E MASSA MOLAR DOS PRODUTOS
tmoles=sum(nu);% Numero total de moles dos produtos / moles de ar
Y=nu/tmoles; % Fração molar de cada produto divido pelo número total de moles dos produtos
Mres=sum(Y.*M); % Peso molecular dos produtos Kg/kmol 
X=(M/Mres).*Y ;% Fração em Massa de cada Especie


% FRAÇÃO MOLAR DOS REAGENTES
% Massa total dos reagentes por mol de Oxigênio.
Fs= eps*(12.01*alpha+1.008*beta+16*gamma+14.01*delta)/28.85;% Razão Combustivel-ar estequimetrico;
F= phi*Fs; % Razão Combustivel / Ar - Real
Lambdaest=1/Fs; % Lambda Razão Ar/Combustivel
F= phi*Fs; % Razão Combustivel / Ar - Real

fuel=eps*phi/(1+eps*phi); % FRAÇÃO MOLAR DE COMBUSTIVEL
o2=0.21/(1+eps*phi); % FRAÇÃO MOLAR DE OXIGÊNIO
n2=0.79/(1+eps*phi); % FRAÇÃO MOLAR DE NITROGÊNIO
Mfa=fuel*(12.01*alpha+1.008*beta+16*gamma+14.01*delta)+ 32*o2+28.02*n2; % Massa dos reagentes Kg/kmol
Mar=32*o2+28.02*n2; % Kg/kmol
fuelest=eps/(1+eps);% Fração de Combustivel Estequimetrico molar
fuelestm=Fs/(1+Fs); % Fração de Combustivel Estequimetrico massa
% Fração Molar de gás Combustivel-Ar-residuais
%Yres=f/(f+Mres/Mfa*(1-f));

%%%%CALCULO DA FRAÇÃO RESIDUAL DE UMA MISTURA - AR RESIDUAL + AR FRESCO + COMBUSTIVEL %%%%
Yres=(1+ Mres/Mfa*(1/f -1))^-1; % Fração residual Molar. EQUAÇÃO 3.37 do FERGUSON.
%Yres1=(1+ Mres/Mar*(1/f -1))^-1 % Fração de Produtos Residuais e Massa de ar fresco.
Y=Y*Yres;
Yfuel=fuel*(1-Yres);% Fração molar do combustivel
Y(3)=Y(3)+n2*(1-Yres); % fração de N2
Y(4)=Y(4)+o2*(1-Yres); % fração de O2

%CALCULO DAS PROPRIEDADES TERMODINAMICAS DOS PRODUTOS
Tcp0=[1 T T^2 T^3 T^4]'; % EQUAÇÂO 3.22 
Th0=[1 T/2 T^2/3 T^3/4 T^4/5 1/T]';% EQUAÇÂO 3.23
Ts0=[log(T) T T^2/2 T^3/3 T^4/4 1]';% EQUAÇÂO 3.24
cp0=A(1:6,1:5)*Tcp0;
h0=A(1:6,1:6)*Th0;% 
s0=A(1:6,[1:5 7])*Ts0; %

%CALCULO DAS PROPRIEDADES TERMODINAMICAS DO COMBUSTIVEL
Mfuel=12.01*alpha+1.008*beta+16.000*gamma+14.01*delta; % MASSA DE COMBUSTIVEL
a0=Afuel(1);b0=Afuel(2);c0=Afuel(3);d0=Afuel(6);e0=Afuel(7); % INDICIES PARA DIESEL FERGUSON
cpfuel=Afuel(1:5)*[1 T T^2 T^3 1/T^2]'; % CALOR ESPECIFICO DO COMBUSTIVEL
hfuel=Afuel(1:6)*[1 T/2 T^2/3 T^3/4 -1/T^2 1/T]';% ENTALPIA DO COMBUSTIVEL
s0fuel=Afuel([1:5, 7])*[log(T) T  T^2/2 T^3/3 -1/T^2/2 1]';% ENTROPIA DO COMBUSTIVEL

%DEFINE VALOR MINIMO OARA A COMPOSIÇÃO DAS ESPECIES
%CALCULO DO LOG
if Yfuel<MinMol
Yfuel=MinMol;
end
i=find(Y<MinMol); % Retorna os indices da Matrix Y<MinMol não nulos
Y(i)=ones(length(i),1)*MinMol;

%PROPRIEDADES DA MISTURA DE AR RESIDUAL + VAPOR DE COMBUSTIVEL PARA FASE
%ANTES DA COMBUSTÃO

h=hfuel*Yfuel+sum(h0.*Y);
s=(s0fuel-log(Yfuel))*Yfuel+sum((s0-log(Y)).*Y);
cp=cpfuel*Yfuel+sum(cp0.*Y)+sum(h0.*table*T*dcdT*Yres/tmoles);
MW=Mfuel*Yfuel+sum(Y.*M);
Rmix=Ru/MW;
hmix=Rmix*T*h;
umix=h-Rmix*T;
vmix=Rmix*T/p;
smix=Rmix*(-log(p/101.325e3)+s);
cpmix=Rmix*cp; % calculo o cp da mistura ar+ar residual + vapor de combustivel
cvmix =cpmix-Rmix;


% CALCULO DOS PRODUTOS DA COMBUSTÃO(PÓS QUEIMA) CONSIDERANDO APENAS
%%% CO2 H20 N2 O2

 % FRAÇÃO MOLAR E MASSA MOLAR DOS PRODUTOS
tmoles=sum(nu); % Numero total de moles dos produtos / moles de ar
Y1=nu/tmoles ;% Fração molar de cada produto divido pelo número total de moles dos produtos
Mres1=sum(Y1.*M); % Peso molecular dos produtos Kg/kmol 
X=(M/Mres1).*Y1; % Fração em Massa de cada Especie

%PROPRIEDADES DOS PRODUTOS DA COMBUSTÃO
Rp=Ru/sum(Y1.*M);% CONSTANTE UNIVERSAL PARA PRODUTOS DA COMBUSTÃO
hp=sum(h0.*Y1);
cpp=sum(cp0.*Y1);
Y1=log(Y1(1:4));
Y1(5)=0;
Y1(6)=0;
sp=sum((s0-Y1).*Y1);
hp=Rp*T*hp;
up=hp-Rp*T;
vp=Rp*T/p;
sp=Rp*(-log(p/101.325e3)+sp);
cpp=Rp*cpp;
cvp=cpp-Rp;


% %PROPRIEDADES DO COMBUSTIVEL
% Mfuel=Mfuel*Yfuel
% R3=Ru/Mfuel; % CONSTANTE UNIVERSAL PARA COMBUSTIVEL
% hf=hfuel*Yfuel;
% cpf=cpfuel*Yfuel;
% sf=(s0fuel-log(Yfuel))*Yfuel;
% hf=R3*T*hf;
% uf=hf-R3*T;
% vf=R3*T/p;
% sf=R3*(-log(p/101.325e3)+sf);
% cpf=R3*cpf

% %PROPRIEDADES DOS PRODUTOS DA COMBUSTÃO
% R2=Ru/sum(Y.*M)% CONSTANTE UNIVERSAL PARA PRODUTOS DA COMBUSTÃO
% hp=sum(h0.*Y)
% cpp=sum(cp0.*Y)+ sum(h0.*table*T*dcdT*Yres/tmoles)
% sp=sum((s0-log(Y)).*Y)
% hp=R2*T*hp
% up=hp-R2*T
% vp=R2*T/p
% sp=R2*(-log(p/101.325e3)+sp)
% cpp=R2*cpp
% 
%PROPRIEDADES DO AR COMO GÁS IDEAL [ FAIXA DE 273 - 1800 K] PAPER LIBEREC
Tar0=[1 T T^2 T^3]';
Rar=Ru/Mar;% J/kg.K
C=[28.11 0.1967e-2 0.4802e-5 -1.966e-9]';
cpar=sum(C.*Tar0);
cpar=(cpar/Mar)*1000;% J/kg.K
har=cpar*T;
uar=har-Rar*T;
var=Rar*T/p;
cvar=cpar-Rar;


end
