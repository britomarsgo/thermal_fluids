% Full Thermo-Fluid in-house code developed to simulate internal combustion engine
% Code written from scratch by Alisson Vinicius Brito Lopes as part of M.Sc. in
% Thermal and Fluid Science during the years of 2013-2014 at UNICAMP-Brazil

% crackshaft delta degree used for numerical integration
deltathetap=0.01; 

% Definition of geometrical and initial model parameters
RPM=1800; % Engine speed
r=17;% Compression ratio
phi=0.8;% fuel/air equivalence ratio
Tw=450;% Engine wall (surface) temperature . 
thetaij=350;% Angle of injection start
thetab=75 ; % combustion duration. 15° pre-mixtre phase and 60° difusive phase.

%Dimensionless lift/diameter ratio
YDmaxe=0.30; % Maximum YD exhasut valve
YDmaxi=0.25; % Maximum YD intake valve

rve=-4; % Exhaust valve accelaration coefficient
rvi=-4; % Intake valve accelaration coefficient

De=38.76/1000; % exhaust valve diameter estimated based on TAYLOR for a Flame cylinder head
Di=48.88/1000; % intake valve diameter estimated based on TAYLOR for a Flame cylinder head

%Valve timing angles
% Intake
thetafva=930; %Intake valve closing angle - 30° after BDC
thetaava=700; %Intake valve opening - 30° before TDC

% Exhaust
thetaave=510; %Exhaust valve opening - 30° before BDC
thetafve=725; %Exhaust valve closing - 10 ° after TDC

% Model initial guessings
f=0.05;% Residual gas fraction 
p1=0.9*10^5 ;% Initial pressure [N/m²]
T1=330; % Initial temperature in K

% indexes  
ir=1;
je=1;  
ji=1;  
jee=1; 
jei=1; 
jo=1;  % index for valve overlapping 
joo=1;
jfo=1;
jfi=1;
jej=1;

%Memory allocation procedure

%Politropic compression
pf=[];  Vfcp=[];    thetacp=[]; z1=[];  nep=[]%
% Combustion
pc=[];  Tc=[];  Vc=[];  dQdtheta=[];    dQdthetaw=[];   dPdtheta=[];    dTdtheta=[];    dWdtheta=[] 
yc=[];  Awc=[]; thetac=[];  nx=[];  
% Politropic expansion
pfe=[];      Vep=[]; ne=[]%
% Exhaust process 
p_cyle=[];   Ve=[]; %
% Valve overlapping process
p_cylo=[];   Vo=[];%
%Intake process
p_cyla=[];   Vi=[];% 

%Numerical tolerances and initialisation
ER1=10;
ER2=10;
ER3=10;
tol1=1e-6; 
tol2=1e-6;
tol3=1e-6;

%%%% Main loop that tests the convergence of the entire model %%%%

while ( ER1>=tol1) && (ER2>=tol2 )&& (ER3>=tol3)
RPM=1800; 
phi=0.82;
Tw=450;
thetaij=350;
thetab=75 ; 
%Geometry of the MWM D220-4 Diesel Engine
b=0.102;% Bore (m)
curso=0.120; 
stroke=curso;% Stroke
l=0.207; %connecting rod (m)
eps=curso/(2*l); % Curso/2*l
RR=2*l/curso;
a=curso/2;
Vd=pi/4*b^2*curso; % displaced volume
Vpms=Vd/(r-1); %Volume at TDC
Vtdc=Vpms;%
Vpmi=pi/4*b^2*curso+Vpms;% Volume at BDC
Vbdc=Vpmi;
hch=(4*Vpms)/(pi*b^2); % TDC height
Ach=pi*b*hch * pi/4*b^2; % Combustion chamber area
Apc=pi/4*b^2; % Area of piston's top face

% Fuel data 
NC=52; % CEtane number for Diesel
PCI=42940*10^3; % LHV MJ

thetas=thetaij ;% Angle of combustion start
omega=RPM*pi/30; % Engine speed RAD/S
upmean=omega*stroke/pi;% Piston velocity

% Initial condition
theta1=thetafva; % Initial angle, piston 30° up BDC.
V1=Vd/(r-1) + Vd/2*(RR+1-cos(theta1*pi/180)-((RR^2-(sin(theta1*pi/180).^2)).^0.5));

Vthetafva=V1;

[spp,h,cp,cv,R,u,v,Y,Ru,F,Fs]=propres1(p1,T1,phi,f); % this function returns the thermophysical properties of air 

mass1=Vthetafva/v ; % mass1 = air mass + residual gas mass
U1=u*mass1;
S1=spp*mass1;
mf=(1-f)*mass1*F;

% Loop with the beginning of the polytropic compression process
z=(thetaij-(thetafva-720))/deltathetap; % Number of iteration
p=zeros(z+1,1);
T=zeros(z+1,1);

thetacp(1)=thetafva; % the initial angle is equal to the angle of intake valve closing
p(1)=p1;
T(1)=T1;

jj=1;
ii=1;
ll=1;
ss=1;
yyy=1;
uuu=1;

 while(jj<=z+1);
V(jj)=Vd/(r-1) + Vd/2*(RR +1-cos(thetacp(jj)*pi/180)-((RR^2-(sin(thetacp(jj)*pi/180).^2)).^0.5));
thetacp(jj+1)=thetacp(jj)+deltathetap;
jj=jj+1;
 end

while(ii<=z)
nep=1.35; % polytropic coefficient of compression - recommend to use between 1,34 -1,37
p(ii+1)=p(ii)*(V(ii)/V(ii+1))^nep;
T(ii+1)=T(ii)*(V(ii)/V(ii+1))^(nep-1);
thetacp(ii+1)=thetacp(ii)+deltathetap;
ii=ii+1;
end

% Handerberg parameters
pp(1)=100e3 ;% Pressure at the BDC [N/m²] - initial guessing
TT(1)=300; % Temperature at the BDC [K] - initial guessing
thetapms=360;
thetapmi=180;
thetacpp(1)=thetapmi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zzz=(thetapms-thetapmi)/deltathetap;

 while(uuu<=zzz+1);
Vt(uuu)=Vd/(r-1) + Vd/2*(RR +1-cos(thetacpp(uuu)*pi/180)-((RR^2-(sin(thetacpp(uuu)*pi/180).^2)).^0.5));
thetacpp(uuu+1)=thetacpp(uuu)+deltathetap;
uuu=uuu+1;
 end

while(yyy<=zzz)
nepp=1.35; 
pp(yyy+1)=pp(yyy)*(Vt(yyy)/Vt(yyy+1))^nep;
TT(yyy+1)=TT(yyy)*(Vt(yyy)/Vt(yyy+1))^(nep-1);
thetacpp(yyy+1)=thetacpp(yyy)+deltathetap;
yyy=yyy+1;

end

% This function returns the ignition delay time converted in cranckshaft angle

pp2=pp(zzz);% Pressure imidiately before the start of injection for a politropyc process
TT2=TT(zzz);% Temperature imidiately before the start of injection for a politropyc proces

%[IDms, IDtheta, thetaic]=igdelayarai(RPM,phi,p2,T2,thetaij)
[IDtheta, IDms,thetaic]=handerberg(RPM,upmean,TT2,pp2,Ru,NC,thetaij);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z1=(thetaic-(thetafva-720))/deltathetap;
pf=zeros(z1+1,1);
Tf=zeros(z1+1,1);
pf(1)=p1;
Tf(1)=T1;

thetacp(1)=thetafva;

while(ll<=z1+2);
Vfcp(ll)=Vd/(r-1) + Vd/2*(RR +1-cos(thetacp(ll)*pi/180)-((RR^2-(sin(thetacp(ll)*pi/180).^2)).^0.5)) ;
thetacp(ll+1)=thetacp(ll)+deltathetap;
ll=ll+1;
end

while(ss<=z1+1);
pf(ss+1)=pf(ss)*(Vfcp(ss)/Vfcp(ss+1))^nep;
Tf(ss+1)=Tf(ss)*(Vfcp(ss)/Vfcp(ss+1))^(nep-1);
thetacp(ss+1)=thetacp(ss)+deltathetap;
ss=ss+1;
end

%Heat and Work calculation during the reversible adiabatic compression
p2=pf(z1+1);
T2=Tf(z1+1);

[spp,h,cp,cv,R,u,v,Y,Ru,F,Fs]=propres1(p2,T2,phi,f);

theta2=thetacp(z1+1);

V2=Vd/(r-1) + Vd/2*(RR+1-cos(theta2*pi/180)-((RR^2-(sin(theta2*pi/180).^2)).^0.5));
U2=u*mass1;
S2=spp*mass1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p3=pf(z1);% Pressure after the ignition delay period and before the combustion start
T3=Tf(z1);% Temperature after the ignition delay period and before the combustion start

% 2ª part - Heat release using a double Wiebe function - parameters calibrated from experimental data

[x,mf1,dxdtheta,dQdthetac,nx,nx2,thetafq]=wiebeduplateste(mf,thetaic,PCI,deltathetap); % Double Wiebe function

Vc=[]; 
thetac(1)=(thetaic); % Angle of the combustion start (the ignition angle is accounted)
pc(1)=p3; % Pressure at the end of the politropyc compression at the start of the combustion process
Tc(1)=T3; % Temperature at the end of the politropyc compression at the start of the combustion process


for j=1:nx+1
Vc(j)=Vd/(r-1) + Vd/2*(RR +1-cos(thetac(j)*pi/180)-((RR^2-(sin(thetac(j)*pi/180).^2)).^0.5));
dVdthetac(j)=(Vtdc*(r-1)/2*(sin(thetac(j)*pi/180)+ eps/2*sin(2*thetac(j)*pi/180)./sqrt(1-eps^2*sin(thetac(j)*pi/180).^2)));
yc(j)=l+a-(((l^2-a^2*(sin(thetac(j)*pi/180).^2)).^0.5)+a*cos(thetac(j)*pi/180));
Awc(j)=pi*b*(yc(j))+pi/2*b^2;
thetac(j+1)=thetac(j)+deltathetap;% thetac(nx+1) = thetafq , the angle at the end of the combustion process
end

y=1;

while y<=nx
           
[h,cp,cv,R,u,v,Y,Ru,F,Fs]=propres1(pc(y),Tc(y),phi,f);

[MW, Rmix, hmix, umix, vmix, smix, cpmix, cvmix, Rp, hp, up, vp, sp, cpp, cvp]=propres(pc(y),Tc(y),phi,f); % this function returns the thermophysical properties of reacting gases

gammamix(y)= cpmix/cvmix;
gammap(y)= cpp/cvp;
gammaeq(y)= gammamix(y)*(1-x(y))+gammap(y)*(x(y));
cveq(y)=cvmix*(1-x(y))+ cvp*x(y);
Req(y)= Rmix*(1-x(y))+Rp*(x(y));

% Hohenberg correlation during combustion

khoh=130;
a1=1.4;
hh(y)=(khoh*(Vc(y)^(-0.06))*((pc(y)/100000)^(0.8))*(Tc(y)^(-0.4))*(upmean+a1)^(0.8));% 
n1=1/omega;
n3=n1*pi/180;
dQdthetaw(y)=hh(y)*Awc(y)*(Tc(y)-Tw)*n3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dQdtheta(y)=dQdthetac(y)- dQdthetaw(y);
dVdthetan(y+1)=(Vc(y+1)-Vc(y))/deltathetap;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration of the ODE set (dP, dT and dW) for the combustion phase    %% 
%         BUTCHER'S RUNGE KUTTA 5 (FIFTH) ORDER                       %%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k1= ((gammaeq(y)-1)/Vc(y))*(dQdthetac(y)-hh(y)*Awc(y)*(Tc(y)-Tw)*n3)-(gammaeq(y)*pc(y)/Vc(y))*(Vc(y+1)-Vc(y))/deltathetap;
b1=(1/(mass1*Req(y)))*(pc(y)*(Vc(y+1)-Vc(y))/deltathetap+Vc(y)*k1);

k2= ((gammaeq(y)-1)/Vc(y))*(dQdthetac(y)-hh(y)*Awc(y)*((Tc(y)+b1/4*deltathetap)-Tw)*n3)- (gammaeq(y)*(pc(y)+k1/4*deltathetap)/Vc(y))*(Vc(y+1)-Vc(y))/deltathetap;
b2=(1/(mass1*Req(y)))*((pc(y)+k1/4*deltathetap)*(Vc(y+1)-Vc(y))/deltathetap+Vc(y)*k1);

k3=((gammaeq(y)-1)/Vc(y))*(dQdthetac(y)-hh(y)*Awc(y)*((Tc(y)+b1/8*deltathetap + b2/8*deltathetap)-Tw)*n3)-(gammaeq(y)*(pc(y)+ k1/8*deltathetap + k2/8*deltathetap)/Vc(y))*(Vc(y+1)-Vc(y))/deltathetap;
b3=(1/(mass1*Req(y)))*((pc(y)+k1/8*deltathetap + k2/8*deltathetap)*(Vc(y+1)-Vc(y))/deltathetap+Vc(y)*k1);

k4=((gammaeq(y)-1)/Vc(y))*(dQdthetac(y)-hh(y)*Awc(y)*((Tc(y)-b2/2*deltathetap + b3*deltathetap)-Tw)*n3)-(gammaeq(y)*(pc(y)-k2/2*deltathetap + k3*deltathetap)/Vc(y))*(Vc(y+1)-Vc(y))/deltathetap;
b4=(1/(mass1*Req(y)))*((pc(y)-k2/2*deltathetap+k3*deltathetap)*(Vc(y+1)-Vc(y))/deltathetap+Vc(y)*k1);

k5=((gammaeq(y)-1)/Vc(y))*(dQdthetac(y)-hh(y)*Awc(y)*((Tc(y)+3/16*b1*deltathetap+9/16*b4*deltathetap)-Tw)*n3)-(gammaeq(y)*(pc(y)+3/16*k1*deltathetap + 9/16*k4*deltathetap)/Vc(y))*(Vc(y+1)-Vc(y))/deltathetap;
b5=(1/(mass1*Req(y)))*((pc(y)+3/16*k1*deltathetap + 9/16*k4*deltathetap)*(Vc(y+1)-Vc(y))/deltathetap+Vc(y)*k1);

k6=((gammaeq(y)-1)/Vc(y))*(dQdthetac(y)-hh(y)*Awc(y)*((Tc(y)-3/7*b1*deltathetap+2/7*b2*deltathetap+12/7*b3*deltathetap -12/7*b4*deltathetap + 8/7*b5*deltathetap)-Tw)*n3)-(gammaeq(y)*(pc(y)-3/7*k1*deltathetap + 2/7*k2*deltathetap - 12/7*k4*deltathetap +8/7*k5*deltathetap)/Vc(y))*(Vc(y+1)-Vc(y))/deltathetap;
b6=(1/(mass1*Req(y)))*((pc(y)-3/7*k1*deltathetap + 2/7*k2*deltathetap - 12/7*k4*deltathetap +8/7*k5*deltathetap)*(Vc(y+1)-Vc(y))/deltathetap+Vc(y)*k1);

pc(y+1,1)=pc(y)+ 1/90*(7*k1+ 32*k2+ 12*k4+ 32*k5 +7*k6)*deltathetap;
Tc(y+1,1)=Tc(y)+ 1/90*(7*b1+ 32*b2+ 12*b4+ 32*b5 +7*b6)*deltathetap;

dPdtheta(y)= 1/90*(7*k1+ 32*k2+ 12*k4+ 32*k5 +7*k6);
dTdtheta(y)= 1/90*(7*b1+ 32*b2+ 12*b4+ 32*b5 +7*b6);
dWdtheta(y)= pc(y)*dVdthetan(y);

y=y+1;

end

% Total mass at the end of the combustion process = MASS1(Mass air + Mass residuals) + fuel mass (mf)                
mt=mf+mass1;

% Polytropic expansion process
ne=gammaeq(nx); % Polytropic coefficient
oo=1;
uu=1;
thetaep(1)=thetafq;

% hx is equal to the number of iteration between the angle of the end of
% the combustion process(TEHTAFQ) until the opening angle of the exhaust valve (THETAAVE) 

hx=(thetaave-thetafq)/deltathetap;
pfe(1)=pc(nx+1);
Tfe(1)=Tc(nx+1);

while(oo<=hx+2)
Vep(oo)=Vd/(r-1) + Vd/2*(RR +1-cos(thetaep(oo)*pi/180)-((RR^2-(sin(thetaep(oo)*pi/180).^2)).^0.5)) ;
thetaep(oo+1)=thetaep(oo)+deltathetap;
oo=oo+1;
end

while(uu<=hx+1);
pfe(uu+1)=pfe(uu)*(Vep(uu)/Vep(uu+1))^ ne;
Tfe(uu+1)=Tfe(uu)*(Vep(uu)/Vep(uu+1))^(ne-1);
theta(uu+1)=thetaep(uu)+deltathetap;
uu=uu+1;
end

% Calculation of the expansion work
p4=pfe(1);
T4=Tfe(1);

[MW, Rmix, hmix, umix, vmix, smix, cpmix, cvmix, Rp, hp, up, vp, sp, cpp, cvp]=propres(p4,T4,phi,f);

theta4=thetaep(1);
V4=Vd/(r-1) + Vd/2*(RR+1-cos(theta4*pi/180)-((RR^2-(sin(theta4*pi/180).^2)).^0.5));
U4=up*mt;
S4=sp*mt;
p5=pfe(hx+1);
T5=Tfe(hx+1);

[MW, Rmix, hmix, umix, vmix, smix, cpmix, cvmix, Rp, hp, up, vp, sp, cpp, cvp]=propres(p5,T5,phi,f);

theta5=theta(hx+1);
V5=Vd/(r-1) + Vd/2*(RR+1-cos(theta5*pi/180)-((RR^2-(sin(theta5*pi/180).^2)).^0.5));
U5=up*mt;
S5=sp*mt

% Calculation of the exhaust process

 [YDe,YDi]= valvulageometria(YDmaxe,YDmaxi,rve,rvi,thetafve,thetaave,thetafva,thetaava,deltathetap); % returns valve geometric features

% [YDe]=valvulaescap(thetaave,thetafve,deltathetap);
% [YDi]=valvulaescap(thetaava,thetafva,deltathetap);

% Intake and exhaust valves lift 
Le=YDe*De % exhaust lift
Lmaxe=YDmaxe*De;% max exhaust lift
Li=YDi*Di % intake lift
Lmaxi=YDmaxi*Di; % max intake lift

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zze=length(YDe); % Number of elements in the vector of YD exhaust 
zzi=length(YDi); % Number of elements in the vector of YD intake 
ddde=1;
dddi=1;

thetaio=thetaave;  % thetaee = start opening exhaust valve
thetaeo=thetaave;

% Loops to compute the effective area of the exhaust valve according to Gallo(1990)and Kastner (197X)
while(ddde<=zze)
    
  if YDe(ddde)<=0.125
Aeffe(ddde)= 2.221441*Le(ddde)*(De+Le(ddde)/2); % Minimum flow area occurs for Y/D <=0.125
  end
  if YDe(ddde)>0.125 && YDe(ddde)<=0.274049;
Aeffe(ddde)=3.337942*De*(Le(ddde)^2 -Le(ddde)*De/8 +De^2/128)^0.5;
  end
 if YDe(ddde)>0.274049
Aeffe(ddde)=0.736311*De^2 ;
 end
 if Aeffe(ddde)<=0
     Aeffe(ddde)=0;
 end
thetaeo(ddde+1)=thetaeo(ddde) +deltathetap;
ddde=ddde+1;
end

% Loops to compute the effective area of the intake valve according to Gallo(1990)and Kastner (197X)
while(dddi<=zzi)
    
  if YDi(dddi)<=0.125
Aeffi(dddi)= 2.221441*Li(dddi)*(Di+Li(dddi)/2); % Minimum flow area occurs for Y/D <=0.125
  end
  if YDi(dddi)>0.125 && YDi(dddi)<=0.274049;
Aeffi(dddi)=3.337942*Di*(Li(dddi)^2 - Li(dddi)*Di/8 +Di^2/128)^0.5;
  end
 if YDi(dddi)>0.274049
Aeffi(dddi)=0.736311*Di^2 ;
 end
 if Aeffi(dddi)<=0
     Aeffi(dddi)=0;
 end
 
thetaio(dddi+1)=thetaio(dddi) +deltathetap;
dddi=dddi+1;
end

[YDe, Cde, Cdbe, YDi, Cdi, Cdbi]=coefdescarga(YDi,YDe) % returns the Cd values 

% The flow rate through the valve
nx4=(thetaava-thetaave)/deltathetap; %Number of iteration between OEV and OVI (Without valve overlapping - only expansion)
nx5=(thetafva-thetafve)/deltathetap; %Number of iteration between EVC and IVO (Without valve overlapping - only intake)
nx6=(thetafve-thetaava)/deltathetap; % (Overlapping period)- there is possibility of backflow
nx7=(thetafva-thetaave)/deltathetap; %Number of iterations of the entire open phase process - this is used to compute the cylinder volume during the open phase.

me(1)= mt;
T_cyle(1)=Tfe(hx+1);% Temperature at end of expansion process
p_cyle(1)=pfe(hx+1);% Pressure at end of expansion process

thetae(1)=thetaave; % Angle of exhaust valve opening start
thetai(1)=thetafve; % Angle of exhaust valve closeing start
thetao(1)=thetaava; % Angle used to compute the overlappig process

Tadm=300; % from Velasquez
Te=750; % Estimated Gas temperature inside the exhaust manifold

pe=1.05*10^5; % Intake manifold pressure [PA].
padm=0.85*10^5; % Exhasut manifold pressure [PA].

if thetaava == thetafve % this verifies whether there is overlapping or not
while je<=nx4+1
Ve(je)=Vd/(r-1) + Vd/2*(RR +1-cos(thetae(je)*pi/180)-((RR^2-(sin(thetae(je)*pi/180).^2)).^0.5));
dVdthetae(je)=(Vtdc*(r-1)/2*(sin(thetae(je)*pi/180)+ eps/2*sin(2*thetae(je)*pi/180)./sqrt(1-eps^2*sin(thetae(je)*pi/180).^2)))*pi()/180;
ye(je)=l+a-(((l^2-a^2*(sin(thetae(je)*pi/180).^2)).^0.5)+a*cos(thetae(je)*pi/180));
Awe(je)=pi*b*(ye(je))+pi/2*b^2;
thetae(je+1)=thetae(je)+deltathetap;
je=je+1;
end

% Pressure variation during the expansion period
while jee<=nx4+1
    
[MW, Rmix, hmix, umix, vmix, smix, cpmix, cvmix, Rp, hp, up, vp, sp, cpp, cvp]=propres(p_cyle(jee),T_cyle(jee),phi,f);

k_exh(jee)=1.28; % Polytropic cofficient of gases inside the exhaust manifold
cveqe(jee)=cvp;
cpeqe(jee)=cpp;
ueqe(jee)=up;
k_cyl(jee)=cpp/cvp;
heqe(jee)=hp;
Reqe(jee)= Rp;
pho(jee)=1/vp;
co(jee)=sqrt(k_cyl(jee)*Reqe(jee)*T_cyle(jee));

% Hohenberg correlation for the exhaust process
khohe=130;
a1e=1.4;
hhee(jee)=(khohe*(Ve(jee)^(-0.06))*((p_cyle(jee)/100000)^(0.8))*(T_cyle(jee)^(-0.4))*(upmean+a1e)^(0.8));

%Nishikwaki correlation for the exhaust period   
hhe(jee)=(25.2*(b^(-0.422))*((upmean*((p_cyle(jee)/100000)/0.9807))^(0.578))*(T_cyle(jee))^(-0.131))*1.161111; % Tranforming KCAL/H PARA J/S

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dQdthetawe(jee)=(hhe(jee)*Awe(jee)*(T_cyle(jee)-Tw))*n3;
dQdthetae(jee)= - dQdthetawe(jee);
dVdt(jee)=dVdthetae(jee);
Cdb(jee)=-1.5*YDe(jee)+ 1;


% Integration
if p_cyle(jee)>pe
    
    
    if (pe/p_cyle(jee))<=((2/(k_cyl(jee)+1))^(k_cyl(jee)/(k_cyl(jee)-1))) % sonic flow
    
        dmdte(jee)=-((Cde(jee)*Aeffe(jee)*p_cyle(jee)/sqrt(Reqe(jee)*T_cyle(jee)))*sqrt(k_cyl(jee))*((2/(k_cyl(jee)+1))^((k_cyl(jee)+1)/(2*(k_cyl(jee)-1)))))*n3;
        
    else % subsonic flow
      
        dmdte(jee)=-((Cde(jee)*Aeffe(jee)*p_cyle(jee)/sqrt(Reqe(jee)*T_cyle(jee)))*((pe/p_cyle(jee))^(1/k_cyl(jee)))*sqrt(2*k_cyl(jee)/(k_cyl(jee)-1))*sqrt(1-(pe/p_cyle(jee))^((k_cyl(jee)-1)/k_cyl(jee))))*n3;
        
    end
    
    else  % p_cyl is less than p_exh
  
       if (p_cyle(jee)/pe)<=((2/(k_cyl(jee)+1))^(k_cyl(jee)/(k_cyl(jee)-1)))
           
        dmdte(jee)=((Cde(jee)*Aeffe(jee)*pe/sqrt(Reqe(jee)*Te))*sqrt(k_cyl(jee))*((2/(k_cyl(jee)+1))^((k_cyl(jee)+1)/(2*(k_cyl(jee)-1)))))*n3;
        
       else
           
        dmdte(jee)=((Cde(jee)*Aeffe(jee)*pe/sqrt(Reqe(jee)*Te))*((p_cyle(jee)/pe)^(1/k_cyl(jee)))*sqrt(2*k_cyl(jee)/(k_cyl(jee)-1))*sqrt(1-(p_cyle(jee)/pe)^((k_cyl(jee)-1)/k_cyl(jee))))*n3;
           
       end
       
end

      if p_cyle(jee)>pe
          
      dtdte(jee)=(dQdthetae(jee) - p_cyle(jee)*dVdt(jee)  - dmdte(jee)*T_cyle(jee)*(cveqe(jee)-cpeqe(jee)))/(me(jee)*cveqe(jee)) ;
      dpdte(jee)=p_cyle(jee)*((-1/Ve(jee))*dVdt(jee) + (1/me(jee))*dmdte(jee) + (1/T_cyle(jee))*dtdte(jee));
     
      else         
     
      dtdte(jee)=(dQdthetae(jee) - p_cyle(jee)*dVdt(jee) - dmdte(jee)*(T_cyle(jee)*(cpeqe(jee)-Reqe(jee))- cpeqe(jee)*Te))/(me(jee)*cveqe(jee));
      dpdte(jee)=p_cyle(jee)*((-1/Ve(jee))*dVdt(jee) + (1/me(jee))*dmdte(jee) + (1/T_cyle(jee))*dtdte(jee));
     
      end   
   
      T_cyle(jee+1)= T_cyle(jee) + dtdte(jee)*deltathetap;
      p_cyle(jee+1)= p_cyle(jee) + dpdte(jee)*deltathetap;
      me(jee+1)=me(jee) + dmdte(jee)*deltathetap;

      dWdthetae(jee+1)= p_cyle(jee)*dVdt(jee);
     
      Vpss=(1/Aeffe(jee))*dVdt(jee);
      Vps(jee)=Vpss;
      jee=jee+1;

end

% figure;plot(thetae(1:(nx4+1)),p_cyle(1:(nx4+1)),'blue','linewidth',1.6);title('PRESSÃO APÓS AVE'),grid on,hold on
% figure;plot(thetae(1:(nx4+1)),T_cyle(1:(nx4+1)),'red','linewidth',1.6);title('TEMPERATURA APÓS AVE'),grid on,hold on
% figure;plot(thetae(1:(nx4+1)),-dmdte(1:(nx4+1)),'blue','linewidth',1.6);title('DMDT'),grid on,hold on
% figure;plot(thetae(1:(nx4+1)),me(1:(nx4+1)),'red','linewidth',1.6);title('VARIAÇÃO DA MASSA NO CILINDRO'),grid on,hold on

while ji<=nx5+1
Vi(ji)=Vd/(r-1) + Vd/2*(RR +1-cos(thetai(ji)*pi/180)-((RR^2-(sin(thetai(ji)*pi/180).^2)).^0.5));
dVdthetai(ji)=(Vtdc*(r-1)/2*(sin(thetai(ji)*pi/180)+ eps/2*sin(2*thetai(ji)*pi/180)./sqrt(1-eps^2*sin(thetai(ji)*pi/180).^2)))*pi()/180;
yi(ji)=l+a-(((l^2-a^2*(sin(thetai(ji)*pi/180).^2)).^0.5)+a*cos(thetai(ji)*pi/180));
Awi(ji)=pi*b*(yi(ji))+pi/2*b^2;
thetai(ji+1)=thetai(ji)+deltathetap;
ji=ji+1;
end

% Intake period
mi(1)=me(nx4+1);
T_cyla(1)=T_cyle(nx4+1);
p_cyla(1)=p_cyle(nx4+1);

while jei<=nx5+1

[MW, Rmix, hmix, umix, vmix, smix, cpmix, cvmix, Rp, hp, up, vp, sp, cpp, cvp]=propres(p_cyla(jei),T_cyla(jei),phi,f);
 
Reqei(jei)=Rp;
cvqei(jei)=cvp;
cpqei(jei)=cpp;
heqei(jei)=hp;
kcili(jei)=cpqei(jei)/cvqei(jei); % k gases inside the cylinder
Rar= 0.2879*10^3; % SHAPIRO MORAN
har= 300.19*10^3; 
uar= 214.07*10^3; 
cpar= 1.005*10^3; 
cvar= 0.718*10^3; 
kint(jei)=cpar/cvar; % Kint gases inside the intake manifold
phoar=1.184; % 
var=1/phoar  ; % SHAPIRO MORAN

% Hohenberg correlation for intake phase
khohi=100;
a1i=1.4;
hhei(jei)=(khohi*(Vi(jei)^(-0.06))*((p_cyla(jei)/100000)^(0.8))*(T_cyla(jei)^(-0.4))*(upmean+a1i)^(0.8));
 
% Nishiwaki correlation for intake phase
hhi(jei)=(584*(b^(-0.193))*((upmean*((p_cyla(jei)/100000)/0.9807))^(0.807))*(T_cyla(jei))^(-0.534))*1.161111; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dQdthetawi(jei)=(hhi(jei)*Awi(jei)*(T_cyla(jei)-Tw))*n3;
dQdthetaa(jei)= - dQdthetawi(jei);
dVdti(jei)=dVdthetai(jei);

if padm>p_cyla(jei)
     
    if (p_cyla(jei)/padm)<=((2/(kint(jei)+1))^(kint(jei)/(kint(jei)-1))) % ESCOAMENTO SÔNICO / SUPERSONICO
    
     dmdti(jei)=((Cdi(jei)*Aeffi(jei)*padm/sqrt(Rar*Tadm))*sqrt(kint(jei))*((2/(kint(jei)+1))^((kint(jei)+1)/(2*(kint(jei)-1)))))*n3;
        
    else % ESCOAMENTO NORMAL
      
      dmdti(jei)=((Cdi(jei)*Aeffi(jei)*padm/sqrt(Rar*Tadm))*((p_cyla(jei)/padm)^(1/kint(jei)))*sqrt(2*kint(jei)/(kint(jei)-1))*sqrt(abs(1-(p_cyla(jei)/padm)^((kint(jei)-1)/kint(jei)))))*n3;
        
    end
    
    else  % p_cyl é MAIOR que padm
            
       if (padm/p_cyla(jei))<=((2/(kcili(jei)+1))^(kcili(jei)/(kcili(jei)-1)))
           
       dmdti(jei)=-((Cdi(jei)*Aeffi(jei)*p_cyla(jei)/sqrt(Reqei(jei)*T_cyla(jei)))*sqrt(kcili(jei))*((2/(kcili(jei)+1))^((kcili(jei)+1)/(2*(kcili(jei)-1)))))*n3;
        
       else
           
           dmdti(jei)=-((Cdi(jei)*Aeffi(jei)*p_cyla(jei)/sqrt(Reqei(jei)*T_cyla(jei)))*((padm/p_cyla(jei))^(1/kcili(jei)))*sqrt(2*kcili(jei)/(kcili(jei)-1))*sqrt(abs(1-(padm/p_cyla(jei))^((kcili(jei)-1)/kcili(jei)))))*n3;
         
       end
       
end
  
    if padm>p_cyla(jei)
        
        dtdta(jei)=(dQdthetaa(jei) - p_cyla(jei)*dVdti(jei) + cpar*Tadm*dmdti(jei) - cvqei(jei)*dmdti(jei)* T_cyla(jei))/(mi(jei)*cvqei(jei));
        dpdta(jei)=p_cyla(jei)*((-1/Vi(jei))*dVdti(jei) + (1/mi(jei))*dmdti(jei) + (1/T_cyla(jei))*dtdta(jei));     
                
         else          
        
        dtdta(jei)=(dQdthetaa(jei) - p_cyla(jei)*dVdti(jei) + cpqei(jei)*T_cyla(jei)*dmdti(jei) - cvqei(jei)*T_cyla(jei)*dmdti(jei) )/(mi(jei)*cvqei(jei));
        dpdta(jei)=p_cyla(jei)*((-1/Vi(jei))*dVdti(jei)+(1/mi(jei))*dmdti(jei) + (1/T_cyla(jei))*dtdta(jei));     
        
    end

        T_cyla(jei+1)= T_cyla(jei) + dtdta(jei)*deltathetap;
        p_cyla(jei+1)= p_cyla(jei) + dpdta(jei)*deltathetap;
        mi(jei+1)=  mi(jei) + dmdti(jei)*deltathetap;
     
        dWdthetaa(jei+1)= p_cyla(jei)*dVdti(jei);
 
        jei=jei+1;
        jej=jej+1;
        
end

% Otherwise Compute the Overlapping phase %

else

while je<=nx4+1
Ve(je)=Vd/(r-1) + Vd/2*(RR +1-cos(thetae(je)*pi/180)-((RR^2-(sin(thetae(je)*pi/180).^2)).^0.5));
dVdthetae(je)=(Vtdc*(r-1)/2*(sin(thetae(je)*pi/180)+ eps/2*sin(2*thetae(je)*pi/180)./sqrt(1-eps^2*sin(thetae(je)*pi/180).^2)))*pi()/180;
ye(je)=l+a-(((l^2-a^2*(sin(thetae(je)*pi/180).^2)).^0.5)+a*cos(thetae(je)*pi/180));
Awe(je)=pi*b*(ye(je))+pi/2*b^2;
thetae(je+1)=thetae(je)+deltathetap;
je=je+1;
end

% Pressure variation during the exhaust period

while jee<=nx4+1
    
[MW, Rmix, hmix, umix, vmix, smix, cpmix, cvmix, Rp, hp, up, vp, sp, cpp, cvp]=propres(p_cyle(jee),T_cyle(jee),phi,f);

k_exh(jee)=1.28; % Ke dos gases no coletor de escape;
cveqe(jee)=cvp;
cpeqe(jee)=cpp;
ueqe(jee)=up;
k_cyl(jee)=cpp/cvp;
heqe(jee)=hp;
Reqe(jee)= Rp;
pho(jee)=1/vp;
co(jee)=sqrt(k_cyl(jee)*Reqe(jee)*T_cyle(jee));

% Hohenberg
khohe=130;
a1e=1.4;
hhee(jee)=(khohe*(Ve(jee)^(-0.06))*((p_cyle(jee)/100000)^(0.8))*(T_cyle(jee)^(-0.4))*(upmean+a1e)^(0.8));

%Nishiwaki

hhe(jee)=(45.2*(b^(-0.422))*((upmean*((p_cyle(jee)/100000)/0.9807))^(0.578))*(T_cyle(jee))^(-0.131))*1.161111; % TRANSFORMAÇÃO DE KCAL/H PARA J/S

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dQdthetawe(jee)=(hhe(jee)*Awe(jee)*(T_cyle(jee)-Tw))*n3;
dQdthetae(jee)= - dQdthetawe(jee);
dVdt(jee)=dVdthetae(jee);
Cdb(jee)=-1.5*YDe(jee)+ 1;

% Integration
 
if p_cyle(jee)>pe
    
     if (pe/p_cyle(jee))<=((2/(k_cyl(jee)+1))^(k_cyl(jee)/(k_cyl(jee)-1))) % ESCOAMENTO SÔNICO / SUPERSONICO
    
        dmdte(jee)=-((Cde(jee)*Aeffe(jee)*p_cyle(jee)/sqrt(Reqe(jee)*T_cyle(jee)))*sqrt(k_cyl(jee))*((2/(k_cyl(jee)+1))^((k_cyl(jee)+1)/(2*(k_cyl(jee)-1)))))*n3;
        
    else % ESCOAMENTO NORMAL
      
        dmdte(jee)=-((Cde(jee)*Aeffe(jee)*p_cyle(jee)/sqrt(Reqe(jee)*T_cyle(jee)))*((pe/p_cyle(jee))^(1/k_cyl(jee)))*sqrt(2*k_cyl(jee)/(k_cyl(jee)-1))*sqrt(1-(pe/p_cyle(jee))^((k_cyl(jee)-1)/k_cyl(jee))))*n3;
        
    end
    
    else  % p_cyl é MENOR que p_exh
  
       if (p_cyle(jee)/pe)<=((2/(k_cyl(jee)+1))^(k_cyl(jee)/(k_cyl(jee)-1)))
           
        dmdte(jee)=((Cde(jee)*Aeffe(jee)*pe/sqrt(Reqe(jee)*Te))*sqrt(k_cyl(jee))*((2/(k_cyl(jee)+1))^((k_cyl(jee)+1)/(2*(k_cyl(jee)-1)))))*n3;
        
       else
           
        dmdte(jee)=((Cde(jee)*Aeffe(jee)*pe/sqrt(Reqe(jee)*Te))*((p_cyle(jee)/pe)^(1/k_cyl(jee)))*sqrt(2*k_cyl(jee)/(k_cyl(jee)-1))*sqrt(1-(p_cyle(jee)/pe)^((k_cyl(jee)-1)/k_cyl(jee))))*n3;
           
       end
       
end

      if p_cyle(jee)>pe
          
      dtdte(jee)=(dQdthetae(jee) - p_cyle(jee)*dVdt(jee)  - dmdte(jee)*T_cyle(jee)*(cveqe(jee)-cpeqe(jee)))/(me(jee)*cveqe(jee)) ;
      dpdte(jee)=p_cyle(jee)*((-1/Ve(jee))*dVdt(jee) + (1/me(jee))*dmdte(jee) + (1/T_cyle(jee))*dtdte(jee));
     
      else         
     
      dtdte(jee)=(dQdthetae(jee) - p_cyle(jee)*dVdt(jee) - dmdte(jee)*(T_cyle(jee)*(cpeqe(jee)-Reqe(jee))- cpeqe(jee)*Te))/(me(jee)*cveqe(jee));
      dpdte(jee)=p_cyle(jee)*((-1/Ve(jee))*dVdt(jee) + (1/me(jee))*dmdte(jee) + (1/T_cyle(jee))*dtdte(jee));
     
      end   
   
      T_cyle(jee+1)= T_cyle(jee) + dtdte(jee)*deltathetap;
      p_cyle(jee+1)= p_cyle(jee) + dpdte(jee)*deltathetap;
      me(jee+1)=me(jee) + dmdte(jee)*deltathetap;

      dWdthetae(jee+1)= p_cyle(jee)*dVdt(jee);
      jee=jee+1;

end

% OVERLAPPING PERIOD

mee(1)=0;
mei(1)=0;

mo(1)=me(nx4+1);
T_cylo(1)=T_cyle(nx4+1);
p_cylo(1)=p_cyle(nx4+1);
jfo=nx4+1;
Tm=(Tadm+Te)/2;

% Volume during overlapping period
while joo<=nx6+1
    
Vo(joo)=Vd/(r-1) + Vd/2*(RR +1-cos(thetao(joo)*pi/180)-((RR^2-(sin(thetao(joo)*pi/180).^2)).^0.5));
dVdthetao(joo)=(Vtdc*(r-1)/2*(sin(thetao(joo)*pi/180)+ eps/2*sin(2*thetao(joo)*pi/180)./sqrt(1-eps^2*sin(thetao(joo)*pi/180).^2)))*pi()/180;
yo(joo)=l+a-(((l^2-a^2*(sin(thetao(joo)*pi/180).^2)).^0.5)+a*cos(thetao(joo)*pi/180));
Awo(joo)=pi*b*(yo(joo))+pi/2*b^2;
thetao(joo+1)=thetao(joo)+deltathetap;
joo=joo+1;
end

while jo<=nx6+1
 
[MW, Rmox, hmox, umox, vmox, smox, cpmox, cvmox, Rp, hp, up, vp, sp, cpp, cvp]=propres(p_cylo(jo),T_cylo(jo),phi,f);
 
Reqeo(jo)=Rp;
cvqeo(jo)=cvp;
cpqeo(jo)=cpp;
heqeo(jo)=hp;
kcilo(jo)=cpqeo(jo)/cvqeo(jo); % k dos gases no interior do Cilindro;

cpe=1.205*10^3;
Rar= 0.2879*10^3; % SHAPIRO MORAN
har= 300.19*10^3; 
uar= 214.07*10^3; 
cpar= 1.005*10^3; 
cvar= 0.718*10^3; 
kint(jo)=cpar/cvar; % Kint dos gases no coletor de adimissão;
phoar=1.184; % W
var=1/phoar; % SHAPIRO MORAN

%HOHENBERG
khoho=100;
a1o=1.4;
hho(jo)=(khoho*(Vo(jo)^(-0.06))*((p_cylo(jo)/100000)^(0.8))*(T_cylo(jo)^(-0.4))*(upmean+a1o)^(0.8));
 
%NISHIWAKI
hhio(jo)=(584*(b^(-0.193))*((upmean*((p_cylo(jo)/100000)/0.9807))^(0.807))*(T_cylo(jo))^(-0.534))*1.161111; % TRANSFORMAÇÃO DE KCAL/H PARA J/S

dQdthetawo(jo)=(hho(jo)*Awo(jo)*(T_cylo(jo)-Tw))*n3;
dQdthetao(jo)= - dQdthetawo(jo);
dVdto(jo)=dVdthetao(jo); 

%Evaluating the mass flow rate derivative sign during the exhaust process
 
if padm>p_cylo(jo)

    if (p_cylo(jo)/padm)<=((2/(kint(jo)+1))^(kint(jo)/(kint(jo)-1))) % ESCOAMENTO SÔNICO / SUPERSONICO
    
     dmdtio(jo)=((Cdi(jo)*Aeffi(jo)*padm/sqrt(Rar*Tadm))*sqrt(kint(jo))*((2/(kint(jo)+1))^((kint(jo)+1)/(2*(kint(jo)-1)))))*n3;
        
    else % ESCOAMENTO NORMAL
      
      dmdtio(jo)=((Cdi(jo)*Aeffi(jo)*padm/sqrt(Rar*Tadm))*((p_cylo(jo)/padm)^(1/kint(jo)))*sqrt(2*kint(jo)/(kint(jo)-1))*sqrt(abs(1-(p_cylo(jo)/padm)^((kint(jo)-1)/kint(jo)))))*n3;
     
    end  
    
    else  % p_cyl é MAIOR que padm
        
       if (padm/p_cylo(jo))<=((2/(kcilo(jo)+1))^(kcilo(jo)/(kcilo(jo)-1)))           
      
       dmdtio(jo)=-((Cdi(jo)*Aeffi(jo)*p_cylo(jo)/sqrt(Reqeo(jo)*T_cylo(jo)))*sqrt(kcilo(jo))*((2/(kcilo(jo)+1))^((kcilo(jo)+1)/(2*(kcilo(jo)-1)))))*n3;
      
       else
           
       dmdtio(jo)=-((Cdi(jo)*Aeffi(jo)*p_cylo(jo)/sqrt(Reqeo(jo)*T_cylo(jo)))*((padm/p_cylo(jo))^(1/kcilo(jo)))*sqrt(2*kcilo(jo)/(kcilo(jo)-1))*sqrt(abs(1-(padm/p_cylo(jo))^((kcilo(jo)-1)/kcilo(jo)))))*n3;
          
       end
       
end

%Evaluating the mass flow rate derivative sign during the exhaust process

if p_cylo(jo)>pe    
   
    if (pe/p_cylo(jo))<=((2/(kcilo(jo)+1))^(kcilo(jo)/(kcilo(jo)-1))) % ESCOAMENTO SÔNICO / SUPERSONICO
    
        dmdteo(jo)=-((Cde(jfo)*Aeffe(jfo)*p_cylo(jo)/sqrt(Reqeo(jo)*T_cylo(jo)))*sqrt(kcilo(jo))*((2/(kcilo(jo)+1))^((kcilo(jo)+1)/(2*(kcilo(jo)-1)))))*n3;
        
    else % ESCOAMENTO NORMAL
      
        dmdteo(jo)=-((Cde(jfo)*Aeffe(jfo)*p_cylo(jo)/sqrt(Reqeo(jo)*T_cylo(jo)))*((pe/p_cylo(jo))^(1/kcilo(jo)))*sqrt(2*kcilo(jo)/(kcilo(jo)-1))*sqrt(1-(pe/p_cylo(jo))^((kcilo(jo)-1)/kcilo(jo))))*n3;
        
    end
    
    else  % p_cyl é MENOR que p_exh
            
        if (p_cylo(jo)/pe)<=((2/(kcilo(jo)+1))^(kcilo(jo)/(kcilo(jo)-1)))
           
        dmdteo(jo)=((Cde(jfo)*Aeffe(jfo)*pe/sqrt(Reqeo(jo)*Te))*sqrt(kcilo(jo))*((2/(kcilo(jo)+1))^((kcilo(jo)+1)/(2*(kcilo(jo)-1)))))*n3;
        
            else % ESCOAMENTO NORMAL  
                
        dmdteo(jo)=((Cde(jfo)*Aeffe(jfo)*pe/sqrt(Reqeo(jo)*Te))*((p_cylo(jo)/pe)^(1/kcilo(jo)))*sqrt(2*kcilo(jo)/(kcilo(jo)-1))*sqrt(1-(p_cylo(jo)/pe)^((kcilo(jo)-1)/kcilo(jo))))*n3;
            
         end
       
end

    if padm>p_cylo(jo) && p_cylo(jo)>pe  
        
        dmdto(jo)=dmdtio(jo)+ dmdteo(jo);         
        dtdto(jo)=(-dmdto(jo)*cvqeo(jo)* T_cylo(jo)+ dQdthetao(jo) - p_cylo(jo)*dVdto(jo) + dmdtio(jo)*cpar*Tadm + dmdteo(jo)*cpqeo(jo)*T_cylo(jo) )/(mo(jo)*cvqeo(jo));
        dpdto(jo)=p_cylo(jo)*((-1/Vo(jo))*dVdto(jo) + (1/mo(jo))*dmdto(jo) + (1/T_cylo(jo))*dtdto(jo));     
                
    end    
            
    if   padm>p_cylo(jo) && pe > p_cylo(jo)
         
        dmdto(jo)=dmdtio(jo)+ dmdteo(jo);             
        dtdto(jo)=(-dmdto(jo)*cvqeo(jo)* T_cylo(jo)+ dQdthetao(jo) - p_cylo(jo)*dVdto(jo) + dmdteo(jo)*cpe*Te + dmdtio(jo)*cpar*Tadm  )/((mo(jo))*cvqeo(jo));
        dpdto(jo)=p_cylo(jo)*((-1/Vo(jo))*dVdto(jo)+(1/(mo(jo)))*dmdto(jo) + (1/T_cylo(jo))*dtdto(jo));   
        
    end    
        
   if   p_cylo(jo)>padm && p_cylo(jo)>= pe 
         
        dmdto(jo)=dmdtio(jo) + dmdteo(jo);
        
        dtdto(jo)=(-dmdto(jo)*cvqeo(jo)* T_cylo(jo)+ dQdthetao(jo) - p_cylo(jo)*dVdto(jo) + dmdteo(jo)*cpqeo(jo)*T_cylo(jo) + dmdtio(jo)*cpqeo(jo)*T_cylo(jo)  )/(mo(jo)*cvqeo(jo));
        
        dpdto(jo)=p_cylo(jo)*((-1/Vo(jo))*dVdto(jo)+(1/mo(jo))*dmdto(jo) + (1/T_cylo(jo))*dtdto(jo)); 
        
           end             
            
   if   p_cylo(jo)>padm &&  pe>p_cylo(jo)
            
        dmdto(jo)=dmdteo(jo)+ dmdtio(jo);
        
        dtdto(jo)=(-dmdto(jo)*cvqeo(jo)*T_cylo(jo) + dQdthetao(jo) - p_cylo(jo)*dVdto(jo) + dmdtio(jo)*cpqeo(jo)*T_cylo(jo) + dmdteo(jo)*cpe*Te)/((mo(jo))*cvqeo(jo));
        
        dpdto(jo)=p_cylo(jo)*((-1/Vo(jo))*dVdto(jo)+(1/(mo(jo)))*dmdto(jo) + (1/T_cylo(jo))*dtdto(jo));  
        
   end
             
        mee(jo+1)= mee(jo)- dmdteo(jo)*deltathetap;
        mei(jo+1)= mei(jo)- dmdtio(jo)*deltathetap;
                       
        T_cylo(jo+1)= T_cylo(jo) + dtdto(jo)*deltathetap;
        p_cylo(jo+1)= p_cylo(jo) + dpdto(jo)*deltathetap;
        
        mo(jo+1)=  mo(jo) + dmdto(jo)*deltathetap;
                       
        dWdthetaa(jo+1)= p_cylo(jo)*dVdto(jo); 
               
        jo=jo+1;
        jfo=jfo+1;
            
end

% Intake period

jei=1;
   
while ji<=nx5+1
    
Vi(ji)=Vd/(r-1) + Vd/2*(RR +1-cos(thetai(ji)*pi/180)-((RR^2-(sin(thetai(ji)*pi/180).^2)).^0.5));
dVdthetai(ji)=(Vtdc*(r-1)/2*(sin(thetai(ji)*pi/180)+ eps/2*sin(2*thetai(ji)*pi/180)./sqrt(1-eps^2*sin(thetai(ji)*pi/180).^2)))*pi()/180;
yi(ji)=l+a-(((l^2-a^2*(sin(thetai(ji)*pi/180).^2)).^0.5)+a*cos(thetai(ji)*pi/180));
Awi(ji)=pi*b*(yi(ji))+pi/2*b^2;
thetai(ji+1)=thetai(ji)+deltathetap;
ji=ji+1;
end
    
mi(1)=mo(nx6+1);
T_cyla(1)=T_cylo(nx6+1);
p_cyla(1)=p_cylo(nx6+1);
jfi=nx6+1;

while jei<=nx5+1

[MW, Rmix, hmix, umix, vmix, smix, cpmix, cvmix, Rp, hp, up, vp, sp, cpp, cvp]=propres(p_cyla(jei),T_cyla(jei),phi,f);
 
Reqei(jei)=Rp;
cvqei(jei)=cvp;
cpqei(jei)=cpp;
heqei(jei)=hp;
kcili(jei)=cpqei(jei)/cvqei(jei); % k dos gases no interior do Cilindro;

Rar= 0.2879*10^3; % SHAPIRO MORAN
har= 300.19*10^3; 
uar= 214.07*10^3; 
cpar= 1.005*10^3; 
cvar= 0.718*10^3; 
kint(jei)=cpar/cvar; % Kint dos gases no coletor de admissão;
phoar=1.184; % 
var=1/phoar  ; % SHAPIRO MORAN

% HOHENBERG

khohi=130;
a1i=1.4;
hhei(jei)=(khohi*(Vi(jei)^(-0.06))*((p_cyla(jei)/100000)^(0.8))*(T_cyla(jei)^(-0.4))*(upmean+a1i)^(0.8));
 
%NISHIWAKI
 
hhi(jei)=(584*(b^(-0.193))*((upmean*((p_cyla(jei)/100000)/0.9807))^(0.807))*(T_cyla(jei))^(-0.534))*1.161111; % TRANSFORMAÇÃO DE KCAL/H PARA J/S
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
dQdthetawi(jei)=(hhei(jei)*Awi(jei)*(T_cyla(jei)-Tw))*n3;
dQdthetaa(jei)= - dQdthetawi(jei);
dVdti(jei)=dVdthetai(jei);

if padm>p_cyla(jei)
      
    if (p_cyla(jei)/padm)<=((2/(kint(jei)+1))^(kint(jei)/(kint(jei)-1))) % ESCOAMENTO SÔNICO / SUPERSONICO
%       if (p_cyla(jei)/padm)<=((2/(kcili(jei)+1))^(kcili(jei)/(kcili(jei)-1))) % ESCOAMENTO SÔNICO / SUPERSONICO
  
     dmdti(jei)=((Cdi(jfi)*Aeffi(jfi)*padm/sqrt(Rar*Tadm))*sqrt(kint(jei))*((2/(kint(jei)+1))^((kint(jei)+1)/(2*(kint(jei)-1)))))*n3;
        
        else % ESCOAMENTO NORMAL
      
     dmdti(jei)=((Cdi(jfi)*Aeffi(jfi)*padm/sqrt(Rar*Tadm))*((p_cyla(jei)/padm)^(1/kint(jei)))*sqrt(2*kint(jei)/(kint(jei)-1))*sqrt(abs(1-(p_cyla(jei)/padm)^((kint(jei)-1)/kint(jei)))))*n3;
        
%         dmdti(jei)=((Cdi(jfi)*Aeffi(jfi)*padm/sqrt(Rar*Tadm))*((p_cyla(jei)/padm)^(1/kcili(jei)))*sqrt(2*kcili(jei)/(kcili(jei)-1))*sqrt(abs(1-(p_cyla(jei)/padm)^((kcili(jei)-1)/kcili(jei)))))*n3;
     
    end
    
        else  % p_cyl é MAIOR que padm
   
            
                if (padm/p_cyla(jei))<=((2/(kcili(jei)+1))^(kcili(jei)/(kcili(jei)-1)))
           
       dmdti(jei)=-((Cdi(jfi)*Aeffi(jfi)*p_cyla(jei)/sqrt(Reqei(jei)*T_cyla(jei)))*sqrt(kcili(jei))*((2/(kcili(jei)+1))^((kcili(jei)+1)/(2*(kcili(jei)-1)))))*n3;
        
                 else
         
       dmdti(jei)=-((Cdi(jfi)*Aeffi(jfi)*p_cyla(jei)/sqrt(Reqei(jei)*T_cyla(jei)))*((padm/p_cyla(jei))^(1/kcili(jei)))*sqrt(2*kcili(jei)/(kcili(jei)-1))*sqrt(abs(1-(padm/p_cyla(jei))^((kcili(jei)-1)/kcili(jei)))))*n3;
           
       end
       
end  
             if padm>p_cyla(jei)
        
        dtdta(jei)=(- cvqei(jei)*dmdti(jei)* T_cyla(jei)+ dQdthetaa(jei) - p_cyla(jei)*dVdti(jei) + cpar*Tadm*dmdti(jei) )/(mi(jei)*cvqei(jei));
        dpdta(jei)=p_cyla(jei)*((-1/Vi(jei))*dVdti(jei) + (1/mi(jei))*dmdti(jei) + (1/T_cyla(jei))*dtdta(jei));     
                
                 else          
        
        dtdta(jei)=(dQdthetaa(jei) - p_cyla(jei)*dVdti(jei) + cpqei(jei)*T_cyla(jei)*dmdti(jei)- cvqei(jei)*T_cyla(jei)*dmdti(jei)  )/(mi(jei)*cvqei(jei));
        dpdta(jei)=p_cyla(jei)*((-1/Vi(jei))*dVdti(jei)+(1/mi(jei))*dmdti(jei) + (1/T_cyla(jei))*dtdta(jei));     
        
    end

        T_cyla(jei+1)= T_cyla(jei) + dtdta(jei)*deltathetap;
        p_cyla(jei+1)= p_cyla(jei) + dpdta(jei)*deltathetap;
        mi(jei+1)=  mi(jei) + dmdti(jei)*deltathetap;
     
        dWdthetaa(jei+1)= p_cyla(jei)*dVdti(jei);
 
        jei=jei+1;
        jfi=jfi+1;        
end

end

if thetaava ==thetafve 
    
ddd=find(mi==min(mi)); 
ff(ir)=mi(1)/mi(nx4+1); 
pff(ir)=p_cyla(nx5+1); 
Tff(ir)=T_cyla(nx5+1);
Vf=Vi(nx5+1);   
mff=mi(nx5+1); 

IntervaloCruz= thetafve-thetaava;
    
else
    
% ff(ir)=mo(nx6+1)/mi(nx5+1) ;   % mo(1) >> É a massa no cilindro no instante em que vai Abrir valvula de adimissão AVA, ou seja, vai começar o Overlap
ff(ir)=(mo(nx6+1)+mei(nx6+1))/(mi(nx5+1))
pff(ir)=p_cyla(nx5+1); 
Tff(ir)=T_cyla(nx5+1); 
IntervaloCruz= thetafve-thetaava;

end

ER1= (abs(ff(ir)-f))/ff(ir);
ER2= (abs(pff(ir)-p1))/pff(ir);
ER3= (abs(Tff(ir)-T1))/Tff(ir);

f=ff(ir);
p1=pff(ir);
T1=Tff(ir);

ir=ir+1;

end

dQdthetae=dQdthetae*(180/pi());
dQdthetao=dQdthetao*(180/pi());
dQdthetaa=dQdthetaa*(180/pi());
dWdtheta=dWdtheta*(180/pi());

% Engine Performance parameters

% QC
WC=U1-U2 % Compression work for polytropic process
% WC1=(p2*V2-p1*V1)/(1-nep)
QC=trapz(thetac(1:nx+1),dQdthetaa)
% QE
WE=U4-U5 % Expansion work for polytropic process
% WE3=(p5*V5-p4*V4)/(1-ne)
neta_v=(mi(nx5+1)-mei(nx6+1))/(phoar*Vd)
% neta_i = Wt/(mf*PCI)

% PLOT OF GRAPHICS%

if thetaava ==thetafve 
    
% figure;plot(thetai(1:(nx5+1)),p_cyla(1:(nx5+1)),'blue','linewidth',1.6);title('PRESSÃO APÓS AVA'),xlabel('\theta (°CA)'),ylabel('xxxx (xxx)'),grid on,hold on
% figure;plot(thetai(1:(nx5+1)),T_cyla(1:(nx5+1)),'green','linewidth',1.6);title('TEMPERATURA APÓS AVA'),xlabel('\theta (°CA)'),ylabel('xxxx (xxx)'),grid on,hold on
% figure;plot(thetai(1:(nx5+1)),dmdti(1:(nx5+1)),'blue','linewidth',1.6);title('DMDT'),grid on,hold on
% figure;plot(thetai(1:(nx5+1)),mi(1:(nx5+1)),'blue','linewidth',1.6);title('VARIAÇÃO DA MASSA NO CILINDRO'),xlabel('\theta (°CA)'),ylabel('xxxx (xxx)'),grid on,hold on
% 
% figure;plot(thetae(1:(nx4+1)),p_cyle(1:(nx4+1)),'blue','linewidth',1.6);title('PRESSÃO APÓS AVE'),xlabel('\theta (°CA)'),ylabel('xxxx (xxx)'),grid on,hold on
% figure;plot(thetae(1:(nx4+1)),T_cyle(1:(nx4+1)),'red','linewidth',1.6);title('TEMPERATURA APÓS AVE'),xlabel('\theta (°CA)'),ylabel('xxxx (xxx)'),grid on,hold on
% figure;plot(thetae(1:(nx4+1)),dmdte(1:(nx4+1)),'blue','linewidth',1.6);title('DMDT'),xlabel('\theta (°CA)'),ylabel('xxxx (xxx)'),grid on,hold on
% figure;plot(thetae(1:(nx4+1)),me(1:(nx4+1)),'red','linewidth',1.6);title('VARIAÇÃO DA MASSA NO CILINDRO'),xlabel('\theta (°CA)'),ylabel('xxxx (xxx)'),grid on,hold on
% 
% figure;plot(thetae(1:(nx4+1)),dmdte(1:(nx4+1)),'red',thetai(1:(nx6+1)),dmdti(1:(nx6+1)),'red','linewidth',1.8);title('TAXA COMPLETA DA VAZÃO DE MASSA NO ESCAPE '),xlabel('\theta (°CA)'),ylabel('xxxx (xxx)'),grid on,hold on

else
    
% figure;plot(mee,'red','linewidth',1.8);title('MASSA ACUMULADA NO COLETOR DE ESCAPE ( Backflow)'),xlabel('\theta (°CA)'),ylabel('Massa (Kg)'),grid on,hold on
 figure;plot(mei,'blue','linewidth',1.8);title('MASSA ACUMULADA NO COLETOR DE ADMISSÃO ( Backflow)'),xlabel('\theta (°CA)'),ylabel('Massa (Mg)'),grid on,hold on
% 
% figure;plot(thetae(1:(nx4+1)),p_cyle(1:(nx4+1)),'red','linewidth',1.8);title('PRESSÃO APÓS AVE ATÉ A AVA (INICIO DO OVERLLAPING)'),xlabel('\theta (°CA)'),ylabel('Pressão (Pa)'),grid on,hold on
% figure;plot(thetae(1:(nx4+1)),T_cyle(1:(nx4+1)),'red','linewidth',1.8);title('TEMPERATURA APÓS AVE ATÉ A AVA (INICIO DO OVERLLAPING)'),xlabel('\theta (°CA)'),ylabel('Temperatura (K)'),grid on,hold on
% figure;plot(thetae(1:(nx4+1)),dmdte(1:(nx4+1)),'red','linewidth',1.8);title('TAXA MASSICA NO CILINDRO ATÉ A AVA (INICIO DO OVERLLAPING)'),xlabel('\theta (°CA)'),ylabel('Fluxo de Massa (Kg/°CA)'),grid on,hold on
 figure;plot(thetae(1:(nx4+1)),me(1:(nx4+1)),'red','linewidth',1.8);title('VARIAÇÃO DA MASSA NO CILINDRO ATÉ A AVA (INICIO DO OVERLLAPING)'),xlabel('\theta (°CA)'),ylabel('Massa (Kg)'),grid on,hold on
% 
% figure;plot(thetai(1:(nx5+1)),p_cyla(1:(nx5+1)),'blue','linewidth',1.8);title('PRESSÃO NO CILINDRO APÓS TÉRMINO DO OVERLAPPING'),xlabel('\theta (°CA)'),ylabel('Pressão (Pa)'),grid on,hold on
% figure;plot(thetai(1:(nx5+1)),T_cyla(1:(nx5+1)),'green','linewidth',1.8);title('TEMPERATURA NO CILINDRO APÓS TÉRMINO DO OVERLAPPING'),xlabel('\theta (°CA)'),ylabel('Temperatura (K)'),grid on,hold on
 figure;plot(thetai(1:(nx5+1)),dmdti(1:(nx5+1)),'blue','linewidth',1.8);title('VAZÃO MÁSSICA DE AR APÓS TÉRMINO DO OVERLAPPING'),xlabel('\theta (°CA)'),ylabel('Fluxo de Massa (Kg/°CA)'),grid on,hold on
% figure;plot(thetai(1:(nx5+1)),mi(1:(nx5+1)),'blue','linewidth',1.8);title('VARIAÇÃO DA MASSA NO CILINDRO APÓS TÉRMINO DO OVERLAPPING'),xlabel('\theta (°CA)'),ylabel('Massa (Kg)'),grid on,hold on
% 
% figure;plot(thetao(1:(nx6+1)),p_cylo(1:(nx6+1)),'black','linewidth',1.8);title('PRESSÃO GASES DURANTE O OVERLLAPING'),xlabel('\theta (°CA)'),ylabel('Pressão (Pa)'),grid on,hold on
% figure;plot(thetao(1:(nx6+1)),T_cylo(1:(nx6+1)),'black','linewidth',1.8);title('TEMPERATURA GASES DURANTE O OVERLLAPING'),xlabel('\theta (°CA)'),ylabel('Temperatura (K)'),grid on,hold on
 figure;plot(thetao(1:(nx6+1)),mo(1:(nx6+1)),'black','linewidth',1.8);title('VARIAÇÃO DA MASSA DOS GASES NO CILINDRO DURANTE O OVERLLAPING'),xlabel('\theta (°CA)'),ylabel('Massa (Kg)'),grid on,hold on

thetacp=thetacp-720;
thetai=thetai-720;
% figure;plot(thetae(1:(nx4+1)),p_cyle(1:(nx4+1)),'green',thetao(1:(nx6+1)),p_cylo(1:(nx6+1)),'green',thetai(1:(nx5+1)),p_cyla(1:(nx5+1)),'green','linewidth',1.8);title('VARIAÇÃO DA PRESSÃO DO GASES NO CILINDRO DURANTE TODA A FASE ABERTA DO CICLO'),xlabel('\theta (°CA)'),ylabel('Pressão (Pa)'),grid on,hold on
 figure;plot(thetai(1:(nx5+1)),T_cyla(1:(nx5+1)),'red',thetacp(1:(z1+1)),Tf(1:(z1+1)),'red',thetac(1:nx+1),Tc(1:nx+1),'red',thetaep(1:(hx+1)),Tfe(1:(hx+1)),'red',thetae(1:(nx4+1)),T_cyle(1:(nx4+1)),'red',thetao(1:(nx6+1)),T_cylo(1:(nx6+1)),'red','linewidth',1.8);title('VARIAÇÃO DA TEMPERATURA DO GASES NO CILINDRO DURANTE TODA A FASE ABERTA DO CICLO'),xlabel('\theta (°CA)'),ylabel('Temperatura (K)'),grid on,hold on
 
 figure;plot(thetai(1:(nx5+1)),p_cyla(1:(nx5+1)),'red',thetacp(1:(z1+1)),pf(1:(z1+1)),'red',thetac(1:nx+1),pc(1:nx+1),'red',thetaep(1:(hx+1)),pfe(1:(hx+1)),'red',thetae(1:(nx4+1)),p_cyle(1:(nx4+1)),'red',thetao(1:(nx6+1)),p_cylo(1:(nx6+1)),'red','linewidth',1.8);title('VARIAÇÃO DA TEMPERATURA DO GASES NO CILINDRO DURANTE TODA A FASE ABERTA DO CICLO'),xlabel('\theta (°CA)'),ylabel('Temperatura (K)'),grid on,hold on
 
 figure;plot(Vi(1:(nx5+1)),p_cyla(1:(nx5+1)),'red',Vfcp(1:(z1+1)),pf(1:(z1+1)),'red',Vc(1:nx+1),pc(1:nx+1),'red',Vep(1:(hx+1)),pfe(1:(hx+1)),'red',Ve(1:(nx4+1)),p_cyle(1:(nx4+1)),'red',Vo(1:(nx6+1)),p_cylo(1:(nx6+1)),'red','linewidth',1.2);title('VARIAÇÃO DA PRESSÃO DO GASES NO CILINDRO DURANTE TODA A FASE ABERTA DO CICLO'),xlabel('\theta (°CA)'),ylabel('Pressão (K)'),grid on,hold on

 figure;plot(Vi(1:(nx5+1)),p_cyla(1:(nx5+1)),'red',Vfcp(1:(z1+1)),pf(1:(z1+1)),'red',Vep(1:(hx+1)),pfe(1:(hx+1)),'red',Ve(1:(nx4+1)),p_cyle(1:(nx4+1)),'red',Vo(1:(nx6+1)),p_cylo(1:(nx6+1)),'red','linewidth',1.2);title('VARIAÇÃO DA PRESSÃO DO GASES NO CILINDRO DURANTE TODA A FASE ABERTA DO CICLO'),xlabel('\theta (°CA)'),ylabel('Pressão (K)'),grid on,hold on
 figure;plot(log(Vi(1:(nx5+1))),log(p_cyla(1:(nx5+1))),'red',log(Vfcp(1:(z1+1))),log(pf(1:(z1+1))),'red',log(Vc(1:nx+1)),log(pc(1:nx+1)),'red',log(Vep(1:(hx+1))),log(pfe(1:(hx+1))),'red',log(Ve(1:(nx4+1))),log(p_cyle(1:(nx4+1))),'red',log(Vo(1:(nx6+1))),log(p_cylo(1:(nx6+1))),'red','linewidth',1.2);title('VARIAÇÃO DA PRESSÃO DO GASES NO CILINDRO DURANTE TODA A FASE ABERTA DO CICLO'),xlabel('\theta (°CA)'),ylabel('Pressão (K)'),grid on,hold on
 figure;plot((Vi(1:(nx5+1))),log(p_cyla(1:(nx5+1))),(Vfcp(1:(z1+1))),log(pf(1:(z1+1))),'black',(Ve(1:(nx4+1))),log(p_cyle(1:(nx4+1))),'black',(Vo(1:(nx6+1))),log(p_cylo(1:(nx6+1))),'black','linewidth',1.6);title('VARIAÇÃO DA PRESSÃO DO GASES NO CILINDRO DURANTE TODA A FASE ABERTA DO CICLO'),xlabel('\theta (°CA)'),ylabel('Pressão (K)'),grid on,hold on

 figure;plot(Vi(1:(nx5+1)),p_cyla(1:(nx5+1)),'red',Vfcp(1:(z1+1)),pf(1:(z1+1)),'red',Vep(1:(hx+1)),pfe(1:(hx+1)),'red',Ve(1:(nx4+1)),p_cyle(1:(nx4+1)),'red',Vo(1:(nx6+1)),p_cylo(1:(nx6+1)),'red','linewidth',1.2);title('VARIAÇÃO DA PRESSÃO DO GASES NO CILINDRO DURANTE TODA A FASE ABERTA DO CICLO'),xlabel('\theta (°CA)'),ylabel('Pressão (K)'),grid on,hold on

 
 figure;plot(thetae(1:(nx4+1)),me(1:(nx4+1)),'green',thetao(1:(nx6+1)),mo(1:(nx6+1)),'green',thetai(1:(nx5+1)),mi(1:(nx5+1)),'green','linewidth',1.8);title('VARIAÇÃO DA MASSA DE GASES NO CILINDRO DURANTE TODA A FASE ABERTA DO CICLO'),xlabel('\theta (°CA)'),ylabel('Massa (Kg)'),grid on,hold on
% 
 figure;plot(thetae(1:(nx4+1)),dmdte(1:(nx4+1)),'red',thetao(1:(nx6+1)),dmdteo(1:(nx6+1)),'red','linewidth',1.8);title('TAXA COMPLETA DA VAZÃO DE MASSA NO ESCAPE '),xlabel('\theta (°CA)'),ylabel('Fluxo de Massa (Kg/°CA)'),grid on,hold on
 figure;plot(thetao(1:(nx6+1)),dmdtio(1:(nx6+1)),'blue',thetai(1:(nx5+1)),dmdti(1:(nx5+1)),'blue','linewidth',1.8);title('TAXA COMPLETA DA VAZÃO DE MASSA NA ADMISSÃO'),xlabel('\theta (°CA)'),ylabel('Fluxo de Massa (Kg/°CA)'),grid on,hold on
 figure;plot(thetae(1:(nx4+1)),dmdte(1:(nx4+1)),'red',thetao(1:(nx6+1)),dmdteo(1:(nx6+1)),'red',thetao(1:(nx6+1)),dmdtio(1:(nx6+1)),'blue',thetai(1:(nx5+1)),dmdti(1:(nx5+1)),'blue','linewidth',1.8);title('TAXA COMPLETA DOS FLUXOS DE MASSA ADM E ESCAPE DURANTE TODA A FASE ABERTA DO CICLO'),xlabel('\theta (°CA)'),ylabel('Fluxo de Massa (Kg/°CA)'),grid on,hold on
 figure;plot(thetao(1:(nx6+1)),dmdtio(1:(nx6+1)),'blue',thetao(1:(nx6+1)),dmdteo(1:(nx6+1)),'red','linewidth',1.8);title('TAXA DE MASSA ADMISSÃO E TAXA DE MASSA ESCAPE ( DURANTE OVERLAPPING)'),xlabel('\theta (°CA)'),ylabel('Fluxo de Massa (Kg/°CA)'),grid on,hold on

jii=(nx5+1);
jiii=1;

joi=(nx6+1);
joii=1;

while joi>0

thetaoi(joii)=thetao(joi);

joi=joi-1;

joii=joii+1;

end

while jii>0

thetaii(jiii)=thetai(jii);

jii=jii-1;

jiii=jiii+1;

end

% figure;plot(thetae(1:(nx4+1)),YDe(1:(nx4+1)),'red','linewidth',1.8);title('Lift Escape '),ylabel('Y/D (Adimensional)'),grid on,hold on
% figure;plot(thetaii(1:(nx5+1)),YDi(1:(nx5+1)),'blue','linewidth',1.8);title('Lift Admissão '),ylabel('Y/D (Adimensional)'),grid on,hold on
% figure;plot(thetae(1:(nx4+1)),YDe(1:(nx4+1)),'red',thetaoi(1:(nx6+1)),YDe(1:(nx6+1)),'red',thetao(1:(nx6+1)),YDi(1:(nx6+1)),'blue',thetaii(1:(nx5+1)),YDi(1:(nx5+1)),'blue','linewidth',1.8);title('Lift dos Cruzamentos de Válvula '),xlabel('\theta (°CA)'),ylabel('Y/D (Adimensional)'),grid on,hold on
% figure;plot(thetao(1:(nx6+1)),YDi(1:(nx6+1)),'blue',thetaoi(1:(nx6+1)),YDe(1:(nx6+1)),'red','linewidth',1.8);title('Lift dos Cruzamentos de Válvula '),xlabel('\theta (°CA)'),ylabel('Y/D (Adimensional)'),grid on,hold on

% figure;plot(thetac(1:(nx)),dQdthetaw(1:(nx)),'red','linewidth',1.8),grid on,hold on
% figure;plot(thetae(1:(nx4+1)),dQdthetae(1:(nx4+1)),'red','linewidth',1.8),grid on,hold on
% figure;plot(thetao(1:(nx6+1)),dQdthetao(1:(nx6+1)),'yellow','linewidth',1.8),grid on,hold on
% figure;plot(thetai(1:(nx5+1)),dQdthetaa(1:(nx5+1)),'blue','linewidth',1.8),grid on,hold on
% 
 
% 
 figure;plot(thetac(1:nx),dWdtheta(1:nx),'red','linewidth',1.8);title(' TAXA DO TRABALHO FASE FECHADA (COMBUSTÃO) '),grid on,hold on
% figure;plot(thetac(1:nx),dQdthetaw(1:nx),'red','linewidth',1.8);title('TAXA DE TROCA DE CALOR COM AS PAREDES (COMBUSTÃO) '),grid on,hold on
% figure;plot(thetac(1:nx),dQdtheta(1:nx),'red','linewidth',1.8);title('TAXA TOTAL DE TROCA DE CALOR (QCOMBUSTÃO - QPAREDE'),grid on,hold on
% 
% figure;plot(thetac(1:nx),dPdtheta(1:nx),'red','linewidth',1.8);title(' TAXA DE VARIAÇÃO DA PRESSÃO ( COMBUSTÃO)'),grid on,hold on
% figure;plot(thetac(1:nx),dTdtheta(1:nx),'red','linewidth',1.8);title(' TAXA DE VARIAÇÃO DA TEMPERATURA (COMBUSTÃO) '),grid on,hold on
% 
figure;plot(thetac(1:nx),pc(1:nx),'red','linewidth',1.8);title(' PRESSÃO '),grid on,hold on
figure;plot(thetac(1:nx),Tc(1:nx),'red','linewidth',1.8);title(' TEMPERATURA '),grid on,hold on
figure;plot(thetac(1:nx),hh(1:nx),'red','linewidth',1.8);title(' COEFICIENTE DE PELICULA '),grid on,hold on
% 

Peexp=[1.179301665
1.109761662
1.078892946
1.02095006
0.978441532
0.939791593
0.889544539
0.858547736
0.796703566
0.796468741
0.834926551
0.776855579
0.826803766
0.884554523
0.926956313
0.96929406
0.945993089
0.972811087
1.007453002
1.026617864
0.98792523
0.960765669
0.987583667
0.937229874
0.90625442
0.902246396
0.871121508
0.913480603
0.951981107
0.994404245
1.032926097
1.07144795
1.113871087
1.140753128
1.152286201
1.12898523
1.128835796
1.113273352
1.116982508
1.105278654
1.082041726
1.054860818
1.016189532
];


thetaeexp=[185.577343017
189.618455052
189.644072244
192.679709494
195.702538148
198.722164653
202.747265943
208.748093165
213.77866924
224.733220962
230.676409501
239.687256781
249.604312227
255.531490022
257.48800306
262.432121112
269.422412375
278.362812376
282.317466389
288.276665673
293.288028855
300.281522266
309.221922268
318.22636525
323.231324133
330.20560465
342.181641902
346.129891616
350.081343479
351.041988179
353.997571704
356.953155229
357.913799928
363.866594915
365.848725144
372.839016407
379.810094775
385.798113401
392.76598962
398.750806097
402.753492344
410.742854093
414.758348936
];

end
