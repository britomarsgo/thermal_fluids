



function [YDn]=valvulaescap(thetaav,thetafv,deltathetap)
%%%%%%% DIAGRAMA DE ABERTURA E FECHAMENTO DE VÁLVULAS %%%%%%%%%%%%%%%%%

deltatheta=deltathetap;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

YDmax=0.25;

De=38.76/1000;  % DIÂMETRO DA VALVULA DE ESCAPE
Di=48.88/1000;  % DIÂMETRO DA VALVULA DE ADMISSÃO

Lmaxe=YDmax*De; % MÁXIMO LIFT no ESCAPE
Lmaxi=YDmax*Di; % MÁXIMO LIFT na ADMISSÃO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a1=9;
m=3;

thetad=63;  % Corresponde ao valor onde a partir do inicio de abertura a válvuala atinge máxima abertura
thetad2=63; % Corresponde ao valor onde a partir do inicio de abertura a válvuala atinge máxima abertura

thetaw(1)=thetaav;
thetaw1(1)=thetaav+thetad2;

ta=0;
nxx=(thetafv-thetaav)/deltatheta;
nxxx=((thetafv)-(thetaav-ta))/deltatheta;
ii=1;


while(ii<=nxx+1)
    
 xx1(ii)=1-(exp(-a1*(((thetaw(ii)-thetaav)./thetad)).^(m+1)));
 
 dxx1(ii)=a1*(m+1)/thetad*((thetaw(ii)-thetaav)./thetad).^m*(exp(-a1*(((thetaw(ii)-thetaav)./thetad)).^(m+1)));
 
 ddxx1(ii)=m*a1*(m+1)/(thetad^2)*((thetaw(ii)-thetaav)./thetad).^(m-1)*(exp(-a1*(((thetaw(ii)-thetaav)./thetad)).^(m+1))) - (a1^2)*((m+1)^2)/(thetad^2)*((thetaw(ii)-thetaav)./thetad).^m*((thetaw(ii)-thetaav)./thetad).^m*(exp(-a1*(((thetaw(ii)-thetaav)./thetad)).^(m+1)));
 
 
 if xx1(ii)>=0.9999
     xx1(ii)=1;
     xx1(ii)=xx1(ii);  
  end 
         if xx1(ii)<=1*10^-5
            xx1(ii)=0;
            xx1(ii)=xx1(ii);
         end
                     if xx1(ii)==1 
                       break          
                     end                     
                 
 thetaw(ii+1)=thetaw(ii)+ deltatheta;
 
 ii=ii+1;
end

xx1=xx1';
dxx1=dxx1';
ddxx1=ddxx1';

nxx2=length(xx1+1);
thetaww(1)=thetaw(nxx2);

thetadx1=thetaww(1)-thetaav;

ta=(thetafv-thetadx1);
nit=(ta -thetaww(1))/deltatheta
mm=1;

while(mm<=nit)
xx2(mm)=1;
dxx2(mm)=0;
ddxx2(mm)=0;

thetaww(mm+1)=thetaww(mm)+deltatheta;
mm=mm+1;

end

nxx3=length(xx2);
thetawww(1)=thetaww(nxx3+1);

nxxx3=(thetafv-thetawww(1))/deltatheta;
nxxx3=ceil(nxxx3);

hhh=1;

while(hhh<=nxxx3+1)
   
  xx3(hhh)=1-(exp(-a1*(((thetafv-thetawww(hhh))./thetad)).^(m+1)));
   
  dxx3(hhh)=a1*(m+1)/thetad*((-thetafv + thetawww(hhh))./thetad).^m*(exp(-a1*(((-thetafv + thetawww(hhh))./thetad)).^(m+1)));
  
  ddxx3(hhh)=m*a1*(m+1)/(thetad^2)*((-thetafv + thetawww(hhh))./thetad).^(m-1)*(exp(-a1*(((-thetafv + thetawww(hhh))./thetad)).^(m+1))) - (a1^2)*((m+1)^2)/(thetad^2)*((-thetafv + thetawww(hhh))./thetad).^m*((-thetafv + thetawww(hhh))./thetad).^m*(exp(-a1*(((-thetafv + thetawww(hhh))./thetad)).^(m+1)));

  
  if xx3(hhh)>=0.9999
     xx3(hhh)=1;
     xx3(hhh)=xx3(hhh) ;
  end 
          if xx3(hhh)<=1*10^-5
            xx3(hhh)=0;
            xx3(hhh)=xx3(hhh);
          end
         
 thetawww(hhh+1)=thetawww(hhh)+deltatheta;
 hhh=hhh+1;
end

razacel=min(ddxx3)/max(ddxx3)

xx1=xx1*YDmax;
xx2=xx2'*YDmax;
xx3=xx3'*YDmax;
dxx1=dxx1*YDmax;
dxx3=dxx3'*YDmax;
ddxx1=ddxx1*YDmax;
ddxx3=ddxx3'*YDmax;

iii=1;
jji=1;
jii=1;

aaa=(nxx2+nit+nxxx3+1)

YDn=zeros(aaa,1);

for iii=1:(nxx2)
   YDn(iii,1)=xx1(iii);
   iii=iii+1;
end

    
   
   for jji=1:nit
   YDn((nxx2+jji),1)=xx2(jji);
   jji=jji+1;
  
   end
   
   
   for jii=1:length(xx3)
    YDn((nxx2+floor(nit))+jii,1)=xx3(jii)  ;
    jii=jii+1 ;
   end
   
thetaw=thetaw';
thetaww=thetaww';
thetawww=thetawww';

figure;plot(YDn)
figure; plot(thetaw(1:nxx2),xx1(1:nxx2),'black','linewidth',1.8);title('Elevação xx1 x \theta'),grid on,hold on
plot(thetawww(1:nxxx3),xx3(1:nxxx3),'black','linewidth',1.8);title('Elevação xx1 x \theta'),grid on,hold on
plot(thetaww(1:nit),xx2(1:nit),'red','linewidth',1.8);title('Elevação xx1 x \theta'),grid on,hold on

figure; plot(thetaw(1:nxx2),dxx1(1:nxx2),'red','linewidth',1.8);title('Elevação xx2 x \theta'),grid on,hold on
plot(thetawww(1:nxxx3),dxx3(1:nxxx3),'red','linewidth',1.8);title('Elevação xx2 x \theta'),grid on,hold on
% plot(thetaww(1:nit),dxx2(1:nit),'green','linewidth',1.8);title('Elevação xx1 x \theta'),grid on,hold on

figure; plot(thetaw(1:nxx2),ddxx1(1:nxx2),'green','linewidth',1.8);title('Elevação xx2 x \theta'),grid on,hold on
plot(thetawww(1:nxxx3),ddxx3(1:nxxx3),'green','linewidth',1.8);title('Elevação xx2 x \theta'),grid on,hold on
% plot(thetaww(1:nit),ddxx2(1:nit),'black','linewidth',1.8);title('Elevação xx1 x \theta'),grid on,hold on


end
