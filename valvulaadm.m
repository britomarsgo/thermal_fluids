
global thetafve thetaave deltatheta


deltatheta=0.1

%%%%%%% DIAGRAMA DE ABERTURA E FECHAMENTO DE VÁLVULAS %%%%%%%%%%%%%%%%%

thetafva=210;
thetaava=0;
thetaave=510;
thetafve=720;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

YDmax=0.30;

De=38.76/1000; % DIAMETRO DA VALVULA DE ESCAPE
Di=48.88/1000; % DIAMETRO DA VALVULA DE ADMISSÃO

Lmaxe=YDmax*De; % MÁXIMO LIFT no ESCAPE
Lmaxi=YDmax*Di; % MÁXIMO LIFT na ADMISSÃO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a1=6.0;
m=3;
thetad=63;
thetad2=53;% Corresponde ao valor onde a partir do inicio de abertura a válvuala atinge máxima abertura
thetaw(1)=thetaave;
thetaw1(1)=thetaave+thetad;

nxx=(thetafve-thetaave)/deltatheta;
% YD=linspace(0,YDmax,nxx);
ii=1;



while(ii<=nxx+1)
    
 xx1(ii)=1-(exp(-a1*(((thetaw(ii)-thetaave)./thetad)).^(m+1)))
 xx2(ii)=1-(exp(-a1*(((thetaw(ii)-(thetafve)./thetad2)).^(m+1))))
 
  if xx1(ii)>=0.9999 
     xx1(ii)=1;
     xx1(ii)=xx1(ii) ; 
  end 
  
     if xx2(ii)>=xx1(ii)
        xx2(ii)=0;
     end 
         if xx1(ii)<=1*10^-9
            xx1(ii)=0;
            xx1(ii)=xx1(ii);
         end
              if xx2(ii)<=0.9999
                xx2(ii)=0;
              end  
                 if xx2(ii)>=0.9999
                    xx2(ii)=1;
                 end 
                      if xx1(ii)==1
                         xx2(ii)=1-(exp(-a1*(((thetaw(ii)-thetafve)./thetad)).^(m+1)));
                      end
                      
                      xx3(ii)= xx1(ii)+xx2(ii);
                      
 thetaw(ii+1)=thetaw(ii)+deltatheta;
 ii=ii+1;

 end
 
figure; plot(thetaw(1:nxx+1),xx1(1:nxx+1),'black','linewidth',1.8);title('Elevação xx1 x \theta'),grid on,hold on
plot(thetaw(1:nxx+1),xx2(1:nxx+1),'green','linewidth',1.8);title('Elevação xx2 x \theta'),grid on,hold on
legend('Curva de Elevação Durante Abertura','Curva de Retorno para Sede')
% 
figure; plot(thetaw(1:nxx+1),xx3(1:nxx+1),'black','linewidth',1.8);title('Elevação xx1 x \theta'),grid on,hold on
legend('Função XX1 + XX3','Curva de Retorno para Sede')

figure; plot(thetaw(1:nxx),xx11(1:nxx),'black','linewidth',1.8);title('Elevação Admensional x \theta'), grid on,hold on
plot(thetaw1(1:nxx),xx22(1:nxx),'green','linewidth',1.8);title('Elevação Admensional x \theta'),grid on,hold on
legend('Curva Admensionalisada (L/D) Abertura ','Curva Admensionalisada (L/D)Fechamento')
 
