
%global mf thetaic PCI

 function [x,mf1,dxdtheta,dQdthetac,nxw,nx2,thetafq]=wiebeduplateste(mf,thetaic,PCI,deltathetap)

a1=6.9078;
mp=2.5;%Fator de Forma de Câmara fase PRÉ-MISTURADA
yp=0.45; % Fração que Queima na fase de Pré-mistura
thetaps=thetaic; % INICIO DA FASE DE PRÉ-MISTURADA CONHCIDE COM INICIO COMBUSTÃO
thetapd=18; % DURAÇÃO da PRÉ-MISTURADA

md=1.5; %Fator de Forma de Câmara fase difusiva
yd=1-yp;
thetadss =thetapd+thetaps; % Inicio da fase Difusiva
thetads=thetadss;% Inicio da Fase difusiva = inicio da combustao mais a duração da fase pré-misturada 
thetadd =55; % Duração da Combustão Difusiva.

thetattc=thetapd+thetadd; % Duração do periodo total de liberação de calor
thetafq=thetattc+thetaic; % theta de termino da combustao = soma da duração da fase pré-mistura + difusiva + o angulo de inicio da combustao

nxw=(thetafq-thetaic)/deltathetap; % NUMERO TOTAL DE ITERAÇÕES FINAL COMBUSTAO - INICIO COMBUSTÃO
nx2=(thetads-thetaps)/deltathetap; % NUMERO DE ITERAÇÕES INICIO DIFUSIVA - PRÉ-MISTURADA.
%nx2=(thetads-339.4)/deltatheta; % TESTE TESTE

thetaw(1)=thetaps;
thetaw2(1)=thetads;% INICIO DA FASE DE PRÉ-MISTURADA CONHCIDE COM INICIO COMBUSTÃO
mf1(1)=mf;

j=1;
ii=1;

%%%%%%%%% WIEBE DIFUSIVA INICIANDO JUNTO COM PRÉ-MISTURADA %%%%%%%%%%%%%%

while(ii<=nxw+1)
    
 x1(ii)=1-(exp(-a1*(abs((thetaw(ii)-thetaps)./thetapd)).^(mp+1)));
 x2(ii)=1-(exp(-a1*(abs((thetaw2(ii)-thetads)./thetadd)).^(md+1))); 
 
 x(ii)=yp*x1(ii)+ yd*x2(ii);
 
 if x(ii) <0.0001
     x(ii)=0.0001;
 end
     if x(ii)>0.9999
         x(ii)=0.9999;
     end
     
 thetaw(ii+1)=thetaw(ii)+deltathetap;
 thetaw2(ii+1)=thetaw2(ii)+deltathetap;
 
 ii=ii+1;
end
  
kk=1;
% dxdtheta(1)=0;
% tb=length(x);
% x(tb)=0;

while (kk<=nxw);
    
 dxdtheta(kk)=(x(kk+1)-x(kk))/deltathetap;
 mf1(kk+1)=mf1(kk)*(1-x(kk));
 kk=kk+1;

 
end

netac=0.98; % EFICIÊCIA DA COMBUSTÃO
dQdthetac=mf*PCI*dxdtheta*netac; % TAXA DE CALOR LIBERADO FUNÇÃO DE WIEBE


% figure
% plot(thetaw(1:nx),dxdtheta(1:nx),'black'),xlabel('\theta (°CA)'),ylabel(' Taxa de Queima do Combustível(Adimensional)')
% title (' dxdtheta x \theta')
% % axis([350,5,440]);
% grid on
% hold on
% 
% figure
% plot(thetaw(1:nxw),dQdthetac(1:nxw),'black','linewidth',1.8),xlabel('\theta (°CA)'),ylabel('Taxa de Liberação de Calor (J/°CA)')
% title (' dQdtheta x \theta')
% grid on
% hold on
% 
% figure
% plot(thetaw(1:nx),x1(1:nx),'black','linewidth',1.8),xlabel('\theta (°CA)'),ylabel('Primeira Taxa de Queima (Adimensional)')
% title (' x1 x \theta')
% grid on
% hold on
% 
% figure
% plot(thetaw(1:nx),x2(1:nx),'green','linewidth',1.8),xlabel('\theta (°CA)'),ylabel('Segunda Taxa de Queima (Adimensional)')
% title (' x2 x \theta')
% grid on
% hold on
% 
% %%%%%%%%%%%%%%%%%%PLOT DUPLA DE WIEBE DIFUSIVA %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure
% plot(thetaw(1:nx),x(1:nx),'red','linewidth',1.8),xlabel('\theta (°CA)'),ylabel('Dupla Função de Wiebe (Adimensional)')
% title (' x x \theta')
% grid on
% hold on
% 
% %%%%%%%%%%%%%%%%%%TODAS JUNTAS %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure
% plot(thetaw(1:nx),x1(1:nx),'black',thetaw(1:nx),x2(1:nx),'green',thetaw(1:nx),x(1:nx),'red','linewidth',1.8),xlabel('\theta (°CA)'),ylabel(' Função dupla de Wiebe(Adimensional)')
% grid on
% hold on
%  end
 
 CalorLC=trapz(thetaw(1:nxw),dQdthetac) % Quantidade de Calor total liberado durante a Combustão
 
%   CalorLCN=splinetx(thetaw(1:nx),dQdthetac) % Quantidade de Calor total liberado durante a Combustão

%%%%%%%%%%% FRAÇÃO de QUEIMA x d(ÂNGULO DE MANIVELA) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(1)
% subplot(1,2,1)
% pos = get(gcf, 'Position');
% plot(thetaw(1:nxw),x1(1:nxw),'k',thetaw(1:nxw),x2(1:nxw),'green',thetaw(1:nxw),x(1:nxw),'red','linewidth',1.8)
% ylabel('Fração de Combustível Queimado [Adimensional]','FontSize',12)
% xlabel('Ângulo de Manivela [grau]','FontSize',12);
% legend('x1','x2','xt','Location','SouthEast')
% set(gca,'XTick',[350,360,370,380,400,430])
% grid on;
% 
% 
% subplot
% subplot(1,2,2)
% pos = get(gcf, 'Position');
% plot(thetaw(1:nxw),dxdtheta(1:nxw),'black','linewidth',1.8)
% ylabel('Taxa de Combustível Queimado [Adimensional]','FontSize',12)
% xlabel('Ângulo de Manivela [grau]','FontSize',12);
% set(gca,'XTick',[350,360,370,380,400,430])
% grid on;
% saveas(gcf,'fraçãoetaxaqueimawiebe.pdf');
