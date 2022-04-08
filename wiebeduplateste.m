
%global mf thetaic PCI

 function [x,mf1,dxdtheta,dQdthetac,nxw,nx2,thetafq]=wiebeduplateste(mf,thetaic,PCI,deltathetap)

a1=6.9078;
mp=2.5;%Fator de Forma de C�mara fase PR�-MISTURADA
yp=0.45; % Fra��o que Queima na fase de Pr�-mistura
thetaps=thetaic; % INICIO DA FASE DE PR�-MISTURADA CONHCIDE COM INICIO COMBUST�O
thetapd=18; % DURA��O da PR�-MISTURADA

md=1.5; %Fator de Forma de C�mara fase difusiva
yd=1-yp;
thetadss =thetapd+thetaps; % Inicio da fase Difusiva
thetads=thetadss;% Inicio da Fase difusiva = inicio da combustao mais a dura��o da fase pr�-misturada 
thetadd =55; % Dura��o da Combust�o Difusiva.

thetattc=thetapd+thetadd; % Dura��o do periodo total de libera��o de calor
thetafq=thetattc+thetaic; % theta de termino da combustao = soma da dura��o da fase pr�-mistura + difusiva + o angulo de inicio da combustao

nxw=(thetafq-thetaic)/deltathetap; % NUMERO TOTAL DE ITERA��ES FINAL COMBUSTAO - INICIO COMBUST�O
nx2=(thetads-thetaps)/deltathetap; % NUMERO DE ITERA��ES INICIO DIFUSIVA - PR�-MISTURADA.
%nx2=(thetads-339.4)/deltatheta; % TESTE TESTE

thetaw(1)=thetaps;
thetaw2(1)=thetads;% INICIO DA FASE DE PR�-MISTURADA CONHCIDE COM INICIO COMBUST�O
mf1(1)=mf;

j=1;
ii=1;

%%%%%%%%% WIEBE DIFUSIVA INICIANDO JUNTO COM PR�-MISTURADA %%%%%%%%%%%%%%

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

netac=0.98; % EFICI�CIA DA COMBUST�O
dQdthetac=mf*PCI*dxdtheta*netac; % TAXA DE CALOR LIBERADO FUN��O DE WIEBE


% figure
% plot(thetaw(1:nx),dxdtheta(1:nx),'black'),xlabel('\theta (�CA)'),ylabel(' Taxa de Queima do Combust�vel(Adimensional)')
% title (' dxdtheta x \theta')
% % axis([350,5,440]);
% grid on
% hold on
% 
% figure
% plot(thetaw(1:nxw),dQdthetac(1:nxw),'black','linewidth',1.8),xlabel('\theta (�CA)'),ylabel('Taxa de Libera��o de Calor (J/�CA)')
% title (' dQdtheta x \theta')
% grid on
% hold on
% 
% figure
% plot(thetaw(1:nx),x1(1:nx),'black','linewidth',1.8),xlabel('\theta (�CA)'),ylabel('Primeira Taxa de Queima (Adimensional)')
% title (' x1 x \theta')
% grid on
% hold on
% 
% figure
% plot(thetaw(1:nx),x2(1:nx),'green','linewidth',1.8),xlabel('\theta (�CA)'),ylabel('Segunda Taxa de Queima (Adimensional)')
% title (' x2 x \theta')
% grid on
% hold on
% 
% %%%%%%%%%%%%%%%%%%PLOT DUPLA DE WIEBE DIFUSIVA %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure
% plot(thetaw(1:nx),x(1:nx),'red','linewidth',1.8),xlabel('\theta (�CA)'),ylabel('Dupla Fun��o de Wiebe (Adimensional)')
% title (' x x \theta')
% grid on
% hold on
% 
% %%%%%%%%%%%%%%%%%%TODAS JUNTAS %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure
% plot(thetaw(1:nx),x1(1:nx),'black',thetaw(1:nx),x2(1:nx),'green',thetaw(1:nx),x(1:nx),'red','linewidth',1.8),xlabel('\theta (�CA)'),ylabel(' Fun��o dupla de Wiebe(Adimensional)')
% grid on
% hold on
%  end
 
 CalorLC=trapz(thetaw(1:nxw),dQdthetac) % Quantidade de Calor total liberado durante a Combust�o
 
%   CalorLCN=splinetx(thetaw(1:nx),dQdthetac) % Quantidade de Calor total liberado durante a Combust�o

%%%%%%%%%%% FRA��O de QUEIMA x d(�NGULO DE MANIVELA) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(1)
% subplot(1,2,1)
% pos = get(gcf, 'Position');
% plot(thetaw(1:nxw),x1(1:nxw),'k',thetaw(1:nxw),x2(1:nxw),'green',thetaw(1:nxw),x(1:nxw),'red','linewidth',1.8)
% ylabel('Fra��o de Combust�vel Queimado [Adimensional]','FontSize',12)
% xlabel('�ngulo de Manivela [grau]','FontSize',12);
% legend('x1','x2','xt','Location','SouthEast')
% set(gca,'XTick',[350,360,370,380,400,430])
% grid on;
% 
% 
% subplot
% subplot(1,2,2)
% pos = get(gcf, 'Position');
% plot(thetaw(1:nxw),dxdtheta(1:nxw),'black','linewidth',1.8)
% ylabel('Taxa de Combust�vel Queimado [Adimensional]','FontSize',12)
% xlabel('�ngulo de Manivela [grau]','FontSize',12);
% set(gca,'XTick',[350,360,370,380,400,430])
% grid on;
% saveas(gcf,'fra��oetaxaqueimawiebe.pdf');
