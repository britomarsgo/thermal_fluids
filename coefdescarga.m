
  function [YDe Cde Cdbe YDi Cdi Cdbi]=coefdescarga(YDi,YDe)

 
%   global YDe YDi 
% 
%  deltatheta=0.1;
%   
%  
% %%%%%%% DIAGRAMA DE ABERTURA E FECHAMENTO DE VÁLVULAS %%%%%%%%%%%%%%%%%
% 
% thetafva=210;
% thetaava=0;
% thetaave=510;
% thetafve=720;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% YDmaxe=0.32;
% YDmaxi=0.32;
%  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% De=38.76/1000 % DIÂMETRO DA VÁLVULA DE ESCAPE [m]
% Di=48.88/1000 % DIÂMETRO DA VÁLVULA DE ADMISSÃO [m]
% 
% Lmaxe=YDmaxe*De % MÁXIMO LIFT no ESCAPE [m]
% Lmaxi=YDmaxi*Di % MÁXIMO LIFT na ADMISSÃO [m]




  %%%%%%%%%%%%%%%%%%%%% CÁLCULO DO COEFICIENTE DE ARRASTO - GALLO %%%%%%%%%%%%%%%%
   
  AA0=0.9999876;
  AA(1)=7.633573*10^-1;
  AA(2)=-4.089484*10^2;
  AA(3)=1.885862*10^4;
  AA(4)=-4.016319*10^5;
  AA(5)=4.720187*10^6;
  AA(6)=-3.295265*10^7;
  AA(7)=1.401485*10^8;
  AA(8)=-3.567916*10^8;
  AA(9)=5.00400*10^8;
  AA(10)=-2.977371*10^8;

 
  jjj=1; % CONTADOR
  
  n=10;  % CONTADOR PARA OS COEFICIENTES AA DO POLINÔMIO
  
  nn1e=length(YDe);
  nn1i=length(YDi);
    
  for jjj=1:nn1e
  for iii=1:n       %%%%%%%% SOMATORIO PARA CALCULO DO CD %%%%%%%%%%%%%%
 
  AA0=AA0+AA(iii)*YDe(jjj)^(iii);
  Cdbe(jjj)=-1.5*YDe(jjj)+1;
  Cde(jjj)=AA0;
  iii=iii+1;
  end
   jjj=jjj+1;
  AA0=0.9999876; % Coeficiente A0 do GALLO
  end
    
  
  for jjj=1:nn1i
  for iii=1:n       %%%%%%%% SOMATORIO PARA CALCULO DO CD %%%%%%%%%%%%%%
 
  AA0=AA0+AA(iii)*YDi(jjj)^(iii);
  Cdbi(jjj)=-1.5*YDi(jjj)+1;
  Cdi(jjj)=AA0;
  iii=iii+1;
  end
   jjj=jjj+1;
  AA0=0.9999876; % Coeficiente A0 do GALLO
  end
  
  
% figure;plot(YDe,Cdbe(1:nn1e),'r-','linewidth',1.8);title('Coeficiente de Descarga (Cd) x YD'),grid on,hold on
% plot(YDe,Cde(1:nn1e),'Green-','linewidth',1.8);title('Coeficiente de Descarga x YD'),grid on,hold on
% legend('ESCAPE - Cd - Polinomio Grau 1 Blumberg', 'Cd - Polinomio Grau 10')
%  
% figure;plot(YDi,Cdbi(1:nn1i),'k-','linewidth',1.8);title('Coeficiente de Descarga (Cd) x YD'),grid on,hold on
% plot(YDi,Cdi(1:nn1i),'Green-','linewidth',1.8);title('Coeficiente de Descarga x YD'),grid on,hold on
% legend('Cd - Polinomio Grau 1 Blumberg','ADMISSÃO - Cd - Polinomio Grau 10')
   end
