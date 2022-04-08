
function [IDtheta,IDms,thetaic]=handerberg(RPM,upmean,T,p,Ru,NC,thetaij)
Ea=(6.18840*10^5)/(NC+25);% Energia de Ativa��o
IDtheta=(0.36 + 0.22*upmean)*exp(Ea*((1/((Ru/10^3)*T)) - (1/17190))+((21.2/((p/10^5)-12.4)).^0.63));% P e T  s�o os mesmo obtidos atrav�s de um processo politropico de compress�o, obtidos no PMS se n�o houvesse combust�o
IDms=IDtheta/(0.006*RPM);
IDtheta=round(IDtheta);
thetaic=thetaij+IDtheta; % �ngulo do �nicio da combust�o
end