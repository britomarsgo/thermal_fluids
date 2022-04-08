
function [IDtheta,IDms,thetaic]=handerberg(RPM,upmean,T,p,Ru,NC,thetaij)
Ea=(6.18840*10^5)/(NC+25);% Energia de Ativação
IDtheta=(0.36 + 0.22*upmean)*exp(Ea*((1/((Ru/10^3)*T)) - (1/17190))+((21.2/((p/10^5)-12.4)).^0.63));% P e T  são os mesmo obtidos através de um processo politropico de compressão, obtidos no PMS se não houvesse combustão
IDms=IDtheta/(0.006*RPM);
IDtheta=round(IDtheta);
thetaic=thetaij+IDtheta; % Ângulo do ínicio da combustão
end