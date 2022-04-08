
%%%%%%%%%%%%%% CÁLCULO DE GEOMETRIA DE VÁLVULA - MODELO PARABÓLICO %%%%%%%%%%%%%%%%%%%%%%%%%


 function [YDe,YDi]=valvulageometria(YDmaxe,YDmaxi,rve,rvi,thetafve,thetaave,thetafva,thetaava,deltatheta)

thetave=thetafve-thetaave; % VÁLVULA DE ESCAPE
thetavi=thetafva-thetaava  % VÁLVULA DE ADMISSÃO

nne=2-2*rve;  % VÁLVULA DE ESCAPE
nni=2-2*rvi;  % VÁLVULA DE ADMISSÃO


bb1e=0;
bb1i=0;


cc1e=bb1e;
cc1i=bb1i;


aa1e=(4*YDmaxe*(1-rve))/(thetave^2);
aa1i=(4*YDmaxi*(1-rvi))/(thetavi^2);


aa3e=aa1e;
aa3i=aa1i;


bb2e=(-aa1e*thetave)/rve;
bb2i=(-aa1i*thetavi)/rvi;


bb3e=-2*aa1e*thetave;
bb3i=-2*aa1i*thetavi;


cc3e=aa1e*(thetave^2);
cc3i=aa1i*(thetavi^2);


cc2e=aa1e*(thetave^2)*(1/(4*rve*(1-rve)));
cc2i=aa1i*(thetavi^2)*(1/(4*rvi*(1-rvi)));


aa2e=aa1e/rve;
aa2i=aa1i/rvi;


nnne=(nne-1);
nnni=(nni-1);


gge=1;
ggi=1;


for theta1e=0:deltatheta:thetave
    if theta1e<=thetave/nne
Le(gge)=aa1e*theta1e^2 +bb1e*theta1e+ cc1e ;
    end
    
    if theta1e>=thetave/nne & theta1e<=nnne*(thetave)/nne
Le(gge)=aa2e*theta1e^2+bb2e*theta1e+cc2e;
    end
    
    if theta1e>=nnne*(thetave)/nne & theta1e<=thetave
Le(gge)=aa1e*theta1e^2 +bb3e*theta1e+ cc3e;
    end
    gge=gge+1;
end

te=(0:deltatheta:thetave);
YDe=Le';


for theta1i=0:deltatheta:thetavi
    if theta1i<=thetavi/nni
Li(ggi)=aa1i*theta1i^2 +bb1i*theta1i+ cc1i ;
    end
    
    if theta1i>=thetavi/nni & theta1i<=nnni*(thetavi)/nni
Li(ggi)=aa2i*theta1i^2+bb2i*theta1i+cc2i;
    end
    
    if theta1i>=nnni*(thetavi)/nni & theta1i<=thetavi
Li(ggi)=aa1i*theta1i^2 +bb3i*theta1i+ cc3i;
    end
    ggi=ggi+1;
end

ti=(0:deltatheta:thetavi);
YDi=Li';

% figure;plot(te,Le,'red','linewidth',1.8);title('ESCAPE - Elevação Modelo Parabólico (L/D) x \theta');legend('Elevação Admensional ESCAPE');grid on;hold on
% plot(ti,Li,'blue','linewidth',1.8);title('ADMISSÃO - Elevação Modelo Parabólico (L/D) x \theta');legend('Elevação Admensional ADMISSÃO');grid on;hold on

 end