clear 
load To_IRFs_SP_rob.mat


LengthIRF  = 20;        %length of the IRFs


X=0:1:(LengthIRF);

subplot(1,2,1);
plot(X,Piw_v(1,:),'m--',X,Piw_v(2,:),'b--',X,Piw_v(3,:),'k-',X,Piw_v(4,:),'g--',X,Piw_v(5,:),'r--','LineWidth',1.5);
%   ylim([-1.5 .1]);
title('1. Wage Inflation  (\Pi^W, %)');


subplot(1,2,2);
plot(X,Pip_v(1,:),'m--',X,Pip_v(2,:),'b--',X,Pip_v(3,:),'k-',X,Pip_v(4,:),'g--',X,Pip_v(5,:),'r--','LineWidth',1.5);
  ylim([-.05 .04]);
title('2. Price Inflation (\Pi^P, %)');
legend('0.0062', '0.019', '0.05', '0.085', '0.24', 'Location', 'southeast');

% slope_p = [0.0062, 0.019, 0.05 ,0.085, 0.24]
% slope_p = [0.0062, 0.019, 0.05 ,0.085, 0.23] ;

set(gcf, 'Position',[100, 100, 800, 400])
saveas(gcf,'../../Graph/Robustness_Slope_SP.png')