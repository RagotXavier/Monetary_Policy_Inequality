

load To_IRFs_SP_unequal 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% IRFs Eco 1_2  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure;
%%% Plotting figures
X=0:1:(LengthIRF);

subplot(5,2,1);
plot(X,Yv(1,:),'k-',X,Yv(2,:),'b--','LineWidth',1.5);
ylim([-0.6 0.1]);
title('1. GDP (Y, %)');

subplot(5,2,2);
plot(X,Cv(1,:),'k-',X,Cv(2,:),'b--','LineWidth',1.5);
ylim([-.3 .1]);
title('2. Consumption (C, %)');


subplot(5,2,3);
plot(X,pipv(1,:),'k-',X,pipv(2,:),'b--','LineWidth',1.5);
ylim([-0.1 0.2]);
title('3. Capital (K, %)');


subplot(5,2,4);
plot(X,Lv(1,:),'k-',X,Lv(2,:),'b--','LineWidth',1.5);
ylim([-1. 1.]);
title('4. Labor (L, %)');

subplot(5,2,5);
plot(X,rnv(1,:),'k-',X,rnv(2,:),'b--','LineWidth',1.5);
ylim([-0.2 0.2]);
title('5. Nominal interest rate (R^{N}, %)');

subplot(5,2,6);
plot(X,rv(1,:),'k-',X,rv(2,:),'b--','LineWidth',1.5);
ylim([-0.05 0.05]);
title('6. Real Rate (r, %)');


subplot(5,2,7);
plot(X,Bv(1,:),'k-',X,Bv(2,:),'b--','LineWidth',1.5);
ylim([-.1 0.8]);
title('7. Debt over GDP (B/Y, %)');


subplot(5,2,8);
plot(X,Tv(1,:),'k-',X,Tv(2,:),'b--','LineWidth',1.5);
ylim([-6 6.]);
% ylim([-.1 2]);
title('8. Transfer over GDP (T/Y, %)');
        
subplot(5,2,9);
plot(X,pipv(1,:),'k-',X,pipv(2,:),'b--','LineWidth',1.5);
ylim([-0.1 0.2]);
title('9. Price Inflation (\Pi^P, %)');


subplot(5,2,10);
plot(X,piwv(1,:),'k-',X,piwv(2,:),'b--','LineWidth',1.5);
ylim([-1 1]);
title('10. Wage Inflation (\Pi^W, %)');

pos = get(gcf, 'Position');
set(gcf, 'Position',[ 616   600   800   680])

saveas(gcf,'IRFs_SP_uneq_Eco_1_2.png')





