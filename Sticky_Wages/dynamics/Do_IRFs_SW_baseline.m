load To_IRFs_SW

% %%% Black line = benchmark (Optimal Monetary Policy)
% %%% Blue Dashed line = Eco 2 (Pi^w = 0)
% %%% Red Dashed line = Taylor 


calib = 'taylor'

fileName = ['To_IRFs_SW_',calib,'.mat' ];
load(fileName)
         
figure;
%%% Plotting figures
X=0:1:(LengthIRF);

subplot(5,2,1);
plot(X,Yv(1,:),'k-',X,Yv(4,:),'b--',X,Yv_taylor,'r--','LineWidth',1.5);
ylim([-.8 0.1]);
title('1. GDP (Y, %)');

subplot(5,2,2);
plot(X,Cv(1,:),'k-',X,Cv(4,:),'b--',X,Cv_taylor,'r--','LineWidth',1.5);
% ylim([-.3 .1]);
ylim([-.4 .1]);
title('2. Consumption (C, %)');

subplot(5,2,3);
plot(X,Cv(1,:),'k-',X,Cv(4,:),'b--',X,Cv_taylor,'r--','LineWidth',1.5);
% ylim([-.3 .1]);
ylim([-0.4 .1]);
title('3. Capital (K, %)');


subplot(5,2,4);
plot(X,Lv(1,:),'k-',X,Lv(4,:),'b--',X,Lv_taylor,'r--','LineWidth',1.5);
ylim([-.7 .1]);
title('4. Labor (L, %)');


subplot(5,2,5);
plot(X,RBv(1,:),'k-',X,RBv(4,:),'b--',X,RBNv_taylor,'r--','LineWidth',1.5);
ylim([-0.05 0.09]);
title('5. Nominal interest rate (R^{N}, %)');

subplot(5,2,6);
plot(X,rv(1,:),'k-',X,rv(4,:),'b--',X,rv_taylor,'r--','LineWidth',1.5);
ylim([-0.05 0.05]);
title('6. Real Rate (r, %)');


subplot(5,2,7);
plot(X,Bv(1,:),'k-',X,Bv(4,:),'b--',X,Bv_taylor,'r--','LineWidth',1.5);
ylim([-.1 1.]);
title('7. Debt over GDP (B/Y, %)');


subplot(5,2,8);
plot(X,Tv(1,:),'k-',X,Tv(4,:),'b--',X,Tv_taylor,'r--','LineWidth',1.5);
ylim([-2 9.]);
title('8. Transfer over GDP (T/Y, %)');

subplot(5,2,9);
plot(X,PIpv(1,:),'k-',X,PIpv(4,:),'b--',X,PIpv_taylor,'r--','LineWidth',1.5);
ylim([-0.05 0.4]);
title('9. Price Inflation (\Pi^P, %)');


subplot(5,2,10);
plot(X,PIwv(1,:),'k-',X,PIwv(4,:),'b--',X,PIwv_taylor,'r--','LineWidth',1.5);
ylim([-0.5 0.5]);
title('10. Wage Inflation (\Pi^W, %)');


pos = get(gcf, 'Position');
set(gcf, 'Position',[ 616   600   800   680])
graphName = ['IRFs_SW_','Eco_1_4_taylor','.png' ];
saveas(gcf,graphName)

close




