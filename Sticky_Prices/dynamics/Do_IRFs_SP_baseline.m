
%%% Generate the Figure 1 with IRFS in the Model Sticky Price 

%%% Load Baseline Results 
calib = 'baseline' ;
fileName = ['To_IRFs_SP_',calib,'.mat' ];
load(fileName)

%%% Load Taylor Rule Results 
calib = 'taylor' ;
fileName = ['To_IRFs_SP_',calib];
load(fileName)


figure;
%%% Plotting figures
X=0:1:(LengthIRF);

subplot(5,2,1);
plot(X,Yv(1,:),'k-',X,Yv(2,:),'b--',X,Yv_taylor,'r--','LineWidth',1.5);
ylim([-0.6 0.1]);
title('1. GDP (Y, %)');


subplot(5,2,2);
plot(X,Cv(1,:),'k-',X,Cv(2,:),'b--',X,Cv_taylor,'r--','LineWidth',1.5);
ylim([-.4 .1]);
title('2. Consumption (C, %)');

subplot(5,2,3);
plot(X,Kv(1,:),'k-',X,Kv(2,:),'b--',X,Cv_taylor,'r--','LineWidth',1.5);
ylim([-0.4 .1]);
title('3. Capital (K, %)');

subplot(5,2,4);
plot(X,Lv(1,:),'k-',X,Lv(2,:),'b--',X,Lv_taylor,'r--','LineWidth',1.5);
ylim([-.2 .2]);
title('4. Labor (L, %)');

subplot(5,2,5);
plot(X,rnv(1,:),'k-',X,rnv(2,:),'b--',X,rnv_taylor,'r--','LineWidth',1.5);
ylim([-0.05 0.05]);
title('5. Nominal interest rate (R^{N}, %)');


subplot(5,2,6);
plot(X,rv(1,:),'k-',X,rv(2,:),'b--',X,rv_taylor,'r--','LineWidth',1.5);
ylim([-0.05 0.05]);
title('6. Real Rate (r, %)');


subplot(5,2,7);
plot(X,Bv(1,:),'k-',X,Bv(2,:),'b--',X,Bv_taylor,'r--','LineWidth',1.5);
ylim([-.1 0.8]);
title('7. Debt over GDP (B/Y, %)');

subplot(5,2,8);
plot(X,Tv(1,:),'k-',X,Tv(2,:),'b--',X,Tv_taylor,'r--','LineWidth',1.5);
ylim([-2 6.]);
title('8. Transfer over GDP (T/Y, %)');


subplot(5,2,9);
plot(X,pipv(1,:),'k-',X,pipv(2,:),'b--',X,pipv_taylor,'r--','LineWidth',1.5);
ylim([-0.05 0.05]);
title('9. Price Inflation (\Pi^P, %)');


subplot(5,2,10);
plot(X,piwv(1,:),'k-',X,piwv(2,:),'b--',X,piwv_taylor,'r--','LineWidth',1.5);
ylim([-1 0.7]);
title('10. Wage Inflation (\Pi^W, %)');



pos = get(gcf, 'Position');
set(gcf, 'Position',[ 616   600   800   680])
graphName = ['IRFs_SP_','Eco_1_2_taylor','.png' ];
saveas(gcf,graphName)

close


