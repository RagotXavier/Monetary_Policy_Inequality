clear 
load For_Welfare_SW


figure(1);
plot(state_idio,CEy1,'k-','LineWidth',1.5);
title('Welfare gains (Consumption Equivalent, %)');
xlabel('Productivity Levels');
saveas(gcf,'../../Graph/Welfare_SW_baseline_Eco_1_4.png')

