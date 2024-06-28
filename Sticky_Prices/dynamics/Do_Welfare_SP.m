load For_Welfare_SP


figure(1);
plot(states,CEy,'k-','LineWidth',1.5);
title('Welfare gains (Consumption Equivalent, %)');
xlabel('Productivity Levels');
saveas(gcf,'Welfare_SP_baseline_Eco_1_2.png')

% w

close 
