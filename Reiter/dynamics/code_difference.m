% Computing average and maximum difference


load to_code_difference_refined 


% load moments_Comp_Refined
index_Y = find(list_table == "GDP");
index_C = find(list_table == "Ctot");
index_L = find(list_table == "Ltot");
index_w = find(list_table == "w");
index_K = find(list_table == "K");




Ymean = mean_refined(index_Y);
Cmean = mean_refined(index_C);
Lmean = mean_refined(index_L);
Kmean = mean_refined(index_K);
wmean = mean_refined(index_w);



load to_code_difference_reiter 

Nperiod = length(GDP1);
dY_av = sum(abs(GDP1 - GDP2))/Nperiod/Ymean;
dY_max = max(abs(GDP1 - GDP2))/Ymean;

dC_av = sum(abs(C1 - C2))/Nperiod/Cmean;
dC_max = max(abs(C1 - C2))/Cmean;

dK_av = sum(abs(K1 - K2))/Nperiod/Kmean;
dK_max = max(abs(K1 - K2))/Kmean;

dw_av = sum(abs(w1 - w2))/Nperiod/wmean;
dw_max = max(abs(w1 - w2))/wmean;

dr_av = sum(abs(r1 - r2))/Nperiod;
dr_max = max(abs(r1 - r2));

dL_av = sum(abs(L1 - L2))/Nperiod/Lmean;
dL_max = max(abs(L1 - L2))/Lmean;

disp({'Y mean',dY_av;'C mean',dC_av;'K mean ',dK_av;'L mean',dL_av;'r mean',dr_av;'w mean',dw_av;})

disp({'Y max',dY_max;'C max',dC_max;'K max ',dK_max;'L max',dL_max;'r max',dr_max;'w max',dw_max;})


list = ["Y" "C"  "K"  "L"  "r" 'w' ] ;
list_var = cellstr(list)

mean_diff = [dY_av dC_av dK_av  dL_av  dr_av dw_av ]
max_diff = [dY_max dC_max dK_max  dL_max  dr_max dw_max ]


file = 'diff_Reiter.mat'
save(file,"list_var","mean_diff", "max_diff") 

