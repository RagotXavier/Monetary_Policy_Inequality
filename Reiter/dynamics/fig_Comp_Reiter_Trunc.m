

load todiff_Reiter

u1   = ut;
GDP1 = GDPt;
C1   = Ct;
C1b   = Ctott;

w1   = wt;
r1   = rt;
K1   = Kt;
L1 = Lt;


load todiff_Comp_Refined


u2   = ut;
GDP2 = GDPt;
C2   = Ct;
w2   = wt;
r2   = rt;
K2   = Kt;
L2 = Lt;


%LL = 100;%length(u1);%length(B_eps);
LL = length(u1); %length(B_eps);

%
X=1:1:LL;

figure;
 
   
subplot(2,3,1);
plot(X,GDP1(1:LL),'k-',X,GDP2(1:LL),'b--','LineWidth',1.5);
title('1. $Y  ~ (\%)$','Interpreter','latex');


subplot(2,3,2);
plot(X,C1(1:LL),'k-',X,C2(1:LL),'b--','LineWidth',1.5);
title('2. $C ~ (\%)$','Interpreter','latex');
%ylim([0. .2]);

subplot(2,3,3);
plot(X,K1(1:LL),'k-',X,K2(1:LL),'b--','LineWidth',1.5);
title('3. $K  ~ (\%)$','Interpreter','latex');
%ylim([-2 1]);
 
%ylim([0. .3]);
legend('Reiter','Truncation')
legend boxoff 

subplot(2,3,4);
plot(X,w1(1:LL),'k-',X,w2(1:LL),'b--','LineWidth',1.5);
title('4. $w  ~ (\%)$','Interpreter','latex');
%ylim([-2 1]);

subplot(2,3,5);
plot(X,r1(1:LL),'k-',X,r2(1:LL),'b--','LineWidth',1.5);
title('5. $r  ~ (\%)$','Interpreter','latex');
%ylim([0 6]);

subplot(2,3,6);
plot(X,L1(1:LL),'k-',X,L2(1:LL),'b--','LineWidth',1.5);
title('6. $L  ~ (\%)$','Interpreter','latex');



saveas(figure(1),'../../Graph/Comp_Reiter_Trunc.png')



% 
% 
% clear 
% %load tofigcomptrunc
% 
% load tofigcomptrunc_5
% %Cv_trunc = [c1_eps c2_eps c3_eps c4_eps c5_eps c6_eps c7_eps c8_eps c9_eps c10_eps c11_eps c12_eps c13_eps c14_eps c15_eps c16_eps c17_eps c18_eps c19_eps c20_eps c21_eps c22_eps c23_eps c24_eps c25_eps c26_eps c27_eps c28_eps c29_eps c30_eps c31_eps c32_eps c33_eps c34_eps c35_eps c36_eps c37_eps c38_eps c39_eps c40_eps c41_eps c42_eps c43_eps c44_eps c45_eps c46_eps c47_eps c48_eps c49_eps c50_eps c51_eps c52_eps c53_eps c54_eps c55_eps c56_eps c57_eps c58_eps c59_eps c60_eps c61_eps c62_eps c63_eps c64_eps c65_eps c66_eps c67_eps c68_eps c69_eps c70_eps c71_eps c72_eps c73_eps c74_eps c75_eps c76_eps c77_eps c78_eps c79_eps c80_eps c81_eps c82_eps c83_eps c84_eps c85_eps c86_eps c87_eps c88_eps c89_eps c90_eps c91_eps c92_eps c93_eps c94_eps c95_eps c96_eps c97_eps c98_eps c99_eps c100_eps c101_eps c102_eps c103_eps c104_eps c105_eps c106_eps c107_eps c108_eps c109_eps ] % equation 454
% %Cv_trunc= [c1_eps c2_eps c3_eps c4_eps c5_eps c6_eps c7_eps c8_eps c9_eps c10_eps c11_eps c12_eps c13_eps c14_eps c15_eps c16_eps c17_eps c18_eps c19_eps c20_eps c21_eps c22_eps c23_eps c24_eps c25_eps c26_eps c27_eps c28_eps c29_eps c30_eps c31_eps c32_eps c33_eps c34_eps c35_eps c36_eps c37_eps c38_eps c39_eps c40_eps c41_eps c42_eps c43_eps c44_eps c45_eps c46_eps c47_eps c48_eps c49_eps c50_eps c51_eps c52_eps c53_eps c54_eps c55_eps c56_eps c57_eps c58_eps c59_eps c60_eps c61_eps c62_eps c63_eps c64_eps c65_eps c66_eps c67_eps c68_eps c69_eps c70_eps c71_eps c72_eps c73_eps c74_eps c75_eps c76_eps c77_eps c78_eps c79_eps c80_eps c81_eps c82_eps c83_eps c84_eps c85_eps c86_eps c87_eps c88_eps c89_eps c90_eps c91_eps c92_eps c93_eps c94_eps c95_eps c96_eps c97_eps c98_eps c99_eps c100_eps c101_eps c102_eps c103_eps c104_eps c105_eps c106_eps c107_eps c108_eps c109_eps c110_eps c111_eps c112_eps c113_eps c114_eps c115_eps c116_eps c117_eps c118_eps c119_eps c120_eps c121_eps c122_eps c123_eps c124_eps c125_eps c126_eps c127_eps c128_eps c129_eps c130_eps c131_eps c132_eps c133_eps c134_eps c135_eps c136_eps c137_eps c138_eps c139_eps c140_eps c141_eps c142_eps c143_eps c144_eps c145_eps c146_eps c147_eps c148_eps c149_eps c150_eps c151_eps c152_eps c153_eps c154_eps c155_eps c156_eps c157_eps c158_eps c159_eps c160_eps c161_eps c162_eps c163_eps c164_eps c165_eps c166_eps c167_eps c168_eps c169_eps c170_eps c171_eps c172_eps c173_eps c174_eps c175_eps c176_eps c177_eps c178_eps c179_eps c180_eps c181_eps c182_eps c183_eps c184_eps c185_eps c186_eps c187_eps c188_eps c189_eps c190_eps c191_eps c192_eps c193_eps c194_eps c195_eps c196_eps c197_eps c198_eps c199_eps c200_eps c201_eps c202_eps c203_eps c204_eps c205_eps c206_eps c207_eps c208_eps c209_eps c210_eps c211_eps c212_eps c213_eps c214_eps c215_eps c216_eps c217_eps c218_eps c219_eps c220_eps c221_eps c222_eps c223_eps c224_eps c225_eps c226_eps c227_eps c228_eps c229_eps c230_eps c231_eps c232_eps c233_eps c234_eps c235_eps c236_eps c237_eps c238_eps c239_eps c240_eps c241_eps c242_eps c243_eps c244_eps c245_eps c246_eps c247_eps c248_eps c249_eps c250_eps c251_eps c252_eps c253_eps c254_eps c255_eps c256_eps c257_eps c258_eps c259_eps c260_eps c261_eps c262_eps c263_eps c264_eps c265_eps c266_eps c267_eps c268_eps c269_eps c270_eps c271_eps c272_eps c273_eps c274_eps c275_eps c276_eps c277_eps c278_eps c279_eps c280_eps c281_eps c282_eps c283_eps c284_eps c285_eps c286_eps c287_eps c288_eps c289_eps c290_eps c291_eps c292_eps c293_eps c294_eps c295_eps c296_eps c297_eps c298_eps c299_eps c300_eps c301_eps c302_eps c303_eps c304_eps c305_eps c306_eps c307_eps c308_eps c309_eps c310_eps c311_eps c312_eps c313_eps c314_eps c315_eps c316_eps c317_eps c318_eps c319_eps c320_eps c321_eps c322_eps c323_eps c324_eps c325_eps c326_eps c327_eps c328_eps c329_eps c330_eps c331_eps c332_eps c333_eps c334_eps c335_eps c336_eps c337_eps c338_eps c339_eps c340_eps c341_eps c342_eps c343_eps c344_eps c345_eps c346_eps c347_eps c348_eps c349_eps c350_eps c351_eps c352_eps c353_eps c354_eps c355_eps c356_eps c357_eps c358_eps c359_eps c360_eps c361_eps c362_eps c363_eps c364_eps c365_eps c366_eps c367_eps c368_eps c369_eps c370_eps c371_eps c372_eps c373_eps c374_eps c375_eps c376_eps c377_eps c378_eps c379_eps c380_eps c381_eps c382_eps c383_eps c384_eps c385_eps c386_eps c387_eps c388_eps c389_eps c390_eps c391_eps c392_eps c393_eps c394_eps c395_eps c396_eps c397_eps c398_eps c399_eps c400_eps c401_eps c402_eps c403_eps c404_eps c405_eps c406_eps c407_eps c408_eps c409_eps c410_eps c411_eps c412_eps c413_eps c414_eps c415_eps c416_eps c417_eps c418_eps c419_eps c420_eps c421_eps c422_eps c423_eps c424_eps c425_eps c426_eps c427_eps c428_eps c429_eps c430_eps c431_eps c432_eps c433_eps c434_eps c435_eps c436_eps c437_eps c438_eps c439_eps c440_eps c441_eps c442_eps c443_eps c444_eps c445_eps c446_eps c447_eps c448_eps c449_eps c450_eps c451_eps c452_eps c453_eps c454_eps c455_eps c456_eps c457_eps c458_eps c459_eps c460_eps c461_eps c462_eps c463_eps c464_eps c465_eps c466_eps c467_eps c468_eps c469_eps c470_eps c471_eps c472_eps c473_eps c474_eps c475_eps c476_eps c477_eps c478_eps c479_eps c480_eps c481_eps c482_eps c483_eps c484_eps c485_eps c486_eps c487_eps c488_eps c489_eps c490_eps c491_eps c492_eps c493_eps c494_eps c495_eps c496_eps c497_eps c498_eps c499_eps c500_eps c501_eps c502_eps c503_eps c504_eps c505_eps c506_eps c507_eps c508_eps c509_eps c510_eps c511_eps c512_eps c513_eps c514_eps c515_eps c516_eps c517_eps c518_eps c519_eps c520_eps c521_eps c522_eps c523_eps c524_eps c525_eps c526_eps c527_eps c528_eps c529_eps c530_eps c531_eps c532_eps c533_eps c534_eps c535_eps c536_eps c537_eps c538_eps c539_eps c540_eps c541_eps c542_eps c543_eps c544_eps c545_eps c546_eps c547_eps c548_eps c549_eps c550_eps c551_eps c552_eps c553_eps c554_eps c555_eps c556_eps c557_eps c558_eps c559_eps c560_eps c561_eps c562_eps c563_eps c564_eps c565_eps c566_eps c567_eps c568_eps c569_eps c570_eps c571_eps c572_eps c573_eps c574_eps c575_eps c576_eps c577_eps c578_eps c579_eps c580_eps c581_eps c582_eps c583_eps c584_eps c585_eps c586_eps c587_eps c588_eps c589_eps c590_eps c591_eps c592_eps c593_eps c594_eps c595_eps c596_eps c597_eps c598_eps c599_eps c600_eps c601_eps c602_eps c603_eps c604_eps c605_eps c606_eps c607_eps c608_eps c609_eps c610_eps c611_eps c612_eps c613_eps c614_eps c615_eps c616_eps c617_eps c618_eps c619_eps c620_eps c621_eps c622_eps c623_eps c624_eps c625_eps c626_eps c627_eps c628_eps c629_eps c630_eps c631_eps c632_eps c633_eps c634_eps c635_eps c636_eps c637_eps c638_eps c639_eps c640_eps c641_eps c642_eps c643_eps c644_eps c645_eps c646_eps c647_eps c648_eps c649_eps c650_eps c651_eps c652_eps c653_eps c654_eps c655_eps c656_eps c657_eps c658_eps c659_eps c660_eps c661_eps c662_eps c663_eps c664_eps c665_eps c666_eps c667_eps c668_eps c669_eps c670_eps c671_eps c672_eps c673_eps c674_eps c675_eps c676_eps c677_eps c678_eps c679_eps c680_eps c681_eps c682_eps c683_eps c684_eps c685_eps c686_eps c687_eps c688_eps c689_eps c690_eps c691_eps c692_eps c693_eps c694_eps c695_eps c696_eps c697_eps c698_eps c699_eps c700_eps c701_eps c702_eps c703_eps c704_eps c705_eps c706_eps c707_eps c708_eps c709_eps c710_eps c711_eps c712_eps c713_eps c714_eps c715_eps c716_eps c717_eps c718_eps c719_eps c720_eps c721_eps c722_eps c723_eps c724_eps c725_eps c726_eps c727_eps ] ;% equation 2926
% str = ['Cv',' = [']; fprintf(fid, str);
%         for h = 1:Nbin
%            str = ['c',num2str(h),'_eps',' ']; fprintf(fid, str);
%         end
% str = [']']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1; 
% 
% 
% 
% 
% % Truncation
% Cv_trunc = [c1_eps c2_eps c3_eps c4_eps c5_eps c6_eps c7_eps c8_eps c9_eps c10_eps c11_eps c12_eps c13_eps c14_eps c15_eps c16_eps c17_eps c18_eps c19_eps c20_eps c21_eps c22_eps c23_eps c24_eps c25_eps c26_eps c27_eps c28_eps c29_eps c30_eps c31_eps c32_eps c33_eps c34_eps c35_eps c36_eps c37_eps c38_eps c39_eps c40_eps c41_eps c42_eps c43_eps c44_eps c45_eps c46_eps c47_eps c48_eps c49_eps c50_eps c51_eps c52_eps c53_eps c54_eps c55_eps c56_eps c57_eps c58_eps c59_eps ] ; % equation 289
% 
% C_ss = zeros(Nbin,1) ;
% 
% for h = 1:Nbin
%     C_ss(h) = eco.cp(h)/Sp(h) ; % This is the value per capita at ss 
% end
% 
% Css = repmat(C_ss', Length,1) ; %'
% 
% Cdev = Cv_trunc + Css ;  % Here I recover the variable in level and no deviation from steady state 
% 
% [out_ss idx_ss] = sort(C_ss,1) ;
% Csort_ss = C_ss(idx_ss) ;
% SPsort_ss = Sp(idx_ss);
% 
% 
% Cumdis_css = cumsum(SPsort_ss) ;
% ind_1_css = find(Cumdis_css <= 0.10);
% ind_2_css = find(Cumdis_css > 0.10 & Cumdis_css < 0.35 );
% ind_3_css = find(Cumdis_css > 0.35 & Cumdis_css <= 0.70 );
% ind_4_css = find(Cumdis_css > 0.70 & Cumdis_css <= 0.93555 );
% ind_5_css = find(Cumdis_css > 0.93555 );
% 
% 
% 
% 
% index_1_css = idx_ss(ind_1_css); % find the point zhere 10 and 90 are located at ss
% index_2_css = idx_ss(ind_2_css); % find the point zhere 10 and 90 are located at ss
% index_3_css = idx_ss(ind_3_css);
% index_4_css = idx_ss(ind_4_css);
% index_5_css = idx_ss(ind_5_css);
% 
% 
% 
% for i=1:Length
%     D = Cdev(i,:);
%     E1 = Sp(index_1_css)'* D(index_1_css)' ; %
%     E2 = Sp(index_2_css)'* D(index_2_css)' ; %
%     E3 = Sp(index_3_css)'* D(index_3_css)' ; %
%     E4 = Sp(index_4_css)'* D(index_4_css)' ; %
%     E5 = Sp(index_5_css)'* D(index_5_css)' ; %
%     C1_c(i,:) = E1 / sum(Sp(index_1_css)); % conso moyenne par tete des 10% qui ont le moins
%     C2_c(i,:) = E2 / sum(Sp(index_2_css)); % conso moyenne par tete des 10% qui ont le moins
%     C3_c(i,:) = E3 / sum(Sp(index_3_css)); % conso moyenne par tete des 10% qui ont le plus plus
%     C4_c(i,:) = E4 / sum(Sp(index_4_css)); % conso moyenne par tete des 10% qui ont le moins
%     C5_c(i,:) = E5 / sum(Sp(index_5_css)); % conso moyenne par tete des 10% qui ont le plus plus
% end    
% 
% 
% 
% E1 = Sp(index_1_css)'* C_ss(index_1_css) ; %'
% E2 = Sp(index_2_css)'* C_ss(index_2_css) ; %'
% E3 = Sp(index_3_css)'* C_ss(index_3_css);  %'
% E4 = Sp(index_4_css)'* C_ss(index_4_css) ; %'
% E5 = Sp(index_5_css)'* C_ss(index_5_css);  %'
% 
% 
% 
% C1_css = E1 / sum(Sp(index_1_css)) ;
% C2_css = E2 / sum(Sp(index_2_css)) ;
% C3_css = E3 / sum(Sp(index_3_css));
% C4_css = E4 / sum(Sp(index_4_css)) ;
% C5_css = E5 / sum(Sp(index_5_css));
% 
% 
% 
% percentile_1cv_trunc = 100*(C1_c -C1_css)/ C1_css
% percentile_2cv_trunc = 100*(C2_c -  C2_css) / C2_css
% percentile_3cv_trunc = 100*(C3_c -  C3_css) / C3_css
% percentile_4cv_trunc = 100*(C4_c -  C4_css) / C4_css
% percentile_5cv_trunc = 100*(C5_c -  C5_css) / C5_css
% 
% percentile_1cl_trunc = (C1_c -C1_css)
% percentile_2cl_trunc = (C2_c -  C2_css) 
% percentile_3cl_trunc = (C3_c -  C3_css) 
% percentile_4cl_trunc = (C4_c -  C4_css) 
% percentile_5cl_trunc = (C5_c -  C5_css) 
% 
% 
% 
% 
% X=1:1:(Length);
% figure(9)
% plot(X,percentile_1cv_trunc,X,percentile_2cv_trunc,X, percentile_3cv_trunc,X, percentile_4cv_trunc,X, percentile_5cv_trunc);
% legend('bottom 7 percent','2','3','4','top')
% title('1. Evolution of consumption after an MP shock with respect to ss value');
% saveas(figure(9),'percent_consum_decile_trun.png')
% 
% 
% X=1:1:(Length);
% figure(99)
% plot(X,percentile_1cl_trunc,X,percentile_2cl_trunc,X, percentile_3cl_trunc,X, percentile_4cl_trunc,X, percentile_5cl_trunc);
% legend('bottom 7 percent','2','3','4','top')
% title('1. Evolution level  consumption after an MP shock with respect to ss value');
% saveas(figure(99),'level_consum_decile_trun.png')
% 
% 
% 
% 
% X=1:1:(Length);
% figure(100)
% plot(X,percentile_1cl_trunc);
% legend('bottom 7 percent','2','3','4','top')
% title('1. Evolution level  consumption after an MP shock with respect to ss value');
% saveas(figure(100),'level1_consum_decile_trun.png')
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% save tofigcomptrunc
% 
% % Reiter
% clear 
% %load tofigreiter 
% 
% load tofigreiter_5
% Cv = [c1_eps c2_eps c3_eps c4_eps c5_eps c6_eps c7_eps c8_eps c9_eps c10_eps c11_eps c12_eps c13_eps c14_eps c15_eps c16_eps c17_eps c18_eps c19_eps c20_eps c21_eps c22_eps c23_eps c24_eps c25_eps c26_eps c27_eps c28_eps c29_eps c30_eps c31_eps c32_eps c33_eps c34_eps c35_eps c36_eps c37_eps c38_eps c39_eps c40_eps c41_eps c42_eps c43_eps c44_eps c45_eps c46_eps c47_eps c48_eps c49_eps c50_eps c51_eps c52_eps c53_eps c54_eps c55_eps c56_eps c57_eps c58_eps c59_eps c60_eps c61_eps c62_eps c63_eps c64_eps c65_eps c66_eps c67_eps c68_eps c69_eps c70_eps c71_eps c72_eps c73_eps c74_eps c75_eps c76_eps c77_eps c78_eps c79_eps c80_eps c81_eps c82_eps c83_eps c84_eps c85_eps c86_eps c87_eps c88_eps c89_eps c90_eps c91_eps c92_eps c93_eps c94_eps c95_eps c96_eps c97_eps c98_eps c99_eps c100_eps c101_eps c102_eps c103_eps c104_eps c105_eps c106_eps c107_eps c108_eps c109_eps c110_eps c111_eps c112_eps c113_eps c114_eps c115_eps c116_eps c117_eps c118_eps c119_eps c120_eps c121_eps c122_eps c123_eps c124_eps c125_eps c126_eps c127_eps c128_eps c129_eps c130_eps c131_eps c132_eps c133_eps c134_eps c135_eps c136_eps c137_eps c138_eps c139_eps c140_eps c141_eps c142_eps c143_eps c144_eps c145_eps c146_eps c147_eps c148_eps c149_eps c150_eps c151_eps c152_eps c153_eps c154_eps c155_eps c156_eps c157_eps c158_eps c159_eps c160_eps c161_eps c162_eps c163_eps c164_eps c165_eps c166_eps c167_eps c168_eps c169_eps c170_eps c171_eps c172_eps c173_eps c174_eps c175_eps c176_eps c177_eps c178_eps c179_eps c180_eps c181_eps c182_eps c183_eps c184_eps c185_eps c186_eps c187_eps c188_eps c189_eps c190_eps c191_eps c192_eps c193_eps c194_eps c195_eps c196_eps c197_eps c198_eps c199_eps c200_eps c201_eps c202_eps c203_eps c204_eps c205_eps c206_eps c207_eps c208_eps c209_eps c210_eps c211_eps c212_eps c213_eps c214_eps c215_eps c216_eps c217_eps c218_eps c219_eps c220_eps c221_eps c222_eps c223_eps c224_eps c225_eps c226_eps c227_eps c228_eps c229_eps c230_eps c231_eps c232_eps c233_eps c234_eps c235_eps c236_eps c237_eps c238_eps c239_eps c240_eps c241_eps c242_eps c243_eps c244_eps c245_eps c246_eps c247_eps c248_eps c249_eps c250_eps c251_eps c252_eps c253_eps c254_eps c255_eps c256_eps c257_eps c258_eps c259_eps c260_eps c261_eps c262_eps c263_eps c264_eps c265_eps c266_eps c267_eps c268_eps c269_eps c270_eps c271_eps c272_eps c273_eps c274_eps c275_eps c276_eps c277_eps c278_eps c279_eps c280_eps c281_eps c282_eps c283_eps c284_eps c285_eps c286_eps c287_eps c288_eps c289_eps c290_eps c291_eps c292_eps c293_eps c294_eps c295_eps c296_eps c297_eps c298_eps c299_eps c300_eps c301_eps c302_eps c303_eps c304_eps c305_eps c306_eps c307_eps c308_eps c309_eps c310_eps c311_eps c312_eps c313_eps c314_eps c315_eps c316_eps c317_eps c318_eps c319_eps c320_eps c321_eps c322_eps c323_eps c324_eps c325_eps c326_eps c327_eps c328_eps c329_eps c330_eps c331_eps c332_eps c333_eps c334_eps c335_eps c336_eps c337_eps c338_eps c339_eps c340_eps c341_eps c342_eps c343_eps c344_eps c345_eps c346_eps c347_eps c348_eps c349_eps c350_eps c351_eps c352_eps c353_eps c354_eps c355_eps c356_eps c357_eps c358_eps c359_eps c360_eps c361_eps c362_eps c363_eps c364_eps c365_eps c366_eps c367_eps c368_eps c369_eps c370_eps c371_eps c372_eps c373_eps c374_eps c375_eps c376_eps c377_eps c378_eps c379_eps c380_eps c381_eps c382_eps c383_eps c384_eps c385_eps c386_eps c387_eps c388_eps c389_eps c390_eps c391_eps c392_eps c393_eps c394_eps c395_eps c396_eps c397_eps c398_eps c399_eps c400_eps c401_eps c402_eps c403_eps c404_eps c405_eps c406_eps c407_eps c408_eps c409_eps c410_eps c411_eps c412_eps c413_eps c414_eps c415_eps c416_eps c417_eps c418_eps c419_eps c420_eps c421_eps c422_eps c423_eps c424_eps c425_eps c426_eps c427_eps c428_eps c429_eps c430_eps c431_eps c432_eps c433_eps c434_eps c435_eps c436_eps c437_eps c438_eps c439_eps c440_eps c441_eps c442_eps c443_eps c444_eps c445_eps c446_eps c447_eps c448_eps c449_eps c450_eps c451_eps c452_eps c453_eps c454_eps c455_eps c456_eps c457_eps c458_eps c459_eps c460_eps c461_eps c462_eps c463_eps c464_eps c465_eps c466_eps c467_eps c468_eps c469_eps c470_eps c471_eps c472_eps c473_eps c474_eps c475_eps c476_eps c477_eps c478_eps c479_eps c480_eps c481_eps c482_eps c483_eps c484_eps c485_eps c486_eps c487_eps c488_eps c489_eps c490_eps c491_eps c492_eps c493_eps c494_eps c495_eps c496_eps c497_eps c498_eps c499_eps c500_eps c501_eps c502_eps c503_eps c504_eps c505_eps c506_eps c507_eps c508_eps c509_eps c510_eps c511_eps c512_eps c513_eps c514_eps c515_eps c516_eps c517_eps c518_eps c519_eps c520_eps c521_eps c522_eps c523_eps c524_eps c525_eps c526_eps c527_eps c528_eps c529_eps c530_eps c531_eps c532_eps c533_eps c534_eps c535_eps c536_eps c537_eps c538_eps c539_eps c540_eps c541_eps c542_eps c543_eps c544_eps c545_eps c546_eps c547_eps c548_eps c549_eps c550_eps c551_eps c552_eps c553_eps c554_eps c555_eps c556_eps c557_eps c558_eps c559_eps c560_eps c561_eps c562_eps c563_eps c564_eps c565_eps c566_eps c567_eps c568_eps c569_eps c570_eps c571_eps c572_eps c573_eps c574_eps c575_eps c576_eps c577_eps c578_eps c579_eps c580_eps c581_eps c582_eps c583_eps c584_eps c585_eps c586_eps c587_eps c588_eps c589_eps c590_eps c591_eps c592_eps c593_eps c594_eps c595_eps c596_eps c597_eps c598_eps c599_eps c600_eps c601_eps c602_eps c603_eps c604_eps c605_eps c606_eps c607_eps c608_eps c609_eps c610_eps c611_eps c612_eps c613_eps c614_eps c615_eps c616_eps c617_eps c618_eps c619_eps c620_eps c621_eps c622_eps c623_eps c624_eps c625_eps c626_eps c627_eps c628_eps c629_eps c630_eps c631_eps c632_eps c633_eps c634_eps c635_eps c636_eps c637_eps c638_eps c639_eps c640_eps c641_eps c642_eps c643_eps c644_eps c645_eps c646_eps c647_eps c648_eps c649_eps c650_eps c651_eps c652_eps c653_eps c654_eps c655_eps c656_eps c657_eps c658_eps c659_eps c660_eps c661_eps c662_eps c663_eps c664_eps c665_eps c666_eps c667_eps c668_eps c669_eps c670_eps c671_eps c672_eps c673_eps c674_eps c675_eps c676_eps c677_eps c678_eps c679_eps c680_eps c681_eps c682_eps c683_eps c684_eps c685_eps c686_eps c687_eps c688_eps c689_eps c690_eps c691_eps c692_eps c693_eps c694_eps c695_eps c696_eps c697_eps c698_eps c699_eps c700_eps ] ;
% file = ['code_dynare_Comp_Reiter','.mod'];
% fid = fopen(file, 'w') ;
% 
% formatSpec  = '%16.12f';
% 
% str = ['Cv',' = [']; fprintf(fid, str);
%         for h = 1:Ntot
%            str = ['c',num2str(h),'_eps',' ']; fprintf(fid, str);
%         end
% str = [']']; fprintf(fid, str); 
% 
% 
% Cv_reit = [c1_eps c2_eps c3_eps c4_eps c5_eps c6_eps c7_eps c8_eps c9_eps c10_eps c11_eps c12_eps c13_eps c14_eps c15_eps c16_eps c17_eps c18_eps c19_eps c20_eps c21_eps c22_eps c23_eps c24_eps c25_eps c26_eps c27_eps c28_eps c29_eps c30_eps c31_eps c32_eps c33_eps c34_eps c35_eps c36_eps c37_eps c38_eps c39_eps c40_eps c41_eps c42_eps c43_eps c44_eps c45_eps c46_eps c47_eps c48_eps c49_eps c50_eps c51_eps c52_eps c53_eps c54_eps c55_eps c56_eps c57_eps c58_eps c59_eps c60_eps c61_eps c62_eps c63_eps c64_eps c65_eps c66_eps c67_eps c68_eps c69_eps c70_eps c71_eps c72_eps c73_eps c74_eps c75_eps c76_eps c77_eps c78_eps c79_eps c80_eps c81_eps c82_eps c83_eps c84_eps c85_eps c86_eps c87_eps c88_eps c89_eps c90_eps c91_eps c92_eps c93_eps c94_eps c95_eps c96_eps c97_eps c98_eps c99_eps c100_eps c101_eps c102_eps c103_eps c104_eps c105_eps c106_eps c107_eps c108_eps c109_eps c110_eps c111_eps c112_eps c113_eps c114_eps c115_eps c116_eps c117_eps c118_eps c119_eps c120_eps c121_eps c122_eps c123_eps c124_eps c125_eps c126_eps c127_eps c128_eps c129_eps c130_eps c131_eps c132_eps c133_eps c134_eps c135_eps c136_eps c137_eps c138_eps c139_eps c140_eps c141_eps c142_eps c143_eps c144_eps c145_eps c146_eps c147_eps c148_eps c149_eps c150_eps c151_eps c152_eps c153_eps c154_eps c155_eps c156_eps c157_eps c158_eps c159_eps c160_eps c161_eps c162_eps c163_eps c164_eps c165_eps c166_eps c167_eps c168_eps c169_eps c170_eps c171_eps c172_eps c173_eps c174_eps c175_eps c176_eps c177_eps c178_eps c179_eps c180_eps c181_eps c182_eps c183_eps c184_eps c185_eps c186_eps c187_eps c188_eps c189_eps c190_eps c191_eps c192_eps c193_eps c194_eps c195_eps c196_eps c197_eps c198_eps c199_eps c200_eps c201_eps c202_eps c203_eps c204_eps c205_eps c206_eps c207_eps c208_eps c209_eps c210_eps c211_eps c212_eps c213_eps c214_eps c215_eps c216_eps c217_eps c218_eps c219_eps c220_eps c221_eps c222_eps c223_eps c224_eps c225_eps c226_eps c227_eps c228_eps c229_eps c230_eps c231_eps c232_eps c233_eps c234_eps c235_eps c236_eps c237_eps c238_eps c239_eps c240_eps c241_eps c242_eps c243_eps c244_eps c245_eps c246_eps c247_eps c248_eps c249_eps c250_eps c251_eps c252_eps c253_eps c254_eps c255_eps c256_eps c257_eps c258_eps c259_eps c260_eps c261_eps c262_eps c263_eps c264_eps c265_eps c266_eps c267_eps c268_eps c269_eps c270_eps c271_eps c272_eps c273_eps c274_eps c275_eps c276_eps c277_eps c278_eps c279_eps c280_eps c281_eps c282_eps c283_eps c284_eps c285_eps c286_eps c287_eps c288_eps c289_eps c290_eps c291_eps c292_eps c293_eps c294_eps c295_eps c296_eps c297_eps c298_eps c299_eps c300_eps c301_eps c302_eps c303_eps c304_eps c305_eps c306_eps c307_eps c308_eps c309_eps c310_eps c311_eps c312_eps c313_eps c314_eps c315_eps c316_eps c317_eps c318_eps c319_eps c320_eps c321_eps c322_eps c323_eps c324_eps c325_eps c326_eps c327_eps c328_eps c329_eps c330_eps c331_eps c332_eps c333_eps c334_eps c335_eps c336_eps c337_eps c338_eps c339_eps c340_eps c341_eps c342_eps c343_eps c344_eps c345_eps c346_eps c347_eps c348_eps c349_eps c350_eps c351_eps c352_eps c353_eps c354_eps c355_eps c356_eps c357_eps c358_eps c359_eps c360_eps c361_eps c362_eps c363_eps c364_eps c365_eps c366_eps c367_eps c368_eps c369_eps c370_eps c371_eps c372_eps c373_eps c374_eps c375_eps c376_eps c377_eps c378_eps c379_eps c380_eps c381_eps c382_eps c383_eps c384_eps c385_eps c386_eps c387_eps c388_eps c389_eps c390_eps c391_eps c392_eps c393_eps c394_eps c395_eps c396_eps c397_eps c398_eps c399_eps c400_eps c401_eps c402_eps c403_eps c404_eps c405_eps c406_eps c407_eps c408_eps c409_eps c410_eps c411_eps c412_eps c413_eps c414_eps c415_eps c416_eps c417_eps c418_eps c419_eps c420_eps c421_eps c422_eps c423_eps c424_eps c425_eps c426_eps c427_eps c428_eps c429_eps c430_eps c431_eps c432_eps c433_eps c434_eps c435_eps c436_eps c437_eps c438_eps c439_eps c440_eps c441_eps c442_eps c443_eps c444_eps c445_eps c446_eps c447_eps c448_eps c449_eps c450_eps c451_eps c452_eps c453_eps c454_eps c455_eps c456_eps c457_eps c458_eps c459_eps c460_eps c461_eps c462_eps c463_eps c464_eps c465_eps c466_eps c467_eps c468_eps c469_eps c470_eps c471_eps c472_eps c473_eps c474_eps c475_eps c476_eps c477_eps c478_eps c479_eps c480_eps c481_eps c482_eps c483_eps c484_eps c485_eps c486_eps c487_eps c488_eps c489_eps c490_eps c491_eps c492_eps c493_eps c494_eps c495_eps c496_eps c497_eps c498_eps c499_eps c500_eps ];
% 
% 
% 
% for h = 1:Ntot
%     C_ss(h) = Cpol(h) ; % This is the value per capita at ss 
% end
% 
% 
% 
% Css = repmat(C_ss, Length,1) ; %
% 
% Cdev_reit = Cv_reit + Css ;  % Here I recover the variable in level and no deviation from steady state 
% 
% [out_ss idx_ss] = sort(C_ss,2) ;
% Csort_ss = C_ss(idx_ss) ;
% SPsort_ss = D_ss(idx_ss);
% 
% 
% Cumdis_css = cumsum(SPsort_ss) ;
% ind_1_css = find(Cumdis_css <= 0.10);
% ind_2_css = find(Cumdis_css > 0.10 & Cumdis_css < 0.35 );
% ind_3_css = find(Cumdis_css > 0.35 & Cumdis_css <= 0.70 );
% ind_4_css = find(Cumdis_css > 0.70 & Cumdis_css <= 0.93555 );
% ind_5_css = find(Cumdis_css > 0.93555 );
% 
% 
% 
% 
% index_1_css = idx_ss(ind_1_css); % find the point zhere 10 and 90 are located at ss
% index_2_css = idx_ss(ind_2_css); % find the point zhere 10 and 90 are located at ss
% index_3_css = idx_ss(ind_3_css);
% index_4_css = idx_ss(ind_4_css);
% index_5_css = idx_ss(ind_5_css);
% 
% 
% sum(D_ss(index_1_css)) + sum(D_ss(index_2_css)) + sum(D_ss(index_3_css)) + sum(D_ss(index_4_css)) +sum(D_ss(index_5_css)) 
% 
% 
% for i=1:Length
%     D = Cdev_reit(i,:);
%     E1 = D_ss(index_1_css)'* D(index_1_css)' ; %
%     E2 = D_ss(index_2_css)'* D(index_2_css)' ; %
%     E3 = D_ss(index_3_css)'* D(index_3_css)' ; %
%     E4 = D_ss(index_4_css)'* D(index_4_css)' ; %
%     E5 = D_ss(index_5_css)'* D(index_5_css)' ; %
%     C1_c(i,:) = E1 / sum(D_ss(index_1_css)); % conso moyenne par tete des 10% qui ont le moins
%     C2_c(i,:) = E2 / sum(D_ss(index_2_css)); % conso moyenne par tete des 10% qui ont le moins
%     C3_c(i,:) = E3 / sum(D_ss(index_3_css)); % conso moyenne par tete des 10% qui ont le plus plus
%     C4_c(i,:) = E4 / sum(D_ss(index_4_css)); % conso moyenne par tete des 10% qui ont le moins
%     C5_c(i,:) = E5 / sum(D_ss(index_5_css)); % conso moyenne par tete des 10% qui ont le plus plus
% end    
% 
% 
% E1 = D_ss(index_1_css)'* C_ss(index_1_css)' ; %'
% E2 = D_ss(index_2_css)'* C_ss(index_2_css)' ; %'
% E3 = D_ss(index_3_css)'* C_ss(index_3_css)';  %'
% E4 = D_ss(index_4_css)'* C_ss(index_4_css)' ; %'
% E5 = D_ss(index_5_css)'* C_ss(index_5_css)';  %
% 
% 
% C1_css = E1 / sum(D_ss(index_1_css)) ;
% C2_css = E2 / sum(D_ss(index_2_css)) ;
% C3_css = E3 / sum(D_ss(index_3_css));
% C4_css = E4 / sum(D_ss(index_4_css)) ;
% C5_css = E5 / sum(D_ss(index_5_css));
% 
% 
% percentile_1cv_reit = 100*(C1_c -C1_css)/ C1_css    ;
% percentile_2cv_reit = 100*(C2_c -  C2_css) / C2_css;
% percentile_3cv_reit = 100*(C3_c -  C3_css) / C3_css;
% percentile_4cv_reit = 100*(C4_c -  C4_css) / C4_css;
% percentile_5cv_reit = 100*(C5_c -  C5_css) / C5_css;
% 
% percentile_1cl_reit = (C1_c -C1_css)    ;
% percentile_2cl_reit = (C2_c -  C2_css) ;
% percentile_3cl_reit = (C3_c -  C3_css) ;
% percentile_4cl_reit = (C4_c -  C4_css) ;
% percentile_5cl_reit = (C5_c -  C5_css) ;
% 
% 
% 
% X=1:1:(Length);
% figure(33)
% plot(X,percentile_1cv_reit,X,percentile_2cv_reit,X, percentile_3cv_reit,X, percentile_4cv_reit,X, percentile_5cv_reit);
% legend('bottom 7 percent','2','3','4','top')
% title('1. Percent evo of consumption after an MP shock with respect to ss value');
% saveas(figure(33),'percent_consum_decile_reiter.png')
% 
% 
% 
% X=1:1:(Length);
% figure(35)
% plot(X,percentile_1cl_reit,X,percentile_2cl_reit,X, percentile_3cl_reit,X, percentile_4cl_reit,X, percentile_5cl_reit);
% legend('bottom 7 percent','2','3','4','top')
% title('1. Level ev consumption after an MP shock with respect to ss value');
% saveas(figure(35),'cl_decile_reiter.png')
% 
% 
% X=1:1:(Length);
% figure(36)
% plot(X,percentile_1cl_reit);
% legend('bottom 7 percent')
% title('1. Level evo consumption after an MP shock with respect to ss value');
% saveas(figure(36),'cl_decile_reiter.png')
% 
% save tofigreiter
% 
% 
% load tofigcomptrunc
% 
% X=1:1:(Length);
% figure(44)
% plot(X,percentile_1cv_trunc,X,percentile_2cv_trunc,X, percentile_3cv_trunc,X, percentile_4cv_trunc,X, percentile_5cv_trunc);
% legend('bottom 7 percent','2','3','4','top')
% title('1. Evolution of log consumption after an MP shock with respect to ss value');
% saveas(figure(44), 'cons_decile_cc.png')
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% formatSpec = '%20.12f';
% 
% file = ['code_dynare_Trunc_Comp','.mod'];
% fid = fopen(file, 'w') ;
% 
% 
% str = ['lv',' = [']; fprintf(fid, str);
%         for h = 1:Nbin
%            str = ['l',num2str(h),'_eps',' ']; fprintf(fid, str);
%         end
% str = ['];']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1; 
% 
% 
% lv = [l1_eps l2_eps l3_eps l4_eps l5_eps l6_eps l7_eps l8_eps l9_eps l10_eps l11_eps l12_eps l13_eps l14_eps l15_eps l16_eps l17_eps l18_eps l19_eps l20_eps l21_eps l22_eps l23_eps l24_eps l25_eps l26_eps l27_eps l28_eps l29_eps l30_eps l31_eps l32_eps l33_eps l34_eps l35_eps l36_eps l37_eps l38_eps l39_eps l40_eps l41_eps l42_eps l43_eps l44_eps l45_eps l46_eps l47_eps l48_eps l49_eps l50_eps l51_eps l52_eps l53_eps l54_eps l55_eps l56_eps l57_eps l58_eps l59_eps ]; % equation 256
% 
% 
% l_ss = zeros(Nbin,1) ;
% 
% for h = 1:Nbin
%     l_ss(h) = eco.lp(h)/Sp(h) ; % This is the value per capita at ss 
% end
% 
% lss = repmat(l_ss', Length,1) ; %' I need it at each period 
% ldev = lv + lss ;  % Here I recover variables in level 
% 
% % steady state 
% 
% 
% %%%%%%%%%%%%%%%% MEASURE 3 : log difference between the 90th and 10th percentile of the cross sectonal distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % steady state
% 
% Here. I look at labor supply with decile of consumption at steady stqte 
% 
% 
% 
% index_1_css = idx_ss(ind_1_css); % find the point zhere 10 and 90 are located at ss
% index_2_css = idx_ss(ind_2_css); % find the point zhere 10 and 90 are located at ss
% index_3_css = idx_ss(ind_3_css);
% index_4_css = idx_ss(ind_4_css);
% index_5_css = idx_ss(ind_5_css);
% 
% 
% 
% for i=1:Length
%     D = ldev(i,:);
%     E1 = Sp(index_1_css)'* D(index_1_css)' ; %
%     E2 = Sp(index_2_css)'* D(index_2_css)' ; %
%     E3 = Sp(index_3_css)'* D(index_3_css)' ; %
%     E4 = Sp(index_4_css)'* D(index_4_css)' ; %
%     E5 = Sp(index_5_css)'* D(index_5_css)' ; %
%     C1_l(i,:) = E1 / sum(Sp(index_1_css)); % conso moyenne par tete des 10% qui ont le moins
%     C2_l(i,:) = E2 / sum(Sp(index_2_css)); % conso moyenne par tete des 10% qui ont le moins
%     C3_l(i,:) = E3 / sum(Sp(index_3_css)); % conso moyenne par tete des 10% qui ont le plus plus
%     C4_l(i,:) = E4 / sum(Sp(index_4_css)); % conso moyenne par tete des 10% qui ont le moins
%     C5_l(i,:) = E5 / sum(Sp(index_5_css)); % conso moyenne par tete des 10% qui ont le plus plus
% end    
% 
% 
% 
% E1 = Sp(index_1_css)'* l_ss(index_1_css) ; %'
% E2 = Sp(index_2_css)'* l_ss(index_2_css) ; %'
% E3 = Sp(index_3_css)'* l_ss(index_3_css);  %'
% E4 = Sp(index_4_css)'* l_ss(index_4_css) ; %'
% E5 = Sp(index_5_css)'* l_ss(index_5_css);  %'
% 
% 
% 
% C1_lss = E1 / sum(Sp(index_1_css)) ;
% C2_lss = E2 / sum(Sp(index_2_css)) ;
% C3_lss = E3 / sum(Sp(index_3_css));
% C4_lss = E4 / sum(Sp(index_4_css)) ;
% C5_lss = E5 / sum(Sp(index_5_css));
% 
% 
% percentile_1lv = log(C1_l) -  log(C1_lss) 
% percentile_2lv = log(C2_l) -  log(C2_lss) 
% percentile_3lv = log(C3_l) -  log(C3_lss) 
% percentile_4lv = log(C4_l) -  log(C4_lss) 
% percentile_5lv = log(C5_l) -  log(C5_lss) 
% 
% 
% 
% 
% X=1:1:(Length);
% figure(18)
% plot(X,percentile_1lv,X,percentile_2lv,X, percentile_3lv,X, percentile_4lv,X, percentile_5lv);
% legend('bottom 7 percent','2','3','4','top')
% title('1. Evolution of log labor earnings after an MP shock with respect to ss value');
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%
% lv = [l1_eps l2_eps l3_eps l4_eps l5_eps l6_eps l7_eps l8_eps l9_eps l10_eps l11_eps l12_eps l13_eps l14_eps l15_eps l16_eps l17_eps l18_eps l19_eps l20_eps l21_eps l22_eps l23_eps l24_eps l25_eps l26_eps l27_eps l28_eps l29_eps l30_eps l31_eps l32_eps l33_eps l34_eps l35_eps l36_eps l37_eps l38_eps l39_eps l40_eps l41_eps l42_eps l43_eps l44_eps l45_eps l46_eps l47_eps l48_eps l49_eps l50_eps l51_eps l52_eps l53_eps l54_eps l55_eps l56_eps l57_eps l58_eps l59_eps ]; % equation 256
% 
% 
% l_ss = zeros(Nbin,1) ;
% 
% for h = 1:Nbin
%     l_ss(h) = eco.lp(h)/Sp(h) ; % This is the value per capita at ss 
% end
% 
% lss = repmat(l_ss', Length,1) ; %' I need it at each period 
% ldev = lv + lss ;  % Here I recover variables in level 
% 
% 
% for h=1:Nbin
%     ldev_y(:,h) =  ldev(:,h)*eco.ytype(h) ;  % Here I recover the labor supply times productivity by history 
% end
% 
% wdev = w_eps + eco.w ;  % Here I recover the variable in level 
% 
% for i=1:Length
%     ldev_yw(i,:) =  ldev_y(i,:)*wdev(i) ;  % Here I recover the labor supply times productivity times wage by t 
% end    
% 
% 
% 
% 
% 
% index_1_css = idx_ss(ind_1_css); % find the point zhere 10 and 90 are located at ss
% index_2_css = idx_ss(ind_2_css); % find the point zhere 10 and 90 are located at ss
% index_3_css = idx_ss(ind_3_css);
% index_4_css = idx_ss(ind_4_css);
% index_5_css = idx_ss(ind_5_css);
% 
% 
% 
% 
% for i=1:Length
%     D = ldev_yw(i,:);
%     E1 = Sp(index_1_css)'* D(index_1_css)' ; %
%     E2 = Sp(index_2_css)'* D(index_2_css)' ; %
%     E3 = Sp(index_3_css)'* D(index_3_css)' ; %
%     E4 = Sp(index_4_css)'* D(index_4_css)' ; %
%     E5 = Sp(index_5_css)'* D(index_5_css)' ; %
%     C1_l(i,:) = E1 / sum(Sp(index_1_css)); % conso moyenne par tete des 10% qui ont le moins
%     C2_l(i,:) = E2 / sum(Sp(index_2_css)); % conso moyenne par tete des 10% qui ont le moins
%     C3_l(i,:) = E3 / sum(Sp(index_3_css)); % conso moyenne par tete des 10% qui ont le plus plus
%     C4_l(i,:) = E4 / sum(Sp(index_4_css)); % conso moyenne par tete des 10% qui ont le moins
%     C5_l(i,:) = E5 / sum(Sp(index_5_css)); % conso moyenne par tete des 10% qui ont le plus plus
% end    
% 
% 
% 
% % steady state 
% for h=1:Nbin
%     llss(h) =  l_ss(h)*eco.ytype(h)*eco.w ; 
% end
% 
% llss = llss' ; %'
% 
% 
% 
% E1 = Sp(index_1_css)'* llss(index_1_css) ; %'
% E2 = Sp(index_2_css)'* llss(index_2_css) ; %'
% E3 = Sp(index_3_css)'* llss(index_3_css);  %'
% E4 = Sp(index_4_css)'* llss(index_4_css) ; %'
% E5 = Sp(index_5_css)'* llss(index_5_css);  %'
% 
% 
% 
% C1_lss = E1 / sum(Sp(index_1_css)) ;
% C2_lss = E2 / sum(Sp(index_2_css)) ;
% C3_lss = E3 / sum(Sp(index_3_css));
% C4_lss = E4 / sum(Sp(index_4_css)) ;
% C5_lss = E5 / sum(Sp(index_5_css));
% 
% 
% percentile_1lv = 100*(C1_l -  C1_lss)/C1_lss
% percentile_2lv = 100*(C2_l -  C2_lss)/ C2_lss
% percentile_3lv = 100*(C3_l -  C3_lss)/ C3_lss
% percentile_4lv = 100*(C4_l -  C4_lss)/ C4_lss
% percentile_5lv = 100*(C5_l -  C5_lss)/ C5_lss
% 
% percentile_1ll = (C1_l -  C1_lss)
% percentile_2ll = (C2_l -  C2_lss)
% percentile_3ll = (C3_l -  C3_lss)
% percentile_4ll = (C4_l -  C4_lss)
% percentile_5ll = (C5_l -  C5_lss)
% 
% 
% 
% X=1:1:(Length);
% figure(18)
% plot(X,percentile_1lv,X,percentile_2lv,X, percentile_3lv,X, percentile_4lv,X, percentile_5lv);
% legend('bottom 7 percent','2','3','4','top')
% title('1. Evolution percentage of wyl after an MP shock decile of consum');
% saveas(figure(18), 'learning_decile_percent_trunc.png')
% 
% 
% 
% X=1:1:(Length);
% figure(77)
% plot(X,percentile_1ll,X,percentile_2ll,X, percentile_3ll,X, percentile_4ll,X, percentile_5ll);
% legend('bottom 7 percent','2','3','4','top')
% title('1. Evolution level of wyl after an MP shock decile of consum');
% saveas(figure(77), 'learning_decile_level_trunc.png')
% 
% X=1:1:(Length);
% figure(7)
% plot(X,percentile_1ll);
% legend('bottom 7 percent')
% title('1. Evolution of wyl after an MP shock decile of consum');
% saveas(figure(7), 'learning_decile_level1_trunc.png')
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% formatSpec = '%20.12f';
% 
% file = ['code_dynare_Trunc_Comp','.mod'];
% fid = fopen(file, 'w') ;
% 
% 
% 
% str = ['at',' = [']; fprintf(fid, str);
%         for h = 1:Nbin
%            str = ['at',num2str(h),'_eps',' ']; fprintf(fid, str);
%         end
% str = ['];']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1; 
% 
% 
% at = [at1_eps at2_eps at3_eps at4_eps at5_eps at6_eps at7_eps at8_eps at9_eps at10_eps at11_eps at12_eps at13_eps at14_eps at15_eps at16_eps at17_eps at18_eps at19_eps at20_eps at21_eps at22_eps at23_eps at24_eps at25_eps at26_eps at27_eps at28_eps at29_eps at30_eps at31_eps at32_eps at33_eps at34_eps at35_eps at36_eps at37_eps at38_eps at39_eps at40_eps at41_eps at42_eps at43_eps at44_eps at45_eps at46_eps at47_eps at48_eps at49_eps at50_eps at51_eps at52_eps at53_eps at54_eps at55_eps at56_eps at57_eps at58_eps at59_eps ]; % equation 255
% % c = wly + (1+r)at + TT - a
% 
% at_ss = zeros(Nbin,1) ;
% 
% for h = 1:Nbin
%     at_ss(h) = eco.abp(h)/Sp(h) ; % This is the value per capita at ss 
% end
% 
% 
% atss = repmat(at_ss', Length,1) ; %'
% atdev = at + atss ;  % Here I recover the variable in level 
% 
% 
% rdev = r_eps + (eco.R-1) ;  % Here I recover the variable in level 
% 
% 
% for i=1:Length
%     at_dev(i,:) =  atdev(i,:)*(1+rdev(i)) ;  % Here I recover the return 
% end    
% 
% 
% 
% for h=1:Nbin
%     a_t_ss(h) =  at_ss(h)*((1+eco.R-1)) ; 
% end
% 
% 
% 
% 
% for i=1:Length
%     D = at_dev(i,:);
%     E1 = Sp(index_1_css)'* D(index_1_css)' ; %
%     E2 = Sp(index_2_css)'* D(index_2_css)' ; %
%     E3 = Sp(index_3_css)'* D(index_3_css)' ; %
%     E4 = Sp(index_4_css)'* D(index_4_css)' ; %
%     E5 = Sp(index_5_css)'* D(index_5_css)' ; %
%     C1_at(i,:) = E1 / sum(Sp(index_1_css)); % conso moyenne par tete des 10% qui ont le moins
%     C2_at(i,:) = E2 / sum(Sp(index_2_css)); % conso moyenne par tete des 10% qui ont le moins
%     C3_at(i,:) = E3 / sum(Sp(index_3_css)); % conso moyenne par tete des 10% qui ont le plus plus
%     C4_at(i,:) = E4 / sum(Sp(index_4_css)); % conso moyenne par tete des 10% qui ont le moins
%     C5_at(i,:) = E5 / sum(Sp(index_5_css)); % conso moyenne par tete des 10% qui ont le plus plus
% end    
% 
% 
% 
% a_t_ss = a_t_ss' %'
% 
% 
% 
% E11 = Sp(index_1_css)'* a_t_ss(index_1_css) ; %'
% E22 = Sp(index_2_css)'* a_t_ss(index_2_css) ; %'
% E33 = Sp(index_3_css)'* a_t_ss(index_3_css);  %'
% E44 = Sp(index_4_css)'* a_t_ss(index_4_css) ; %'
% E55 = Sp(index_5_css)'* a_t_ss(index_5_css);  %'
% 
% 
% A1_atss = E11 / sum(Sp(index_1_css)) ;
% A2_atss = E22 / sum(Sp(index_2_css)) ;
% A3_atss = E33 / sum(Sp(index_3_css));
% A4_atss = E44 / sum(Sp(index_4_css)) ;
% A5_atss = E55 / sum(Sp(index_5_css));
% 
% 
% percentile_1atv_trunc = 100*(C1_at -A1_atss)/ A1_atss
% percentile_2atv_trunc = 100*(C2_at -  A2_atss) / A2_atss
% percentile_3atv_trunc = 100*(C3_at -  A3_atss) / A3_atss
% percentile_4atv_trunc = 100*(C4_at -  A4_atss) / A4_atss
% percentile_5atv_trunc = 100*(C5_at -  A5_atss) / A5_atss
% 
% 
% percentile_1atl_trunc = (C1_at -A1_atss)
% percentile_2atl_trunc = (C2_at -  A2_atss) 
% percentile_3atl_trunc = (C3_at -  A3_atss) 
% percentile_4atl_trunc = (C4_at -  A4_atss) 
% percentile_5atl_trunc = (C5_at -  A5_atss) 
% 
% 
% X=1:1:(Length);
% figure(19)
% plot(X,percentile_1atv,X,percentile_2atv,X, percentile_3atv,X, percentile_4atv,X, percentile_5atv);
% legend('bottom 7 percent','2','3','4','top')
% title('1. Evolution of (1+r)at an MP shock with respect to ss value');
% saveas(figure(19), 'atearning_decile_percent_trunc.png')
% 
% X=1:1:(Length);
% figure(81)
% plot(X,percentile_1atl_trunc,X,percentile_2atl_trunc,X, percentile_3atl_trunc,X, percentile_4atl_trunc,X, percentile_5atl_trunc);
% legend('bottom 7 percent','2','3','4','top')
% title('1. Evolution Level of (1+r)at');
% saveas(figure(81), 'atearning_decile_percent_trunc.png')
% 
% X=1:1:(Length);
% figure(81)
% plot(X,percentile_1atl_trunc);
% legend('bottom 7 percent','2','3','4','top')
% title('1. Evolution Level1 of (1+r)at');
% saveas(figure(81), 'atearning_decile1_percent_trunc.png')
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% earn_dev = at_dev + ldev_yw
% 
% earn_ss =  a_t_ss + llss
% 
% 
% for i=1:Length
%     D = earn_dev(i,:);
%     E1 = Sp(index_1_css)'* D(index_1_css)' ; %
%     E2 = Sp(index_2_css)'* D(index_2_css)' ; %
%     E3 = Sp(index_3_css)'* D(index_3_css)' ; %
%     E4 = Sp(index_4_css)'* D(index_4_css)' ; %
%     E5 = Sp(index_5_css)'* D(index_5_css)' ; %
%     EA1_at(i,:) = E1 / sum(Sp(index_1_css)); % conso moyenne par tete des 10% qui ont le moins
%     EA2_at(i,:) = E2 / sum(Sp(index_2_css)); % conso moyenne par tete des 10% qui ont le moins
%     EA3_at(i,:) = E3 / sum(Sp(index_3_css)); % conso moyenne par tete des 10% qui ont le plus plus
%     EA4_at(i,:) = E4 / sum(Sp(index_4_css)); % conso moyenne par tete des 10% qui ont le moins
%     EA5_at(i,:) = E5 / sum(Sp(index_5_css)); % conso moyenne par tete des 10% qui ont le plus plus
% end    
% 
% 
% 
% 
% EA11 = Sp(index_1_css)'* earn_ss(index_1_css) ; %'
% EA22 = Sp(index_2_css)'* earn_ss(index_2_css) ; %'
% EA33 = Sp(index_3_css)'* earn_ss(index_3_css);  %'
% EA44 = Sp(index_4_css)'* earn_ss(index_4_css) ; %'
% EA55 = Sp(index_5_css)'* earn_ss(index_5_css);  %'
% 
% 
% 
% 
% 
% EA1_atss = EA11 / sum(Sp(index_1_css)) ;
% EA2_atss = EA22 / sum(Sp(index_2_css)) ;
% EA3_atss = EA33 / sum(Sp(index_3_css));
% EA4_atss = EA44 / sum(Sp(index_4_css)) ;
% EA5_atss = EA55 / sum(Sp(index_5_css));
% 
% percentile_1ea_trunc = 100*(EA1_at -EA1_atss)/ EA1_atss
% percentile_2ea_trunc = 100*(EA2_at -  EA2_atss) / EA2_atss
% percentile_3ea_trunc = 100*(EA3_at -  EA3_atss) / EA3_atss
% percentile_4ea_trunc = 100*(EA4_at -  EA4_atss) / EA4_atss
% percentile_5ea_trunc = 100*(EA5_at -  EA5_atss) / EA5_atss
% 
% percentile_1eal_trunc = EA1_at -EA1_atss
% percentile_2eal_trunc = (EA2_at -  EA2_atss) 
% percentile_3eal_trunc = (EA3_at -  EA3_atss) 
% percentile_4eal_trunc = (EA4_at -  EA4_atss) 
% percentile_5eal_trunc = (EA5_at -  EA5_atss) 
% 
% 
% X=1:1:(Length);
% figure(23)
% plot(X,percentile_1ea_trunc,X,percentile_2ea_trunc,X, percentile_3ea_trunc,X, percentile_4ea_trunc,X, percentile_5ea_trunc);
% legend('bottom 7 percent','2','3','4','top')
% title('1. earning percent_trunc');
% saveas(figure(23), 'earning_decile_level_trunc.png')
% 
% X=1:1:(Length);
% figure(21)
% plot(X,percentile_1eal_trunc,X,percentile_2eal_trunc,X, percentile_3eal_trunc,X, percentile_4eal_trunc,X, percentile_5eal_trunc);
% legend('bottom 7 percent','2','3','4','top')
% title('1. earning level  after an MP shock with respect to ss value');
% saveas(figure(21), 'earning_decile_level_trunc.png')
% 
% 
% X=1:1:(Length);
% figure(24)
% plot(X,percentile_1eal_trunc);
% legend('bottom 7 percent')
% title('1. earning level  after an MP shock with respect to ss value');
% saveas(figure(24), 'earning_decile_level1_trunc.png')
% 
% 
% ###############
% 
% 
% percentile_1cv_trunc = 100*(C1_c -C1_css)/ C1_css
% percentile_2cv_trunc = 100*(C2_c -  C2_css) / C2_css
% percentile_3cv_trunc = 100*(C3_c -  C3_css) / C3_css
% percentile_4cv_trunc = 100*(C4_c -  C4_css) / C4_css
% percentile_5cv_trunc = 100*(C5_c -  C5_css) / C5_css
% 
% 
% 
% X=1:1:(Length);
% figure(94)
% plot(X,C1_c -C1_css);
% legend('bottom 7 percent','2','3','4','top')
% title('1. Evolution of consumption after an MP shock with respect to ss value');
% 
% 
% 
% X=1:1:(Length);
% figure(944)
% plot(X,EA1_at -EA1_atss);
% legend('bottom 7 percent','2','3','4','top')
% title('1. Evolution of earn after an MP shock with respect to ss value');
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % REITER
% %%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%
% % wt*lt*y
% %%%%%%%%%%%%%%%
% 
% 
% 
% formatSpec  = '%16.12f';
% file = ['code_dynare_Comp_Reiter','.mod'];
% fid = fopen(file, 'w') ;
% 
% str = ['lv',' = [']; fprintf(fid, str);
%         for h = 1:Ntot
%            str = ['l',num2str(h),'_eps',' ']; fprintf(fid, str);
%         end
% str = [']']; fprintf(fid, str); 
% 
% lv = [l1_eps l2_eps l3_eps l4_eps l5_eps l6_eps l7_eps l8_eps l9_eps l10_eps l11_eps l12_eps l13_eps l14_eps l15_eps l16_eps l17_eps l18_eps l19_eps l20_eps l21_eps l22_eps l23_eps l24_eps l25_eps l26_eps l27_eps l28_eps l29_eps l30_eps l31_eps l32_eps l33_eps l34_eps l35_eps l36_eps l37_eps l38_eps l39_eps l40_eps l41_eps l42_eps l43_eps l44_eps l45_eps l46_eps l47_eps l48_eps l49_eps l50_eps l51_eps l52_eps l53_eps l54_eps l55_eps l56_eps l57_eps l58_eps l59_eps l60_eps l61_eps l62_eps l63_eps l64_eps l65_eps l66_eps l67_eps l68_eps l69_eps l70_eps l71_eps l72_eps l73_eps l74_eps l75_eps l76_eps l77_eps l78_eps l79_eps l80_eps l81_eps l82_eps l83_eps l84_eps l85_eps l86_eps l87_eps l88_eps l89_eps l90_eps l91_eps l92_eps l93_eps l94_eps l95_eps l96_eps l97_eps l98_eps l99_eps l100_eps l101_eps l102_eps l103_eps l104_eps l105_eps l106_eps l107_eps l108_eps l109_eps l110_eps l111_eps l112_eps l113_eps l114_eps l115_eps l116_eps l117_eps l118_eps l119_eps l120_eps l121_eps l122_eps l123_eps l124_eps l125_eps l126_eps l127_eps l128_eps l129_eps l130_eps l131_eps l132_eps l133_eps l134_eps l135_eps l136_eps l137_eps l138_eps l139_eps l140_eps l141_eps l142_eps l143_eps l144_eps l145_eps l146_eps l147_eps l148_eps l149_eps l150_eps l151_eps l152_eps l153_eps l154_eps l155_eps l156_eps l157_eps l158_eps l159_eps l160_eps l161_eps l162_eps l163_eps l164_eps l165_eps l166_eps l167_eps l168_eps l169_eps l170_eps l171_eps l172_eps l173_eps l174_eps l175_eps l176_eps l177_eps l178_eps l179_eps l180_eps l181_eps l182_eps l183_eps l184_eps l185_eps l186_eps l187_eps l188_eps l189_eps l190_eps l191_eps l192_eps l193_eps l194_eps l195_eps l196_eps l197_eps l198_eps l199_eps l200_eps l201_eps l202_eps l203_eps l204_eps l205_eps l206_eps l207_eps l208_eps l209_eps l210_eps l211_eps l212_eps l213_eps l214_eps l215_eps l216_eps l217_eps l218_eps l219_eps l220_eps l221_eps l222_eps l223_eps l224_eps l225_eps l226_eps l227_eps l228_eps l229_eps l230_eps l231_eps l232_eps l233_eps l234_eps l235_eps l236_eps l237_eps l238_eps l239_eps l240_eps l241_eps l242_eps l243_eps l244_eps l245_eps l246_eps l247_eps l248_eps l249_eps l250_eps l251_eps l252_eps l253_eps l254_eps l255_eps l256_eps l257_eps l258_eps l259_eps l260_eps l261_eps l262_eps l263_eps l264_eps l265_eps l266_eps l267_eps l268_eps l269_eps l270_eps l271_eps l272_eps l273_eps l274_eps l275_eps l276_eps l277_eps l278_eps l279_eps l280_eps l281_eps l282_eps l283_eps l284_eps l285_eps l286_eps l287_eps l288_eps l289_eps l290_eps l291_eps l292_eps l293_eps l294_eps l295_eps l296_eps l297_eps l298_eps l299_eps l300_eps l301_eps l302_eps l303_eps l304_eps l305_eps l306_eps l307_eps l308_eps l309_eps l310_eps l311_eps l312_eps l313_eps l314_eps l315_eps l316_eps l317_eps l318_eps l319_eps l320_eps l321_eps l322_eps l323_eps l324_eps l325_eps l326_eps l327_eps l328_eps l329_eps l330_eps l331_eps l332_eps l333_eps l334_eps l335_eps l336_eps l337_eps l338_eps l339_eps l340_eps l341_eps l342_eps l343_eps l344_eps l345_eps l346_eps l347_eps l348_eps l349_eps l350_eps l351_eps l352_eps l353_eps l354_eps l355_eps l356_eps l357_eps l358_eps l359_eps l360_eps l361_eps l362_eps l363_eps l364_eps l365_eps l366_eps l367_eps l368_eps l369_eps l370_eps l371_eps l372_eps l373_eps l374_eps l375_eps l376_eps l377_eps l378_eps l379_eps l380_eps l381_eps l382_eps l383_eps l384_eps l385_eps l386_eps l387_eps l388_eps l389_eps l390_eps l391_eps l392_eps l393_eps l394_eps l395_eps l396_eps l397_eps l398_eps l399_eps l400_eps l401_eps l402_eps l403_eps l404_eps l405_eps l406_eps l407_eps l408_eps l409_eps l410_eps l411_eps l412_eps l413_eps l414_eps l415_eps l416_eps l417_eps l418_eps l419_eps l420_eps l421_eps l422_eps l423_eps l424_eps l425_eps l426_eps l427_eps l428_eps l429_eps l430_eps l431_eps l432_eps l433_eps l434_eps l435_eps l436_eps l437_eps l438_eps l439_eps l440_eps l441_eps l442_eps l443_eps l444_eps l445_eps l446_eps l447_eps l448_eps l449_eps l450_eps l451_eps l452_eps l453_eps l454_eps l455_eps l456_eps l457_eps l458_eps l459_eps l460_eps l461_eps l462_eps l463_eps l464_eps l465_eps l466_eps l467_eps l468_eps l469_eps l470_eps l471_eps l472_eps l473_eps l474_eps l475_eps l476_eps l477_eps l478_eps l479_eps l480_eps l481_eps l482_eps l483_eps l484_eps l485_eps l486_eps l487_eps l488_eps l489_eps l490_eps l491_eps l492_eps l493_eps l494_eps l495_eps l496_eps l497_eps l498_eps l499_eps l500_eps ];
% 
% 
% for h = 1:Ntot
%     L_ss(h) = Lpol(h) ; % This is the value per capita at ss 
% end
% 
% lss = repmat(L_ss, Length,1) ; % I need it at each period 
% ldev = lv + lss ;  % Here I recover variables in level 
% 
% 
% 
% for h=1:Ntot
%     l_y_reit(:,h) =  ldev(:,h)*ys(yv(h)) ;  % Here I recover the labor supply times productivity by history 
% end
% 
% wdev = w_eps + w ;  % Here I recover the variable in level 
% 
% for i=1:Length
%     l_yw(i,:) =  l_y_reit(i,:)*wdev(i) ;  % Here I recover the labor supply times productivity times wage by t 
% end    
% 
% % steady state 
% for h=1:Ntot
%     llss(h) =  L_ss(h)*ys(yv(h))*w ; 
% end
% 
% llss = llss' %'
% 
% 
% for i=1:Length
%     D = l_yw(i,:);
%     L1 = D_ss(index_1_css)'* D(index_1_css)' ; %
%     L2 = D_ss(index_2_css)'* D(index_2_css)' ; %
%     L3 = D_ss(index_3_css)'* D(index_3_css)' ; %
%     L4 = D_ss(index_4_css)'* D(index_4_css)' ; %
%     L5 = D_ss(index_5_css)'* D(index_5_css)' ; %
%     L1_reit(i,:) = L1 / sum(D_ss(index_1_css)); % conso moyenne par tete des 10% qui ont le moins
%     L2_reit(i,:) = L2 / sum(D_ss(index_2_css)); % conso moyenne par tete des 10% qui ont le moins
%     L3_reit(i,:) = L3 / sum(D_ss(index_3_css)); % conso moyenne par tete des 10% qui ont le plus plus
%     L4_reit(i,:) = L4 / sum(D_ss(index_4_css)); % conso moyenne par tete des 10% qui ont le moins
%     L5_reit(i,:) = L5 / sum(D_ss(index_5_css)); % conso moyenne par tete des 10% qui ont le plus plus
% end    
% 
% 
% L11 = D_ss(index_1_css)'* llss(index_1_css) ; %'
% L22 = D_ss(index_2_css)'* llss(index_2_css) ; %'
% L33 = D_ss(index_3_css)'* llss(index_3_css);  %'
% L44 = D_ss(index_4_css)'* llss(index_4_css) ; %'
% L55 = D_ss(index_5_css)'* llss(index_5_css);  %'
% 
% 
% L1_ss_reit = L11 / sum(D_ss(index_1_css)) ;
% L2_ss_reit = L22 / sum(D_ss(index_2_css)) ;
% L3_ss_reit = L33 / sum(D_ss(index_3_css));
% L4_ss_reit = L44 / sum(D_ss(index_4_css)) ;
% L5_ss_reit = L55 / sum(D_ss(index_5_css));
% 
% 
% percentile_1lv_reit = 100*(L1_reit -L1_ss_reit)/ L1_ss_reit   ;
% percentile_2lv_reit = 100*(L2_reit -  L2_ss_reit) / L2_ss_reit;
% percentile_3lv_reit = 100*(L3_reit -  L3_ss_reit) / L3_ss_reit;
% percentile_4lv_reit = 100*(L4_reit -  L4_ss_reit) / L4_ss_reit;
% percentile_5lv_reit = 100*(L5_reit -  L5_ss_reit) / L5_ss_reit;
% 
% percentile_1ll_reit = (L1_reit -L1_ss_reit)   ;
% percentile_2ll_reit = (L2_reit -  L2_ss_reit) ;
% percentile_3ll_reit = (L3_reit -  L3_ss_reit) ;
% percentile_4ll_reit = (L4_reit -  L4_ss_reit) ;
% percentile_5ll_reit = (L5_reit -  L5_ss_reit) ;
% 
% 
% X=1:1:(Length);
% figure(90)
% plot(X,percentile_1lv_reit,X,percentile_2lv_reit,X, percentile_3lv_reit,X, percentile_4lv_reit,X, percentile_5lv_reit);
% legend('bottom 7 percent','2','3','4','top')
% title('1. deviation w*l*y percent ');
% saveas(figure(90), 'learning_deviation_percen_reit.png')
% 
% X=1:1:(Length);
% figure(91)
% plot(X,percentile_1ll_reit,X,percentile_2ll_reit,X, percentile_3ll_reit,X, percentile_4ll_reit,X, percentile_5ll_reit);
% legend('bottom 7 percent','2','3','4','top')
% title('1. deviation w*l*y level ');
% saveas(figure(91), 'learning_deviation_level_reit.png')
% 
% X=1:1:(Length);
% figure(91)
% plot(X,percentile_1ll_reit);
% legend('bottom 7 percent')
% title('1. deviation w*l*y level ');
% saveas(figure(91), 'learning_deviation_level1_reit.png')
% 
% 
% %%%%%%%%%%%%%%%
% % (1+r)*at
% %%%%%%%%%%%%%%%
% 
% aGrid = aGrid(1:500)
% 
% a = repmat(aGrid', Length,1 ) ; %'
% a= a'   %'
% 
% at_ss_reit = zeros(Ntot,1) ;
% 
% for h = 1:Ntot
%     at_ss_reit(h) = aGrid(h) ; % This is the value per capita at ss 
% end
% 
% 
% atss = repmat(at_ss_reit', Length,1) ; %'
% 
% 
% rdev_reit = r_eps + (R-1) ;  % Here I recover the variable in level 
% 
% 
% 
% 
% at_dev_reit = zeros(Length,Ntot); 
% for i=1:Length
%     at_dev_reit(i,:) =  a(i,:).*(1+rdev_reit(i)) ;  % Here I recover the return 
% end    
% 
% 
% for h=1:Ntot
%     a_t_ss_reit(h) =  aGrid(h)*((1+R-1)) ; 
% end
% 
% 
% 
% for i=1:Length
%     D = at_dev_reit(i,:);
%     AT1 = D_ss(index_1_css)'* D(index_1_css)' ; %
%     AT2 = D_ss(index_2_css)'* D(index_2_css)' ; %
%     AT3 = D_ss(index_3_css)'* D(index_3_css)' ; %
%     AT4 = D_ss(index_4_css)'* D(index_4_css)' ; %
%     AT5 = D_ss(index_5_css)'* D(index_5_css)' ; %
%     AT1_reit(i,:) = AT1 / sum(D_ss(index_1_css)); % conso moyenne par tete des 10% qui ont le moins
%     AT2_reit(i,:) = AT2 / sum(D_ss(index_2_css)); % conso moyenne par tete des 10% qui ont le moins
%     AT3_reit(i,:) = AT3 / sum(D_ss(index_3_css)); % conso moyenne par tete des 10% qui ont le plus plus
%     AT4_reit(i,:) = AT4 / sum(D_ss(index_4_css)); % conso moyenne par tete des 10% qui ont le moins
%     AT5_reit(i,:) = AT5 / sum(D_ss(index_5_css)); % conso moyenne par tete des 10% qui ont le plus plus
% end    
% 
% 
% a_t_ss_reit  = a_t_ss_reit' %'
% 
% AT11 = D_ss(index_1_css)'* a_t_ss_reit(index_1_css) ; %'
% AT22 = D_ss(index_2_css)'* a_t_ss_reit(index_2_css) ; %'
% AT33 = D_ss(index_3_css)'* a_t_ss_reit(index_3_css);  %'
% AT44 = D_ss(index_4_css)'* a_t_ss_reit(index_4_css) ; %'
% AT55 = D_ss(index_5_css)'* a_t_ss_reit(index_5_css);  %'
% 
% 
% AT1_ss_reit = AT11 / sum(D_ss(index_1_css)) ;
% AT2_ss_reit = AT22 / sum(D_ss(index_2_css)) ;
% AT3_ss_reit = AT33 / sum(D_ss(index_3_css));
% AT4_ss_reit = AT44 / sum(D_ss(index_4_css)) ;
% AT5_ss_reit = AT55 / sum(D_ss(index_5_css));
% 
% 
% percentile_1atv_reit = 100*(AT1_reit -   AT1_ss_reit)/ AT1_ss_reit   ;
% percentile_2atv_reit = 100*(AT2_reit -  AT2_ss_reit) / AT2_ss_reit;
% percentile_3atv_reit = 100*(AT3_reit -  AT3_ss_reit) / AT3_ss_reit;
% percentile_4atv_reit = 100*(AT4_reit -  AT4_ss_reit) / AT4_ss_reit;
% percentile_5atv_reit = 100*(AT5_reit -  AT5_ss_reit) / AT5_ss_reit;
% 
% percentile_1atl_reit = (AT1_reit -  AT1_ss_reit)   ;
% percentile_2atl_reit = (AT2_reit -  AT2_ss_reit) ;
% percentile_3atl_reit = (AT3_reit -  AT3_ss_reit) ;
% percentile_4atl_reit = (AT4_reit -  AT4_ss_reit) ;
% percentile_5atl_reit = (AT5_reit -  AT5_ss_reit) ;
% 
% 
% X=1:1:(Length);
% figure(90)
% plot(X,percentile_1atv_reit,X,percentile_2atv_reit,X, percentile_3atv_reit,X, percentile_4atv_reit,X, percentile_5atv_reit);
% legend('bottom 7 percent','2','3','4','top')
% title('1. deviation at(1+r) percent ');
% saveas(figure(90), 'atearning_deviation_percen_reit.png')
% 
% X=1:1:(Length);
% figure(91)
% plot(X,percentile_1atl_reit,X,percentile_2atl_reit,X, percentile_3atl_reit,X, percentile_4atl_reit,X, percentile_5atl_reit);
% legend('bottom 7 percent','2','3','4','top')
% title('1. deviation at(1+r) level ');
% saveas(figure(91), 'atearning_deviation_level_reit.png')
% 
% X=1:1:(Length);
% figure(91)
% plot(X,percentile_1atl_reit);
% legend('bottom 7 percent')
% title('1. deviation at(1+r) level ');
% saveas(figure(91), 'atearning_deviation_level1_reit.png')
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% earn_rei_dev = at_dev_reit + l_yw ; 
% earn_rei_ss = (a_t_ss_reit + llss)' %'
% 
% 
% 
% for i=1:Length
%     D = earn_rei_dev(i,:);
%     EA1 = D_ss(index_1_css)'* D(index_1_css)' ; %
%     EA2 = D_ss(index_2_css)'* D(index_2_css)' ; %
%     EA3 = D_ss(index_3_css)'* D(index_3_css)' ; %
%     EA4 = D_ss(index_4_css)'* D(index_4_css)' ; %
%     EA5 = D_ss(index_5_css)'* D(index_5_css)' ; %
%     EA1_at_reit(i,:) = EA1 / sum(D_ss(index_1_css)); % conso moyenne par tete des 10% qui ont le moins
%     EA2_at_reit(i,:) = EA2 / sum(D_ss(index_2_css)); % conso moyenne par tete des 10% qui ont le moins
%     EA3_at_reit(i,:) = EA3 / sum(D_ss(index_3_css)); % conso moyenne par tete des 10% qui ont le plus plus
%     EA4_at_reit(i,:) = EA4 / sum(D_ss(index_4_css)); % conso moyenne par tete des 10% qui ont le moins
%     EA5_at_reit(i,:) = EA5 / sum(D_ss(index_5_css)); % conso moyenne par tete des 10% qui ont le plus plus
% end    
% 
% 
% EA11 = D_ss(index_1_css)'* earn_rei_ss(index_1_css) ; %'
% EA22 = D_ss(index_2_css)'* earn_rei_ss(index_2_css) ; %'
% EA33 = D_ss(index_3_css)'* earn_rei_ss(index_3_css);  %'
% EA44 = D_ss(index_4_css)'* earn_rei_ss(index_4_css) ; %'
% EA55 = D_ss(index_5_css)'* earn_rei_ss(index_5_css);  %'
% 
% 
% EA1_atss_reit = EA11 / sum(D_ss(index_1_css)) ;
% EA2_atss_reit = EA22 / sum(D_ss(index_2_css)) ;
% EA3_atss_reit = EA33 / sum(D_ss(index_3_css));
% EA4_atss_reit = EA44 / sum(D_ss(index_4_css)) ;
% EA5_atss_reit = EA55 / sum(D_ss(index_5_css));
% 
% 
% percentile_1cv_reit = 100*(EA1_at_reit -EA1_atss_reit)/ EA1_atss_reit   ;
% percentile_2cv_reit = 100*(EA2_at_reit -  EA2_atss_reit) / EA2_atss_reit;
% percentile_3cv_reit = 100*(EA3_at_reit -  EA3_atss_reit) / EA3_atss_reit;
% percentile_4cv_reit = 100*(EA4_at_reit -  EA4_atss_reit) / EA4_atss_reit;
% percentile_5cv_reit = 100*(EA5_at_reit -  EA5_atss_reit) / EA5_atss_reit;
% 
% percentile_1i_reit = (EA1_at_reit -EA1_atss_reit)   ;
% percentile_2i_reit = (EA2_at_reit -  EA2_atss_reit) ;
% percentile_3i_reit = (EA3_at_reit -  EA3_atss_reit) ;
% percentile_4i_reit = (EA4_at_reit -  EA4_atss_reit) ;
% percentile_5i_reit = (EA5_at_reit -  EA5_atss_reit) ;
% 
% 
% X=1:1:(Length);
% figure(90)
% plot(X,percentile_1cv_reit,X,percentile_2cv_reit,X, percentile_3cv_reit,X, percentile_4cv_reit,X, percentile_5cv_reit);
% legend('bottom 7 percent','2','3','4','top')
% title('1. deviation (1+r)at + w*l*y percent ');
% saveas(figure(90), 'income_deviation_percen_reit.png')
% 
% X=1:1:(Length);
% figure(91)
% plot(X,percentile_1i_reit,X,percentile_2i_reit,X, percentile_3i_reit,X, percentile_4i_reit,X, percentile_5i_reit);
% legend('bottom 7 percent','2','3','4','top')
% title('1. deviation (1+r)at + w*l*y level ');
% saveas(figure(91), 'income_deviation_level_reit.png')
% 
% X=1:1:(Length);
% figure(91)
% plot(X,percentile_1i_reit);
% legend('bottom 7 percent','2','3','4','top')
% title('1. deviation (1+r)at + w*l*y level ');
% saveas(figure(91), 'income_deviation_level1_reit.png')
% 
% 
% 
% 
% E1 = D_ss(index_1_css)'* C_ss(index_1_css)' ; %'
% E2 = D_ss(index_2_css)'* C_ss(index_2_css)' ; %'
% E3 = D_ss(index_3_css)'* C_ss(index_3_css)';  %'
% E4 = D_ss(index_4_css)'* C_ss(index_4_css)' ; %'
% E5 = D_ss(index_5_css)'* C_ss(index_5_css)';  %
% 
% 
% 
% C1_css = E1 / sum(D_ss(index_1_css)) ;
% C2_css = E2 / sum(D_ss(index_2_css)) ;
% C3_css = E3 / sum(D_ss(index_3_css));
% C4_css = E4 / sum(D_ss(index_4_css)) ;
% C5_css = E5 / sum(D_ss(index_5_css));
% 
% 
% percentile_1cv = log(C1_c) -  log(C1_css) 
% percentile_2cv = log(C2_c) -  log(C2_css) 
% percentile_3cv = log(C3_c) -  log(C3_css) 
% percentile_4cv = log(C4_c) -  log(C4_css) 
% percentile_5cv = log(C5_c) -  log(C5_css) 
% 
% 
% percentile_1cv_reit = 100*(C1_c -C1_css)/ C1_css
% percentile_2cv_reit = 100*(C2_c -  C2_css) / C2_css
% percentile_3cv_reit = 100*(C3_c -  C3_css) / C3_css
% percentile_4cv_reit = 100*(C4_c -  C4_css) / C4_css
% percentile_5cv_reit = 100*(C5_c -  C5_css) / C5_css
% 
% 
% percentile_1i_reit = 100*(C1_c -C1_css) 
% percentile_2i_reit = 100*(C2_c -  C2_css) 
% percentile_3i_reit = 100*(C3_c -  C3_css) 
% percentile_4i_reit = 100*(C4_c -  C4_css) 
% percentile_5i_reit = 100*(C5_c -  C5_css) 
% 
% 
% 
% 
% 
% X=1:1:(Length);
% figure(33)
% plot(X,percentile_1cv_reit,X,percentile_2cv_reit,X, percentile_3cv_reit,X, percentile_4cv_reit,X, percentile_5cv_reit);
% legend('bottom 7 percent','2','3','4','top')
% title('1. Evolution of log consumption after an MP shock with respect to ss value');
% saveas(figure(33),'decile_reiter_7.png')
% 
% 
% 
% 
