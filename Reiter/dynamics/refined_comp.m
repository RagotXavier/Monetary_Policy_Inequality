
clear 
warning('off','all');

calib = 'Refined' 

filename_root = 'todynare_Comp_' ;
fileExtension = '.mat';
FileNameTot = [filename_root, calib, fileExtension];
load (FileNameTot)


Nbin  = eco.Nbin;    %total number of bins
PI    = eco.Matab;
Sp    = eco.Sp;
K     = eco.A;
Y     = K^alpha*Ltot^(1-alpha);
rho_u        = 0.95;

P= ones(Nbin,1);
P(eco.indcc) = 0;

Length = 100; %length of the simulation


formatSpec = '%20.12f';

        
file = ['code_dynare_Comp_', calib ,'.mod'];
fid = fopen(file, 'w') ;



str=['var\n'] ; fprintf(fid, str);

str=['u r w  A K TFP Ltot Ctot GDP\n'];fprintf(fid, str);
str = ['ut Kt Ct wt rt GDPt Lt Zt \n'] ; fprintf(fid, str); 


for p = 1:Nbin
    str = ['a',num2str(p),' '] ; fprintf(fid,str);     %  end-of-period
    str = ['at',num2str(p),' '] ; fprintf(fid, str);   % beginning-of-period  
    str = ['c',num2str(p),' '] ; fprintf(fid, str);    
    str = ['l',num2str(p),' '] ; fprintf(fid, str);   % beginning-of-period  
    if mod(p,10)==0; %multiple de 10
        str = ['\n'] ;         fprintf(fid, str);
    end;
end;


str = ['; \n\n'] ; fprintf(fid, str);


str = ['varexo eps; \n\n'] ; fprintf(fid, str);
str = ['parameters\n'] ; fprintf(fid, str);

str = ['beta alpha phi abar delta gamma chi rho_u G ;\n\n'] ; fprintf(fid, str);


str = ['alpha','   = ',num2str(alpha,formatSpec),';\n']; fprintf(fid, str); 
str = ['beta','    = ',num2str(beta,formatSpec),';\n']; fprintf(fid, str); 
str = ['phi','     = ',num2str(phi,formatSpec),';\n']; fprintf(fid, str); 
str = ['abar','   = ',num2str(min(eco.aep)),';\n']; fprintf(fid, str); 
str = ['delta','   = ',num2str(delta,formatSpec),';\n']; fprintf(fid, str);
str = ['gamma','   = ',num2str(gamma,formatSpec),';\n']; fprintf(fid, str);
str = ['rho_u','   = ',num2str(rho_u,formatSpec),';\n']; fprintf(fid, str); %HERE
str = ['G','       = ',num2str(G,formatSpec),';\n']; fprintf(fid, str);
str = ['chi','     = ',num2str(chi,formatSpec),';\n']; fprintf(fid, str);



%%
equa = 1; % Numbering equations
tol = 10^-6;  % threshold for transition, to avoid to consider very small transitions

%~ str = ['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'] ; fprintf(fid, str); 
str = ['\n\n'] ; fprintf(fid, str);
str = ['model;\n\n']; fprintf(fid, str);

for h = 1:Nbin  
    nsi = eco.ytype(h);


    %%%%%%%% Budget constraint  
                str = ['c',num2str(h),' =' 'w*l',num2str(h),'*',num2str(nsi,formatSpec),' + (1+r)*','at',num2str(h),'-a', num2str(h),';']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    
    %%%%%%%% Definition of at   
        str = ['at',num2str(h),' = 10^-10 ']; fprintf(fid, str);
        for hi = 1:Nbin
            if (PI(h,hi))>tol   % hi to h
              str = ['+ ',num2str(Sp(hi)*PI(h,hi)/Sp(h),formatSpec),'*a',num2str(hi),'(-1)']; fprintf(fid, str);
            end
        end
        str = [';']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    

    %%%%%%%% Euler equation
        if P(h)==1 % not constrained
            str = [num2str(eco.xsiue(h),formatSpec),'*(c',num2str(h),')^-gamma =', num2str(eco.Resh(h),formatSpec),' + beta*(1+r(+1))*(10^-10']; fprintf(fid, str);
            for hp = 1:Nbin                
                    if (PI(hp,h)*eco.xsiue(hp))>tol % h to hp
                        str = ['+ ',num2str((PI(hp,h))*eco.xsiue(hp),formatSpec),'*(c',num2str(hp),'(+1))^-gamma']; fprintf(fid, str);
                    end
                    % if mod(hi,10)==0;str = ['\n'] ;fprintf(fid, str); end;                
            end
            str = [');']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        else
            str = ['a',num2str(h),' = ',num2str(eco.aep(h)/Sp(h),formatSpec),';']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        end
        
        
%%%%%%%% Labor supply

str = [num2str(eco.xsiv1(h),formatSpec) '*(1/chi)*(l',num2str(h),'^(1/phi)) = w','*',num2str(eco.xsiu1(h),formatSpec),'*',num2str(nsi,formatSpec),'* c',num2str(h),'^(-gamma);']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
end;


%%%%%%%%%%%%%%%% Aggregate Constraints %%%%%%%%%%%%%%%%
    
%%%%%%%% Financial market eq
str = ['A = K;'] ; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;


%%%%%%%%%%%%%%%% Interest Rate Constraints:

%%%%%%%% Factor prices : 

str = ['1 = (1/(alpha*TFP))*(r + delta)*(K(-1)/Ltot)^(1-alpha);'] ; fprintf(fid, str);      str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['w = (1 - alpha)*TFP*(K(-1)/Ltot)^alpha;'] ; fprintf(fid, str);   str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1; 


%%%%%%%% Production side 
str = ['GDP = TFP*K(-1)^alpha*Ltot^(1-alpha);'] ; fprintf(fid, str);      str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

% Shock
str = ['TFP = 1 - u;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['u = rho_u*u(-1) + eps;'] ; fprintf(fid, str);    str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;


%%%%%%%% Total consumption
str = ['\n'] ;         fprintf(fid, str);
str = ['Ctot ='] ; fprintf(fid, str);    
for h = 1:Nbin
        str = ['+',num2str(Sp(h),formatSpec),'*c',num2str(h)] ; fprintf(fid, str);  
         if mod(h,5)==0; %multiple de 10
        str = ['\n'] ;         fprintf(fid, str);
    end; 
end;
str = [';']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1; 



%%%%%%%% Total labor
    str = ['Ltot = '] ; fprintf(fid, str);    
for h = 1:Nbin
        str = ['+', num2str(Sp(h),formatSpec),'*l',num2str(h),'*',num2str(eco.ytype(h),formatSpec)] ; fprintf(fid, str);    
    if mod(p,10)==10; tr = ['\n'] ;fprintf(fid, str);
    end;
end;
str = [';']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;



%%%%%%%% Saving agents
    str = ['A = '] ; fprintf(fid, str);    
for h = 1:Nbin
        str = ['+', num2str(Sp(h),formatSpec),'*a',num2str(h)] ; fprintf(fid, str);    
    if mod(p,10)==10; tr = ['\n'] ;fprintf(fid, str);
    end;
end;
str = [';']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

%%%%%%%% variables in percent deviations 


str = ['ut = 100*u;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Kt = 100*(K/',num2str(K,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Ct = 100*(Ctot/steady_state(Ctot)-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['wt = 100*(w/',num2str(eco.w,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['rt = 100*(r - ',num2str(eco.R-1,formatSpec),');'] ; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['GDPt = 100*(GDP/',num2str(Y,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Zt = -ut;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Lt = 100*(Ltot/steady_state(Ltot)-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;


str = ['\n'] ;         fprintf(fid, str);
str = ['end;\n\n'] ; fprintf(fid, str);  



%~ str = ['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'] ; fprintf(fid, str);    
str = ['\n\n'] ; fprintf(fid, str);    
%%


str = ['steady_state_model;\n\n'] ; fprintf(fid, str);   
for h = 1:Nbin
    str = ['a',num2str(h),' = ',num2str(eco.aep(h)/Sp(h),formatSpec),';\n'] ; fprintf(fid, str);
    str = ['at',num2str(h),' = ',num2str(eco.abp(h)/Sp(h),formatSpec),';\n'] ; fprintf(fid, str);
    str = ['c',num2str(h),' = ',num2str(eco.cp(h)/Sp(h),formatSpec),';\n'] ; fprintf(fid, str);    
    str = ['l',num2str(h),' = ',num2str(eco.lp(h)/Sp(h),formatSpec),';\n'] ; fprintf(fid, str);    
end;
 
str = ['r = ', num2str(eco.R-1,formatSpec),';\n'] ; fprintf(fid, str);
str = ['w = ', num2str(eco.w,formatSpec),';\n'] ; fprintf(fid, str);
str = ['K = ', num2str(K,formatSpec),';\n'] ; fprintf(fid, str);
str = ['Ltot =' num2str(Ltot,formatSpec),' ;\n'] ; fprintf(fid, str); 
str = ['TFP = ', num2str(1),';\n'] ; fprintf(fid, str);
str = ['A = ', num2str(eco.A,formatSpec),';\n'] ; fprintf(fid, str);
str = ['u = ', num2str(0),';\n'] ; fprintf(fid, str);

str = ['Ctot = ', num2str(Ctot,formatSpec),';\n'] ; fprintf(fid, str);
str = ['GDP = ', num2str(Y,formatSpec),';\n'] ; fprintf(fid, str);

str = ['ut = 0;\n'] ; fprintf(fid, str);  
str = ['Kt = 0;\n'] ; fprintf(fid, str);  
str = ['Ct  = 0;\n'] ; fprintf(fid, str);  
str = ['wt = 0;\n'] ; fprintf(fid, str);  
str = ['rt = 0;\n'] ; fprintf(fid, str); 
str = ['Zt = 0;\n'] ; fprintf(fid, str);
str = ['Lt = 0;\n'] ; fprintf(fid, str); 
str = ['GDPt = 0;\n'] ; fprintf(fid, str);

str = ['end;\n\n'] ; fprintf(fid, str);


% str = ['resid;\n\n'] ; fprintf(fid, str);




%  
%  str = ['steady;\n\n'] ; fprintf(fid, str);
%  str = ['check;\n\n'] ; fprintf(fid, str);


%%%%%%%%% for the simulation
str = ['shocks;\n var eps = 0.00001;\n end;\n'];fprintf(fid, str);
str = ['stoch_simul (order=1,irf=',num2str(Length),',nograph,noprint,periods=10000) GDP Ctot K Ltot r w  u Zt GDPt Ct Kt rt wt Lt ut;']; fprintf(fid, str);



command = ['dynare ',file,' noclearall']
eval(command);

    


irf = oo_.irfs;
ut =irf.ut_eps;
Zt =irf.Zt_eps;
GDPt=irf.GDPt_eps;
Ct =irf.Ct_eps;
Kt=irf.Kt_eps;
rt =irf.rt_eps;
wt=irf.wt_eps;
Lt=irf.Lt_eps;
    
filename_graph = sprintf('todiff_Comp_%s.mat', calib); 
save(filename_graph, "ut","Zt", "GDPt", "Ct", "Kt", "rt", "wt", "Lt") 

    

list = oo_.var_list;
var_simu = diag(oo_.var);
std_simu = var_simu.^(1/2);
mean_simu     = oo_.mean ;
std_norm = (std_simu./mean_simu)*100;
index_u= find(list == "u");
std_norm(index_u) = std_simu(index_u) *100 ;   
index_r = find(list == "r");
std_norm(index_r) = std_simu(index_r) *100 ;
std_norm = round(std_norm,4);
mean_var = round(oo_.mean,4);

index_Y = find(list == "GDP");
index_K= find(list == "K");
index_L= find(list == "Ltot");
index_C= find(list == "Ctot");
index_r= find(list == "r");
index_w= find(list == "w");
 
 

corr_c_y = oo_.var(index_Y,index_C)./(std_simu(index_Y)*std_simu(index_C) );% Y, c
auto_corr_y_y1 =  oo_.autocorr{1, 1}(index_Y,index_Y)  ; % Y, Y(-1)
corr_list = ["corr(C,Y)" "corr(Y,Y(-1))"] ;
corr = round([corr_c_y auto_corr_y_y1],4) ;
corr_list = cellstr(corr_list) 
list_table     = [list(index_Y) ; list(index_C) ; list(index_L); list(index_K); list(index_r) ; list(index_w) ; list(index_u)] ;
 mean_simu = [mean_var(index_Y) ; mean_var(index_C) ; mean_var(index_L); mean_var(index_K); mean_var(index_r) ; mean_var(index_w) ; mean_var(index_u)] ;
 std_norm_table = [std_norm(index_Y) ; std_norm(index_C) ; std_norm(index_L); std_norm(index_K); std_norm(index_r) ; std_norm(index_w) ; std_norm(index_u)] ;
  




filename_momemt  = sprintf('moments_Comp_%s.mat', calib); 
save(filename_momemt,"list_table","mean_simu", "std_norm_table", "corr_list", "corr") 


list_endo_name =   M_.endo_names 
index_Y_endo = find(list_endo_name == "GDP");
index_C_endo = find(list_endo_name == "Ctot");
index_w_endo = find(list_endo_name == "w");
index_r_endo = find(list_endo_name == "r");
index_K_endo = find(list_endo_name == "K");
index_L_endo = find(list_endo_name == "Ltot");
GDP1 = oo_.endo_simul(index_Y_endo,:) ; 
C1 = oo_.endo_simul(index_C_endo,:) ; 
w1 = oo_.endo_simul(index_w_endo,:) ; 
r1 = oo_.endo_simul(index_r_endo,:) ; 
K1 = oo_.endo_simul(index_K_endo,:) ; 
L1 = oo_.endo_simul(index_L_endo,:) ; 
    
    

mean_refined = mean_simu ; 
save('to_code_difference_refined.mat',"list_table",'C1',"GDP1","w1" ,"r1" ,"K1","L1", 'mean_refined') ; 
    






