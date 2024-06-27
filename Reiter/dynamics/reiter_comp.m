
clear
clc
warning('off','all');
load todynare_Comp_Reiter.mat 


Ntot = na*ny; 
Ny = ny; % idiosyncratic shocks
thres       = 0; % threshold to consider 0 (not to have too many equations)
Y = K^alpha*Ltot^(1-alpha); % GDP
I = delta*K; % investment

%% AGGREGATE SHOCKS IN THE ECONOMY
rho_u       = 0.95;
sigma_u     = 0.00312;

file = ['code_dynare_Comp_Reiter','.mod'];
fid = fopen(file, 'w') ;

formatSpec  = '%16.12f';

Length = 100; %length of the IRFs

%% DEFINING THE ENDOGENOUS VARIABLES 

str = ['var\n'] ; fprintf(fid, str);
str = ['r w Z K GDP u Pop I C Ctot Ltot \n'] ; fprintf(fid, str);  % variables
str = ['rt wt Zt Kt GDPt ut It Ct Ctott Lt \n'] ; fprintf(fid, str);

for i = 1:Ntot
    str     = ['ap',num2str(i),' '] ; fprintf(fid,str);        
    str     = ['c',num2str(i),' '] ; fprintf(fid, str);  
    str     = ['S',num2str(i),' '] ; fprintf(fid, str); 
    str     = ['ww',num2str(i),' '] ; fprintf(fid, str); 
    str     = ['l',num2str(i),' '] ; fprintf(fid, str); 
    if mod(i,10)==0; % multiple de 10
        str = ['\n'] ;         fprintf(fid, str);
    end;    
end;

str         = ['; \n\n'] ; fprintf(fid, str);

str         = ['varexo eps; \n\n'] ; fprintf(fid, str);
str         = ['parameters\n'] ; fprintf(fid, str);

str         = ['beta alpha abar delta gamma rho_u chi phi  \n'] ; fprintf(fid, str);
str = ['; \n\n'] ; fprintf(fid, str);

abar        = abar;
 
str         = ['alpha','   = ',num2str(alpha),';\n']; fprintf(fid, str); 
str         = ['beta ','   = ',num2str(beta),';\n']; fprintf(fid, str); 
str         = ['abar ','   = ',num2str(abar),';\n']; fprintf(fid, str); 
str         = ['gamma','   = ',num2str(gamma),';\n']; fprintf(fid, str); 
str         = ['rho_u','   = ',num2str(rho_u),';\n']; fprintf(fid, str); 
str         = ['delta','   = ',num2str(delta),';\n']; fprintf(fid, str);
str         = ['phi','   = ',num2str(phi),';\n']; fprintf(fid, str);
str         = ['chi','   = ',num2str(chi,formatSpec),';\n']; fprintf(fid, str);

str = ['\n\n'] ; fprintf(fid, str);
yv = ytype; % vector of productivity type from 1 to ns*na (for budget constraint)
yvv = [1:1:ny]; % vector of productivity type from 1 to ns 
Th = 10^-8;
indc    = find(Resid(:) > 0.0001); % credit constraind agents
indnc = find(Resid(:) <= 0.0001); 



%% EQUATIONS

str     = ['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'] ; fprintf(fid, str); 
str     = ['\n\n'] ; fprintf(fid, str);    

str     = ['model;\n\n']; fprintf(fid, str);
equa    = 1;

%%%%%%%%%%%%%%%%% LABOR SUPPLY %%%%%%%%%%%%%%%%%%%%
for i = 1:Ntot
    str = ['(1/chi)*l',num2str(i),'^(1/phi)',' =  w','*',num2str(ys(yv(i)),formatSpec),'*c',num2str(i),'^(-gamma)',';']; fprintf(fid, str);
    str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1; 
end;

%%%%%%%%%%%%%%%%% BUDGET CONSTRAINT %%%%%%%%%%%%%%%
for i = 1:Ntot
    str = ['c',num2str(i),'+ap',num2str(i), ' = w*',num2str(ys(yv(i)),formatSpec),'*l',num2str(i),'+ (1+r)*',num2str(aGrid(i),formatSpec) ,';']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;     
end;
str = ['\n\n']; fprintf(fid, str);

%%%%%%%%%%%%%%%%  WEALTH TRANSITION  %%%%%%%%%%%%%%
for i = 1:Ntot 
     str = ['ap',num2str(i),' =  ww',num2str(i),'*',num2str(aGrid(Vind(i)),formatSpec),'+ (1 -ww',num2str(i),')*',num2str(aGrid(Vind(i)+1),formatSpec),';']; fprintf(fid, str);
     str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;  
end;    
str = ['\n\n']; fprintf(fid, str);


%%%%%%%%%% EULER EQUATIONS %%%%%%%%%%%%%%%%%%%%%%%%
str     = ['%%%%%%%% equation ',num2str(equa),': Euler\n\n']; fprintf(fid, str);
 
for i = 1:length(indc)    
     str = ['ap',num2str(indc(i)),' = ',num2str(abar,formatSpec) ,';']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1; 
end;

for k = 1:length(indnc)
    i = indnc(k);
    yi = ytype(i); % current productivity status
    threse = 0;
          str = ['(c',num2str(i),' - ','0)^-gamma  = beta*(1+r(+1))*('];fprintf(fid, str);
          for j = 1:ny % next period productivity status
                  if (Trans(j,yi)>threse) % does agent i reach agent j ? effective transition toward the bin j
                     str = ['+',num2str(Trans(j,yi),formatSpec),'*( ww',num2str(i),'*c',num2str(Vind(i)+(j-1)*na ),'(+1)+(1-ww',num2str(i),')*c',num2str(Vind(i)+1+(j-1)*na),'(+1) )^-gamma']; fprintf(fid, str);                               
                  end
          end     
            str = [')+',num2str(resE(i),formatSpec),';']; fprintf(fid, str);
            str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;     
end     

%%%%%%%%%%%%%%%%%%%%%%%%%  SHARES  %%%%%%%%%%%%%%%%
 str    = [' %%%%%%%% equation ',num2str(equa),': Shares \n\n']; fprintf(fid, str);
for i = 1:na
    for j=1:ny   
        pass=0;
        ind =   i + (j-1)*na;
        str    = ['S',num2str(ind),' = ']; fprintf(fid, str); 
        
        for ip = 1:na
            for jp=1:ny
             indp =   ip + (jp-1)*na ;
               if (Vind(indp)==i)&(Trans(j,jp)>0) % goes from indp to ip (i+1 appears later)
                        str    = ['+S',num2str(indp),'(-1)*',num2str(Trans(j,jp),formatSpec),'*(',num2str(aGrid(Vind(indp)+1),formatSpec) ,'-ap',num2str(indp),'(-1))/',num2str(aGrid(Vind(indp)+1)-aGrid(Vind(indp)),formatSpec)]; fprintf(fid, str);pass=1;
               end
               if (Vind(indp)==i-1)&(Trans(j,jp)>0) % goes from indp to ip (i+1 appears later)
                        str    = ['+S',num2str(indp),'(-1)*',num2str(Trans(j,jp),formatSpec),'*(ap',num2str(indp),'(-1)-',num2str(aGrid(Vind(indp)),formatSpec) ,')/',num2str(aGrid(Vind(indp)+1)-aGrid(Vind(indp)),formatSpec)]; fprintf(fid, str);pass=1;
              end
            end
       end
         if pass==0
          str    = [num2str(D_ss(ind),formatSpec)]; fprintf(fid, str); 
        end
       str    = [';']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1; 
    end
end
  
%% AGGREGATIONS
str = [' %%%%%%%% equation ',num2str(equa),': Aggretations \n\n']; fprintf(fid, str);
str = ['Pop = ']; fprintf(fid, str);
for i = 1:Ntot
     str    = ['+S',num2str(i)]; fprintf(fid, str);
      if mod(i,10)==0; %multiple de 10
         str = ['\n'] ;         fprintf(fid, str);
     end;  
end;
str = [';'] ;         fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n\n']; fprintf(fid, str); equa = equa+1;  

str = ['K = ']; fprintf(fid, str);
for i = 1:Ntot
     str    = ['+S',num2str(i),'(-1)*ap',num2str(i),'(-1)']; fprintf(fid, str);
      if mod(i,10)==0; % multiple de 10
         str = ['\n'] ;         fprintf(fid, str);
     end;  
end;
str = [';'] ;         fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n\n']; fprintf(fid, str); equa = equa+1;  
   

str = ['Ltot = ']; fprintf(fid, str);
for i = 1:Ntot
     str    = ['+S',num2str(i),'*l',num2str(i),'*',num2str(ys(yv(i)),formatSpec)]; fprintf(fid, str);
      if mod(i,10)==0; % multiple de 10
         str = ['\n'] ;         fprintf(fid, str);
     end;  
end;
str = [';'] ;         fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n\n']; fprintf(fid, str); equa = equa+1;  


str = ['w = (1 - alpha)*Z*(K/Ltot)^(alpha);']; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   
str = ['(r + delta)/(alpha*Z) = (K/Ltot)^(alpha-1);']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   


str = ['Z = 1-u;']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   
str = ['u = rho_u*u(-1)+eps;']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   
str = ['GDP = Z*K^alpha*Ltot^(1-alpha);'] ; fprintf(fid, str);      str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

%%%%%%%% total consumption
str = ['\n'] ;         fprintf(fid, str);
str = ['Ctot ='] ; fprintf(fid, str);    
for h = 1:Ntot
        str = ['+S',num2str(h),'*c',num2str(h)] ; fprintf(fid, str);  
         if mod(h,5)==0; %multiple de 10
        str = ['\n'] ;         fprintf(fid, str);
    end; 
end;
str = [';']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;


str = ['I =  K(+1) - ( 1 - delta)*K;']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   
str = ['C  =  GDP - I;']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;   


%%%%%%%%%% VARIABLES IN PERCENT DEVIATIONS %%%%%%%%%%%%%%%%%%%%%%%%
str = ['ut = 100*u;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Kt = 100*(K/',num2str(K,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Lt = 100*(Ltot/',num2str(Ltot,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

str = ['Ctott = 100*(Ctot/',num2str(Ctot,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;



str = ['Ct = 100*(C/steady_state(C)-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['wt = 100*(w/',num2str(w,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['rt = 100*(r - ',num2str(R-1,formatSpec),');'] ; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Zt = -ut;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['GDPt = 100*(GDP/',num2str(Y,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['It = 100*(I/',num2str(I,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;



str = ['end;\n\n'] ; fprintf(fid, str);  

%% STEADY STATE

str = ['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'] ; fprintf(fid, str);   
str = ['\n\n'] ; fprintf(fid, str);    
%str = ['initval;\n\n'] ; fprintf(fid, str);  

str = ['steady_state_model;\n\n'] ; fprintf(fid, str);   

for i = 1:Ntot
    str = ['ap',num2str(i),' = ',num2str(apol(i),formatSpec),';\n'] ; fprintf(fid, str);       
    str = ['c',num2str(i), ' = ',num2str(Cpol(i),formatSpec),';\n'] ; fprintf(fid, str); 
    str = ['S',num2str(i), ' = ',num2str(D_ss(i),formatSpec),';\n'] ; fprintf(fid, str);    
     str = ['ww',num2str(i), ' = ',num2str(Wind(i),formatSpec),';\n'] ; fprintf(fid, str);    
    str = ['l',num2str(i), ' = ',num2str(Lpol(i),formatSpec),';\n'] ; fprintf(fid, str);    
end;


str = ['r = ', num2str(R-1,formatSpec),';\n'] ; fprintf(fid, str);
str = ['K = ', num2str(K,formatSpec),';\n'] ; fprintf(fid, str);
str = ['Ltot = ', num2str(Ltot,formatSpec),';\n'] ; fprintf(fid, str);

str = ['w = ', num2str(w,formatSpec),';\n'] ; fprintf(fid, str);
str = ['Z  = 1;\n\n']; fprintf(fid, str);
str = ['u  = 0;\n\n']; fprintf(fid, str);
str = ['GDP = Z*K^alpha*Ltot^(1-alpha);\n']; fprintf(fid, str);
str = ['Pop = 1;\n']; fprintf(fid, str);


str = ['I =  delta*K;\n']; fprintf(fid, str); 
str = ['C =  GDP - I;']; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Ctot = ', num2str(Ctot,formatSpec),';\n'] ; fprintf(fid, str);


str = ['ut = 0;\n'] ; fprintf(fid, str);
str = ['Kt = 0;\n'] ; fprintf(fid, str);
str = ['Lt = 0;\n'] ; fprintf(fid, str);

str = ['Ct  = 0;\n'] ; fprintf(fid, str);
str = ['Ctott  = 0;\n'] ; fprintf(fid, str);
str = ['wt = 0;\n'] ; fprintf(fid, str);
str = ['rt = 0;\n'] ; fprintf(fid, str);
str = ['Zt = 0;\n'] ; fprintf(fid, str);
str = ['GDPt = 0;\n'] ; fprintf(fid, str);
str = ['It = 0;\n'] ; fprintf(fid, str);

str = ['end;\n\n'] ; fprintf(fid, str);     


%%%%%%%%% for the simulation
str = ['shocks;\n var eps = 0.00001;\n end;\n'];fprintf(fid, str);
str = ['stoch_simul (order=1,irf=',num2str(Length),',noprint,nograph,periods=10000)  GDP Ctot K Ltot r w u ut Zt GDPt Ct Ctott Kt rt wt Lt;']; fprintf(fid, str);


%% CALLING DYNARE
command = ['dynare ',file,' noclearall'];
eval(command);


irf = oo_.irfs;
ut =irf.ut_eps;
Zt =irf.Zt_eps;
GDPt=irf.GDPt_eps;
Ct =irf.Ct_eps;
Ctott =irf.Ctott_eps;
Kt=irf.Kt_eps;
rt =irf.rt_eps;
wt=irf.wt_eps;
Lt=irf.Lt_eps;

save todiff_Reiter ut Zt GDPt Ct Ctott Kt rt wt Lt


%%%%%%%%% for the simulation

list = oo_.var_list ;
var_simu = diag(oo_.var) ;
std_simu = var_simu.^(1/2) ;
mean_simu     = oo_.mean ;
std_norm = (std_simu./mean_simu)*100 ; % normalization std by mean 
index_u= find(list == "u");
std_norm(index_u) = std_simu(index_u) *100 ; 
index_r = find(list == "r");
std_norm(index_r) = std_simu(index_r) *100 ; 
std_norm = round(std_norm,4) ;
mean_var     = round(oo_.mean,4) ; 

    

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
list_table     = [list(index_Y) ; list(index_C) ; list(index_L); list(index_K); list(index_r) ; list(index_w) ; list(index_u)] ;
mean_simu = [mean_var(index_Y) ; mean_var(index_C) ; mean_var(index_L); mean_var(index_K); mean_var(index_r) ; mean_var(index_w) ; mean_var(index_u)] ;
std_norm_table = [std_norm(index_Y) ; std_norm(index_C) ; std_norm(index_L); std_norm(index_K); std_norm(index_r) ; std_norm(index_w) ; std_norm(index_u)] ;
  
corr_list = cellstr(corr_list) ;
filename = sprintf('moments_Comp_Reiter.mat');
save(filename,"list_table","mean_simu", "std_norm_table", "corr_list", "corr") ;


list_endo_name =   M_.endo_names ;
index_Y_endo = find(list_endo_name == "GDP");
index_C_endo = find(list_endo_name == "Ctot");
index_w_endo = find(list_endo_name == "w");
index_r_endo = find(list_endo_name == "r");
index_K_endo = find(list_endo_name == "K");
index_L_endo = find(list_endo_name == "Ltot");
C2   = oo_.endo_simul(index_C_endo,:) ; 
GDP2 = oo_.endo_simul(index_Y_endo,:) ; 
w2   = oo_.endo_simul(index_w_endo,:) ; 
r2   = oo_.endo_simul(index_r_endo,:) ; 
K2   = oo_.endo_simul(index_K_endo,:) ; 
L2   = oo_.endo_simul(index_L_endo,:) ; 


mean_reiter = mean_simu; 
save('to_code_difference_reiter.mat',"list",'C2',"GDP2","w2" ,"r2" ,"K2","L2","mean_reiter") ; 

