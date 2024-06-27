clear

calib = 'taylor' ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%         Output parameters         %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LengthIRF  = 20;        %length of the IRFs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%       Beginning of the code       %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


warning('off','all');
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
    warning('off', 'Octave:array-as-logical');
end

% loading steady-state values and parameters
% load steady_state_dynare.mat

load todynare_SP_baseline.mat


alpha = alphaa ; 
beta = betaa ;
gamma = gammaa ; 

% initialization
Nbin  = eco.Nbin;
PI    = eco.Matab;
Sp    = eco.Sp;
K     = eco.A-eco.B;
Y     = K^alpha*Ltot^(1-alpha);
rho_u = 0.95;
phi_taylor   = 1.5 ;  %Taylor Rule 

P = ones(Nbin,1);
P(eco.indcc) = 0;

formatSpec = '%25.20f'; %print format

% Initialization of return results
EcoM = 1; 
Economy = 1

rv_taylor = zeros(EcoM,LengthIRF+1);  % path of nominal rate
rnv_taylor = zeros(EcoM,LengthIRF+1); % path of real rate
Yv_taylor = zeros(EcoM,LengthIRF+1);
pipv_taylor = zeros(EcoM,LengthIRF+1);
Bv_taylor = zeros(EcoM,LengthIRF+1);
Cv_taylor = zeros(EcoM,LengthIRF+1);
taulv_taylor = zeros(EcoM,LengthIRF+1);
taukv_taylor = zeros(EcoM,LengthIRF+1);
Lv_taylor = zeros(EcoM,LengthIRF+1);
piwv_taylor = zeros(EcoM,LengthIRF+1);
Uv_taylor = zeros(EcoM,LengthIRF+1);
Tv_taylor = zeros(EcoM,LengthIRF+1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Specifications of the 3 economies %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% common fiscal rule parameters
coefBL = 4.0;
coefBu = 8.5;

coefl  = 0.;   % reaction of labor tax to debt
coefk  = 0.;
coefkL = 0.;
coeflL = 0.;


   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%     opening the mod file     %%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
file = ['code_dynare_', calib,'.mod'];
fid = fopen(file, 'w') ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     variables and params     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str = ['var\n'] ; fprintf(fid, str);
str = ['u r w FK FL A B K TFP Ltot tauk taul Ctot GDP TT BsYt TsYt piw piwt \n'];fprintf(fid, str);
str = ['ut Kt Ct Bt wt rt rnt GDPt PIt rtildeK wtilde Tt tkt tlt Lt Zt zeta M PI RBN \n'] ; fprintf(fid, str);


for p = 1:Nbin
    str = ['a',num2str(p),' '] ; fprintf(fid,str);
    str = ['at',num2str(p),' '] ; fprintf(fid, str);
    str = ['c',num2str(p),' '] ; fprintf(fid, str);
    str = ['l',num2str(p),' '] ; fprintf(fid, str);
    if mod(p,10)==0;
        str = ['\n'] ;         fprintf(fid, str);
    end
end


str = ['; \n\n'] ; fprintf(fid, str);
str = ['varexo eps; \n\n'] ; fprintf(fid, str);
str = ['parameters\n'] ; fprintf(fid, str);
str = ['coefl coefk coefBu coefBL beta alpha phi abar delta gamma chi rho_u G kappa epss coefkL coeflL phi_taylor;\n\n'] ; fprintf(fid, str);
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%        initialisation        %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str = ['alpha   = ',num2str(alpha,formatSpec),';\n']; fprintf(fid, str);
str = ['beta    = ',num2str(beta,formatSpec),';\n']; fprintf(fid, str);
str = ['phi     = ',num2str(phi,formatSpec),';\n']; fprintf(fid, str);
str = ['abar    = ',num2str(min(eco.aep)),';\n']; fprintf(fid, str);
str = ['delta   = ',num2str(delta,formatSpec),';\n']; fprintf(fid, str);
str = ['gamma   = ',num2str(gamma,formatSpec),';\n']; fprintf(fid, str);
str = ['rho_u   = ',num2str(rho_u,formatSpec),';\n']; fprintf(fid, str);
str = ['G       = ',num2str(G,formatSpec),';\n']; fprintf(fid, str);
str = ['chi     = ',num2str(chi,formatSpec),';\n']; fprintf(fid, str);
str = ['kappa   = ',num2str(kappa,formatSpec),';\n']; fprintf(fid, str);
str = ['epss   = ',num2str(epss,formatSpec),';\n']; fprintf(fid, str);
str = ['coefl  = ',num2str(coefl,formatSpec),';\n']; fprintf(fid, str);
str = ['coefk  = ',num2str(coefk,formatSpec),';\n']; fprintf(fid, str);
str = ['coefBu = ',num2str(coefBu,formatSpec),';\n']; fprintf(fid, str);
str = ['coefBL = ',num2str(coefBL,formatSpec),';\n']; fprintf(fid, str);
str = ['coefkL = ',num2str(coefkL,formatSpec),';\n']; fprintf(fid, str);
str = ['coeflL = ',num2str(coeflL,formatSpec),';\n']; fprintf(fid, str);
str = ['phi_taylor','   = ',num2str(phi_taylor,formatSpec),';\n']; fprintf(fid, str); 




equa = 1; % Numbering equations
tol = 10^(-16);  % threshold for transition, to avoid to consider very small transitions
str = ['\n\n'] ; fprintf(fid, str);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%       model definition       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str = ['model;\n\n']; fprintf(fid, str);
for h = 1:Nbin
    nsi = eco.ytype(h);
    %%%%%%%% Budget constraint
    str = ['c',num2str(h),' =' 'w*l',num2str(h),'*',num2str(nsi,formatSpec),' + (1+r)*at',num2str(h),'-a', num2str(h),' +TT; %% equation ', num2str(equa),'\n'];
    fprintf(fid, str); equa = equa+1;

    %%%%%%%% Definition of at
    str = ['at',num2str(h),' = 10^-10 ']; fprintf(fid, str);
    for hi = 1:Nbin
        if (PI(h,hi))>tol   % hi to h
          str = ['+ ',num2str(Sp(hi)*PI(h,hi)/Sp(h),formatSpec),'*a',num2str(hi),'(-1)']; fprintf(fid, str);
        end
    end
    str = ['; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    %%%%%%%% Euler equation

    if P(h)==1 % not constrained
        str = [num2str(eco.xsiue(h),formatSpec),'*(c',num2str(h),')^-gamma =', num2str(eco.Resh(h),formatSpec),' + beta*(1+r(+1))*(10^-10']; fprintf(fid, str);
        for hp = 1:Nbin
                if (PI(hp,h)*eco.xsiue(hp))>tol % h to hp
                    str = ['+ ',num2str((PI(hp,h))*eco.xsiue(hp),formatSpec),'*(c',num2str(hp),'(+1))^-gamma']; fprintf(fid, str);
                end
                % if mod(hi,10)==0;str = ['\n'] ;fprintf(fid, str); end;
        end
        str = [');  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    else
        str = ['a',num2str(h),' = ',num2str(eco.aep(h)/Sp(h),formatSpec),';  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    end

    %%%%%%%%% Labor supply
    str = [num2str(eco.xsiv1(h),formatSpec) '*(1/chi)*(l',num2str(h),'^(1/phi)) = w*',num2str(eco.xsiu1(h),formatSpec),'*',num2str(nsi,formatSpec),...
                '* c',num2str(h),'^(-gamma); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
end


 %%%%%%%% Definition Pricing kernel M
 str = ['M = 1; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
 %%%%%%%%% Taylor Rule 
 str = ['RBN/steady_state(RBN) = PI^phi_taylor ;'] ; fprintf(fid, str);      str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
%%%%%%%%%%%%%%%% Aggregate Constraints %%%%%%%%%%%%%%%%
%%%%%%%% Definition piw
str = ['(1+piw)*w(-1) = PI*w; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
%%%%%%%% Budget contraint of the govt
str = ['B + GDP*(1-kappa*(PI-1)*(PI-1)/2) -delta*K(-1) = G + TT + B(-1) +r*A(-1) + w*Ltot; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
%%%%%%%% NK Philipps curve
str = ['GDP*M*PI*(PI-1) = ((epss-1)/kappa)*(zeta-1)*GDP*M + beta*PI(+1)*(PI(+1)-1)*GDP(+1)*M(+1); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
%%%%%%%% Definition of zeta
str = ['zeta = (1/(alpha*TFP))*(rtildeK + delta)*(K(-1)/Ltot)^(1-alpha); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
%%%%%%%% Financial market eq
str = ['A = B + K; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
%%%%%%%%%%%%%%%% Interest Rate Constraints:
str = ['(r - (1-tauk)*rtildeK)*A(-1) = (1-tauk)*(RBN(-1)/PI - 1- rtildeK)*B(-1); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['RBN/PI(+1) = (1+rtildeK(+1)); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
%%%%%%%% Factor prices :
str = ['K(-1)/Ltot = alpha*wtilde/((1-alpha)*(rtildeK+delta)); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
%%%%%%%% Production side
str = ['GDP = TFP*K(-1)^alpha*Ltot^(1-alpha); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['FL = (1 - alpha)*TFP*(K(-1)/Ltot)^alpha; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['FK = alpha*TFP*(K(-1)/Ltot)^(alpha-1)-delta; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
%%%%%%%% Aggregate shock
str = ['TFP = 1 - u;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['u = rho_u*u(-1) + eps;'] ; fprintf(fid, str);    str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
%%%%%%%% Total consumption
str = ['Ctot = 0'] ; fprintf(fid, str);
for h = 1:Nbin
    str = ['+',num2str(Sp(h),formatSpec),'*c',num2str(h)] ; fprintf(fid, str);
    if mod(h,5)==0; str = ['\n']; fprintf(fid, str);end
end;
str = ['; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
%%%%%%%% Aggregate labor supply
str = ['Ltot = 0'] ; fprintf(fid, str);
for h = 1:Nbin
        str = ['+', num2str(Sp(h),formatSpec),'*l',num2str(h),'*',num2str(eco.ytype(h),formatSpec)] ; fprintf(fid, str);
    if mod(p,10)==10; tr = ['\n'] ;fprintf(fid, str); end;
end;
str = ['; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
%%%%%%%% Aggregate savings
str = ['A = 0'] ; fprintf(fid, str);
for h = 1:Nbin
    str = ['+', num2str(Sp(h),formatSpec),'*a',num2str(h)] ; fprintf(fid, str);
    if mod(p,10)==10; tr = ['\n'] ;fprintf(fid, str);
    end;
end;
str = ['; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

%%%%%%%% Definition of w
str = ['w = (1-taul)*wtilde;   %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
%%%%%%%%%%%%%%%% Fiscal rules %%%%%%%%%%%%%%%%
str = ['taul = ',num2str(eco.tl,formatSpec),' + coefl*u + coeflL*u(-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
str = ['tauk = ',num2str(eco.tk,formatSpec),' + coefk*u + coefkL*u(-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
str = ['TT = ',num2str(eco.TT,formatSpec),'+ coefBu*u - coefBL*(B - ',num2str(eco.B,formatSpec),'); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
%%%%%%%%%%%%%%%% log-deviation variables %%%%%%%%%%%%%%%%
% There are denoted with 't': for variable x, xt = (x -x_ss)/x_ss, where x_ss is the steady state value of x.
str = ['piwt = 100*piw;  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['ut = 100*u;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Kt = 100*(K/',num2str(K,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Ct = 100*(Ctot/steady_state(Ctot)-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Bt = 100*((B)/',num2str(eco.B,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['wt = 100*(w/',num2str(eco.w,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['rt = 100*(r - ',num2str(eco.R-1,formatSpec),');'] ; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['rnt = 100*(RBN - ',num2str(1+(eco.R-1)/(1-eco.tk),formatSpec),');'] ; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Tt = 100*(TT/',num2str(eco.TT,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['tkt = 100*(tauk -',num2str(eco.tk,formatSpec),');'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['tlt = 100*(taul -',num2str(eco.tl,formatSpec),');'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['PIt = 100*(PI-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['GDPt = 100*(GDP/',num2str(Y,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Zt = -ut;'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Lt = 100*(Ltot/steady_state(Ltot)-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['TsYt = Tt - GDPt ;'] ; fprintf(fid, str);   str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['BsYt = Bt - GDPt ;'] ; fprintf(fid, str);   str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

%%%%%%%%%%%%%%%% end of the model part
str = ['end;\n\n'] ; fprintf(fid, str);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     Steady-state model       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str = ['steady_state_model;\n\n'] ; fprintf(fid, str);
for h = 1:Nbin
    str = ['a',num2str(h),' = ',num2str(eco.aep(h)/Sp(h),formatSpec),';\n'] ; fprintf(fid, str);
    str = ['at',num2str(h),' = ',num2str(eco.abp(h)/Sp(h),formatSpec),';\n'] ; fprintf(fid, str);
    str = ['c',num2str(h),' = ',num2str(eco.cp(h)/Sp(h),formatSpec),';\n'] ; fprintf(fid, str);
    str = ['l',num2str(h),' = ',num2str(eco.lp(h)/Sp(h),formatSpec),';\n'] ; fprintf(fid, str);
end;
str = ['r = ', num2str(eco.R-1,formatSpec),';\n'] ; fprintf(fid, str);
str = ['w = ', num2str(eco.w,formatSpec),';\n'] ; fprintf(fid, str);
str = ['FK = ', num2str((1/beta-1),formatSpec),';\n'] ; fprintf(fid, str);
str = ['FL = ', num2str(1*eco.w/(1-eco.tl),formatSpec),';\n'] ; fprintf(fid, str);
str = ['K = ', num2str(K,formatSpec),';\n'] ; fprintf(fid, str);
str = ['Ltot =' num2str(Ltot,formatSpec),' ;\n'] ; fprintf(fid, str);
str = ['TFP = ', num2str(1),';\n'] ; fprintf(fid, str);
str = ['A = ', num2str(eco.A,formatSpec),';\n'] ; fprintf(fid, str);
str = ['B = ', num2str(eco.B,formatSpec),';\n'] ; fprintf(fid, str);
str = ['u = ', num2str(0),';\n'] ; fprintf(fid, str);
str = ['tauk = ', num2str(eco.tk,formatSpec),';\n'] ; fprintf(fid, str);
str = ['taul = ', num2str(eco.tl,formatSpec),';\n'] ; fprintf(fid, str);
str = ['Ctot = ', num2str(Ctot,formatSpec),';\n'] ; fprintf(fid, str);
str = ['GDP = ', num2str(Y,formatSpec),';\n'] ; fprintf(fid, str);
str = ['TT = ', num2str(eco.TT,formatSpec),';\n'] ; fprintf(fid, str);
str = ['PI   = 1;\n'] ; fprintf(fid, str);
str = ['zeta = 1;\n'] ; fprintf(fid, str);
str = ['M = ', num2str(1,formatSpec),';\n'] ; fprintf(fid, str);
str = ['rtildeK = r/(1-tauk);\n'] ; fprintf(fid, str);
str = ['RBN = 1+rtildeK;\n'] ; fprintf(fid, str);
str = ['wtilde =FL;\n'] ; fprintf(fid, str);
str = ['ut = 0;\n'] ; fprintf(fid, str);
str = ['Kt = 0;\n'] ; fprintf(fid, str);
str = ['Ct  = 0;\n'] ; fprintf(fid, str);
str = ['Bt = 0;\n'] ; fprintf(fid, str);
str = ['wt = 0;\n'] ; fprintf(fid, str);
str = ['rt = 0;\n'] ; fprintf(fid, str);
str = ['rnt = 0;\n'] ; fprintf(fid, str);
str = ['Tt = 0;\n'] ; fprintf(fid, str);
str = ['tkt = 0;\n'] ; fprintf(fid, str);
str = ['tlt = 0;\n'] ; fprintf(fid, str);
str = ['Zt = 0;\n'] ; fprintf(fid, str);
str = ['PIt = 0;\n'] ; fprintf(fid, str);
str = ['Lt = 0;\n'] ; fprintf(fid, str);
str = ['GDPt = 0;\n'] ; fprintf(fid, str);
str = ['TsYt = 0 ;'] ; fprintf(fid, str);
str = ['BsYt = 0 ;'] ; fprintf(fid, str);
str = ['piw   = 0;\n'] ; fprintf(fid, str);
str = ['piwt  = 0;\n'] ; fprintf(fid, str);
 %%%%%%%%%%% end of the steady state part
 str = ['end;\n\n'] ; fprintf(fid, str);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%   Steady-state computation    %%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %~ str = ['resid;\n'] ; fprintf(fid, str);
 % str = ['steady(nocheck);\n\n']; fprintf(fid, str);
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%       Aggregate shock        %%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 str = ['shocks;\n var eps = 0.00001;\nend;\n'];fprintf(fid, str);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%         Stoch simul          %%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
str = ['stoch_simul (order=1,irf=',num2str(LengthIRF),',noprint,nograph,nofunctions) TFP GDPt Tt Bt PIt rt rnt Lt Ct wt taul tauk Kt Ctot K GDP B Ltot TT BsYt r RBN PI u ut piw TsYt;']; 
fprintf(fid, str);

if isOctave
     fflush(fid);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%        Launch Dynare         %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

command = ['dynare ',file, ' noclearall'];
eval(command);

irfs  =oo_.irfs ;
Yv_taylor(Economy,:)   = [0 irfs.GDPt_eps];
Cv_taylor(Economy,:)   = [0 irfs.Ct_eps];
Kv_taylor(Economy,:)   = [0 irfs.Kt_eps];
Lv_taylor(Economy,:)   = [0 irfs.Lt_eps];
rnv_taylor(Economy,:)  = [0 irfs.rnt_eps];
rv_taylor(Economy,:)   = [0 irfs.rt_eps];
Bv_taylor(Economy,:)   = [0 irfs.BsYt_eps];
Tv_taylor(Economy,:)   = [0 irfs.TsYt_eps];
pipv_taylor(Economy,:) = [0 irfs.PIt_eps];
piwv_taylor(Economy,:) = [0 irfs.piwt_eps];
taukv_taylor(Economy,:) = [0 100*irfs.tauk_eps];
taulv_taylor(Economy,:) = [0 100*irfs.taul_eps];
Uv_taylor(Economy,:)    = [0 irfs.ut_eps];

list = oo_.var_list ;
mean_var = oo_.mean ;
var = oo_.var ; 
var = diag(var);
std = var.^(1/2);
std_norm = (std./mean_var)*100;
std_norm = round(std_norm,3); 
mean_var = round(mean_var,3); 

index_Y = find(list == "GDP");
index_C= find(list == "Ctot");
index_L= find(list == "Ltot");
index_K= find(list == "K");
index_RBN= find(list == "RBN");
index_r= find(list == "r");
index_B = find(list == "B");
index_T = find(list == "TT");
index_pip= find(list == "PI");
index_piw= find(list == "piw");
index_u= find(list == "u");
index_tauk = find(list == "tauk");
index_taul = find(list == "taul");

% Those std are not normalized by the mean (rates)
std_norm(index_u) = round(std(index_u) *100,3);
std_norm(index_piw) = round(std(index_piw) *100,3);
std_norm(index_pip) = round(std(index_pip) *100,3);
std_norm(index_RBN) = round(std(index_RBN) *100,3);
std_norm(index_r) = round(std(index_r) *100,3);

var_init = evalin('base', 'oo_.var');
autocorr = evalin('base', 'oo_.autocorr');
list_table     = [list(index_Y) ; list(index_C) ; list(index_L); list(index_K); list(index_RBN) ; list(index_r) ; list(index_B) ; list(index_T) ; list(index_pip); list(index_piw) ; list(index_taul) ; list(index_tauk) ; list(index_u)] ;
mean_var_table = [mean_var(index_Y) ; mean_var(index_C) ; mean_var(index_L); mean_var(index_K); mean_var(index_RBN) ; mean_var(index_r) ; mean_var(index_B) ; mean_var(index_T) ; mean_var(index_pip); mean_var(index_piw) ; mean_var(index_taul) ; mean_var(index_tauk) ; mean_var(index_u)] ;
std_norm_table = [std_norm(index_Y) ; std_norm(index_C) ; std_norm(index_L); std_norm(index_K); std_norm(index_RBN) ; std_norm(index_r) ; std_norm(index_B) ; std_norm(index_T) ; std_norm(index_pip); std_norm(index_piw) ; std_norm(index_taul) ; std_norm(index_tauk) ; std_norm(index_u)] ;
corr_pi_y = var_init(index_Y,index_pip)./( std(index_Y)*std(index_pip) ) ;% Y, pip
corr_piw_y = var_init(index_Y,index_piw)./( std(index_Y)*std(index_piw) ) ;% Y, piw
corr_tauk_y = var_init(index_Y,index_tauk)./(std(index_Y)*std(index_tauk) ); % Y, tauk
corr_taul_y = var_init(index_Y,index_taul)./(std(index_Y)*std(index_taul) ) ;% Y, taul
corr_b_y = var_init(index_Y,index_B)./(std(index_Y)*std(index_B) ) ;% Y, b
corr_c_y = var_init(index_Y,index_C)./(std(index_Y)*std(index_C) );% Y, c
auto_corr_y_y1 =  autocorr{1, 1}(index_Y,index_Y)  ; % Y, Y(-1)
auto_corr_b_b1 =  autocorr{1, 1}(index_B,index_B)  ; % B, B(-1)

corr_list = ["corr(Pip,Y)" "corr(Piw,Y)" "corr(tau^k, Y)"  "corr(tau^L,Y)" "corr(B,Y)" "corr(C,Y)" "corr(Y,Y(-1))" "corr(B,B(-1))"] ;
corr = round([corr_pi_y  corr_piw_y      corr_tauk_y       corr_taul_y     corr_b_y    corr_c_y     auto_corr_y_y1  auto_corr_b_b1],3) ;
corr_list = cellstr(corr_list) ;


baseFileName_momemt  = sprintf('moments'); 
fullFileName_momemt = sprintf('%s_%s.mat', baseFileName_momemt, calib);
save(fullFileName_momemt,"list_table","mean_var_table", "std_norm_table", "corr_list", "corr") 


fileName = ['To_IRFs_SP_',calib,'.mat' ];
save(fileName, "rnv_taylor", "pipv_taylor", "piwv_taylor", "Yv_taylor", "Bv_taylor", "Tv_taylor", "Cv_taylor", "rv_taylor", "taulv_taylor", "taukv_taylor", "Kv_taylor", "Lv_taylor", "Uv_taylor", "LengthIRF") 