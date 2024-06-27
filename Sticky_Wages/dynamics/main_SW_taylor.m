clear
warning('off','all');

calib = 'taylor' ;

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
isOctave && warning('off', 'Octave:array-as-logical');

filename_root = 'todynare_SW_' ;
fileExtension = '.mat';

load toDynare_SW.mat 
FileNameTot = [filename_root, calib, fileExtension];


trunc  = Ramsey.truncatedModel.truncatedAllocation;
Ntot  = Ramsey.truncatedModel.Ntot(1);    %total number of bins
Pi_h  = trunc.Pi_h';
S_h   = trunc.S_h;
K     = solution.K(1);
alpha = eco.alpha(1);
gamma = eco.gamma(1);
chi   = eco.chi(1);
phi   = eco.phi(1);
psiw  = eco.psiw(1);
beta  = eco.beta(1);
Y     = K^alpha*solution.L(1)^(1-alpha);
xsis  = Ramsey.truncatedModel.xsis;
ys    = eco.ys;
rho_u = 0.95;

phi_taylor   = 1.5 ;  %Taylor Rule 

P     = ones(Ntot,1);
P(trunc.ind_cc_h) = 0;


 
LengthIRF = 20 ; 

EcoM = 1;
Economy = EcoM ;


Yv_taylor = zeros(EcoM,LengthIRF+1);
Cv_taylor = zeros(EcoM,LengthIRF+1);
Lv_taylor = zeros(EcoM,LengthIRF+1);
Kv_taylor = zeros(EcoM,LengthIRF+1);
rv_taylor = zeros(EcoM,LengthIRF+1);  
RBNv_taylor = zeros(EcoM,LengthIRF+1); 
Bv_taylor = zeros(EcoM,LengthIRF+1);
Tv_taylor = zeros(EcoM,LengthIRF+1);
PIpv_taylor = zeros(EcoM,LengthIRF+1);
PIwv_taylor = zeros(EcoM,LengthIRF+1);
taulv_taylor = zeros(EcoM,LengthIRF+1);
taukv_taylor = zeros(EcoM,LengthIRF+1);
Utilv_taylor = zeros(EcoM,LengthIRF+1);
Uv_taylor = zeros(EcoM,LengthIRF+1);


coefBL = 4.0;
coefBu = 8.5;
coefl  = 0.;  
coefk  = 0.;
coefkL = 0.;
coeflL = 0.;
    

file = ['code_dynare_SW_',calib,'.mod'];
fid = fopen(file, 'w') ;
formatSpec = '%25.20f';
omegab_h = ones(Ntot,1) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     variables and params     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str = ['var\n'] ; fprintf(fid, str);
str = ['u  r  w  A  L   Ctot Y  piw  pip  T  Util  B  G  Z IUtil \n'] ; fprintf(fid, str);
str = ['ut rt wt At Lt  Ct   Yt piwt pipt Tt Utilt Bt TsYt \n'] ; fprintf(fid, str);
str = ['RBN tauk taul K mkupw  \n'] ; fprintf(fid, str);
str = ['RBNt Kt BsYt \n'] ; fprintf(fid, str);

for h = 1:Ntot
    str = ['a',num2str(h),' '] ; fprintf(fid,str);     %  end-of-period
    str = ['at',num2str(h),' '] ; fprintf(fid, str);   % beginning-of-period
    str = ['c',num2str(h),' '] ; fprintf(fid, str);
    if mod(h,10)==0; %multiple de 10
        str = ['\n'] ;         fprintf(fid, str);
    end;
end;

str = ['; \n\n'] ; fprintf(fid, str);
str = ['varexo eps; \n\n'] ; fprintf(fid, str);
str = ['parameters\n'] ; fprintf(fid, str);
str = ['beta alpha phi abar delta gamma chi psiw epsw rho_u Gss \n\n'] ; fprintf(fid, str);
str = ['phi_taylor\n\n'] ; fprintf(fid, str); 
str = ['coefl coefk coefBu coefBL coefkL coeflL;\n\n'] ; fprintf(fid, str);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%        initialisation        %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


str = ['alpha','   = ',num2str(alpha,formatSpec),';\n']; fprintf(fid, str);
str = ['beta','    = ',num2str(eco.beta,formatSpec),';\n']; fprintf(fid, str);
str = ['phi','     = ',num2str(eco.phi,formatSpec),';\n']; fprintf(fid, str);
str = ['abar','    = ',num2str(eco.a_min,formatSpec),';\n']; fprintf(fid, str);
str = ['delta','   = ',num2str(eco.delta,formatSpec),';\n']; fprintf(fid, str);
str = ['gamma','   = ',num2str(eco.gamma,formatSpec),';\n']; fprintf(fid, str);
str = ['rho_u','   = ',num2str(rho_u,formatSpec),';\n']; fprintf(fid, str); %HERE
str = ['Gss        = ',num2str(solution.G,formatSpec),';\n']; fprintf(fid, str);
str = ['chi','     = ',num2str(eco.chi,formatSpec),';\n']; fprintf(fid, str);
str = ['psiw','    = ',num2str(eco.psiw,formatSpec),';\n']; fprintf(fid, str);
str = ['epsw','    = ',num2str(eco.epsw,formatSpec),';\n']; fprintf(fid, str);
str = ['coefl  = ',num2str(coefl,formatSpec),';\n']; fprintf(fid, str);
str = ['coefk  = ',num2str(coefk,formatSpec),';\n']; fprintf(fid, str);
str = ['coefBu = ',num2str(coefBu,formatSpec),';\n']; fprintf(fid, str);
str = ['coefBL = ',num2str(coefBL,formatSpec),';\n']; fprintf(fid, str);
str = ['coefkL = ',num2str(coefkL,formatSpec),';\n']; fprintf(fid, str);
str = ['coeflL = ',num2str(coeflL,formatSpec),';\n']; fprintf(fid, str);
str = ['phi_taylor','   = ',num2str(phi_taylor,formatSpec),';\n']; fprintf(fid, str); 

equa = 1; % Numbering equations
tol = 10^-16;  % threshold for transition, to avoid to consider very small transitions
str = ['\n\n'] ; fprintf(fid, str);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%       model definition       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


str = ['model;\n\n']; fprintf(fid, str);

for h = 1:Ntot
    
    %%%%%%%% Budget constraint
    str = ['c',num2str(h),' = ',num2str(trunc.y0_h(h),formatSpec),'*w*L + (1+r)*at',num2str(h),'-a', num2str(h),'+T; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    
    %%%%%%%% Definition of atilde (at)
    str = ['at',num2str(h),' = 0 ']; fprintf(fid, str);
    for hi = 1:Ntot
        if (Pi_h(h,hi))>tol   % hi to h
        str = [' + ',num2str(S_h(hi)*Pi_h(h,hi)/S_h(h),formatSpec),'*a',num2str(hi),'(-1)']; fprintf(fid, str);
        end
    end
    str = ['; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

    %%%%%%%% Euler equation
    if P(h)==1 % not constrained
        str = [num2str(xsis.xsiuE(h),formatSpec),'*(c',num2str(h),')^-gamma =', num2str(trunc.resid_E_h(h),formatSpec),' + beta*(1+r(+1))*(0']; fprintf(fid, str);
        for hp = 1:Ntot
            if (Pi_h(hp,h)*xsis.xsiuE(hp))>tol
                str = [' + ',num2str(Pi_h(hp,h)*xsis.xsiuE(hp),formatSpec),'*(c',num2str(hp),'(+1))^-gamma']; fprintf(fid, str);
            end
        end
        str = ['); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    else
        str = ['a',num2str(h),' =  ',num2str(trunc.a_end_h(h),formatSpec),';']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    end

    isOctave && fflush(fid);
end

sumOmegab = 0;
for h = 1:Ntot
    sumOmegab = sumOmegab + omegab_h(h)*S_h(h);
end

%%%%%%%%% Taylor Rule 
str = ['RBN/steady_state(RBN) = (1+pip)^phi_taylor ;'] ; fprintf(fid, str);      str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;


%%%%%%%% price wage inflation relationship
str = ['(1+piw)*w(-1) = (1+pip)*w; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;


%%%%%%%% Phillips curve wage
str = ['mkupw = (1/chi)*L^(1/phi) - (epsw-1)/epsw*w*(0']; fprintf(fid, str);
for hp = 1:Ntot
    str = [' + ',num2str(xsis.xsiu1(hp)*S_h(hp)*trunc.y0_h(hp),formatSpec),'*(c',num2str(hp),')^(-gamma)']; fprintf(fid, str);
end
str = ['); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
mkupw_SS_part = 0 ;
for hp = 1:Ntot
        mkupw_SS_part   = mkupw_SS_part + xsis.xsiu1(hp)*S_h(hp)*trunc.y0_h(hp)*trunc.c_h(hp)^(-eco.gamma); 
end
mkupw_SS =(1/eco.chi)*solution.L^(1/eco.phi)- (eco.epsw-1)/eco.epsw*solution.w*mkupw_SS_part ;
str = ['piw*(piw+1) = (epsw/psiw)*mkupw*L + beta*piw(+1)*(piw(+1)+1); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;


%%%%%%%% Resource contraint
str = ['G = Gss + ',num2str(solution.G,formatSpec),'*0 ;\n'] ; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['G = B + Y -delta*K(-1) - T -  B(-1) - r*A(-1) -  w*L ;'] ; fprintf(fid, str);   str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;


%%%%%%%% Financial market clearing
str = ['B = A - K; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

%%%%%%%% Production side
str = ['u = rho_u*u(-1) + eps;  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Z = 1 - u; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

str = ['Y = Z*K(-1)^alpha*L^(1-alpha);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['rt = alpha*Z*(K(-1)/L)^(alpha-1)-delta;'] ; fprintf(fid, str);   str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1; 
str = ['wt = (1-alpha)*Z*(K(-1)/L)^alpha;'] ; fprintf(fid, str);   str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1; 
str = ['w = wt*(1-taul);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['r = ((1-tauk)*( rt*K(-1) +  B(-1)*(RBN(-1)/(pip+1)-1) ))/A(-1) ;'] ; fprintf(fid, str);      str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;


%%%%%%%%%%%%%%%% Interest Rate Constraint:
str = ['RBN = (pip(+1)+1)*(1+rt(+1));'] ; fprintf(fid, str);      str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;


%%%%%%%% Aggregate consumption
str = ['Ctot = 0'] ; fprintf(fid, str);
for h = 1:Ntot
    str = [' + ',num2str(S_h(h),formatSpec),'*c',num2str(h)] ; fprintf(fid, str);
end;
str = [';  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;


%%%%%%%% Aggregate utility
str = ['Util = -0.5*psiw*piw*piw - 1/(chi*(1+1/phi))*(L^(1+1/phi))'] ; fprintf(fid, str);
for h = 1:Ntot
    if abs(eco.gamma-1) < 1e-6
        str = [' + ',num2str(omegab_h(h)*xsis.xsiu0(h),formatSpec),'*log(c',num2str(h),')'] ; fprintf(fid, str);
    else
        str = [' + ',num2str(omegab_h(h)*xsis.xsiu0(h),formatSpec),'*(c',num2str(h),'^(1-gamma))/(1-gamma)'] ; fprintf(fid, str);
    end
end;
str = [';  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['IUtil = Util + beta*IUtil(+1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;


%%%%%%%% Aggregate savings
str = ['A = 0'] ; fprintf(fid, str);
for h = 1:Ntot
    str = [' + ', num2str(S_h(h),formatSpec),'*a',num2str(h)] ; fprintf(fid, str);
end;
str = [';  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;str = ['At = 0'] ; fprintf(fid, str);
for h = 1:Ntot
    str = [' + ', num2str(S_h(h),formatSpec),'*at',num2str(h)] ; fprintf(fid, str);
end;
str = [';  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;    % %%%%%%%%%%%%%%%%%%%%      Exogenous Taxes

%%%%%%%% Fiscal rules
str = ['taul = ',num2str(eco.taul,formatSpec),' + coefl*u + coeflL*u(-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
str = ['tauk = ',num2str(eco.tauk,formatSpec),' + coefk*u + coefkL*u(-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
str = ['T    = ',num2str(eco.Tt,formatSpec),' + coefBu*u - coefBL*(B - ',num2str(solution.B,formatSpec),'); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    


Util = -1/(eco.chi*(1+1/eco.phi))*(solution.L^(1+1/eco.phi));
for h = 1:Ntot
    if abs(eco.gamma-1) < 1e-6
        Util = Util + omegab_h(h)*xsis.xsiu0(h)*log(trunc.c_h(h));
    else
        Util = Util + omegab_h(h)*xsis.xsiu0(h)*(trunc.c_h(h)^(1-eco.gamma))/(1-eco.gamma);
    end
end


%%%%%%%%%%%%%%%% log-deviation variables %%%%%%%%%%%%%%%%
% There are denoted with 't': for variable x, xt = (x -x_ss)/x_ss, where x_ss is the steady state value of x.
str = ['ut = 100*u;  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Kt = 100*(K/',num2str(solution.K(1),formatSpec),'-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Ct = 100*(Ctot/',num2str(solution.C(1),formatSpec),'-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Bt    = 100*(B/', num2str(solution.B(1),formatSpec),'-1);\n'] ; fprintf(fid, str);
str = ['RBNt = 100*(RBN - ',num2str(1+(solution.R-1)/(1-eco.tauk),formatSpec),');'] ; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Lt = 100*(L/',num2str(solution.L,formatSpec),'-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Yt = 100*(Y/',num2str(solution.Y(1),formatSpec),'-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['piwt = 100*piw;  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['pipt = 100*pip;  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Utilt = 100*(Util-(',num2str(Util,formatSpec),'));  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['Tt = 100*(T/',num2str(eco.Tt,formatSpec),'-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
str = ['BsYt = Bt - Yt ;'] ; fprintf(fid, str);   str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1; 
str = ['TsYt = Tt - Yt ;'] ; fprintf(fid, str);   str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1; 
str = ['end;\n\n'] ; fprintf(fid, str);
str = ['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n'] ; fprintf(fid, str);
%%%%%%%%%%%%%%%% end of the model part


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     Steady-state model       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


str = ['steady_state_model;\n\n'] ; fprintf(fid, str);
for h = 1:Ntot
    str = ['a',num2str(h),'        = ',num2str(trunc.a_end_h(h),formatSpec),';\n'] ; fprintf(fid, str);
    str = ['at',num2str(h),'       = ',num2str(trunc.a_beg_h(h),formatSpec),';\n'] ; fprintf(fid, str);
    str = ['c',num2str(h),'        = ',num2str(trunc.c_h(h),formatSpec),';\n'] ; fprintf(fid, str);
    isOctave && fflush(fid);
end;

str = ['r     = ', num2str(solution.R-1,formatSpec),';\n'] ; fprintf(fid, str);
str = ['w     = ', num2str(solution.w,formatSpec),';\n'] ; fprintf(fid, str);
str = ['L     = ', num2str(solution.L(1),formatSpec),';\n'] ; fprintf(fid, str);
str = ['TFP   = ', num2str(1),';\n'] ; fprintf(fid, str);
str = ['A     = ', num2str(solution.A(1),formatSpec),';\n'] ; fprintf(fid, str);
str = ['B     = ', num2str(solution.B(1),formatSpec),';\n'] ; fprintf(fid, str);
str = ['K     = ', num2str(solution.K(1),formatSpec),';\n'] ; fprintf(fid, str);
str = ['At    = ', num2str(solution.A(1),formatSpec),';\n'] ; fprintf(fid, str);
str = ['u     = ', num2str(0),';\n'] ; fprintf(fid, str);
str = ['Ctot  = ', num2str(solution.C,formatSpec),';\n'] ; fprintf(fid, str);
str = ['Y     = ', num2str(solution.Y,formatSpec),';\n'] ; fprintf(fid, str);
str = ['Z     = 1;\n'] ; fprintf(fid, str);
str = ['rt = r/(1-',num2str(eco.tauk,formatSpec),');\n'] ; fprintf(fid, str);
str = ['wt  = ',num2str(solution.w/(1-eco.taul),formatSpec),';\n'] ; fprintf(fid, str);
str = ['RBN = 1+rt;\n'] ; fprintf(fid, str);
str = ['RBNt = 0;\n'] ; fprintf(fid, str); 
str = ['ut    = 0;\n'] ; fprintf(fid, str);
str = ['BsYt = 0 ;'] ; fprintf(fid, str);  
str = ['TsYt = 0 ;'] ; fprintf(fid, str);  
str = ['Ct    = 0;\n'] ; fprintf(fid, str);
str = ['Bt    = 0;\n'] ; fprintf(fid, str);
str = ['Kt    = 0;\n'] ; fprintf(fid, str);
str = ['wtt    = 0;\n'] ; fprintf(fid, str);
str = ['rtt    = 0;\n'] ; fprintf(fid, str);
str = ['Lt    = 0;\n'] ; fprintf(fid, str);
str = ['Yt    = 0;\n'] ; fprintf(fid, str);
str = ['BoYt  = 0;\n'] ; fprintf(fid, str);
str = ['G     = Gss;\n'] ; fprintf(fid, str);
str = ['piw   = 0;\n'] ; fprintf(fid, str);
str = ['piwt  = 0;\n'] ; fprintf(fid, str);
str = ['mkupw =' num2str(mkupw_SS,formatSpec),';\n'] ; fprintf(fid, str);
str = ['pip   = 0;\n'] ; fprintf(fid, str);
str = ['pipt  = 0;\n'] ; fprintf(fid, str);
str = ['mkupp = 0;\n'] ; fprintf(fid, str);
str = ['taul  = ', num2str(eco.taul,formatSpec),';\n'] ; fprintf(fid, str);
str = ['tauk  = ', num2str(eco.tauk,formatSpec),';\n'] ; fprintf(fid, str);
str = ['T','       = ',num2str(eco.Tt,formatSpec),';\n']; fprintf(fid, str);
str = ['Tt  = 0;\n'] ; fprintf(fid, str);
str = ['Utilt  = 0;\n'] ; fprintf(fid, str);

str = ['Util  = ',num2str(Util,formatSpec),';\n'] ; fprintf(fid, str);
str = ['IUtil  = ',num2str(Util/(1-eco.beta),formatSpec),';\n'] ; fprintf(fid, str);
str = ['end;\n\n'] ; fprintf(fid, str);
%%%%%%%%%%% end of the steady state part

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Steady-state computation    %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%       Aggregate shock        %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%         Stoch simul          %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
str = ['shocks;\n var eps =  0.00001;\nend;\n'];fprintf(fid, str); 
str = ['stoch_simul (order=1,irf=',num2str(LengthIRF), ',nograph, nofunctions,noprint) Y Ctot K L RBN r T B pip piw u ut Yt pipt piwt RBNt rt Lt Ct wt taul tauk Kt BsYt TsYt ;'] ; fprintf(fid, str);

isOctave && fflush(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%        Launch Dynare         %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

command = ['dynare ',file, ' noclearall'];
eval(command);



Yv_taylor(Economy,:)  = [0 oo_.irfs.Yt_eps];
Cv_taylor(Economy,:)  = [0 oo_.irfs.Ct_eps]
Kv_taylor(Economy,:)  = [0 oo_.irfs.Kt_eps];
Lv_taylor(Economy,:) = [0 oo_.irfs.Lt_eps]
RBNv_taylor(Economy,:)=  [0 oo_.irfs.RBNt_eps];
rv_taylor(Economy,:) = [0 oo_.irfs.rt_eps]
Bv_taylor(Economy,:)  = [0 oo_.irfs.BsYt_eps];
Tv_taylor(Economy,:)  = [0 oo_.irfs.TsYt_eps]
PIpv_taylor(Economy,:) = [0  oo_.irfs.pipt_eps];
PIwv_taylor(Economy,:) = [0 oo_.irfs.piwt_eps]
taulv_taylor(Economy,:) = [0 100*oo_.irfs.taul_eps];
taukv_taylor(Economy,:) = [0 100*oo_.irfs.tauk_eps];
Uv_taylor(Economy,:) = [0 oo_.irfs.ut_eps];

     
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% GENERATE Moments_eco1/2/3.MAT %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list = oo_.var_list;
mean_var = oo_.mean
var = diag(oo_.var);
std = var.^(1/2);

index_Y = find(list == "Y");
index_C= find(list == "Ctot")
index_L= find(list == "L");
index_K= find(list == "K")
index_RBN= find(list == "RBN");
index_r= find(list == "r")
index_B = find(list == "B");
index_T = find(list == "T")
index_pip= find(list == "pip");
index_piw= find(list == "piw")
index_u= find(list == "u")
index_tauk = find(list == "tauk");
index_taul = find(list == "taul")

std_norm = (std./mean_var)*100;
std_norm(index_u) = std(index_u) *100;

std_norm = round(std_norm,3);
mean_var = round(mean_var,3)
std_norm(index_piw) = round(std(index_piw) *100,3);
std_norm(index_pip) = round(std(index_pip) *100,3)
std_norm(index_RBN) = round(std(index_RBN) *100,3);
std_norm(index_r) = round(std(index_r) *100,3)
list_table     = [list(index_Y) ; list(index_C) ; list(index_L); list(index_K); list(index_RBN) ; list(index_r) ; list(index_B) ; list(index_T) ; list(index_pip); list(index_piw) ; list(index_taul) ; list(index_tauk) ; list(index_u)] ;
mean_var_table = [mean_var(index_Y) ; mean_var(index_C) ; mean_var(index_L); mean_var(index_K); mean_var(index_RBN) ; mean_var(index_r) ; mean_var(index_B) ; mean_var(index_T) ; mean_var(index_pip); mean_var(index_piw) ; mean_var(index_taul) ; mean_var(index_tauk) ; mean_var(index_u)] ;
std_norm_table = [std_norm(index_Y) ; std_norm(index_C) ; std_norm(index_L); std_norm(index_K); std_norm(index_RBN) ; std_norm(index_r) ; std_norm(index_B) ; std_norm(index_T) ; std_norm(index_pip); std_norm(index_piw) ; std_norm(index_taul) ; std_norm(index_tauk) ; std_norm(index_u)] 

corr_pip_y = oo_.var(index_Y,index_pip)./( std(index_Y)*std(index_pip) ) ;% Y, pip
corr_piw_y = oo_.var(index_Y,index_piw)./( std(index_Y)*std(index_piw) ) ;% Y, pip
corr_tauk_y = oo_.var(index_Y,index_tauk)./(std(index_Y)*std(index_tauk) ); % Y, tauk
corr_taul_y = oo_.var(index_Y,index_taul)./(std(index_Y)*std(index_taul) ) ;% Y, taul
corr_b_y = oo_.var(index_Y,index_B)./(std(index_Y)*std(index_B) ) ;% Y, b
corr_c_y = oo_.var(index_Y,index_C)./(std(index_Y)*std(index_C) );% Y, c
auto_corr_y_y1 =  oo_.autocorr{1, 1}(index_Y,index_Y)  ; % Y, Y(-1)
auto_corr_b_b1 =  oo_.autocorr{1, 1}(index_B,index_B)  ; % B, B(-1
corr_list = ["corr(Pip,Y)" "corr(Piw,Y)"  "corr(tau^k, Y)"  "corr(tau^L,Y)" "corr(B,Y)" "corr(C,Y)" "corr(Y,Y(-1))" "corr(B,B(-1))"] ;
corr = round([corr_pip_y corr_piw_y  corr_tauk_y  corr_taul_y corr_b_y corr_c_y auto_corr_y_y1 auto_corr_b_b1],3) 
corr_list = cellstr(corr_list) 


baseFileName_momemt = sprintf('moments');
fullFileName_momemt = sprintf('%s_%s.mat', baseFileName_momemt, calib);
save(fullFileName_momemt,"list_table","mean_var_table", "std_norm_table", "corr_list", "corr") ;


%%%%%%%%%%%%%%%%% GENERATE For_IRFs_SW.mat %%%%%%%%%%%%%%%%%%%%%%

fileName = ['To_IRFs_SW_',calib,'.mat' ];        
save(fileName, "RBNv_taylor", "PIpv_taylor", "PIwv_taylor", "Yv_taylor", "Bv_taylor", "Cv_taylor", "rv_taylor", "taulv_taylor", "taukv_taylor", "Lv_taylor", "Uv_taylor", "Tv_taylor", "LengthIRF") 

        

