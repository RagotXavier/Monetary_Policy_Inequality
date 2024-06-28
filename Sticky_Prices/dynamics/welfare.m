clear

calib = 'baseline'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%         Output parameters         %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LengthIRF  = 20;        %length of the IRFs
LengthSIM  = 10000;     %length of the simulation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%       Beginning of the code       %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


warning('off','all');
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
    warning('off', 'Octave:array-as-logical');
end



filename_root = 'todynare_SP_' ;
fileExtension = '.mat';
FileNameTot = [filename_root, calib, fileExtension];

if exist(FileNameTot, 'file') == 2
    load(FileNameTot);
else
    error(['File ', FileNameTot, ' does not exist or is not a file.']);
end



alpha = alphaa ; 
beta = betaa ; 
gamma = gammaa ; 

% initialization
Nbin  = eco.Nbin;
Pi_h    = eco.Matab;
Sp    = eco.Sp;
rho_u = 0.95;


Ltotss = Ltot;
Ctotss = Ctot;
Kss     = eco.A-eco.B;
Yss     = Kss^alpha*Ltotss^(1-alpha);



P = ones(Nbin,1);
P(eco.indcc) = 0;

formatSpec = '%25.20f'; %print format

% Initialization of return results
EcoM = 2; % 


Yv = zeros(EcoM,LengthIRF+1);
Cv = zeros(EcoM,LengthIRF+1);
Kv = zeros(EcoM,LengthIRF+1);
Lv = zeros(EcoM,LengthIRF+1);
rv = zeros(EcoM,LengthIRF+1);  % path of nominal rate
rnv = zeros(EcoM,LengthIRF+1); % path of real rate
Bv = zeros(EcoM,LengthIRF+1);
Tv = zeros(EcoM,LengthIRF+1);

PIv = zeros(EcoM,LengthIRF+1);
piwv = zeros(EcoM,LengthIRF+1);

Uv = zeros(EcoM,LengthIRF+1);
taulv = zeros(EcoM,LengthIRF+1);
taukv = zeros(EcoM,LengthIRF+1);



%%
for Economy =1:EcoM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Specifications of the 3 economies %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Economy 1: no time-varying taxes but optimal inflation
    % Economy 2: no time-varying taxes and constant inflation
    % common fiscal rule parameters

    coefBL = 4.0;
    coefBu = 8.5;

    if  Economy == 3
        % time-varying tax
        coefl  = -0.6;
        coeflL =  0.5;
        coefk  =  0.6;   % reaction of labor tax to debt
        coefkL = -0.5;
    else
        % constant tax
        coefl  = 0.;   % reaction of labor tax to debt
        coefk  = 0.;
        coefkL = 0.;
        coeflL = 0.;
    end

    if Economy == 2
        % constant inflation
        opti_PI = 0;
    else
        % optimal inflation
        opti_PI = 1;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%     opening the mod file     %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    file = ['code_dynare_SP_welfare',num2str(Economy),'.mod'];
    fid = fopen(file, 'w') ;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%     variables and params     %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    str = ['var\n'] ; fprintf(fid, str);
    str = ['u r w FK FL A B K mu TFP Ltot tauk taul Ctot GDP TT BsYt TsYt piw piwt \n'];fprintf(fid, str);
    str = ['ut Kt Ct Bt wt rt rnt GDPt PIt rtildeK wtilde Tt tkt tlt mut Lt Zt zeta Lgamma M PI RBN Gamma Upsilon \n'] ; fprintf(fid, str);
    for p = 1:Nbin
        str = ['a',num2str(p),' '] ; fprintf(fid,str);
        str = ['at',num2str(p),' '] ; fprintf(fid, str);
        str = ['c',num2str(p),' '] ; fprintf(fid, str);
        str = ['l',num2str(p),' '] ; fprintf(fid, str);
        str = ['psib',num2str(p),' '] ; fprintf(fid, str);
        str = ['lambdac',num2str(p),' '] ; fprintf(fid, str);
        str = ['lambdatc',num2str(p),' '] ; fprintf(fid, str);
        str = ['lambdal',num2str(p),' '] ; fprintf(fid, str);
        str = ['lambdatl',num2str(p),' '] ; fprintf(fid, str);
        if mod(p,10)==0;
            str = ['\n'] ;         fprintf(fid, str);
        end
    end
    str = ['; \n\n'] ; fprintf(fid, str);

    str = ['varexo eps; \n\n'] ; fprintf(fid, str);
    str = ['parameters\n'] ; fprintf(fid, str);
    str = ['coefl coefk coefBu coefBL beta alpha phi abar delta gamma chi rho_u G kappa epss coefkL coeflL;\n\n'] ; fprintf(fid, str);


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
            if (Pi_h(h,hi))>tol   % hi to h
              str = ['+ ',num2str(Sp(hi)*Pi_h(h,hi)/Sp(h),formatSpec),'*a',num2str(hi),'(-1)']; fprintf(fid, str);
            end
        end
        str = ['; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

        %%%%%%%% Euler equation
        if P(h)==1 % not constrained
            str = [num2str(eco.xsiue(h),formatSpec),'*(c',num2str(h),')^-gamma =', num2str(eco.Resh(h),formatSpec),' + beta*(1+r(+1))*(10^-10']; fprintf(fid, str);
            for hp = 1:Nbin
                    if (Pi_h(hp,h)*eco.xsiue(hp))>tol % h to hp
                        str = ['+ ',num2str((Pi_h(hp,h))*eco.xsiue(hp),formatSpec),'*(c',num2str(hp),'(+1))^-gamma']; fprintf(fid, str);
                    end
                    % if mod(hi,10)==0;str = ['\n'] ;fprintf(fid, str); end;
            end
            str = [');  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        else
            str = ['a',num2str(h),' = ',num2str(eco.aep(h)/Sp(h),formatSpec),';  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        end

        %%%%%%%% Labor supply
        str = [num2str(eco.xsiv1(h),formatSpec) '*(1/chi)*(l',num2str(h),'^(1/phi)) = w*',num2str(eco.xsiu1(h),formatSpec),'*',num2str(nsi,formatSpec),...
                '* c',num2str(h),'^(-gamma); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

        %%%%%%%% Def PSI
        str = ['psib',num2str(h),'= ',num2str(Sp(h)*eco.xsiu0(h),formatSpec ),'*(c',num2str(h),')^-gamma - ( lambdac',num2str(h), '*',num2str(eco.xsiue(h),formatSpec ),...
            '- (1+r)*lambdatc',num2str(h), '*',num2str(eco.xsiue(h),formatSpec ),' - lambdal',num2str(h),'*',num2str(eco.xsiu1(h),formatSpec ),'*w*',num2str(nsi,formatSpec),...
            ')*(-gamma*(c',num2str(h),')^(-gamma-1) );  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;


        %%%%%%%% Foc at : individual savings        : Euler on psib
        if P(h)==1
            str = ['psib',num2str(h),'=', num2str(Sp(h),formatSpec),'*beta*(1+r(+1))*(0']; fprintf(fid, str);
            for hp = 1:Nbin
                if (Pi_h(hp,h))>tol && (Sp(hp))>tol  %be careful transpose from h to hp
                    str = ['+ ',num2str( Pi_h(hp,h),formatSpec ),'*psib',num2str(hp),'(+1)/',num2str(Sp(hp),formatSpec)]; fprintf(fid, str);
                end
            end
            str = [')+',num2str(Sp(h),formatSpec),'*(beta*alpha*w(+1))/K*(']; fprintf(fid, str);
            for hp = 1:Nbin
                if (Sp(hp))>tol   %be careful transpose from h to hp
                  str = ['+ ',num2str(eco.ytype(hp),formatSpec),'*(psib',num2str(hp),'(+1)*l',num2str(hp),'(+1)+lambdatl',num2str(hp),'(+1)*',num2str(eco.xsiu1(hp),formatSpec),'* (c',num2str(hp),'(+1))^(-gamma))']; fprintf(fid, str);
                end
                if mod(hp,10)==0;str = ['\n'] ;fprintf(fid, str);
                end;
            end
            str = [')+',num2str(Sp(h),formatSpec),'*( ']; fprintf(fid, str);
            str = ['+ beta*Gamma(+1)*( r(+1) - (1-tauk(+1))*rtildeK(+1) )']; fprintf(fid, str);
            str = [ '+ beta*mu(+1)/K*(alpha*GDP(+1)*(1-kappa*(PI(+1)-1)*(PI(+1)-1)/2 )- (r(+1)+delta)*K -  alpha*w(+1)*Ltot(+1) )'] ; fprintf(fid, str);
            str = [ '+ beta*alpha*GDP(+1)*M(+1)/K*( (epss-1)/kappa*Lgamma(+1)*(zeta(+1)-1) - ( Lgamma(+1) - Lgamma )*PI(+1)*(PI(+1)-1) )']  ; fprintf(fid, str);
            str = ['- beta*(alpha-1)/K * (rtildeK(+1) + delta)* ( beta^(-1)*Upsilon + Gamma(+1)*(1-tauk(+1))*K ) ']; fprintf(fid, str);
            str = [' ) ; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        else
            str = ['lambdac',num2str(h),' =  0; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        end

        %%%%%%%% lambdatc : lambdac previous period
        str = ['lambdatc',num2str(h),' = 0 ']; fprintf(fid, str);
        for hi = 1:Nbin
            if Pi_h(h,hi)>tol  %be careful transpose
              str = ['+ ',num2str(Pi_h(h,hi),formatSpec),'*lambdac',num2str(hi),'(-1)']; fprintf(fid, str);
            end
        end
        str = ['; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

        %%%%%%%% lambdatl : lambdal previous period
        str = ['lambdatl',num2str(h),' = 0 ']; fprintf(fid, str);
        for hi = 1:Nbin
            if Pi_h(h,hi)>tol  %be careful transpose
              str = ['+ ',num2str(Pi_h(h,hi),formatSpec),'*lambdal',num2str(hi),'(-1)']; fprintf(fid, str);
            end
        end
        str = ['; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

        %%%%%%%% FOC planner l
        str = ['psib',num2str(h), '=', '(1/(w*',num2str(nsi,formatSpec),'))*(',num2str(Sp(h)*eco.xsiv0(h),formatSpec),'*(1/chi)*l',num2str(h),'^(1/phi)+lambdal',num2str(h),'*',num2str(eco.xsiv1(h),formatSpec),'*(1/chi) *(1/phi)*l',num2str(h),'^(1/phi-1))'] ; fprintf(fid, str);
        str = ['+(',num2str(Sp(h),formatSpec),')*('] ; fprintf(fid, str);
        str = [' (1-alpha)*M*GDP/(Ltot*w)*((Lgamma - Lgamma(-1))*PI*(PI-1)-(epss-1)/kappa*Lgamma*(zeta-1))'] ; fprintf(fid, str);
        str = ['-', ' (1-alpha)*mu/(Ltot*w) *(GDP*(1-kappa/2*(PI-1)^2)-w*Ltot)'] ; fprintf(fid, str);
        str = ['+', '(1-alpha)*(rtildeK + delta)/(w*Ltot)*(Gamma*(1-tauk)*K(-1)+ beta^(-1)*Upsilon(-1))'] ; fprintf(fid, str);
        str = ['+alpha*(1/Ltot)*(10^-10']; fprintf(fid, str);
        for hp = 1:Nbin
            if Sp(hp)>tol   %be careful transpose from h to hp
                str = ['+(psib',num2str(hp),'*l',num2str(hp),'+lambdal',num2str(hp),'*',num2str(eco.xsiu1(hp),formatSpec),...
                        '* c',num2str(hp),'^(-gamma))*',num2str(eco.ytype(hp),formatSpec)] ; fprintf(fid, str);
                if mod(hp,10)==0;
                    str = ['\n'] ;fprintf(fid, str);
                end
            end
        end
        str = [')); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        if isOctave
            fflush(fid);
        end
    end

    %%%%%%%% FOC planner r
    str = ['0  = 0']; fprintf(fid, str);
    for hi = 1:Nbin
          if (Sp(hi))>tol
              str = ['+ at',num2str(hi),'*psib',num2str(hi)]; fprintf(fid, str);
          end
         if mod(hi,10)==0;str = ['\n'] ;fprintf(fid, str); end;
    end
    for hp = 1:Nbin
        if Sp(hp)>tol   %be careful transpose from h to hp
            str = ['+ ',num2str(eco.xsiue(hp),formatSpec),'*lambdatc',num2str(hp),'*(c',num2str(hp),'^-gamma)']; fprintf(fid, str);
        end
        if mod(hp,10)==0;str = ['\n'] ;fprintf(fid, str); end;
    end
    str = ['+(Gamma-mu)*(B(-1)+K(-1)); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

    %%%%%%%% FOC planner zeta
    str = ['mu*w*Ltot- (epss-1)/kappa*Lgamma*zeta*GDP + (rtildeK + delta)*(beta^(-1)*Upsilon(-1) + Gamma*(1-tauk)*K(-1) )= w*(10^-10 ']; fprintf(fid, str);
    for h = 1:Nbin
        if Sp(h)>tol   %
            str = ['+',num2str(eco.ytype(h),formatSpec),'*(psib',num2str(h),'*l',num2str(h),'+lambdal',num2str(h),'*',num2str(eco.xsiu1(h),formatSpec),'* c',num2str(h),'^(-gamma))'];
            fprintf(fid, str);
        end
    end
    str = ['); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

    %%%%%%%% FOC planner  RBN
    str = ['beta*(1-tauk(+1))*(Gamma(+1)/PI(+1))*B =  Upsilon/PI(+1); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

    %%%%%%%% FOC planner PI
    if opti_PI
        str = ['0 = mu*kappa*(PI-1) + (Lgamma- Lgamma(-1))*(2*PI-1)*M + (beta^(-1)*Upsilon(-1) - Gamma*(1-tauk)*B(-1)) *RBN(-1)/(PI*PI*GDP); %% equation ', num2str(equa),'\n'];
        fprintf(fid, str); equa = equa+1;
    else
        str = ['PI = 1; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    end
    
    %%%%%%%% FOC planner TT
    str = ['mu = 0'] ; fprintf(fid, str);
    for h = 1:Nbin
            str = ['+ psib',num2str(h)] ; fprintf(fid, str);
        if mod(p,10)==10; str = ['\n'] ;fprintf(fid, str);
        end;
    end;
    str = ['; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

    %%%%%%%% Definition Pricing kernel M
    str = ['M = 1; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;


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
    str = ['Kt = 100*(K/',num2str(Kss,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    str = ['Ct = 100*(Ctot/steady_state(Ctot)-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    str = ['Bt = 100*((B)/',num2str(eco.B,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    str = ['wt = 100*(w/',num2str(eco.w,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    str = ['rt = 100*(r - ',num2str(eco.R-1,formatSpec),');'] ; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    str = ['rnt = 100*(RBN - ',num2str(1+(eco.R-1)/(1-eco.tk),formatSpec),');'] ; fprintf(fid, str); str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    str = ['Tt = 100*(TT/',num2str(eco.TT,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    str = ['tkt = 100*(tauk -',num2str(eco.tk,formatSpec),');'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    str = ['tlt = 100*(taul -',num2str(eco.tl,formatSpec),');'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    str = ['mut = 100*(mu/',num2str(Plann.mu,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    str = ['PIt = 100*(PI-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
    str = ['GDPt = 100*(GDP/',num2str(Yss,formatSpec),'-1);'] ; fprintf(fid, str);  str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
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
        str = ['lambdac',num2str(h),' = ',num2str(Plann.lambdac(h),formatSpec),';\n'] ; fprintf(fid, str);
        str = ['lambdal',num2str(h),' = ',num2str(Plann.lambdal(h),formatSpec),';\n'] ; fprintf(fid, str);
        str = ['lambdatc',num2str(h),' = ',num2str(Plann.lambdact(h),formatSpec),';\n'] ; fprintf(fid, str);
        str = ['lambdatl',num2str(h),' = ',num2str(Plann.lambdalt(h),formatSpec),';\n'] ; fprintf(fid, str);
        str = ['psib',num2str(h),' = ',num2str(Plann.psib(h),formatSpec),';\n'] ; fprintf(fid, str);
    end;

    str = ['r = ', num2str(eco.R-1,formatSpec),';\n'] ; fprintf(fid, str);
    str = ['w = ', num2str(eco.w,formatSpec),';\n'] ; fprintf(fid, str);
    str = ['FK = ', num2str((1/beta-1),formatSpec),';\n'] ; fprintf(fid, str);
    str = ['FL = ', num2str(1*eco.w/(1-eco.tl),formatSpec),';\n'] ; fprintf(fid, str);
    str = ['K = ', num2str(Kss,formatSpec),';\n'] ; fprintf(fid, str);
    str = ['Ltot =' num2str(Ltotss,formatSpec),' ;\n'] ; fprintf(fid, str);
    str = ['mu = ', num2str(Plann.mu,formatSpec),';\n'] ; fprintf(fid, str);
    str = ['TFP = ', num2str(1),';\n'] ; fprintf(fid, str);
    str = ['A = ', num2str(eco.A,formatSpec),';\n'] ; fprintf(fid, str);
    str = ['B = ', num2str(eco.B,formatSpec),';\n'] ; fprintf(fid, str);
    str = ['u = ', num2str(0),';\n'] ; fprintf(fid, str);
    str = ['tauk = ', num2str(eco.tk,formatSpec),';\n'] ; fprintf(fid, str);
    str = ['taul = ', num2str(eco.tl,formatSpec),';\n'] ; fprintf(fid, str);
    str = ['Ctot = ', num2str(Ctotss,formatSpec),';\n'] ; fprintf(fid, str);
    str = ['GDP = ', num2str(Yss,formatSpec),';\n'] ; fprintf(fid, str);
    str = ['TT = ', num2str(eco.TT,formatSpec),';\n'] ; fprintf(fid, str);
    str = ['PI   = 1;\n'] ; fprintf(fid, str);
    str = ['zeta = 1;\n'] ; fprintf(fid, str);
    str = ['Lgamma = ', num2str(Plann.Lgamma,formatSpec),';\n'] ; fprintf(fid, str);
    str = ['Gamma = ', num2str(Plann.Gamma,formatSpec),';\n'] ; fprintf(fid, str);
    str = ['Upsilon = ', num2str(Plann.Upsilon,formatSpec),';\n'] ; fprintf(fid, str);
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
    str = ['mut = 0;\n'] ; fprintf(fid, str);
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
    % if !bprint
        % str = ['options_.noprint=1;\n']; fprintf(fid, str); %1 removes printing of steady
    % end
    % str = ['steady(nocheck);\n\n']; fprintf(fid, str);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%       Aggregate shock        %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    str = ['shocks;\n var eps = 0.00001;\nend;\n'];fprintf(fid, str);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%         Stoch simul          %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        str = ['stoch_simul (order=1,irf=',num2str(LengthIRF),', periods=',num2str(LengthSIM),',nograph,nofunctions,noprint) TFP GDPt Tt Bt PIt rt rnt Lt Ct wt taul tauk Kt Ctot K GDP B Ltot TT BsYt TsYt r RBN PI u ut piw piwt;']; 


        fprintf(fid, str);
    if isOctave
        fflush(fid);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%        Launch Dynare         %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    command = ['dynare ',file, ' noclearall'];
    eval(command);
end




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% WELFARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    close all ;

    isEqual = strcmp(calib, 'baseline');
    if isEqual
        % if WELFARE == 1 

            PIss = Pi_h;

            if gamma==1
                u_c = @(c) log(c);
            else   
                u_c = @(c) (c.^(1-gamma) -1)./(1-gamma);
            end

            v_l =  @(l) (1/chi).*(l.^(1+1/phi))./(1+1/phi);


            Maxseed = 10; 
            Welfare = zeros(Maxseed,2,Nbin); % store the vector of intertemporal utility for each bin

            tic

            for Economy = [1, 2] % 1 is Optimal, 2 is Pi^P = 1


                first_time = 1;
                for i=1:Maxseed +1 % ici + 1 car le premier est pour le if first time 
                    disp('economy');disp(i);
                    toc
                    if first_time
                        file = ['code_dynare_SP_welfare',num2str(Economy),'.mod',' ','noclearall'];
                        command = ['dynare ',file]
                        eval(command)
                        first_time = 0;
                    else
                        % Initialize vectors for consumption and labor by bin and periods
                        conso_simul_vec = zeros(Nbin,LengthSIM); % consumption vector by bin and periods of the simu
                        l_simul_vec = zeros(Nbin,LengthSIM);    % consumption vector by bin and periods of the simu
                        set_dynare_seed(i-1); % Set seed for reproducibility in Dynare simulations

                        % Perform stochastic simulation
                        [info, oo_, options_, M_]  = stoch_simul(M_, options_, oo_, var_list_) ; 
    
                        % Extract consumption and labor position per bin from simulation results
                        % For given seed, results per bin and period 

                        for j=1:Nbin
                            idxc = find(strcmp([M_.endo_names], ['c',num2str(j)])) ;
                            idxl = find(strcmp([M_.endo_names], ['l',num2str(j)])) ;
                            conso_simul_vec(j,:) =  oo_.endo_simul(idxc,:) ;
                            l_simul_vec(j,:) =  oo_.endo_simul(idxl,:) ;
                        end

                        % For given seed, results per bin and period
                        % Calculate utility for each bin and period
                        Ut_simul_vec = zeros(Nbin,LengthSIM); % consumption vector by bin and periods
                        Ut_simul_vec = u_c(conso_simul_vec) - v_l(l_simul_vec);
                        Ut_simul_vec = repmat(Sp,1,LengthSIM).*Ut_simul_vec;
                        W = Ut_simul_vec(:,end);

                        for k=1:(LengthSIM-1)
                            W = Ut_simul_vec(:,LengthSIM-k) +beta*PIss*W;
                        end   
                        Welfare(i-1, Economy,:) = (1-beta)*W./Sp; % percapita

                    end
                end
            end

            % Computing consumption equivalent in percent (I rely here on the log case). 

            CE = zeros(Maxseed,Nbin);
            CE = 100*(exp(Welfare(:,1,:)-Welfare(:,2,:))-1); % per capita per bin per simulation 

            CEav = zeros(Maxseed,1);

            for i =1:Maxseed
                CEav(i,1) =(sum(CE(i,:)'.*Sp)); % (per bin) per simulation over bins 
            end

            meanCE =  zeros(Nbin,1);

            
            for i=1:Nbin
                meanCE(i,1) = mean(CE(:,i)); % Average CE per Bin over simulations
            end

            for i=1:eco.ns
                index = find(eco.ytype == states(i));
                CEy(i) = real(sum(Sp(index).*meanCE(index,1))/sum(Sp(index)));
            end

            w = mean(CEav)

            save("For_Welfare_SP.mat", "CEy", "states",'w') 

            % average aggregate welfare gain in percent over all simulations

        end
