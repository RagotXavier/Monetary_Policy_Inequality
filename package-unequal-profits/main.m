clear 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%         Output parameters         %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bprint  = 0; %set to 1 to get Dynare output (including 2nd-order moments)
Npanels = 3; %set to anyother value to get 8 panels instead of 3. 

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

% loading steady-state values and parameters
load steady_state_dynare.mat

% initialization
Nbin  = eco.Nbin;    
PImat = eco.Matab;
Sp    = eco.Sp;
K     = eco.A-eco.B;
Y     = K^alpha*Ltot^(1-alpha);
rho_u = 0.95;

P = ones(Nbin,1);
P(eco.indcc) = 0;

formatSpec = '%25.20f'; %print format

% Initialization of return results
EcoM = 3; % 3 economies of the papers

rv    = zeros(EcoM,LengthIRF); % path of nominal rate
rnv   = zeros(EcoM,LengthIRF); % path of real rate
Yv    = zeros(EcoM,LengthIRF); % path of GDP
PIv   = zeros(EcoM,LengthIRF); % path of inlfation
Bv    = zeros(EcoM,LengthIRF); % path of public debt
Cv    = zeros(EcoM,LengthIRF); % path of aggregate consmption
taulv = zeros(EcoM,LengthIRF); % path of labor tax
taukv = zeros(EcoM,LengthIRF); % path of capital tax


for Economy =1:EcoM % Loop for the 3 economies

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Specifications of the 3 economies %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Economy 1: no time-varying taxes but optimal inflation
    % Economy 2: no time-varying taxes and constant inflation
    % Economy 3: time-varying taxes and optimal inflation
    
    % common fiscal rule parameters
    coefBL = 4.0;
    coefBu = 8.5;
        
    if  Economy == 3
        % time-varying taxes
        coefl   = -0.6;
        coeflL  = -1.8;
        coefk   =  0.9;
        coefkL  =  0.7;
    else
        % no time-varying taxes
        coefl   = 0.;
        coeflL  = 0.;
        coefk   = 0.;
        coefkL  = 0.;
    end
 
    if Economy == 2
        opti_PI = 0;
    else
        opti_PI = 1;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%     opening the mod file     %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear K
    K     = eco.A-eco.B;
    
    file = ['code_dynare_',num2str(Economy),'.mod'];
    fid = fopen(file, 'w') ;
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%     variables and params     %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    str=['var\n']; fprintf(fid, str);
    str=['u r w FK FL A B K mu TFP Ltot tauk taul Ctot GDP TT BsYt TsYt \n'];fprintf(fid, str);
    str = ['ut Kt Ct Bt wt rt rnt GDPt PIt rtildeK wtilde Tt tkt tlt mut Lt Zt zeta Lgamma M PI RBN Gamma Upsilon \n']; fprintf(fid, str); % Lgamma is the multplier of the NK-Philipps curve

    for p = 1:Nbin
        str = ['a',num2str(p),' ']; fprintf(fid,str);     %  end-of-period
        str = ['at',num2str(p),' ']; fprintf(fid, str);   % beginning-of-period  
        str = ['c',num2str(p),' ']; fprintf(fid, str);    
        str = ['l',num2str(p),' ']; fprintf(fid, str);   % beginning-of-period  
        str = ['psib',num2str(p),' ']; fprintf(fid, str);    
        str = ['lambdac',num2str(p),' ']; fprintf(fid, str);    
        str = ['lambdatc',num2str(p),' ']; fprintf(fid, str);    
        str = ['lambdal',num2str(p),' ']; fprintf(fid, str);    % labor 
        str = ['lambdatl',num2str(p),' ']; fprintf(fid, str);    % labor 
        if mod(p,10)==0; %multiple de 10
            str = ['\n'];         fprintf(fid, str);
        end;
    end;

    for i = 1:length(eco.indcc)
        str = ['LMCC',num2str(eco.indcc(i)),' ']; fprintf(fid,str);     %  Lagrange Multiplier on credit constraint
    end

    str = ['; \n\n']; fprintf(fid, str);

    str = ['varexo eps; \n\n']; fprintf(fid, str);

    str = ['parameters\n']; fprintf(fid, str);
    str = ['coefl coefk coeflL coefkL coefBu coefBL beta alpha phi abar delta gamma chi rho_u G kappa epss nu varPSI;\n\n']; fprintf(fid, str);


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
    str = ['nu      = ',num2str(nu,formatSpec),';\n']; fprintf(fid, str);
    str = ['varPSI  = ',num2str(varPSI,formatSpec),';\n']; fprintf(fid, str);
    str = ['G       = ',num2str(G,formatSpec),';\n']; fprintf(fid, str);
    str = ['chi     = ',num2str(chi,formatSpec),';\n']; fprintf(fid, str);
    str = ['kappa   = ',num2str(kappa,formatSpec),';\n']; fprintf(fid, str); 
    str = ['epss    = ',num2str(epss,formatSpec),';\n']; fprintf(fid, str);
    
    str = ['coefl  = ',num2str(coefl,formatSpec),';\n']; fprintf(fid, str);
    str = ['coeflL = ',num2str(coeflL,formatSpec),';\n']; fprintf(fid, str);
    str = ['coefk  = ',num2str(coefk,formatSpec),';\n']; fprintf(fid, str);
    str = ['coefkL = ',num2str(coefkL,formatSpec),';\n']; fprintf(fid, str);
    str = ['coefBu = ',num2str(coefBu,formatSpec),';\n']; fprintf(fid, str);
    str = ['coefBL = ',num2str(coefBL,formatSpec),';\n']; fprintf(fid, str);

    equa = 1; % Numbering equations
    tol  = 10^(-16);  % threshold for transition, to avoid to consider very small transitions
    str = ['\n\n']; fprintf(fid, str);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%       model definition       %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    str = ['model;\n\n']; fprintf(fid, str);

    for h = 1:Nbin  
        nsi = eco.ytype(h);

        %%%%%%%% Budget constraint  
        str = ['c',num2str(h),' =' 'w*l',num2str(h),'*',num2str(nsi,formatSpec),' + (1+r)*','at',num2str(h),'-a', num2str(h),' + TT + (',...
            num2str(nsi,formatSpec),'^nu/varPSI)*(1-zeta-(kappa/2)*(PI-1)^2)*GDP;  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;

        %%%%%%%% Definition of at   
        str = ['at',num2str(h),' = 10^-10 ']; fprintf(fid, str);
        for hi = 1:Nbin
            if (PImat(h,hi))>tol % hi to h
              str = ['+ ',num2str(Sp(hi)*PImat(h,hi)/Sp(h),formatSpec),'*a',num2str(hi),'(-1)']; fprintf(fid, str);
            end
        end
        str = [';  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;

        %%%%%%%% Euler equation
        if P(h)==1 % not constrained history h
            str = [num2str(eco.xsiue(h),formatSpec),'*(c',num2str(h),')^-gamma =', num2str(eco.Resh(h),formatSpec),' + beta*(1+r(+1))*(0']; fprintf(fid, str);
            for hp = 1:Nbin
                if (PImat(hp,h)*eco.xsiue(hp))>tol % h to hp
                    str = ['+ ',num2str((PImat(hp,h))*eco.xsiue(hp),formatSpec),'*(c',num2str(hp),'(+1))^-gamma']; fprintf(fid, str);
                end
            end
            str = [');  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
        else % constrained history h
            str = ['a',num2str(h),' = ',num2str(eco.aep(h)/Sp(h),formatSpec),';']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;

            % Lagrange multipliers on Euler equations
            str = ['LMCC',num2str(h),'=',num2str(eco.xsiue(h),formatSpec),'*(c',num2str(h),')^-gamma  - beta*(1+r(+1))*(0']; fprintf(fid, str);
            for hp = 1:Nbin
                if (PImat(hp,h)*eco.xsiue(hp))>tol % h to hp
                    str = ['+ ',num2str((PImat(hp,h))*eco.xsiue(hp),formatSpec),'*(c',num2str(hp),'(+1))^-gamma']; fprintf(fid, str);
                end
            end
            str = [');  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
        end

        %%%%%%%% Labor supply
        str = [num2str(eco.xsiv1(h),formatSpec),'*(1/chi)*(l',num2str(h),'^(1/phi)) = w*',num2str(eco.xsiu1(h),formatSpec),'*',num2str(nsi,formatSpec), ...
            '* c',num2str(h),'^(-gamma); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;

        % lambda is by bin here ie lambda*S_yn: it is equivalent to the lambda_bar in the lyx file 
        %%%%%%%% Def PSI
        str = ['psib',num2str(h),'= ',num2str(Sp(h)*eco.xsiu0(h),formatSpec ),'*(c',num2str(h),')^-gamma - (lambdac',num2str(h), '*',num2str(eco.xsiue(h),formatSpec ),...
                '- (1+r)*lambdatc',num2str(h),'*',num2str(eco.xsiue(h),formatSpec ),' - lambdal',num2str(h),'*',num2str(eco.xsiu1(h),formatSpec ),'*w*', ...
                num2str(nsi,formatSpec),')*(-gamma*(c',num2str(h),')^(-gamma-1)); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;

        %%%%%%%% Foc at : individual savings: Euler on psib
        %%%%%%%% remark: psib(h) = S(h)*psi(h) = S(h)*sum_hp(psib(hp)/S(hp))
        if P(h)==1
            str = ['psib',num2str(h),'=', num2str(Sp(h),formatSpec),'*beta*(1+r(+1))*(0']; fprintf(fid, str);
            for hp = 1:Nbin
                if (PImat(hp,h))>tol && (Sp(hp))>tol %be careful transpose from h to hp
                    str = [' + ',num2str( PImat(hp,h),formatSpec ),'*psib',num2str(hp),'(+1)/',num2str(Sp(hp),formatSpec)]; fprintf(fid, str);
                end
            end
            str = [') + ',num2str(Sp(h),formatSpec),'*(beta*alpha*w(+1))/K*(0']; fprintf(fid, str);
            for hp = 1:Nbin
                if (Sp(hp))>tol   %be careful transpose from h to hp
                  str = [' + ',num2str(eco.ytype(hp),formatSpec),'*(psib',num2str(hp),'(+1)*l',num2str(hp),'(+1) + lambdatl',num2str(hp),'(+1)*', ...
                            num2str(eco.xsiu1(hp),formatSpec),' * (c',num2str(hp),'(+1))^(-gamma))']; fprintf(fid, str);
                end
                if mod(hp,10)==0;str = ['\n'];fprintf(fid, str); end;
            end  
            str = [') +',num2str(Sp(h),formatSpec),'*( ']; fprintf(fid, str);
            str = [ ' + beta*Gamma(+1)*( r(+1) - (1-tauk(+1))*rtildeK(+1) )']; fprintf(fid, str);
            str = [ ' + beta*mu(+1)/K*(alpha*GDP(+1)*zeta(+1) - (r(+1)+delta)*K -  alpha*w(+1)*Ltot(+1))']; fprintf(fid, str);
            str = [ ' + beta*alpha*GDP(+1)*M(+1)/K*( (epss-1)/kappa*Lgamma(+1)*(zeta(+1)-1) - ( Lgamma(+1) - Lgamma )*PI(+1)*(PI(+1)-1) )']; fprintf(fid, str);
            str = [ ' + beta*(1-alpha)/K * (rtildeK(+1) + delta)* ( beta^(-1)*Upsilon + Gamma(+1)*(1-tauk(+1))*K ) ']; fprintf(fid, str);
            str = [ ' + beta*alpha*GDP(+1)/K*(1-zeta(+1)-(kappa/2)*(PI(+1)-1)^2)*(0']; fprintf(fid, str);
            for hp = 1:Nbin
                if (Sp(hp))>tol  %be careful transpose from h to hp
                  str = ['+ (',num2str(eco.ytype(hp),formatSpec),'^nu)/varPSI*psib',num2str(hp),'(+1)']; fprintf(fid, str);
                end
                if mod(hp,10)==0;str = ['\n'];fprintf(fid, str);
                end;
            end  
            str = [ ')); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
        else
            str = ['lambdac',num2str(h),' =  0;  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
        end

        %%%%%%%% lambdatc : lambdac previous period
        str = ['lambdatc',num2str(h),' = 0 ']; fprintf(fid, str);
        for hi = 1:Nbin
            if PImat(h,hi)>tol  %be careful transpose
              str = ['+ ',num2str(PImat(h,hi),formatSpec),'*lambdac',num2str(hi),'(-1)']; fprintf(fid, str);
            end
        end
        str = [';  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;

        %%%%%%%% lambdatl : lambdal previous period
        str = ['lambdatl',num2str(h),' = 0 ']; fprintf(fid, str);
        for hi = 1:Nbin
            if PImat(h,hi)>tol  %be careful transpose
              str = ['+ ',num2str(PImat(h,hi),formatSpec),'*lambdal',num2str(hi),'(-1)']; fprintf(fid, str);
            end
        end
        str = [';  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;

        %%%%%%%% FOC planner l
        str = ['psib',num2str(h), '= (1/(w * ',num2str(nsi,formatSpec),'))*(',num2str(Sp(h)*eco.xsiv0(h),formatSpec), ...
                '*(1/chi)*l',num2str(h),'^(1/phi)+lambdal',num2str(h),'*',num2str(eco.xsiv1(h),formatSpec),'*(1/chi) *(1/phi)*l',num2str(h),'^(1/phi-1))']; fprintf(fid, str);
        str = ['+(',num2str(Sp(h),formatSpec),')*(']; fprintf(fid, str);
        str = [' (1-alpha)*M*GDP/(Ltot*w) * ((Lgamma - Lgamma(-1))*PI*(PI-1)-(epss-1)/kappa*Lgamma*(zeta-1))']; fprintf(fid, str);
        str = ['- (1-alpha)*mu/(Ltot*w) * (GDP*zeta - w*Ltot)']; fprintf(fid, str);
        str = ['+ (1-alpha)*(rtildeK + delta)/(w*Ltot)*(Gamma*(1-tauk)*K(-1)+ beta^(-1)*Upsilon(-1))']; fprintf(fid, str);
        str = ['+ alpha/Ltot*(0']; fprintf(fid, str);
        for hp = 1:Nbin
           if Sp(hp)>tol   %be careful transpose from h to hp
                str = ['+(psib',num2str(hp),'*l',num2str(hp),'+lambdal',num2str(hp),'*',num2str(eco.xsiu1(hp),formatSpec),'* c',num2str(hp),'^(-gamma))*', ...
                    num2str(eco.ytype(hp),formatSpec)]; fprintf(fid, str);
                if mod(hp,10)==0; str = ['\n'];fprintf(fid, str); end;
            end
        end
        str = [')'];fprintf(fid, str);
        str = [ '- (1-alpha)*GDP/(w*Ltot)*(1-zeta-(kappa/2)*(PI-1)^2)*(0']; fprintf(fid, str);
        for hp = 1:Nbin
            if (Sp(hp))>tol   %be careful transpose from h to hp
              str = ['+ (',num2str(eco.ytype(hp),formatSpec),'^nu)/varPSI*psib',num2str(hp)]; fprintf(fid, str);
            end
            if mod(hp,10)==0;str = ['\n'];fprintf(fid, str); end;
        end
        str = ['));  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
        
        if isOctave
            fflush(fid);
        end
    end
    
    %%%%%%%% FOC planner r            
    str = ['0  = 0']; fprintf(fid, str); 
    for hi = 1:Nbin
        if (Sp(hi))>tol
            str = [' + at',num2str(hi),'*psib',num2str(hi)]; fprintf(fid, str);
        end
        if mod(hi,10)==0;str = ['\n'];fprintf(fid, str); end;
        end
        for hp = 1:Nbin
        if Sp(hp)>tol   %be careful transpose from h to hp
              str = ['+ ',num2str(eco.xsiue(hp),formatSpec),'*lambdatc',num2str(hp),'*(c',num2str(hp),'^-gamma)']; fprintf(fid, str);
        end
        if mod(hp,10)==0;str = ['\n'];fprintf(fid, str); end; 
    end 
    str = ['+(Gamma-mu)*(B(-1)+K(-1));  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    
    %%%%%%%% FOC planner zeta 
    str = ['mu*w*Ltot - (epss-1)/kappa*Lgamma*zeta*GDP + (rtildeK + delta)*(beta^(-1)*Upsilon(-1) + Gamma*(1-tauk)*K(-1))','= w*(0 ']; fprintf(fid, str);
        for h = 1:Nbin
            if Sp(h)>tol
              str = ['+',num2str(eco.ytype(h),formatSpec),'*(psib',num2str(h),'*l',num2str(h),'+lambdal',num2str(h),'*',num2str(eco.xsiu1(h),formatSpec),'* c',num2str(h),'^(-gamma))']; fprintf(fid, str);
            end
        end
    str = [') + zeta*GDP*(mu - (0 ']; fprintf(fid, str);
    for hp = 1:Nbin
        if (Sp(hp))>tol   %be careful transpose from h to hp
          str = ['+ (',num2str(eco.ytype(hp),formatSpec),'^nu)/varPSI*psib',num2str(hp)]; fprintf(fid, str);
        end
        if mod(hp,10)==0;str = ['\n'];fprintf(fid, str);end;
    end  
    str = ['));  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;

    %%%%%%%% FOC planner RBN 
    str = ['beta*(1-tauk(+1))*(Gamma(+1)/PI(+1))*B =  Upsilon/PI(+1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;

    %%%%%%%% FOC planner PI
    if !opti_PI
        str = ['PI = 1;  %% equation ', num2str(equa),'\n']; fprintf(fid,str); equa=equa+1;
    else
        str = ['0 = (Lgamma- Lgamma(-1))*(2*PI-1)*M + (beta^(-1)*Upsilon(-1) - Gamma*(1-tauk)*B(-1)) *RBN(-1)/(PI*PI*GDP) + kappa*(PI-1)*(0']; fprintf(fid, str);
        for hp = 1:Nbin
            if (Sp(hp))>tol   %be careful transpose from h to hp
              str = ['+ (',num2str(eco.ytype(hp),formatSpec),'^nu)/varPSI*psib',num2str(hp)]; fprintf(fid, str);
            end
            if mod(hp,10)==0;str = ['\n'];fprintf(fid, str);
            end;
        end  
        str = [');  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    end
    
    %%%%%%%% FOC planner TT
    str = ['mu = 0']; fprintf(fid, str);    
    for h = 1:Nbin
        str = ['+ psib',num2str(h)]; fprintf(fid, str);
        if mod(p,10)==10; str = ['\n'];fprintf(fid, str); end;
    end;
    str = ['; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;

    %%%%%%%% Definition Pricing kernel M    
    str = ['M = 1;  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;


    %%%%%%%%%%%%%%%% Aggregate Constraints %%%%%%%%%%%%%%%%

    %%%%%%%% Budget contraint of the govt 
    str = ['B + GDP*zeta -delta*K(-1) = G + TT + B(-1) +r*A(-1) + w*Ltot; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;

    %%%%%%%% NK Philipps curve
    str = ['GDP*M*PI*(PI-1) = ((epss-1)/kappa)*(zeta-1)*GDP*M + beta*PI(+1)*(PI(+1)-1)*GDP(+1)*M(+1); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;

    %%%%%%%% Financial market eq
    str = ['A = B + K;']; fprintf(fid, str);str = [' %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;

    %%%%%%%% Interest rate constraints
    str = ['(r - (1-tauk)*rtildeK)*A(-1) = (1-tauk)*(RBN(-1)/PI - 1- rtildeK)*B(-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['RBN/PI(+1) = (1+rtildeK(+1));  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;

    %%%%%%%% Factor prices : 
    str = ['K(-1)/Ltot = alpha*wtilde/((1-alpha)*(rtildeK+delta)); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;

    %%%%%%%% Production side (here, TPF shock)
    str = ['GDP = TFP*K(-1)^alpha*Ltot^(1-alpha); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['FL = (1 - alpha)*TFP*(K(-1)/Ltot)^alpha;  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1; 
    str = ['FK = alpha*TFP*(K(-1)/Ltot)^(alpha-1)-delta;  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1; 

    %%%%%%%% Aggregate shock 
    str = ['TFP = 1 - u; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['u = rho_u*u(-1) + eps;  %% equation ', num2str(equa),'\n\n']; fprintf(fid, str); equa=equa+1;

    %%%%%%%% Aggregate consumption
    str = ['Ctot = 0']; fprintf(fid, str);    
    for h = 1:Nbin
        str = ['+',num2str(Sp(h),formatSpec),'*c',num2str(h)]; fprintf(fid, str);  
        if mod(h,5)==0; str = ['\n'];fprintf(fid, str);end; 
    end;
    str = [';  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1; 

    %%%%%%%% Aggregate labor supply
    str = ['Ltot = 0 ']; fprintf(fid, str);    
    for h = 1:Nbin
        str = ['+', num2str(Sp(h),formatSpec),'*l',num2str(h),'*',num2str(eco.ytype(h),formatSpec)]; fprintf(fid, str);    
        if mod(p,10)==10; tr = ['\n'];fprintf(fid, str); end;
    end;
    str = [';  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;

    %%%%%%%% Aggregate savings
    str = ['A = ']; fprintf(fid, str);    
    for h = 1:Nbin
        str = ['+', num2str(Sp(h),formatSpec),'*a',num2str(h)]; fprintf(fid, str);    
        if mod(p,10)==10; tr = ['\n'];fprintf(fid, str); end;
    end;
    str = [';  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;

    %%%%%%%% Definition of zeta
    str = ['zeta = (1/(alpha*TFP))*(rtildeK + delta)*(K(-1)/Ltot)^(1-alpha); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;

    %%%%%%%% Definition of w
    str = ['w = (1-taul)*wtilde;   %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;


    %%%%%%%%%%%%%%%% Fiscal rules %%%%%%%%%%%%%%%%
    str = ['taul = ',num2str(eco.tl,formatSpec),' + coefl*u + coeflL*u(-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['tauk = ',num2str(eco.tk,formatSpec),' + coefk*u + coefkL*u(-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['TT = ',num2str(eco.TT,formatSpec),'+ coefBu*u - coefBL*(B - ',num2str(eco.B,formatSpec),'); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1; 

    
    %%%%%%%%%%%%%%%% log-deviation variables %%%%%%%%%%%%%%%%
    % There are denoted with 't': for variable x, xt = (x -x_ss)/x_ss, where x_ss is the steady state value of x.
    
    str = ['ut = 100*u;  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['Kt = 100*(K/',num2str(K,formatSpec),'-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['Ct = 100*(Ctot/',num2str(sum(eco.cp),formatSpec),'-1); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['Bt = 100*(B/',num2str(eco.B,formatSpec),'-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['wt = 100*(w/',num2str(eco.w,formatSpec),'-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['rt = 100*(r - ',num2str(eco.R-1,formatSpec),');  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['rnt = 100*(RBN - ',num2str(1+(eco.R-1)/(1-eco.tk),formatSpec),');  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['Tt = 100*(TT/',num2str(eco.TT,formatSpec),'-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['tkt = 100*(tauk -',num2str(eco.tk,formatSpec),');  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['tlt = 100*(taul -',num2str(eco.tl,formatSpec),');  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['mut = 100*(mu/',num2str(Plann.mu,formatSpec),'-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['PIt = 100*(PI-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['GDPt = 100*(GDP/',num2str(Y,formatSpec),'-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['Zt = -ut;  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['Lt = 100*(Ltot/',num2str(sum(eco.ytype.*eco.lp),formatSpec),'-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
    str = ['TsYt = Tt - GDPt;  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1; 
    str = ['BsYt = Bt - GDPt;  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1; 
    
    %%%%%%%%%%%%%%%% end of the model part
    str = ['end;\n\n\n\n']; fprintf(fid, str);  


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%     Steady-state model       %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    str = ['steady_state_model;\n\n']; fprintf(fid, str);   
    
    for h = 1:Nbin
        str = ['a',num2str(h),' = ',num2str(eco.aep(h)/Sp(h),formatSpec),';\n']; fprintf(fid, str);
        str = ['at',num2str(h),' = ',num2str(eco.abp(h)/Sp(h),formatSpec),';\n']; fprintf(fid, str);
        str = ['c',num2str(h),' = ',num2str(eco.cp(h)/Sp(h),formatSpec),';\n']; fprintf(fid, str);    
        str = ['l',num2str(h),' = ',num2str(eco.lp(h)/Sp(h),formatSpec),';\n']; fprintf(fid, str);    
        str = ['lambdac',num2str(h),' = ',num2str(Plann.lambdac(h),formatSpec),';\n']; fprintf(fid, str);
        str = ['lambdal',num2str(h),' = ',num2str(Plann.lambdal(h),formatSpec),';\n']; fprintf(fid, str);
        str = ['lambdatc',num2str(h),' = ',num2str(Plann.lambdact(h),formatSpec),';\n']; fprintf(fid, str);
        str = ['lambdatl',num2str(h),' = ',num2str(Plann.lambdalt(h),formatSpec),';\n']; fprintf(fid, str);
        str = ['psib',num2str(h),' = ',num2str(Plann.psib(h),formatSpec),';\n']; fprintf(fid, str);
     end;
     
    str = ['r = ', num2str(eco.R-1,formatSpec),';\n']; fprintf(fid, str);
    str = ['w = ', num2str(eco.w,formatSpec),';\n']; fprintf(fid, str);
    str = ['FK = ', num2str((1/beta-1),formatSpec),';\n']; fprintf(fid, str);
    str = ['FL = ', num2str(1*eco.w/(1-eco.tl),formatSpec),';\n']; fprintf(fid, str);
    str = ['K = ', num2str(K,formatSpec),';\n']; fprintf(fid, str);
    str = ['Ltot =' num2str(sum(eco.ytype.*eco.lp),formatSpec),';\n']; fprintf(fid, str); 
    str = ['mu = ', num2str(Plann.mu,formatSpec),';\n']; fprintf(fid, str);
    str = ['TFP = ', num2str(1),';\n']; fprintf(fid, str);
    str = ['A = ', num2str(eco.A,formatSpec),';\n']; fprintf(fid, str);
    str = ['B = ', num2str(eco.B,formatSpec),';\n']; fprintf(fid, str);
    str = ['u = ', num2str(0),';\n']; fprintf(fid, str);

    str = ['tauk = ', num2str(eco.tk,formatSpec),';\n']; fprintf(fid, str);
    str = ['taul = ', num2str(eco.tl,formatSpec),';\n']; fprintf(fid, str);
    str = ['Ctot = ', num2str(sum(eco.cp),formatSpec),';\n']; fprintf(fid, str);
    str = ['GDP = ', num2str(Y,formatSpec),';\n']; fprintf(fid, str);
    str = ['TT = ', num2str(eco.TT,formatSpec),';\n']; fprintf(fid, str);

    str = ['PI   = 1;\n']; fprintf(fid, str);
    str = ['zeta = 1;\n']; fprintf(fid, str);
    str = ['Lgamma = ', num2str(Plann.Lgamma,formatSpec),';\n']; fprintf(fid, str);
    str = ['Gamma = ', num2str(Plann.Gamma,formatSpec),';\n']; fprintf(fid, str);
    str = ['Upsilon = ', num2str(Plann.Upsilon,formatSpec),';\n']; fprintf(fid, str);

    str = ['M = 1 ;\n']; fprintf(fid, str);
    str = ['rtildeK = r/(1-tauk);\n']; fprintf(fid, str);
    str = ['RBN = 1+rtildeK;\n']; fprintf(fid, str);
    str = ['wtilde = FL;\n']; fprintf(fid, str);

    for i = 1:length(eco.indcc)
        str = ['LMCC',num2str(eco.indcc(i)),' = ', num2str(eco.Resh(eco.indcc(i)),formatSpec),';\n']; fprintf(fid, str);
    end;

    str = ['ut   = 0;\n']; fprintf(fid, str);  
    str = ['Kt   = 0;\n']; fprintf(fid, str);  
    str = ['Ct   = 0;\n']; fprintf(fid, str);  
    str = ['Bt   = 0;\n']; fprintf(fid, str);  
    str = ['wt   = 0;\n']; fprintf(fid, str);  
    str = ['rt   = 0;\n']; fprintf(fid, str); 
    str = ['rnt  = 0;\n']; fprintf(fid, str); 
    str = ['Tt   = 0;\n']; fprintf(fid, str);  
    str = ['tkt  = 0;\n']; fprintf(fid, str);  
    str = ['tlt  = 0;\n']; fprintf(fid, str); 
    str = ['mut  = 0;\n']; fprintf(fid, str); 
    str = ['Zt   = 0;\n']; fprintf(fid, str);
    str = ['PIt  = 0;\n']; fprintf(fid, str);
    str = ['Lt   = 0;\n']; fprintf(fid, str); 
    str = ['GDPt = 0;\n']; fprintf(fid, str);
    str = ['TsYt = 0;\n']; fprintf(fid, str);  
    str = ['BsYt = 0;\n']; fprintf(fid, str);  

    
    str = ['end;\n\n']; fprintf(fid, str);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%   Steady-state computation    %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if !bprint
        str = ['options_.noprint=1;\n']; fprintf(fid, str); %1 removes printing of steady
    end
    str = ['steady(nocheck);\n\n']; fprintf(fid, str);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%       Aggregate shock        %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    str = ['shocks;\n var eps = 0.00001;\nend;\n'];fprintf(fid, str);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%         Stoch simul          %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if bprint
        str = ['stoch_simul (order=1,irf=',num2str(LengthIRF),',nograph,nofunctions,periods=',num2str(LengthSIM), ...
            ') TFP GDPt Tt Bt PIt rt rnt Lt Ct wt taul tauk Kt Ctot K GDP B Ltot TT BsYt; \n']; fprintf(fid, str);
    else
        str = ['stoch_simul (order=1,irf=',num2str(LengthIRF),',noprint,nograph,nofunctions,periods=',num2str(LengthSIM), ...
            ') TFP GDPt Tt Bt PIt rt rnt Lt Ct wt taul tauk Kt Ctot K GDP B Ltot TT BsYt; \n']; fprintf(fid, str);
    end
    %%if needed add model diagnostics
    % str = ['model_diagnostics;']; fprintf(fid, str);
    
    if isOctave
        fflush(fid);
    end
        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%        Launch Dynare         %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    command = ['dynare ',file, ' noclearall'];
    eval(command);
    
    rnv(Economy,:)= [rnt_eps];
    PIv(Economy,:) = [PIt_eps];
    Yv(Economy,:)  = [GDPt_eps];
    Bv(Economy,:)  = [BsYt_eps];
    Cv(Economy,:)  = [Ct_eps];
    rv(Economy,:) = [rt_eps];
    taulv(Economy,:) = [100*taul_eps];
    taukv(Economy,:) = [100*tauk_eps];
end
%%



%%% Plotting figures
X=0:1:(LengthIRF-1);

if Npanels == 3

    figure(1);
    subplot(1,3,1);
    plot(X,PIv(1,:),'k-','LineWidth',1.5,X,PIv(2,:),'b--','LineWidth',1.5);
    ylim([-.1 0.3]);
    title('1. Inflation');

    subplot(1,3,2);
    plot(X,Cv(1,:),'k-','LineWidth',1.5,X,Cv(2,:),'b--','LineWidth',1.5);
    ylim([-0.2 .0]);
    title('2. Consumption');

    subplot(1,3,3);
    plot(X,Yv(1,:),'k-','LineWidth',1.5,X,Yv(2,:),'b--','LineWidth',1.5);
    ylim([-0.5 .0]);
    title('3. GDP');

    pos = get(gcf, 'Position');
    set(gcf, 'Position',[ 616   598   560   180])

    saveas(gcf,'IRFs_alternative_calib_3graphs_1.png')

    figure(2);

    subplot(1,3,1);
    plot(X,PIv(3,:),'k-','LineWidth',1.5,X,PIv(2,:),'b--','LineWidth',1.5);
    ylim([-.1 0.3]);
    title('1. Inflation');

    subplot(1,3,2);
    plot(X,Cv(3,:),'k-','LineWidth',1.5,X,Cv(2,:),'b--','LineWidth',1.5);
    ylim([-0.2 .0]);
    title('2. Consumption');

    subplot(1,3,3);
    plot(X,Yv(3,:),'k-','LineWidth',1.5,X,Yv(2,:),'b--','LineWidth',1.5);
    ylim([-0.5 .0]);
    title('3. GDP');

    pos = get(gcf, 'Position');
    set(gcf, 'Position',[ 616   598   560   180])

    saveas(gcf,'IRFs_alternative_calib_3graphs_2.png')

else

    subplot(4,2,1);
    plot(X,Yv(1,:),'k-','LineWidth',1.5,X,Yv(2,:),'r--','LineWidth',1.5,X,Yv(3,:),'b--','LineWidth',1.5);
    ylim([-0.6 0.1]);
    title('1. GDP');

    subplot(4,2,2);
    plot(X,PIv(1,:),'k-','LineWidth',1.5,X,PIv(2,:),'r--','LineWidth',1.5,X,PIv(3,:),'b--','LineWidth',1.5);
    ylim([-0.2 0.2]);
    title('2. Inflation');

    subplot(4,2,3);
    plot(X,rnv(1,:),'k-','LineWidth',1.5,X,rnv(2,:),'r--','LineWidth',1.5,X,rnv(3,:),'b--','LineWidth',1.5);
    ylim([-0.2 0.2]);
    title('3. Nominal interest rate');

    subplot(4,2,4);
    plot(X,rv(1,:),'k-','LineWidth',1.5,X,rv(2,:),'r--','LineWidth',1.5,X,rv(3,:),'b--','LineWidth',1.5);
    ylim([-0.05 0.05]);
    title('4. Real Rate');

    subplot(4,2,5);
    plot(X,Cv(1,:),'k-','LineWidth',1.5,X,Cv(2,:),'r--','LineWidth',1.5,X,Cv(3,:),'b--','LineWidth',1.5);
    ylim([-.4 .2]);
    title('5. Consumption');

    subplot(4,2,6);
    plot(X,Bv(1,:),'k-','LineWidth',1.5,X,Bv(2,:),'r--','LineWidth',1.5,X,Bv(3,:),'b--','LineWidth',1.5);
    ylim([-.1 1]);
    title('6. Debt over GDP');

    subplot(4,2,7);
    plot(X,taulv(1,:),'k-','LineWidth',1.5,X,taulv(2,:),'r--','LineWidth',1.5,X,taulv(3,:),'b--','LineWidth',1.5);
    ylim([-1. 1.0]);
    title('7. Tax on labor');

    subplot(4,2,8);
    plot(X,taukv(1,:),'k-','LineWidth',1.5,X,taukv(2,:),'r--','LineWidth',1.5,X,taukv(3,:),'b--','LineWidth',1.5);
    ylim([-1. 1.]);
    title('8. Tax on capital');

    pos = get(gcf, 'Position');
    set(gcf, 'Position',[ 616   600   800   680])

    saveas(gcf,'IRFs_alternative_calib.png')
end

