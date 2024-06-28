clear

calib = 'rob'


warning('off','all');
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

isOctave && warning('off', 'Octave:array-as-logical');


slope_w_rob = [0.01, 0.03, 0.18, 0.35] ;


LengthIRF = 20 ; 



for Economy = 1:length(slope_w_rob)

    psi_rob = slope_w_rob(Economy)


    psi_floor = floor(psi_rob*100); 
    filename_root = sprintf('todynare_SW_rob_%d', psi_floor);
    fileExtension = '.mat';
    
    FileNameTot = [filename_root, fileExtension];
    load(FileNameTot)

    trunc  = Ramsey.truncatedModel.truncatedAllocation;
    lagMul = Ramsey.lagrangeMult;
    
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
    
    P     = ones(Ntot,1);
    P(trunc.ind_cc_h) = 0;
    omegab_h = ones(Ntot,1);
    formatSpec = '%25.20f';
    



    coefBL = 4.0;
    coefBu = 8.5;
    coefl  = 0.;   % reaction of labor tax to debt
    coefk  = 0.;
    coefkL = 0.;
    coeflL = 0.;




    file = ['code_dynare_SW_',calib,num2str(Economy),'.mod'];
    fid = fopen(file, 'w') ;
    formatSpec = '%25.20f';

    omegab_h = ones(Ntot,1) ;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%     variables and params     %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    str = ['var\n'] ; fprintf(fid, str);

    str = ['u  r  w  A  L  mu  Ctot Y  piw  pip  T  Util  B  G  Z IUtil \n'] ; fprintf(fid, str);
    str = ['ut rt wt At Lt mut Ct Yt piwt pipt Tt Utilt Bt TsYt \n'] ; fprintf(fid, str);
    str = ['RBN tauk taul K mkupw gammaw Lambda LGamma upsilon \n'] ; fprintf(fid, str);
    str = ['RBNt Kt BsYt \n'] ; fprintf(fid, str);



    for h = 1:Ntot
        str = ['a',num2str(h),' '] ; fprintf(fid,str);     %  end-of-period
        str = ['at',num2str(h),' '] ; fprintf(fid, str);   % beginning-of-period
        str = ['c',num2str(h),' '] ; fprintf(fid, str);
        str = ['psib',num2str(h),' '] ; fprintf(fid, str);
        str = ['lambdacb',num2str(h),' '] ; fprintf(fid, str);
        str = ['lambdacbt',num2str(h),' '] ; fprintf(fid, str);    % last period
        if mod(h,10)==0; %multiple de 10
            str = ['\n'] ;         fprintf(fid, str);
        end;

    end;

    str = ['; \n\n'] ; fprintf(fid, str);
    str = ['varexo eps; \n\n'] ; fprintf(fid, str);
    str = ['parameters\n'] ; fprintf(fid, str);
    str = ['beta alpha phi abar delta gamma chi psiw epsw rho_u Gss \n\n'] ; fprintf(fid, str);
    % str = ['coefl coefk coeft coeft2;\n\n'] ; fprintf(fid, str);
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

    



    %%
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


        %%%%%%%% definition psi_hat_bar
        str = ['psib',num2str(h),' = 0 ',' + ',num2str(S_h(h),formatSpec),'*',num2str(omegab_h(h)*xsis.xsiu0(h),formatSpec ),'*(c',num2str(h),')^-gamma - (lambdacb',num2str(h),...
                '*',num2str(xsis.xsiuE(h),formatSpec),' - (1+r)*lambdacbt',num2str(h), '*',num2str(xsis.xsiuE(h),formatSpec ),...
                ' + (epsw-1)/psiw*w*L*gammaw*',num2str(xsis.xsiu1(h)*S_h(h)*trunc.y0_h(h),formatSpec),...
                ')*(-gamma)*(c',num2str(h),')^(-gamma-1); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;



        %%%%%%%% FOC individual savings (at): Euler on psib
        if P(h)==1
            str = ['psib',num2str(h),'  = ', num2str(S_h(h),formatSpec),'*beta*(1+r(+1))*(0']; fprintf(fid, str);
            for hp = 1:Ntot
                if Pi_h(hp,h)>tol   %be careful transpose from h to hp
                    str = ['+ ',num2str(Pi_h(hp,h),formatSpec),'*psib',num2str(hp),'(+1)/',num2str(S_h(hp),formatSpec)]; fprintf(fid, str);
                end
            end
            str = [') + ',num2str(S_h(h),formatSpec),'*( beta*alpha*w(+1)*L(+1)/K*( 0']; fprintf(fid, str);
            % str = [') + ',num2str(S_h(h),formatSpec),'*( beta*alpha*w(+1)*L(+1))/K*( 0']; fprintf(fid, str);
            for hp = 1:Ntot
                if (S_h(hp))>tol   %be careful transpose from h to hp
                str = [' + ',num2str(trunc.y0_h(hp),formatSpec),'*psib',num2str(hp),'(+1)  -  gammaw(+1)*(epsw-1)/psiw*',num2str(S_h(hp)*trunc.y0_h(hp),formatSpec),'*(c',num2str(hp),'(+1))^(-gamma)']; fprintf(fid, str);
                end
                if mod(hp,10)==0;str = ['\n'];fprintf(fid, str); end;
            end  
            str = [') ) +',num2str(S_h(h),formatSpec),'*( ']; fprintf(fid, str);
            str = [ ' + beta*LGamma(+1)*( r(+1) - (1-tauk(+1))*( rt(+1)-(1-alpha)*(rt(+1)+delta) ) )']; fprintf(fid, str);
            str = [ ' + beta*mu(+1)/K*( alpha*Y(+1) - (delta + r(+1))*K  - alpha*w(+1)*L(+1) )']; fprintf(fid, str);
            str = [ ' + beta*alpha*w(+1)/K*( -Lambda(+1)*(pip(+1)+1) + beta*Lambda(+2)*(piw(+2)+1)  )']; fprintf(fid, str);
            str = [ ' + (1-alpha)/K*upsilon*( rt(+1)+delta )']; fprintf(fid, str);
            str = [ '); %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa=equa+1;
        else
            str = ['lambdacb',num2str(h),' =  0; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        end

        %%%%%%%% lambdactc : lambdac previous period
        str = ['lambdacbt',num2str(h),' = 0 ']; fprintf(fid, str);
        for hi = 1:Ntot
            if Pi_h(h,hi)>tol  %be careful transpose
                str = [' + ',num2str(Pi_h(h,hi),formatSpec),'*lambdacb',num2str(hi),'(-1)']; fprintf(fid, str);
            end
        end
        str = ['; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
        

        isOctave && fflush(fid);
    end

    sumOmegab = 0;
    for h = 1:Ntot
        sumOmegab = sumOmegab + omegab_h(h)*S_h(h);
    end


    str = ['(1/chi)*L^(1/phi) = mu*(wt-(1-alpha)*w) + (epsw/psiw)*gammaw*(1/chi)*L^(1/phi)*(1+1/phi) + (1-alpha)*w*(0']; fprintf(fid, str);
    for hp = 1:Ntot 
        str = ['-(epsw-1)/psiw*gammaw*',num2str(S_h(hp)*trunc.y0_h(hp),formatSpec),'*(c',num2str(hp),')^(-gamma)','+',num2str(trunc.y0_h(hp),formatSpec),'*psib',num2str(hp)]; fprintf(fid, str);
    end
    str = [') - alpha*wt*( beta^(-1)*upsilon(-1)/K(-1) + LGamma*(1-tauk) ) -alpha*w/L*( beta*Lambda(+1)*(piw(+1)+1) - Lambda*(pip+1)  )  ; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;


    %%%%%%%% FOC planner r
    str = ['0 = 0 ']; fprintf(fid, str);
    for h = 1:Ntot
        if S_h(h)>0
            str = [' + at',num2str(h),'*psib',num2str(h),' + lambdacbt',num2str(h),'*',num2str(xsis.xsiuE(h),formatSpec),'*(c',num2str(h),'^-gamma)']; fprintf(fid, str);
        end
        %if mod(h,10)==0;str = ['\n'] ;fprintf(fid, str); end
    end
        str = [' + (LGamma - mu)*A(-1) ; %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;


    %%%%%%%% FOC planner RBN : upsilon
    str = ['upsilon = (pip(+1)+1)*beta*B*(1-tauk(+1))*LGamma(+1)/(pip(+1)+1) ; %% equation ', num2str(equa),';\n']; fprintf(fid, str); equa = equa+1;


    %%%%%%%% FOC planner T
        str = ['mu = 0 ']; fprintf(fid, str);
        for h = 1:Ntot
            if S_h(h)>0
                str = [' + psib',num2str(h)]; fprintf(fid, str);
            end
            %if mod(h,10)==0;str = ['\n'] ;fprintf(fid, str); end
        end
        str = [';  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;



    %%%%%%%% FOC planner piw
        str = ['piw*psiw = -(gammaw-gammaw(-1))*(2*piw+1)  + Lambda*w(-1);\n']; fprintf(fid, str); equa = equa+1;

    %%%%%%%% FOC planner pip
        str = [' 0 = Lambda*w +  ( beta^(-1)*upsilon(-1) - LGamma*(1-tauk)*B(-1) )*RBN(-1)/( (pip+1)*(pip+1) ) ;\n']; fprintf(fid, str); equa = equa+1;

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
    str = [';  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;


    str = ['At = 0'] ; fprintf(fid, str);
    for h = 1:Ntot
        str = [' + ', num2str(S_h(h),formatSpec),'*at',num2str(h)] ; fprintf(fid, str);
    end;
    str = [';  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;

    %%%%%%%%%%%%%%%%%%%%      Exogenous Taxes

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
    str = ['mut = 100*(mu/',num2str(Ramsey.lagrangeMult.mu,formatSpec),'-1);  %% equation ', num2str(equa),'\n']; fprintf(fid, str); equa = equa+1;
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
        str = ['lambdacb',num2str(h),' = ',num2str(lagMul.lambdacb(h),formatSpec),';\n'] ; fprintf(fid, str);
        str = ['lambdacbt',num2str(h),' = ',num2str(lagMul.lambdact(h)*trunc.S_h(h),formatSpec),';\n'] ; fprintf(fid, str);
        str = ['psib',num2str(h),'     = ',num2str(lagMul.psib(h),formatSpec),';\n'] ; fprintf(fid, str);   % recall it is psi multiplied by the size of the bin
        isOctave && fflush(fid);
    end;

    str = ['r     = ', num2str(solution.R-1,formatSpec),';\n'] ; fprintf(fid, str);
    str = ['w     = ', num2str(solution.w,formatSpec),';\n'] ; fprintf(fid, str);
    str = ['L     = ', num2str(solution.L(1),formatSpec),';\n'] ; fprintf(fid, str);
    str = ['mu    = ', num2str(lagMul.mu,formatSpec),';\n'] ; fprintf(fid, str);
    str = ['gammaw= ', num2str(lagMul.gammaw,formatSpec),';\n'] ; fprintf(fid, str);
    str = ['upsilon= ', num2str(lagMul.upsilon,formatSpec),';\n'] ; fprintf(fid, str);
    str = ['Lambda= ', num2str(lagMul.Lambda,formatSpec),';\n'] ; fprintf(fid, str);
    str = ['LGamma= ', num2str(lagMul.Gamma,formatSpec),';\n'] ; fprintf(fid, str);
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
    str = ['mut   = 0;\n'] ; fprintf(fid, str);
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%       Aggregate shock        %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%         Stoch simul          %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            str = ['shocks;\n var eps = 0.00001;\nend;\n'];fprintf(fid, str); 
            str = ['stoch_simul (order=1,irf=',num2str(LengthIRF), ',nograph, nofunctions,noprint) Y Ctot K L RBN r T B pip piw u ut Yt pipt piwt RBNt rt Lt Ct wt taul tauk Kt BsYt TsYt;'] ; fprintf(fid, str);
    
        

    isOctave && fflush(fid);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%        Launch Dynare         %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    command = ['dynare ',file, ' noclearall'];
    eval(command);


    irfs = oo_.irfs;

    PIp_v(Economy,:) = [0  irfs.pipt_eps];
    PIw_v(Economy,:) = [0 irfs.piwt_eps];

end

    



fileName = ['To_IRFs_SW_',calib,'.mat' ];
save(fileName, "PIp_v", "PIw_v") 



load To_IRFs_SW_rob.mat
LengthIRF = 20 ; 

X=0:1:(LengthIRF);

subplot(1,2,1);
plot(X,PIw_v(1,:),'m--',X,PIw_v(2,:),'b--',X,PIw_v(3,:),'k-',X,PIw_v(4,:),'g--','LineWidth',1.5);
title('1. Wage Inflation  (\Pi^W, %)');


subplot(1,2,2);
plot(X,PIp_v(1,:),'m--',X,PIp_v(2,:),'b--',X,PIp_v(3,:),'k-',X,PIp_v(4,:),'g--','LineWidth',1.5);
ylim([-.04 .25]);
title('2. Price Inflation (\Pi^P, %)');

legend('0.01', '0.03', '0.18', '0.35', 'Location', 'southeast');
set(gcf, 'Position',[100, 100, 800, 400])
saveas(gcf,'../../Graph/Robustness_Slope_SW.png')

close 

beep; 
