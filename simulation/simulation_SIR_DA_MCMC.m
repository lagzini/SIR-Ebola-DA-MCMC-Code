%% ========================================================================
% simulation_SIR_DA_MCMC.m
% SIR partielly observed : simulation Gillespie + DA-MCMC + MLE
% Diagnostics : ACF, ESS, IACT, MCSE, split-Rhat, taux d'acceptation
% Évaluation optionnelle : biais, MSE, RMSE, couverture et largeur moyenne
%
% Correction principale : la vraisemblance utilise les INCREMENTS
%   dN_I(j) = I(j)-I(j-1) + R(j)-R(j-1)
%   dN_R(j) = R(j)-R(j-1)
% Title =" Bayesian and Frequentist Inference for Partially
%              Observed Stochastic SIR Epidemic Models with
%                      Application to Ebola in Sierra Leone"
%
% authors=" Hamid El Maroufy, Abdelati Lagzini, and Abdelkrim Merbouha"
%% ========================================================================

clear; close all; clc;
rng(12345,'twister');

%% ==================== 1. CONFIGURATION ====================
cfg.Npop       = 100000;
cfg.i0         = 10;
cfg.beta_true  = 0.50;
cfg.gamma_true = 0.490;
cfg.R0_true    = cfg.beta_true/cfg.gamma_true;
cfg.tmax       = 160;
cfg.m          = 4;             % step dt = 1/m day

cfg.H          = 20000;
cfg.burn_in    = 5000;
cfg.r_window   = 5;
cfg.log_freq   = 5000;

% Gamma priors (shape, rate). Moderately low priors around the true
% values only for the simulation study.

cfg.m_beta      = 2;
cfg.lambda_beta = 2/cfg.beta_true;
cfg.m_gamma      = 2;
cfg.lambda_gamma = 2/cfg.gamma_true;

cfg.acf_max_lag    = 100;
cfg.n_split_chains = 4;

% Etude Monte-Carlo : mettre true pour bias/MSE/coverage.
RUN_MONTE_CARLO = false;
B = 50;

% To quickly test the Monte Carlo study, temporarily reduce:
MC_H       = 8000;
MC_burn_in = 2000;

outdir = 'results_simulation_SIR';
if ~exist(outdir,'dir'), mkdir(outdir); end

%% ==================== 2. A COMPLETE SIMULATION ====================
result = run_one_simulation(cfg, 12345, true);

fprintf('\n======================= MAIN RESULTS =======================\n');
fprintf('M==%.1f \n', cfg.m);
fprintf('Valeur vraie : beta=%.4f | gamma=%.4f | R0=%.4f\n', ...
    cfg.beta_true, cfg.gamma_true, cfg.R0_true);
fprintf('Bayes R0     : mean=%.4f | median=%.4f | 95%% CrI=[%.4f, %.4f]\n', ...
    result.BayesMean, result.BayesMedian, result.BayesCI(1), result.BayesCI(2));
fprintf('MLE R0       : mean=%.4f | median=%.4f | interval latent=[%.4f, %.4f]\n', ...
    result.MLEMean, result.MLEMedian, result.MLECI(1), result.MLECI(2));
fprintf('Acceptance   : %.2f%%\n',100*result.AcceptanceRate);
disp(result.DiagnosticsTable);

save(fullfile(outdir,'single_simulation_results.mat'),'result','cfg');
writetable(result.DiagnosticsTable, ...
    fullfile(outdir,'single_simulation_diagnostics.csv'));

%% ==================== 3. ETUDE MONTE-CARLO OPTIONNELLE ====================
if RUN_MONTE_CARLO
    cfgMC = cfg;
    cfgMC.H       = MC_H;
    cfgMC.burn_in = MC_burn_in;

    BayesEst = nan(B,1); BayesLo = nan(B,1); BayesHi = nan(B,1);
    MLEEst   = nan(B,1); MLELo   = nan(B,1); MLEHi   = nan(B,1);
    ESSvec   = nan(B,1); RhatVec = nan(B,1); AccVec = nan(B,1);

    fprintf('\n================  MONTE-CARLO study: B=%d ================\n',B);
    for b = 1:B
        try
            rb = run_one_simulation(cfgMC, 10000+b, false);

            BayesEst(b) = rb.BayesMean;
            BayesLo(b)  = rb.BayesCI(1);
            BayesHi(b)  = rb.BayesCI(2);

            MLEEst(b) = rb.MLEMean;
            MLELo(b)  = rb.MLECI(1);
            MLEHi(b)  = rb.MLECI(2);

            ESSvec(b)   = rb.ESS_R0;
            RhatVec(b)  = rb.Rhat_R0;
            AccVec(b)   = rb.AcceptanceRate;

            fprintf('Replication %3d/%3d | Bayes=%.3f | MLE=%.3f | ESS=%.0f | Rhat=%.3f\n', ...
                b,B,BayesEst(b),MLEEst(b),ESSvec(b),RhatVec(b));
        catch ME
            warning('Réplication %d ignorée : %s',b,ME.message);
        end
    end

    MetricsBayes = evaluation_metrics(BayesEst,BayesLo,BayesHi, ...
        cfg.R0_true,"Bayesian posterior mean");
    MetricsMLE = evaluation_metrics(MLEEst,MLELo,MLEHi, ...
        cfg.R0_true,"Latent-path MLE mean");
    MetricsTable = [MetricsBayes; MetricsMLE];

    ValidReplications = sum(isfinite(BayesEst));
    MeanESS = mean(ESSvec,'omitnan');
    MedianRhat = median(RhatVec,'omitnan');
    MeanAcceptance = mean(AccVec,'omitnan');
    MonteCarloDiagnostics = table(ValidReplications,MeanESS,MedianRhat,MeanAcceptance);

    fprintf('\n================ METRIQUES MONTE-CARLO ================\n');
    disp(MetricsTable);
    disp(MonteCarloDiagnostics);

    writetable(MetricsTable,fullfile(outdir,'simulation_evaluation_metrics.csv'));
    writetable(MonteCarloDiagnostics,fullfile(outdir,'monte_carlo_diagnostics.csv'));
    save(fullfile(outdir,'monte_carlo_results.mat'), ...
        'BayesEst','BayesLo','BayesHi','MLEEst','MLELo','MLEHi', ...
        'ESSvec','RhatVec','AccVec','MetricsTable','MonteCarloDiagnostics','cfgMC');
end

%% ========================================================================
%%                              FuNCTIONS
%% ========================================================================

function result = run_one_simulation(cfg, seed, makeFigures)
    rng(seed,'twister');

    %% A. Données exactes de Gillespie
    Data = Gillespie_SIR(cfg.Npop,cfg.beta_true,cfg.gamma_true,cfg.i0,cfg.tmax);
    Xgrid = interpolate_previous(Data,cfg.m);

    t_grid = Xgrid(:,1);
    I_true = round(Xgrid(:,2));
    R_cont = Xgrid(:,4);
    R_obs  = round(R_cont);
    Nobs   = numel(t_grid);

    if Nobs < 4
        error('Simulated trajectory is too short.');
    end

    % The endpoints are assumed to be known; the interior of I is latent.
    I_curr = initialize_latent_path(I_true,R_obs,cfg.Npop);
    I_curr(1)   = cfg.i0;
    I_curr(end) = I_true(end);

    
    if ~is_valid_path(I_curr,R_obs,cfg.Npop)
        I_curr = I_true;
    end
    if ~is_valid_path(I_curr,R_obs,cfg.Npop)
        error('La trajectoire initiale est incompatible avec le modèle SIR.');
    end

    [ai_vec,bi_vec] = local_ab(I_curr,R_cont,t_grid,cfg.Npop);
    [NI,NR,C1,C2] = sufficient_stats(I_curr,R_obs,ai_vec,bi_vec);

    beta_curr  = gamma_random(NI+cfg.m_beta, C1+cfg.lambda_beta);
    gamma_curr = gamma_random(NR+cfg.m_gamma,C2+cfg.lambda_gamma);
    R0_curr = beta_curr/max(gamma_curr,realmin);

    %% B. Stockage
    H = cfg.H;
    R0_samples   = zeros(H,1);
    beta_samples = zeros(H,1);
    gamma_samples= zeros(H,1);
    R0_MLE_samples = zeros(H,1);
    I_samples = zeros(Nobs,H,'single');

    nProposal = 0;
    nAccept   = 0;

    %% C. DA-MCMC :  MH single-site updates
    for h = 1:H
      
        order = 2:Nobs-1;
        order = order(randperm(numel(order)));

        for jj = order
            currentValue = I_curr(jj);
            lower = max(0,currentValue-cfg.r_window);
            upper = min(cfg.Npop-R_obs(jj),currentValue+cfg.r_window);

            candidates = lower:upper;
            candidates(candidates==currentValue) = [];
            if isempty(candidates), continue; end

            proposedValue = candidates(randi(numel(candidates)));
            I_star = I_curr;
            I_star(jj) = proposedValue;

            % Le point jj affecte les intervalles jj et jj+1.
            intervalIndex = unique([jj,min(jj+1,Nobs)]);
            intervalIndex = intervalIndex(intervalIndex>=2 & intervalIndex<=Nobs);

            logCurrent = local_loglik(I_curr,R_obs,R_cont,t_grid,cfg.Npop, ...
                beta_curr,gamma_curr,intervalIndex);
            logProposed = local_loglik(I_star,R_obs,R_cont,t_grid,cfg.Npop, ...
                beta_curr,gamma_curr,intervalIndex);

            % Correction Hastings due à la troncature aux frontières.
            reverseLower = max(0,proposedValue-cfg.r_window);
            reverseUpper = min(cfg.Npop-R_obs(jj),proposedValue+cfg.r_window);
            reverseCandidates = reverseLower:reverseUpper;
            reverseCandidates(reverseCandidates==proposedValue) = [];

            logHastings = log(numel(candidates)) - log(max(numel(reverseCandidates),1));
            logAlpha = logProposed-logCurrent+logHastings;

            nProposal = nProposal+1;
            if isfinite(logAlpha) && log(rand)<min(0,logAlpha)
                I_curr = I_star;
                nAccept = nAccept+1;
            end
        end

        %% D. Tirages conditionnels corrects de beta et gamma
        [ai_vec,bi_vec] = local_ab(I_curr,R_cont,t_grid,cfg.Npop);
        [NI,NR,C1,C2] = sufficient_stats(I_curr,R_obs,ai_vec,bi_vec);

        beta_curr  = gamma_random(NI+cfg.m_beta, C1+cfg.lambda_beta);
        gamma_curr = gamma_random(NR+cfg.m_gamma,C2+cfg.lambda_gamma);
        R0_curr = beta_curr/max(gamma_curr,realmin);

        % MLE conditionnel au chemin latent courant
        beta_mle  = NI/max(C1,realmin);
        gamma_mle = NR/max(C2,realmin);
        R0_mle = beta_mle/max(gamma_mle,realmin);

        R0_samples(h) = R0_curr;
        beta_samples(h) = beta_curr;
        gamma_samples(h)= gamma_curr;
        R0_MLE_samples(h)=R0_mle;
        I_samples(:,h)=single(I_curr);

        if mod(h,cfg.log_freq)==0
            first = max(cfg.burn_in+1,1);
            if h>=first
                fprintf('It. %5d/%5d | mean R0=%.4f | MLE=%.4f | accept=%.2f%%\n', ...
                    h,H,mean(R0_samples(first:h)),mean(R0_MLE_samples(first:h)), ...
                    100*nAccept/max(nProposal,1));
            end
        end
    end

    %% E. Résumés
    idxPost = (cfg.burn_in+1):H;
    post = R0_samples(idxPost);
    postBeta = beta_samples(idxPost);
    postGamma= gamma_samples(idxPost);
    postMLE = R0_MLE_samples(idxPost);

    result.BayesMean   = mean(post);
    result.BayesMedian = median(post);
    result.BayesSD     = std(post);
    result.BayesCI     = empirical_quantile(post,[0.025 0.975]);

    result.MLEMean   = mean(postMLE);
    result.MLEMedian = median(postMLE);
    result.MLESD     = std(postMLE);
    result.MLECI     = empirical_quantile(postMLE,[0.025 0.975]);

    result.AcceptanceRate = nAccept/max(nProposal,1);

    %% F. Diagnostics sans Econometrics Toolbox
    maxLag = min(cfg.acf_max_lag,numel(post)-1);
    [ESS_R0,IACT_R0,acfR0] = ess_geyer(post,maxLag);
    [ESS_beta,IACT_beta,acfBeta] = ess_geyer(postBeta,maxLag);
    [ESS_gamma,IACT_gamma,acfGamma] = ess_geyer(postGamma,maxLag);

    Rhat_R0 = split_rhat(post,cfg.n_split_chains);
    Rhat_beta = split_rhat(postBeta,cfg.n_split_chains);
    Rhat_gamma= split_rhat(postGamma,cfg.n_split_chains);

    MCSE_R0 = std(post)/sqrt(max(ESS_R0,1));
    MCSE_beta = std(postBeta)/sqrt(max(ESS_beta,1));
    MCSE_gamma= std(postGamma)/sqrt(max(ESS_gamma,1));

    Parameter = ["R0";"beta";"gamma"];
    PosteriorMean = [mean(post);mean(postBeta);mean(postGamma)];
    PosteriorSD = [std(post);std(postBeta);std(postGamma)];
    ESS = [ESS_R0;ESS_beta;ESS_gamma];
    IACT= [IACT_R0;IACT_beta;IACT_gamma];
    MCSE= [MCSE_R0;MCSE_beta;MCSE_gamma];
    SplitRhat=[Rhat_R0;Rhat_beta;Rhat_gamma];
    result.DiagnosticsTable = table(Parameter,PosteriorMean,PosteriorSD,ESS,IACT,MCSE,SplitRhat);

    result.ESS_R0 = ESS_R0;
    result.Rhat_R0= Rhat_R0;
    result.ACF_R0 = acfR0;
    result.R0_samples=R0_samples;
    result.R0_MLE_samples=R0_MLE_samples;
    result.beta_samples=beta_samples;
    result.gamma_samples=gamma_samples;
    result.I_samples=I_samples;
    result.t_grid=t_grid;
    result.I_true=I_true;
    result.R_obs=R_obs;
    result.Data=Data;

    %% G. Figures
    if makeFigures
        figure(1); clf;
        tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

        nexttile;
        histogram(post,50,'Normalization','pdf','EdgeColor','none');
        xline(cfg.R0_true,'r--','LineWidth',1.5);
        xlabel('R_0'); ylabel('Density'); title('(a) Posterior density'); grid on;

        nexttile;
        histogram(postMLE,50,'Normalization','pdf','EdgeColor','none'); hold on;
        xline(cfg.R0_true,'r--','LineWidth',1.5);
        xlabel('R_0'); ylabel('Density'); title('(b) MLE distribution over latent paths'); grid on;

        nexttile;
        plot(R0_samples,'k-','LineWidth',0.7); hold on;
        yline(cfg.R0_true,'r--','LineWidth',1.3);
        xline(cfg.burn_in,'b--'); grid on;
        xlabel('Iteration'); ylabel('R_0'); title('(c) Trace plot');

        nexttile;
        plot(cumsum(post)./(1:numel(post))','LineWidth',1.2); hold on;
        yline(cfg.R0_true,'r--','LineWidth',1.3); grid on;
        xlabel('Post burn-in iteration'); ylabel('Cumulative mean'); title('(d) Cumulative mean');

        nexttile;
        lags=0:(numel(acfR0)-1);
        stem(lags,acfR0,'filled','MarkerSize',3); hold on;
        bound=1.96/sqrt(numel(post));
        yline(bound,'r--'); yline(-bound,'r--'); yline(0,'k-');
        xlim([0 maxLag]); ylim([-1 1]); grid on;
        xlabel('Lag'); ylabel('ACF'); title('(e) ACF');

        nexttile;
        axis off;
        text(0,0.85,sprintf('True R_0 = %.4f',cfg.R0_true),'FontSize',11);
        text(0,0.68,sprintf('Bayes mean = %.4f',result.BayesMean),'FontSize',11);
        text(0,0.51,sprintf('95%% CrI = [%.4f, %.4f]',result.BayesCI),'FontSize',11);
        text(0,0.34,sprintf('ESS = %.0f; split-Rhat = %.4f',ESS_R0,Rhat_R0),'FontSize',11);
        text(0,0.17,sprintf('Acceptance = %.2f%%',100*result.AcceptanceRate),'FontSize',11);
        title('(f) Numerical summary');
        sgtitle('DA-MCMC diagnostics for simulated partially observed SIR model');

        % Trajectoires latentes
        figure(2); clf; hold on;
        nDraw=10;
        drawIndex=idxPost(round(linspace(1,numel(idxPost),nDraw)));
        for k=1:nDraw
            plot(t_grid,double(I_samples(:,drawIndex(k))),'LineWidth',0.7);
        end
        Imean=mean(double(I_samples(:,idxPost)),2);
        Ilo=column_quantile(double(I_samples(:,idxPost)),0.025);
        Ihi=column_quantile(double(I_samples(:,idxPost)),0.975);
        fill([t_grid;flipud(t_grid)],[Ilo;flipud(Ihi)],[0.8 0.8 0.8], ...
            'FaceAlpha',0.35,'EdgeColor','none');
        plot(t_grid,Imean,'k-','LineWidth',2.2);
        plot(t_grid,I_true,'r--','LineWidth',2);
        xlabel('Time'); ylabel('Infectious individuals'); grid on;
        legend({'Posterior draws','','','','','','','','','', ...
            '95% credible band','Posterior mean','Gillespie truth'},'Location','best');
        title('Latent infectious trajectory reconstruction');
    end
    %% ============================================================
% FIGURE SEPAREE 1 : MOYENNE CUMULATIVE DE R0
% ============================================================
hold on;
post_plot = post(:);   % forcer un vecteur colonne

cumulative_mean_R0 = cumsum(post_plot) ./ (1:numel(post_plot))';

figure(10); clf;
plot(1:numel(post_plot), cumulative_mean_R0, ...
    'LineWidth',1.2);
hold on;

yline(cfg.R0_true,'r--','LineWidth',1.3);

grid on;
box on;

xlabel('Post burn-in iteration');
ylabel('Cumulative mean of R_0');
title('Cumulative mean of R_0');

legend({'Cumulative mean','True R_0'}, ...
    'Location','best');
%% ============================================================
% FIGURE SEPAREE 2 : AUTOCORRELATION ACF DE R0
% ============================================================

lags = 0:(numel(acfR0)-1);
confidence_bound = 1.96 / sqrt(numel(post_plot));

figure(11); clf;

stem(lags, acfR0, ...
    'filled', ...
    'MarkerSize',3);
hold on;

yline(confidence_bound, ...
    'r--', ...
    'LineWidth',1.2);

yline(-confidence_bound, ...
    'r--', ...
    'LineWidth',1.2);

yline(0, ...
    'k-', ...
    'LineWidth',0.8);

xlim([0 maxLag]);
ylim([-1 1]);

grid on;
box on;

xlabel('Lag');
ylabel('Autocorrelation');
title('ACF of the post burn-in R_0 chain');

legend({'ACF','95% bounds'}, ...
    'Location','best');
end

function Data = Gillespie_SIR(N,beta,gamma,i0,tmax)
    t=0; S=N-i0; I=i0; R=0;
    Data=[t I S R];

    while I>0 && t<tmax
        rateInf=beta*S*I/N;
        rateRec=gamma*I;
        totalRate=rateInf+rateRec;
        if totalRate<=0, break; end

        tau=-log(max(rand,realmin))/totalRate;
        if t+tau>tmax, break; end

        if rand<rateInf/totalRate
            S=S-1; I=I+1;
        else
            I=I-1; R=R+1;
        end
        t=t+tau;
        Data(end+1,:)=[t I S R]; %#ok<AGROW>
    end

    % Ajouter l'état final à tmax pour une interpolation cohérente.
    if Data(end,1)<tmax
        Data(end+1,:)=[tmax I S R];
    end
end




function X=interpolate_previous(Data,m)
    t=Data(:,1); I=Data(:,2); R=Data(:,4);
    N=Data(1,2)+Data(1,3)+Data(1,4);
    dt=1/m;
    tgrid=(0:dt:t(end))';
    if tgrid(end)<t(end), tgrid(end+1)=t(end); end

    Igrid=interp1(t,I,tgrid,'previous','extrap');
    Rgrid=interp1(t,R,tgrid,'previous','extrap');
    Igrid=max(round(Igrid),0);
    Rgrid=cummax(max(round(Rgrid),0));
    Sgrid=max(N-Igrid-Rgrid,0);
    X=[tgrid Igrid Sgrid Rgrid];
end

function I0=initialize_latent_path(Itrue,Robs,N)
    n=numel(Itrue);
    I0=zeros(n,1);
    I0(1)=Itrue(1);
    I0(end)=Itrue(end);

  
    base=movmean(double(Itrue),5);
    noise=round(0.05*max(base,1).*randn(n,1));
    I0=round(max(base+noise,0));
    I0=min(I0,N-Robs);
    I0(1)=Itrue(1);
    I0(end)=Itrue(end);

  
    for j=2:n
        dR=Robs(j)-Robs(j-1);
        lower=max(0,I0(j-1)-dR);
        I0(j)=max(I0(j),lower);
        I0(j)=min(I0(j),N-Robs(j));
    end
    I0(end)=Itrue(end);
end

function tf=is_valid_path(I,R,N)
    dR=diff(R);
    dInf=diff(I)+dR;
    S=N-I-R;
    tf=all(I>=0) && all(R>=0) && all(S>=0) && ...
       all(dR>=0) && all(dInf>=0) && all(abs(I-round(I))<1e-10);
end

function ll=local_loglik(I,RcontDiscrete,Rcont,t,N,beta,gamma,index)
    %#ok<INUSD> RcontDiscrete is the discrete observed R vector.
    ll=0;
    for j=index(:)'
        dR=RcontDiscrete(j)-RcontDiscrete(j-1);
        dInf=(I(j)-I(j-1))+dR;

        if dR<0 || dInf<0 || I(j)<0 || I(j)>N-RcontDiscrete(j)
            ll=-inf; return;
        end

        dt=t(j)-t(j-1);
        S1=max(N-Rcont(j-1)-I(j-1),0);
        S2=max(N-Rcont(j)-I(j),0);
        ai=(dt/(2*N))*(I(j-1)*S1+I(j)*S2);
        bi=(dt/2)*(I(j-1)+I(j));

        ll=ll+log_poisson(dInf,beta*ai)+log_poisson(dR,gamma*bi);
        if ~isfinite(ll), return; end
    end
end

function lp=log_poisson(k,lambda)
    if k<0 || abs(k-round(k))>1e-10 || lambda<0 || ~isfinite(lambda)
        lp=-inf;
    elseif lambda==0
        if k==0, lp=0; else, lp=-inf; end
    else
        lp=-lambda+k*log(lambda)-gammaln(k+1);
    end
end

function [a,b]=local_ab(I,R,t,N)
    n=numel(t); a=zeros(n,1); b=zeros(n,1);
    for j=2:n
        dt=t(j)-t(j-1);
        S1=max(N-R(j-1)-I(j-1),0);
        S2=max(N-R(j)-I(j),0);
        a(j)=(dt/(2*N))*(I(j-1)*S1+I(j)*S2);
        b(j)=(dt/2)*(I(j-1)+I(j));
    end
end

function [NI,NR,C1,C2]=sufficient_stats(I,R,a,b)
    dR=diff(R);
    dInf=diff(I)+dR;
    if any(dR<0) || any(dInf<0)
        error('Trajectoire SIR non valide : incréments négatifs.');
    end
    NI=sum(dInf);
    NR=sum(dR);
    C1=sum(a(2:end));
    C2=sum(b(2:end));
end

function x=gamma_random(shape,rate)
    if shape<=0 || rate<=0
        error('Paramètres Gamma invalides.');
    end
    x=gamrnd(shape,1/rate);
end

function q=empirical_quantile(x,p)
    x=sort(x(isfinite(x)));
    n=numel(x);
    q=zeros(size(p));
    for k=1:numel(p)
        pos=1+(n-1)*p(k);
        lo=floor(pos); hi=ceil(pos);
        if lo==hi
            q(k)=x(lo);
        else
            q(k)=x(lo)+(pos-lo)*(x(hi)-x(lo));
        end
    end
end

function q=column_quantile(X,p)
    % Quantile ligne par ligne, sans dépendre de quantile(X,2).
    q=zeros(size(X,1),1);
    for i=1:size(X,1)
        q(i)=empirical_quantile(X(i,:),p);
    end
end

function [ess,tau,acfValues]=ess_geyer(x,maxLag)
    x=x(:); x=x(isfinite(x)); n=numel(x);
    acfValues=sample_acf(x,maxLag);
    pairSums=[];
    for k=1:floor(maxLag/2)
        s=acfValues(2*k)+acfValues(2*k+1);
        if ~isfinite(s) || s<=0, break; end
        pairSums(end+1,1)=s; %#ok<AGROW>
    end
    for k=2:numel(pairSums)
        pairSums(k)=min(pairSums(k),pairSums(k-1));
    end
    tau=max(1,1+2*sum(pairSums));
    ess=min(n,n/tau);
end

function acfValues=sample_acf(x,maxLag)
    x=x(:)-mean(x);
    n=numel(x); maxLag=min(maxLag,n-1);
    den=sum(x.^2);
    acfValues=zeros(maxLag+1,1); acfValues(1)=1;
    if den<=0, return; end
    for lag=1:maxLag
        acfValues(lag+1)=sum(x(1:n-lag).*x(1+lag:n))/den;
    end
end

function rhat=split_rhat(x,nSegments)
    % Split-Rhat classique sans tiedrank/norminv.
    x=x(:); x=x(isfinite(x));
    nPer=floor(numel(x)/nSegments);
    if nPer<8, rhat=NaN; return; end
    X=reshape(x(1:nPer*nSegments),nPer,nSegments);
    half=floor(nPer/2);
    X=[X(1:half,:) X(end-half+1:end,:)];
    n=size(X,1);
    W=mean(var(X,0,1));
    B=n*var(mean(X,1),0,2);
    if W<=0, rhat=NaN; return; end
    varPlus=((n-1)/n)*W+B/n;
    rhat=sqrt(varPlus/W);
end

function Metrics=evaluation_metrics(est,lo,hi,trueValue,methodName)
    valid=isfinite(est)&isfinite(lo)&isfinite(hi);
    est=est(valid); lo=lo(valid); hi=hi(valid);
    if isempty(est), error('Aucune réplication valide.'); end

    Method=string(methodName);
    Replications=numel(est);
    TrueValue=trueValue;
    Bias=mean(est-trueValue);
    AbsoluteBias=mean(abs(est-trueValue));
    MSE=mean((est-trueValue).^2);
    RMSE=sqrt(MSE);
    CoverageProbability=mean(lo<=trueValue & trueValue<=hi);
    AverageIntervalWidth=mean(hi-lo);

    Metrics=table(Method,Replications,TrueValue,Bias,AbsoluteBias,MSE,RMSE, ...
        CoverageProbability,AverageIntervalWidth);
end
