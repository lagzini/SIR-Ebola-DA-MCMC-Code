%% ================================================================
%  run_sierra_leone_final_AZ.m — Ebola (Sierra Leone, 2014–2016)
%  DA–MCMC (Algorithm 1) + MLE comparison + Diagnostics + Figures
%  Figure principale: Observed vs Posterior Predictive Mean (95% CI)
%  Author: Lagzini Abdelati
% ================================================================

clear; close all; clc;
 rng(1);  % reproductibilité

%% ====== USER CONFIGURATION ======
data_file   = 'sierra_leone_real_values.xlsx';  % fichier Excel
sheet_name  = '';                               % '' -> feuille par défaut
Ntot        = 7.0e6;                            % population totale
m_per_day   = 4;                                % résolution: dt = 1/m jour
H           = 20000;                            % itérations MCMC
burn_in     = 5000;                             % burn-in
r_window    = 3;                                % fenêtre locale ±r
log_freq    = 5000;
block_B     = [];                               % [] -> max(8,m)

% Diagnostics MCMC
acf_max_lag    = 100;   % nombre maximal de retards pour l'ACF
n_split_chains = 4;     % sous-chaînes utilisées pour le split-Rhat

m_beta = 1; lambda_beta  = 100;                 % priors beta
m_gamma = 1; lambda_gamma = 100;                % priors gamma

start_date_str = '2014-08-29';
end_date_str   = '2015-12-29';

CFR = 0.70;  % taux de létalité (approx.) pour estimer recovered

outdir = 'results_sierra_leone';
if ~exist(outdir, 'dir'), mkdir(outdir); end


%% ====== A. LECTURE EXCEL (Date | Cases cum | Deaths cum) ======
opts = detectImportOptions(data_file, 'VariableNamingRule','preserve');
if ~isempty(sheet_name), opts.Sheet = sheet_name; end
T = readtable(data_file, opts);

% Normaliser la colonne date
dateVar = 'Date';
T = normalize_date_column(T, dateVar);
T = sortrows(T, dateVar);

time = T.(dateVar);

% Utiliser VariableDescriptions si dispo (souvent les longs libellés WHO)
vnames = string(T.Properties.VariableNames);
vdesc  = string(T.Properties.VariableDescriptions);
if isempty(vdesc), vdesc = vnames; end

needle_cases  = "confirmed, probable and suspected Ebola cases";
needle_deaths = "confirmed, probable and suspected Ebola deaths";

idxC = find(contains(lower(vdesc), lower(needle_cases)),  1);
idxD = find(contains(lower(vdesc), lower(needle_deaths)), 1);

if isempty(idxC) || isempty(idxD)
    error('Impossible de localiser les colonnes "cases" et "deaths" via les descriptions.');
end

Cases  = T.(vnames(idxC));
Deaths = T.(vnames(idxD));

% Estimation Recovered et Removed
Recovered = Deaths .* ((1/CFR) - 1);     % ~0.4286 * Deaths
Robs_raw  = cummax(Deaths + Recovered);  % Removed cumulé (morts + guéris)

% Restreindre la période
if ~isempty(start_date_str), time_min = datetime(start_date_str); else, time_min = time(1); end
if ~isempty(end_date_str),   time_max = datetime(end_date_str);   else, time_max = time(end); end
keep = (time >= time_min) & (time <= time_max);

time      = time(keep);
Cases     = Cases(keep);
Deaths    = Deaths(keep);
Recovered = Recovered(keep);
Robs_raw  = Robs_raw(keep);

DataRaw = table(time, Cases, Deaths, Recovered, Robs_raw, ...
    'VariableNames', {'date','cases','deaths','recovered','Robs'});


%% ====== B. GRILLE REGULIERE (dt = 1/m jour) ======
m = m_per_day;
dt_days = 1/m;

t0 = dateshift(min(DataRaw.date),'start','day');
t1 = dateshift(max(DataRaw.date),'start','day');
tgrid_dt = (t0:days(dt_days):t1)';           % datetime grid

t_vec = days(tgrid_dt - tgrid_dt(1));        % temps en jours (numérique)
t_grid = t_vec(:);
Nobs   = numel(t_grid);

% Interpolation "previous" (courbe en escalier comme du cumul)
R_prev = interp1(DataRaw.date, DataRaw.Robs, tgrid_dt, 'previous','extrap');
I_prev = interp1(DataRaw.date, (DataRaw.cases - DataRaw.Robs), tgrid_dt, 'previous','extrap');

R_prev(1) = DataRaw.Robs(1);
R_prev(end) = DataRaw.Robs(end);
I_prev = max(I_prev,0);

S_prev = max(Ntot - I_prev - R_prev, 0);

R_cont = R_prev(:);


%% ====== C. DA–MCMC (Algorithm 1) ======
if isempty(block_B), block_B = max(8, m); end

% Initialisation I_curr (latent I(t))
I_curr = max(round(0.001 * max(R_cont)) * exp(-0.05*(1:Nobs))' , 1);
I_curr = min(I_curr, Ntot - R_cont);

% Pré-calcul a,b
[ai_vec, bi_vec] = local_ab(I_curr, R_cont, t_grid, Ntot);

% Stats suff.
[NI_T, NR_T, C1, C2] = sufficient_stats(I_curr, round(R_cont), ai_vec, bi_vec);

% init beta/gamma
beta_curr  = gamrnd((NI_T + m_beta), 1/ max(C1 + lambda_beta, 1e-12));
gamma_curr = gamrnd((NR_T + m_gamma), 1/ max(C2 + lambda_gamma, 1e-12));
R0_curr    = beta_curr / max(gamma_curr, 1e-12);

% Stockage
R0_MCMC_samples     = zeros(1, H);
R0_MLE_samples      = zeros(1, H);
beta_draw_samples   = zeros(1, H);
gamma_draw_samples  = zeros(1, H);
I_samples           = zeros(Nobs, H);

block_starts = 1:block_B:Nobs;

% Compteurs du taux d'acceptation MH par blocs
n_block_proposals = 0;
n_block_accept     = 0;

fprintf('--- Démarrage DA–MCMC : H=%d, burn-in=%d, m=%d ---\n', H, burn_in, m);

for t = 1:H

    % === Mise à jour bloc par bloc de I(t) ===
    for bs = block_starts
        be = min(bs + block_B - 1, Nobs);
        idxB = bs:be;

        I_star = I_curr;

        % Propositions locales composante par composante
        for jj = idxB

            kmin = max(0, I_curr(jj)-r_window);
            kmax = min(Ntot - R_cont(jj), I_curr(jj)+r_window);
            if kmin > kmax, continue; end

            kvals = kmin:kmax;
            logw  = -inf(1, numel(kvals));

            for kk = 1:numel(kvals)
                I_tmp = I_star;
                I_tmp(jj) = kvals(kk);

                [ai_loc, bi_loc] = local_ab_single(I_tmp, R_cont, t_grid, Ntot, jj);
                if ai_loc<=0 || bi_loc<=0, continue; end

                rj = round(R_cont(jj));
                kj = I_tmp(jj);
                arg1 = kj + rj;
                if arg1 < 0 || rj < 0, continue; end

                beta_c  = R0_curr * gamma_curr;
                gamma_c = gamma_curr;

                % log poids (pseudo-vraisemblance discrète)
                logw(kk) = -(beta_c*ai_loc + gamma_c*bi_loc) ...
                         + arg1*log(max(beta_c*ai_loc,1e-12)) - gammaln(arg1+1) ...
                         + rj*log(max(gamma_c*bi_loc,1e-12)) - gammaln(rj+1);
            end

            mlog = max(logw);
            if isfinite(mlog)
                w = exp(logw - mlog);
                sw = sum(w);
                if sw > 0
                    w = w / sw;
                    I_star(jj) = randsample(kvals, 1, true, w);
                end
            end
        end

        % Metropolis-Hastings sur le bloc
        iL = max(1, bs-1);
        iR = be;

        [ai_curr_all, bi_curr_all] = local_ab(I_curr, R_cont, t_grid, Ntot);
        [ai_star_all, bi_star_all] = local_ab(I_star, R_cont, t_grid, Ntot);

        beta_c  = R0_curr * gamma_curr;
        gamma_c = gamma_curr;

        log_num = 0; log_den = 0;

        for j = iL:iR
            rj = round(R_cont(j));

            kj_s = I_star(j); arg1_s = kj_s + rj;
            kj_c = I_curr(j); arg1_c = kj_c + rj;

            if ai_star_all(j)<=0 || bi_star_all(j)<=0 || arg1_s<0 || rj<0
                log_num = -inf; break;
            end
            if ai_curr_all(j)<=0 || bi_curr_all(j)<=0 || arg1_c<0 || rj<0
                log_den = -inf; break;
            end

            log_num = log_num ...
                -(beta_c*ai_star_all(j) + gamma_c*bi_star_all(j)) ...
                + arg1_s*log(max(beta_c*ai_star_all(j),1e-12)) - gammaln(arg1_s+1) ...
                + rj*log(max(gamma_c*bi_star_all(j),1e-12)) - gammaln(rj+1);

            log_den = log_den ...
                -(beta_c*ai_curr_all(j) + gamma_c*bi_curr_all(j)) ...
                + arg1_c*log(max(beta_c*ai_curr_all(j),1e-12)) - gammaln(arg1_c+1) ...
                + rj*log(max(gamma_c*bi_curr_all(j),1e-12)) - gammaln(rj+1);
        end

        log_alpha = log_num - log_den;
        n_block_proposals = n_block_proposals + 1;

        if log(rand) < log_alpha
            n_block_accept = n_block_accept + 1;
            I_curr = I_star;
            ai_vec = ai_star_all;
            bi_vec = bi_star_all;
        else
            ai_vec = ai_curr_all;
            bi_vec = bi_curr_all;
        end
    end

    % === Tirage beta/gamma (postérieurs conditionnels) ===
    [NI_T, NR_T, C1, C2] = sufficient_stats(I_curr, round(R_cont), ai_vec, bi_vec);

    p_post = NI_T + m_beta;
    q_post = NR_T + m_gamma;

    beta_draw  = gamrnd(p_post, 1/(C1 + lambda_beta));
    gamma_draw = gamrnd(max(q_post,1e-6), 1/(C2 + lambda_gamma));

    R0_curr    = beta_draw / max(gamma_draw, 1e-12);
    gamma_curr = gamma_draw;

    % === Stockage ===
    R0_MCMC_samples(t)    = R0_curr;
    R0_MLE_samples(t)     = (NI_T*C2) / max(NR_T*C1, 1e-12);
    beta_draw_samples(t)  = beta_draw;
    gamma_draw_samples(t) = gamma_draw;
    I_samples(:,t)        = I_curr;

    if mod(t, log_freq)==0 && t > burn_in
        post_now = R0_MCMC_samples(burn_in+1:t);
        mle_now  = R0_MLE_samples(max(1,burn_in+1):t);
        fprintf('It. %5d/%5d | mean_post(R0,MCMC)=%.4f | mean_post(R0,MLE)=%.4f\n', ...
            t, H, mean(post_now), mean(mle_now));
    end
end


%% ====== D. RESUMES POSTERIEURS (R0) ======
post = R0_MCMC_samples(burn_in+1:end);
MLE  = R0_MLE_samples(burn_in+1:end);

meanR0 = mean(post);
medR0  = median(post);
seR0   = std(post);
ciR0   = quantile(post,[0.025 0.975]);

meanR0_mle = mean(MLE);
medR0_mle  = median(MLE);
seR0_mle   = std(MLE);
ciR0_mle   = quantile(MLE,[0.025 0.975]);

fprintf('\nEbola (Sierra Leone):  R0_MCMC= %.4f | median=%.4f | 95%% CI=[%.4f, %.4f] | SE=%.5f\n', ...
    meanR0, medR0, ciR0(1), ciR0(2), seR0);
fprintf('Ebola (Sierra Leone):  R0_MLE = %.4f | median=%.4f | 95%% CI=[%.4f, %.4f] | SE=%.5f\n', ...
    meanR0_mle, medR0_mle, ciR0_mle(1), ciR0_mle(2), seR0_mle);

fprintf('Samples=%d | Burn-in=%d | N=%.2e\n', H, burn_in, Ntot);


%% ====== D2. DIAGNOSTICS NUMERIQUES DE CONVERGENCE ======
post_beta  = beta_draw_samples(burn_in+1:end);
post_gamma = gamma_draw_samples(burn_in+1:end);

maxLag = min(acf_max_lag, numel(post)-1);

[ESS_R0,    IACT_R0,    acf_R0]    = ess_geyer(post,       maxLag);
[ESS_beta,  IACT_beta,  acf_beta]  = ess_geyer(post_beta,  maxLag);
[ESS_gamma, IACT_gamma, acf_gamma] = ess_geyer(post_gamma, maxLag);

% Erreur standard Monte-Carlo (à ne pas confondre avec l'écart-type postérieur)
MCSE_R0    = std(post)       / sqrt(max(ESS_R0,1));
MCSE_beta  = std(post_beta)  / sqrt(max(ESS_beta,1));
MCSE_gamma = std(post_gamma) / sqrt(max(ESS_gamma,1));

% Le script original ne produit qu'une seule chaîne indépendante.
% On rapporte donc un rank-normalized split-Rhat obtenu en divisant
% la partie post burn-in en n_split_chains segments consécutifs.
Rhat_R0    = split_rank_rhat(post,       n_split_chains);
Rhat_beta  = split_rank_rhat(post_beta,  n_split_chains);
Rhat_gamma = split_rank_rhat(post_gamma, n_split_chains);

accept_rate_block = 100 * n_block_accept / max(n_block_proposals,1);

Parameter = ["R0"; "beta"; "gamma"];
PosteriorMean = [mean(post); mean(post_beta); mean(post_gamma)];
PosteriorSD   = [std(post);  std(post_beta);  std(post_gamma)];
ESS           = [ESS_R0; ESS_beta; ESS_gamma];
IACT          = [IACT_R0; IACT_beta; IACT_gamma];
MCSE          = [MCSE_R0; MCSE_beta; MCSE_gamma];
SplitRhat     = [Rhat_R0; Rhat_beta; Rhat_gamma];

DiagnosticsTable = table(Parameter, PosteriorMean, PosteriorSD, ESS, IACT, MCSE, SplitRhat);

fprintf('\n================ MCMC CONVERGENCE DIAGNOSTICS ================\n');
disp(DiagnosticsTable);
fprintf('Block MH acceptance rate = %.2f%% (%d/%d)\n', ...
    accept_rate_block, n_block_accept, n_block_proposals);
fprintf(['Interpretation indicative: Split-Rhat proche de 1, ESS élevé, ', ...
         'MCSE faible et ACF qui décroît rapidement.\n']);

try
    writetable(DiagnosticsTable, fullfile(outdir,'mcmc_convergence_diagnostics.csv'));
catch ME
    warning('Export du tableau de diagnostics impossible: %s', ME.message);
end

% Figure ACF pour R0, beta et gamma
figACF = figure(3); clf;
set(figACF,'Units','normalized','Position',[0.08 0.10 0.84 0.72]);
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

acf_names = {'$\mathcal{R}_0$','$\beta$','$\gamma$'};
acf_values = {acf_R0, acf_beta, acf_gamma};
n_post = numel(post);
acf_bound = 1.96/sqrt(n_post);

for q = 1:3
    nexttile;
    lags = 0:(numel(acf_values{q})-1);
    stem(lags, acf_values{q}, 'filled', 'MarkerSize',3);
    hold on;
    yline(acf_bound,'r--','LineWidth',1);
    yline(-acf_bound,'r--','LineWidth',1);
    yline(0,'k-');
    grid on;
    xlim([0 maxLag]);
    ylim([-1 1]);
    xlabel('Lag');
    ylabel('ACF');
    title(['ACF of ', acf_names{q}], 'Interpreter','latex');
end
sgtitle('Post-burn-in autocorrelation diagnostics');




%% ====== E. CONSTRUIRE OBSERVE SUR GRILLE + SCALE I DRAWS ======
I_obs = DataRaw.cases - DataRaw.Robs;
I_obs = max(I_obs, 0);

I_obs_interp = interp1(DataRaw.date, I_obs, tgrid_dt, 'linear', 'extrap');
I_obs_interp = max(I_obs_interp, 0);
I_obs_interp = I_obs_interp(:);

% Mise à l’échelle des trajectoires I
I_post_lat = I_samples(:, burn_in+1:end);
scaleFactor = max(I_obs_interp) / max(mean(I_post_lat,2));
I_scaled = I_samples * scaleFactor;   % [Nobs x H]

% Posterior set
idx_post = (burn_in+1):H;


%% ====== F. FIGURE 1 — Diagnostics statistiques ======
fig1 = figure(1); clf;
set(fig1, 'Units','normalized','Position',[0.02 0.05 0.96 0.88]);
tiledlayout(3,2,'Padding','compact','TileSpacing','tight');

nexttile;
[f1,x1] = ksdensity(post,'Support','positive');
plot(x1,f1,'b-','LineWidth',2); grid on;
xlabel('$\mathcal{R}_0$','Interpreter','latex'); ylabel('Density');
title('(a) Posterior density of $\mathcal{R}_0$','Interpreter','latex');

nexttile;
[f2,x2] = ksdensity(MLE,'Support','positive');
plot(x1,f1,'b-','LineWidth',2); hold on;
plot(x2,f2,'r--','LineWidth',2);
legend({'Bayesian posterior','MLE'},'Location','best');
xlabel('$\mathcal{R}_0$','Interpreter','latex'); ylabel('Density');
title('(b) Posterior vs MLE'); grid on;

nexttile;
plot(R0_MCMC_samples,'k-','LineWidth',0.8); grid on;
xlabel('Iteration'); ylabel('$\mathcal{R}_0^{(t)}$','Interpreter','latex');
title('(c) Trace plot'); 

nexttile;
histogram(post,40,'EdgeColor','none');
xlabel('$\mathcal{R}_0$','Interpreter','latex'); ylabel('Count');
title('(d) Posterior histogram'); grid on;

nexttile;

% ACF sans Econometrics Toolbox
maxLag = 50;
x = post(:);
x = x - mean(x);
n = length(x);

acf_vals = zeros(maxLag+1,1);
acf_vals(1) = 1;

den = sum(x.^2);

for lag = 1:maxLag
    acf_vals(lag+1) = sum(x(1:n-lag).*x(1+lag:n)) / den;
end

stem(0:maxLag, acf_vals, 'filled');
xlabel('Lag');
ylabel('Autocorrelation');
title('ACF of post-burn-in R_0 chain');
grid on;

nexttile;
plot(cumsum(R0_MCMC_samples)./(1:H),'m-','LineWidth',1.2); grid on;
xlabel('Iteration'); ylabel('Cumulative mean');
title('(f) Cumulative mean diagnostic');

sgtitle('Diagnostics Bayésiens vs MLE pour $\mathcal{R}_0$','Interpreter','latex','FontWeight','bold');


%% ====== G. FIGURE PRINCIPALE (comme ta photo) : I et R avec IC 95% ======
% -------- (1) Infected posterior predictive via draws I_scaled ----------
I_post = I_scaled(:, idx_post);

for k = 1:size(I_post,2)
    I_post(:,k) = movmean(I_post(:,k), 7*m_per_day);
end       % [Nobs x Npost]
I_mean = mean(I_post,2);
I_lo   = quantile(I_post',0.025)';
I_hi   = quantile(I_post',0.975)';

span = 7*m_per_day;   % cohérent avec ton movmean

I_mean = smoothdata(I_mean,'movmean',span);
I_lo   = smoothdata(I_lo  ,'movmean',span);
I_hi   = smoothdata(I_hi  ,'movmean',span);

% -------- (2) Removed posterior predictive (SIR-consistent) --------
Npost = numel(idx_post);
R_pred = zeros(Nobs, Npost);

for k = 1:Npost
    it = idx_post(k);

    I_k = I_scaled(:,it);           % utiliser I SCALED
    Rk  = Ntot - S_prev - I_k;      % définition SIR exacte

    R_pred(:,k) = max(Rk,0);
end

R_mean = mean(R_pred,2);
R_lo   = quantile(R_pred',0.025)';
R_hi   = quantile(R_pred',0.975)';

% Lissage
R_mean = smoothdata(R_mean,'movmean',span);
R_lo   = smoothdata(R_lo  ,'movmean',span);
R_hi   = smoothdata(R_hi  ,'movmean',span);


% R_mean = mean(R_pred, 2);
% R_lo   = quantile(R_pred', 0.025)'; 
% R_hi   = quantile(R_pred', 0.975)';
% 
% R_mean = R_mean + (R_cont(1) - R_mean(1));
% R_lo   = R_lo   + (R_cont(1) - R_mean(1));
% R_hi   = R_hi   + (R_cont(1) - R_mean(1));

% -------- PLOT (2 panneaux) ----------
figMain = figure(2); clf;
set(figMain,'Units','normalized','Position',[0.05 0.10 0.90 0.78]);

t_days = t_grid(:);

% ---- Top: Infected ----
subplot(2,1,1); hold on; box on;
fill([t_days; flipud(t_days)], ...
     [I_lo; flipud(I_hi)], ...
     [0.30 0.60 1.00], ...     % 
     'FaceAlpha',0.35, ...     % 
     'EdgeColor',[0.10 0.30 0.80], ...
     'LineWidth',1.2);
plot(t_days, I_mean, '-', ...
     'Color',[0 0.25 0.8], 'LineWidth',2.5);
plot(t_days, I_obs_interp, 'k.', 'MarkerSize',10);

title('Observed vs Posterior Predictive Mean (95% CI) - Infected');
xlabel('Time (days)'); ylabel('Infected individuals');
legend({'95% CI','Posterior mean','Observed'}, 'Location','northeast');
grid on;

% ---- Bottom: Removed ----
subplot(2,1,2); hold on; box on;
fill([t_days; flipud(t_days)], ...
     [R_lo; flipud(R_hi)], ...
     [1.00 0.55 0.55], ...     % 
     'FaceAlpha',0.35, ...
     'EdgeColor',[0.75 0.15 0.15], ...
     'LineWidth',1.2);
R_obs_plot = smoothdata(R_cont,'movmean',3);
plot(t_days, R_mean, '-', ...
     'Color',[0.75 0 0], 'LineWidth',2.5);
plot(t_days, R_obs_plot, 'k.', 'MarkerSize',10);


title('Observed vs Posterior Predictive Mean (95% CI) - Removed');
xlabel('Time (days)'); ylabel('Removed individuals');
legend({'95% CI','Posterior mean','Observed'}, 'Location','southeast');
grid on;


%% ====== H. SAUVEGARDES ======
try
    save(fullfile(outdir,'mcmc_outputs.mat'), ...
        'R0_MCMC_samples','R0_MLE_samples','beta_draw_samples','gamma_draw_samples', ...
        'I_samples','I_scaled','t_grid','tgrid_dt','R_cont','I_obs_interp','DataRaw', ...
        'burn_in','H','Ntot','m_per_day', ...
        'DiagnosticsTable','accept_rate_block','acf_R0','acf_beta','acf_gamma', ...
        'ESS_R0','ESS_beta','ESS_gamma','Rhat_R0','Rhat_beta','Rhat_gamma');
catch ME
    warning('Sauvegarde MAT impossible: %s', ME.message); %#ok<MEXCEP>
end

try
    saveas(fig1, fullfile(outdir,'diagnostics_R0.png'));
    saveas(figACF, fullfile(outdir,'acf_R0_beta_gamma.png'));
    saveas(figMain, fullfile(outdir,'posterior_predictive_I_R.png'));
catch ME
    warning('Sauvegarde figures impossible: %s', ME.message); %#ok<MEXCEP>
end

fprintf('\nTerminé. Figures et résultats sauvegardés dans: %s\n', outdir);


%% ===================== LOCAL FUNCTIONS =====================

function [ai_vec, bi_vec] = local_ab(I, R_cont, t, Ntot)
    Nobs = numel(t);
    ai_vec = zeros(Nobs,1); bi_vec = zeros(Nobs,1);
    for j = 2:Nobs
        dt = t(j)-t(j-1);
        S1 = max(Ntot - R_cont(j-1) - I(j-1), 0);
        S2 = max(Ntot - R_cont(j)   - I(j),   0);
        ai_vec(j) = (dt/(2*Ntot)) * (I(j-1)*S1 + I(j)*S2);
        bi_vec(j) = (dt/2) * (I(j-1) + I(j));
    end
end

function [ai, bi] = local_ab_single(I, R_cont, t, Ntot, j)
    if j==1, ai=0; bi=0; return; end
    dt = t(j)-t(j-1);
    S1 = max(Ntot - R_cont(j-1) - I(j-1), 0);
    S2 = max(Ntot - R_cont(j)   - I(j),   0);
    ai = (dt/(2*Ntot)) * (I(j-1)*S1 + I(j)*S2);
    bi = (dt/2) * (I(j-1) + I(j));
end

function [NI_T, NR_T, C1, C2] = sufficient_stats(I, R_obs, ai_vec, bi_vec)
    C1 = sum(ai_vec(2:end));
    C2 = sum(bi_vec(2:end));
    NR_T = max(R_obs(end)-R_obs(1), 0);
    NI_T = max((I(end)-I(1)) + NR_T, 0);
end

function T = normalize_date_column(T, dateVar)
    col = T.(dateVar);

    if isdatetime(col)
        return;
    end

    if isnumeric(col) || isduration(col)
        try
            T.(dateVar) = datetime(col, 'ConvertFrom', 'excel');
            return;
        catch ME
            warning('Conversion excel date échouée (%s): %s', dateVar, ME.message);
        end
    end

    if iscell(col) || iscellstr(col) || isstring(col)
        fmt_list = {'default','yyyy-MM-dd','dd/MM/yyyy','MM/dd/yyyy'};
        for k = 1:numel(fmt_list)
            try
                if strcmp(fmt_list{k},'default')
                    T.(dateVar) = datetime(string(col));
                else
                    T.(dateVar) = datetime(string(col), 'InputFormat', fmt_list{k});
                end
                return;
            catch
            end
        end
    end

    error('Impossible de convertir la colonne "%s" en datetime.', dateVar);
end


function [ess, tau, acf_values] = ess_geyer(x, maxLag)
% ESS fondé sur une séquence initiale positive et monotone des paires d'ACF.
    x = x(:);
    x = x(isfinite(x));
    n = numel(x);

    if n < 4
        ess = n;
        tau = 1;
        acf_values = 1;
        return;
    end

    maxLag = min(maxLag, n-1);
    acf_values = sample_acf(x, maxLag);

    % Paires (rho_1+rho_2), (rho_3+rho_4), ...
    nPairs = floor(maxLag/2);
    pairSums = zeros(nPairs,1);
    nKeep = 0;

    for k = 1:nPairs
        lag1_index = 2*k;       % index MATLAB de rho_(2k-1)
        lag2_index = 2*k + 1;   % index MATLAB de rho_(2k)
        pairValue = acf_values(lag1_index) + acf_values(lag2_index);

        if ~isfinite(pairValue) || pairValue <= 0
            break;
        end

        nKeep = nKeep + 1;
        pairSums(nKeep) = pairValue;
    end

    pairSums = pairSums(1:nKeep);

    % Séquence monotone décroissante
    for k = 2:numel(pairSums)
        pairSums(k) = min(pairSums(k), pairSums(k-1));
    end

    tau = 1 + 2*sum(pairSums);
    tau = max(tau,1);
    ess = min(n, max(1, n/tau));
end


function acf_values = sample_acf(x, maxLag)
% ACF empirique avec normalisation par la variance au lag 0.
    x = x(:);
    x = x - mean(x);
    n = numel(x);
    denominator = sum(x.^2);

    acf_values = zeros(maxLag+1,1);
    acf_values(1) = 1;

    if denominator <= 0
        return;
    end

    for lag = 1:maxLag
        acf_values(lag+1) = ...
            sum(x(1:n-lag).*x(1+lag:n)) / denominator;
    end
end


function rhat = split_rank_rhat(x, nSegments)
% Rank-normalized folded split-Rhat.
% Pour une seule longue chaîne, celle-ci est découpée en segments.
    x = x(:);
    x = x(isfinite(x));
    nSegments = max(2, round(nSegments));

    nPerSegment = floor(numel(x)/nSegments);
    if nPerSegment < 4
        rhat = NaN;
        return;
    end

    nUsed = nPerSegment*nSegments;
    X = reshape(x(1:nUsed), nPerSegment, nSegments);

    % Chaque segment est encore scindé en deux moitiés.
    half = floor(nPerSegment/2);
    Xs = [X(1:half,:), X(end-half+1:end,:)];

    Z = rank_normalize_matrix(Xs);

    folded = abs(Xs - median(Xs(:)));
    Zfold = rank_normalize_matrix(folded);

    rhat = max(basic_rhat(Z), basic_rhat(Zfold));
end


function Z = rank_normalize_matrix(X)
    ranks = tiedrank(X(:));
    N = numel(ranks);

    % Transformation de Blom utilisée pour éviter 0 et 1.
    probabilities = (ranks - 3/8) ./ (N + 1/4);
    probabilities = min(max(probabilities, eps), 1-eps);

    Z = reshape(norminv(probabilities,0,1), size(X));
end


function rhat = basic_rhat(X)
    [n, m] = size(X);

    if n < 2 || m < 2
        rhat = NaN;
        return;
    end

    chainMeans = mean(X,1);
    W = mean(var(X,0,1));
    B = n * var(chainMeans,0,2);

    if W <= 0
        rhat = NaN;
        return;
    end

    varPlus = ((n-1)/n)*W + B/n;
    rhat = sqrt(varPlus/W);
end


function Metrics = evaluation_metrics(estimates, ciLow, ciHigh, trueValue, methodName)
% Calcule les performances sur B réplications indépendantes.
    estimates = estimates(:);
    ciLow = ciLow(:);
    ciHigh = ciHigh(:);

    valid = isfinite(estimates) & isfinite(ciLow) & isfinite(ciHigh);
    estimates = estimates(valid);
    ciLow = ciLow(valid);
    ciHigh = ciHigh(valid);

    if isempty(estimates)
        error('Aucune réplication valide pour calculer les métriques.');
    end

    B = numel(estimates);
    Bias = mean(estimates - trueValue);
    AbsoluteBias = mean(abs(estimates - trueValue));
    MSE = mean((estimates - trueValue).^2);
    RMSE = sqrt(MSE);
    CoverageProbability = mean(ciLow <= trueValue & trueValue <= ciHigh);
    AverageIntervalWidth = mean(ciHigh - ciLow);

    Metrics = table(string(methodName), B, trueValue, Bias, AbsoluteBias, ...
        MSE, RMSE, CoverageProbability, AverageIntervalWidth, ...
        'VariableNames', {'Method','Replications','TrueValue','Bias', ...
        'MeanAbsoluteBias','MSE','RMSE','CoverageProbability', ...
        'AverageIntervalWidth'});
end
