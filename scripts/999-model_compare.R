# ============================================================
# Simulation test: TMB (Laplace+MAP) vs Stan (posterior)
# Goal: check whether TMB estimates are in the same ballpark
#       as Stan and whether both recover true parameters.
#
# Notes:
# - Simulation generates B_true from the phylogenetic GP prior:
#     B[,k] ~ MVN(0, Kphy(alpha,rho))  independently for k=1..K
# - We fit with:
#   (A) TMB: integrate out f as random effect (random="f"), keep B as parameters
#       => fast Laplace; returns MAP-like estimates
#   (B) Stan: full posterior sampling; compare posterior means/medians
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
    library(MASS) # mvrnorm
    library(TMB) # Laplace+MAP (proxy for RTMB-style inference)
    library(rstan)
})

rstan_options(auto_write = TRUE)
options(mc.cores = max(1, parallel::detectCores() - 1))

# -----------------------------
# Utility: phylo GP covariance
# -----------------------------
make_Kphy <- function(D_phylo, alpha, rho, jitter = 1e-6) {
    Kphy <- (alpha^2) * exp(-(D_phylo^2) / (2 * rho^2))
    diag(Kphy) <- diag(Kphy) + jitter
    Kphy
}

# -----------------------------
# Simulate from the model
# -----------------------------
simulate_dataset <- function(N = 5000, J = 30, K = 10,
                             alpha = 1.0, rho = 8.0, sigma_f = 0.5,
                             offset_const = 0.0,
                             seed = 1) {
    set.seed(seed)

    species <- sample.int(J, N, replace = TRUE)
    X <- matrix(rnorm(N * K), N, K)
    D_phylo <- as.matrix(dist(1:J))

    Kphy <- make_Kphy(D_phylo, alpha = alpha, rho = rho, jitter = 1e-6)

    B_true <- matrix(NA_real_, J, K)
    for (k in 1:K) {
        B_true[, k] <- as.numeric(MASS::mvrnorm(n = 1, mu = rep(0, J), Sigma = Kphy))
    }

    f_true <- rnorm(N, 0, sigma_f)
    offset <- rep(offset_const, N)

    eta <- rowSums(X * B_true[species, , drop = FALSE]) + f_true + offset
    lambda <- exp(eta)
    y <- rpois(N, lambda)

    list(
        data = list(
            N = N, J = J, K = K, species = species, X = X,
            y = as.numeric(y), D_phylo = D_phylo, offset = offset
        ),
        truth = list(
            alpha = alpha, rho = rho, sigma_f = sigma_f,
            B_true = B_true, f_true = f_true
        )
    )
}

# ============================================================
# TMB fit (Laplace over f, MAP-ish for params including B)
# Assumes OPT2-style params in C++:
#   PARAMETER(alpha), PARAMETER(rho), PARAMETER(sigma_f)
# ============================================================
fit_tmb <- function(sim, cpp_file, model_base = NULL,
                    control_iter = 500, verbose = TRUE) {
    stopifnot(file.exists(cpp_file))
    if (is.null(model_base)) model_base <- sub("\\.cpp$", "", basename(cpp_file))

    so_path <- paste0(dirname(cpp_file), "/", model_base, .Platform$dynlib.ext)
    if (!file.exists(so_path) || file.info(cpp_file)$mtime > file.info(so_path)$mtime) {
        if (verbose) message("Compiling TMB model: ", cpp_file)
        old <- getwd()
        setwd(dirname(cpp_file))
        TMB::compile(basename(cpp_file))
        setwd(old)
    }
    dyn.load(TMB::dynlib(file.path(dirname(cpp_file), model_base)))

    dat <- sim$data

    par <- list(
        B = matrix(0, dat$J, dat$K),
        alpha = 0.8,
        rho = 6.0,
        sigma_f = 0.4,
        f = rep(0, dat$N)
    )

    obj <- TMB::MakeADFun(
        data = dat,
        parameters = par,
        DLL = model_base,
        random = "f",
        silent = TRUE
    )

    t1 <- Sys.time()
    opt <- nlminb(obj$par, obj$fn, obj$gr,
        control = list(iter.max = control_iter, eval.max = control_iter * 2)
    )
    t2 <- Sys.time()

    rep <- TMB::sdreport(obj)

    par_hat <- obj$env$parList(opt$par)

    list(
        opt = opt,
        sdreport = rep,
        estimates = list(
            B = par_hat$B,
            alpha = par_hat$alpha,
            rho = par_hat$rho,
            sigma_f = par_hat$sigma_f
        ),
        runtime_sec = as.numeric(difftime(t2, t1, units = "secs"))
    )
}

# ============================================================
# Stan fit (posterior)
# ============================================================
stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> K;
  int<lower=1> species[N];
  matrix[N, K] X;
  int<lower=0> y[N];
  matrix[J, J] D_phylo;
  vector[N] offset;
}
parameters {
  matrix[J, K] B;
  real<lower=0> alpha;
  real<lower=0> rho;
  real<lower=0> sigma_f;
  vector[N] f;
}
transformed parameters {
  matrix[J, J] K_phylo;
  matrix[J, J] L_B;
  for (i in 1:J) {
    for (j in 1:J) {
      K_phylo[i, j] = alpha^2 * exp(-D_phylo[i, j]^2 / (2 * rho^2));
    }
  }
  K_phylo = K_phylo + diag_matrix(rep_vector(1e-6, J));
  L_B = cholesky_decompose(K_phylo);
}
model {
  alpha ~ normal(0, 1);
  rho ~ normal(0, 1);
  sigma_f ~ normal(0, 1);

  for (k in 1:K) {
    B[, k] ~ multi_normal_cholesky(rep_vector(0, J), L_B);
  }

  f ~ normal(0, sigma_f);

  for (n in 1:N) {
    int s = species[n];
    real log_lambda = dot_product(X[n], B[s, ]) + f[n] + offset[n];
    y[n] ~ poisson_log(log_lambda);
  }
}
"

fit_stan <- function(sim,
                     chains = 1,
                     iter = 1500,
                     warmup = 750,
                     adapt_delta = 0.9,
                     max_treedepth = 10,
                     seed = 1,
                     verbose = TRUE) {
    dat <- sim$data
    standata <- list(
        N = dat$N,
        J = dat$J,
        K = dat$K,
        species = dat$species,
        X = dat$X,
        y = as.integer(dat$y),
        D_phylo = dat$D_phylo,
        offset = dat$offset
    )

    if (verbose) message("Compiling Stan model (once per session)...")
    sm <- rstan::stan_model(model_code = stan_code)

    t1 <- Sys.time()
    fit <- rstan::sampling(
        sm,
        data = standata,
        chains = chains,
        iter = iter,
        warmup = warmup,
        seed = seed,
        control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
        refresh = 0
    )
    t2 <- Sys.time()

    summ <- rstan::summary(fit, pars = c("alpha", "rho", "sigma_f"))$summary
    alpha_hat <- summ["alpha", "mean"]
    rho_hat <- summ["rho", "mean"]
    sigma_f_hat <- summ["sigma_f", "mean"]

    B_draws <- rstan::extract(fit, pars = "B")$B # iters x J x K
    B_mean <- apply(B_draws, c(2, 3), mean)
    B_sd <- apply(B_draws, c(2, 3), sd)

    list(
        fit = fit,
        estimates = list(B_mean = B_mean, B_sd = B_sd, alpha = alpha_hat, rho = rho_hat, sigma_f = sigma_f_hat),
        runtime_sec = as.numeric(difftime(t2, t1, units = "secs"))
    )
}

# ============================================================
# Comparison metrics (existing)
# ============================================================
compare_fits <- function(sim, tmb_fit, stan_fit) {
    B_true <- sim$truth$B_true
    B_tmb <- tmb_fit$estimates$B
    B_stan <- stan_fit$estimates$B_mean

    rmse <- function(a, b) sqrt(mean((a - b)^2))
    cor_vec <- function(a, b) cor(as.numeric(a), as.numeric(b))

    data.frame(
        true_alpha = sim$truth$alpha,
        tmb_alpha = tmb_fit$estimates$alpha,
        stan_alpha = stan_fit$estimates$alpha,
        true_rho = sim$truth$rho,
        tmb_rho = tmb_fit$estimates$rho,
        stan_rho = stan_fit$estimates$rho,
        true_sigma_f = sim$truth$sigma_f,
        tmb_sigma_f = tmb_fit$estimates$sigma_f,
        stan_sigma_f = stan_fit$estimates$sigma_f,
        B_rmse_tmb = rmse(B_true, B_tmb),
        B_rmse_stan = rmse(B_true, B_stan),
        B_cor_tmb = cor_vec(B_true, B_tmb),
        B_cor_stan = cor_vec(B_true, B_stan),
        tmb_time_sec = tmb_fit$runtime_sec,
        stan_time_sec = stan_fit$runtime_sec
    )
}

# ============================================================
# Run a small simulation study
# ============================================================
SIMS <- 3
N <- 2000
J <- 30
K <- 2

TRUE_ALPHA <- 1.0
TRUE_RHO <- 8.0
TRUE_SIGMA_F <- 0.5

cpp_file <- "~/clim_risk_phylosdm/scripts/05_lgcp_corrected.cpp" # <-- EDIT THIS PATH

# Store everything we need for B mean/SD comparison
sim_store <- vector("list", SIMS)
tmb_store <- vector("list", SIMS)
stan_store <- vector("list", SIMS)
results <- vector("list", SIMS)

for (s in 1:SIMS) {
    cat("\n==============================\n")
    cat("SIM", s, "\n")
    cat("==============================\n")

    sim <- simulate_dataset(
        N = N, J = J, K = K,
        alpha = TRUE_ALPHA, rho = TRUE_RHO, sigma_f = TRUE_SIGMA_F,
        offset_const = 0.0,
        seed = 100 + s
    )

    cat("Fitting TMB (Laplace+MAP)...\n")
    tmb_fit <- fit_tmb(sim, cpp_file = cpp_file, control_iter = 400, verbose = TRUE)

    cat("Fitting Stan (posterior)...\n")
    stan_fit <- fit_stan(sim,
        chains = 1,
        iter = 1200,
        warmup = 600,
        adapt_delta = 0.9,
        max_treedepth = 10,
        seed = 200 + s,
        verbose = TRUE
    )

    res <- compare_fits(sim, tmb_fit, stan_fit)
    print(res)

    sim_store[[s]] <- sim
    tmb_store[[s]] <- tmb_fit
    stan_store[[s]] <- stan_fit
    results[[s]] <- res
}

results_df <- do.call(rbind, results)

cat("\n==============================\n")
cat("SUMMARY ACROSS SIMS\n")
cat("==============================\n")
print(results_df)
cat("\nMeans:\n")
print(colMeans(results_df[sapply(results_df, is.numeric)]))

# ============================================================
# NEW: Compare Beta (B) estimates: means and SDs
# ============================================================

# Helper: summarize a JxK matrix across species (rows), per covariate (col)
summ_by_k <- function(M) {
    data.frame(
        k = seq_len(ncol(M)),
        mean_over_species = as.numeric(colMeans(M)),
        sd_over_species = as.numeric(apply(M, 2, sd))
    )
}

beta_comp <- vector("list", SIMS)

for (s in 1:SIMS) {
    B_true <- sim_store[[s]]$truth$B_true
    B_tmb <- tmb_store[[s]]$estimates$B

    B_stan_mean <- stan_store[[s]]$estimates$B_mean
    B_stan_sd <- stan_store[[s]]$estimates$B_sd

    # Summaries across species for:
    # - truth B_true
    # - TMB point estimates B_hat
    # - Stan posterior mean matrix
    # - Stan posterior SD matrix (uncertainty of each B[j,k])
    truth_s <- summ_by_k(B_true)
    tmb_s <- summ_by_k(B_tmb)
    stan_mean_s <- summ_by_k(B_stan_mean)
    stan_sd_s <- summ_by_k(B_stan_sd)

    beta_comp[[s]] <- data.frame(
        sim = s,
        k = truth_s$k,
        true_mean = truth_s$mean_over_species,
        true_sd = truth_s$sd_over_species,
        tmb_mean = tmb_s$mean_over_species,
        tmb_sd = tmb_s$sd_over_species,
        stan_postmean_mean = stan_mean_s$mean_over_species,
        stan_postmean_sd = stan_mean_s$sd_over_species,
        stan_postsd_mean = stan_sd_s$mean_over_species,
        stan_postsd_sd = stan_sd_s$sd_over_species
    )
}

beta_comp_df <- do.call(rbind, beta_comp)

cat("\n==============================\n")
cat("BETA (B) COMPARISON: per-sim, per-k summaries\n")
cat("==============================\n")
print(beta_comp_df)

cat("\n==============================\n")
cat("BETA (B) COMPARISON: averaged across sims (per k)\n")
cat("==============================\n")
beta_comp_mean <- aggregate(. ~ k, data = beta_comp_df[, setdiff(names(beta_comp_df), "sim")], FUN = mean)
print(beta_comp_mean)

# Optional saves
save(results_df, beta_comp_df, beta_comp_mean, file = "sim_tmb_vs_stan_results_with_beta_comparison.RData")
