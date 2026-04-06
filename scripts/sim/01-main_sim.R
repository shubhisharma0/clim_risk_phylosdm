#!/usr/bin/env Rscript
# ============================================================
# Evolutionary-data simulation + RTMB/TMB fit + conditional prediction
#
# What this does:
#  1) Simulates an ultrametric-ish phylogeny (or uses a provided tree)
#  2) Simulates traits along the tree (your traitTrend / simulateAbundanceData1)
#  3) Simulates species-level niche coefficients B from an evolutionary model
#     that is *not* the same as the fitted model (intentional misspecification)
#  4) Simulates observation-level data (long-form) for an LGCP-background-style model
#     y[n] ~ Poisson(exp(x[n,] %*% B[sp[n],] + f[n] + offset[n]))
#  5) Induces (a) fully missing species and (b) data-deficient species
#  6) Fits via TMB/RTMB-style Laplace (fast) and runs conditional prediction
#     of B for missing species using the fitted phylo covariance.
#
# Notes:
# - This script uses TMB. If you use RTMB, the workflow is the same at the R level.
# - You must point `CPP_FILE` to your OPT2 C++ model (alpha,rho,sigma_f as parameters).
# - This is designed to be flexible: vary ntips, #missing species, #data-deficient species,
#   degree of deficiency, N sites, etc.
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
    library(ape)
    library(MASS) # mvrnorm
    library(TMB)
})

# -----------------------------
# USER SETTINGS (edit these)
# -----------------------------
SEED <- 1

# simulation sizes
J_TIPS <- 60 # number of species (tips)
K_COV <- 3 # number of env covariates (include intercept if you want)
N_SITES <- 20000 # number of observation rows (sites x sampled species)

# trait simulation (your existing framework)
NTRAITS <- 2
TREE_TIME <- 1.0
TRAIT_MEAN_ROOT <- rep(0, NTRAITS)
TRAIT_COV_ROOT <- diag(NTRAITS) * 1.0
TRAIT_SIGMA2 <- 1.0
TRAIT_ALPHA <- 0.90 # "rubberband" strength in your traitTrend

# evolutionary beta generation (misspecified vs fitted GP-on-B)
# B_true[j, k] = intercept_k + trait_effects[k,] %*% traits_j + phylo_residual + iid_noise
B_INTERCEPT <- rep(0.0, K_COV)
TRAIT_EFFECTS <- matrix(0.0, nrow = K_COV, ncol = NTRAITS)
TRAIT_EFFECTS[2, 1] <- 1.0 # example: covariate 2 depends on trait-1
TRAIT_EFFECTS[3, 2] <- -0.7 # example: covariate 3 depends on trait-2

# residual phylogenetic structure for B_true (Brownian-ish via vcv)
B_PHYLO_RESID_SD <- 0.4
B_IID_NOISE_SD <- 0.1

# observation-level random effect
SIGMA_F_TRUE <- 0.5

# offsets
OFFSET_CONST <- 0.0

# missingness / deficiency controls
N_MISSING_SPECIES <- 8 # fully missing species (no rows)
N_DEFICIENT_SPECIES <- 6 # species present but heavily thinned to near-zero y
DEF_KEEP_N_PER_SPECIES <- 30 # keep only this many rows for each deficient species

# how to pick missing/deficient species
# options: "alphabetical" or "random"
MISSING_PICK <- "alphabetical"
DEF_PICK <- "alphabetical"

# model to fit (OPT2 C++ for closest equivalence to your Stan)
CPP_FILE <- "~/clim_risk_phylosdm/scripts/05_lgcp_corrected.cpp" # <-- FIX PATH

# TMB optimization
TMB_ITER <- 400

# -----------------------------
# Your trait simulation functions (minimal subset)
# -----------------------------
traitTrend <- function(tr, ntraits,
                       mean_root_vector,
                       covar_traits,
                       sigma2,
                       alpha) {
    tr_table <- as.data.frame(tr$edge)
    colnames(tr_table) <- c("node_from", "node_to")
    tr_table$length <- tr$edge.length

    trait_matrix <- matrix(0, ntraits, nrow = nrow(tr_table))
    root_num <- tr_table[1, 1]
    root_traits <- MASS::mvrnorm(1, mean_root_vector, covar_traits)

    for (i in seq_len(nrow(tr_table))) {
        if (tr_table[i, 1] == root_num) {
            next_traits <- alpha * root_traits +
                MASS::mvrnorm(1,
                    mu = rep(0, ntraits),
                    Sigma = diag(sigma2 * tr_table[i, "length"], nrow = ntraits)
                )
            trait_matrix[i, ] <- next_traits
        } else {
            old_state_num <- which(tr_table[i, 1] == tr_table[, 2])
            old_state <- trait_matrix[old_state_num, ]
            next_traits <- alpha * old_state +
                MASS::mvrnorm(1,
                    mu = rep(0, ntraits),
                    Sigma = diag(sigma2 * tr_table[i, "length"], nrow = ntraits)
                )
            trait_matrix[i, ] <- next_traits
        }
    }

    trait_matrix <- rbind(root_traits, trait_matrix)
    colnames(trait_matrix) <- paste0("trait-", seq_len(ncol(trait_matrix)))
    node_values <- c(root_num, tr_table$node_to)
    trait_df <- as.data.frame(trait_matrix)
    trait_df$node <- node_values

    list(trait_df = trait_df, tr = tr)
}

splitByNodesTips <- function(table, tree) {
    Ntips <- length(tree$tip.label)
    tip_index <- which(table$node %in% 1:Ntips)
    node_index <- setdiff(seq_len(nrow(table)), tip_index)
    list(
        tip_table  = table[tip_index, , drop = FALSE],
        node_table = table[node_index, , drop = FALSE]
    )
}

simulateAbundanceData1 <- function(TREE, simulation_params) {
    tree_time <- simulation_params$tree_time
    ntraits <- simulation_params$ntraits
    mean_root_vector <- simulation_params$mean_root_vector
    covar_traits <- simulation_params$covar_traits
    sigma2 <- simulation_params$sigma2
    alpha <- simulation_params$alpha

    # (tree_time unused in your original; keep for compatibility)
    trait_tree <- traitTrend(
        tr = TREE, ntraits = ntraits,
        mean_root_vector = mean_root_vector,
        covar_traits = covar_traits,
        sigma2 = sigma2, alpha = alpha
    )

    tr <- trait_tree$tr
    trait_tbl <- trait_tree$trait_df

    split_tbl <- splitByNodesTips(trait_tbl, tr)
    tip_traits <- split_tbl$tip_table
    rownames(tip_traits) <- tr$tip.label

    list(tip_traits = tip_traits, trait_tbl = trait_tbl, tr = tr)
}

# -----------------------------
# helper: pick species indices
# -----------------------------
pick_species <- function(tip_labels, n, method = c("alphabetical", "random"), seed = NULL) {
    method <- match.arg(method)
    if (!is.null(seed)) set.seed(seed)
    if (n <= 0) {
        return(character())
    }
    if (method == "alphabetical") {
        return(sort(tip_labels)[seq_len(n)])
    }
    sample(tip_labels, n, replace = FALSE)
}

# -----------------------------
# 1) Simulate phylogeny + traits
# -----------------------------
set.seed(SEED)

# Tree: use birth-death like random tree. If you want ultrametric enforcement,
# you can replace this with a real tree or run chronos().
tr <- ape::rtree(J_TIPS)
tr$tip.label <- sprintf("sp_%03d", seq_len(J_TIPS))

sim_params <- list(
    tree_time = TREE_TIME,
    ntraits = NTRAITS,
    mean_root_vector = TRAIT_MEAN_ROOT,
    covar_traits = TRAIT_COV_ROOT,
    sigma2 = TRAIT_SIGMA2,
    alpha = TRAIT_ALPHA
)

trait_out <- simulateAbundanceData1(tr, sim_params)
tr <- trait_out$tr
tip_traits <- trait_out$tip_traits

# traits matrix (J x NTRAITS)
Z <- as.matrix(tip_traits[, paste0("trait-", seq_len(NTRAITS)), drop = FALSE])
stopifnot(nrow(Z) == length(tr$tip.label))

# phylo distance matrix for the fitted model
D_phylo <- as.matrix(cophenetic(tr))

# Brownian-ish residual covariance for beta simulation
V_phylo <- ape::vcv(tr, corr = TRUE) # J x J

# -----------------------------
# 2) Simulate B_true from evolutionary model
#     (misspecified vs fitted GP-on-B)
# -----------------------------
J <- length(tr$tip.label)
K <- K_COV

B_true <- matrix(NA_real_, nrow = J, ncol = K)
rownames(B_true) <- tr$tip.label
colnames(B_true) <- paste0("x", seq_len(K))

# phylo residual term per covariate k
for (k in seq_len(K)) {
    mu_k <- B_INTERCEPT[k] + as.numeric(Z %*% TRAIT_EFFECTS[k, ])
    phylo_resid_k <- as.numeric(MASS::mvrnorm(1, mu = rep(0, J), Sigma = (B_PHYLO_RESID_SD^2) * V_phylo))
    iid_noise_k <- rnorm(J, mean = 0, sd = B_IID_NOISE_SD)
    B_true[, k] <- mu_k + phylo_resid_k + iid_noise_k
}

# -----------------------------
# 3) Simulate observation-level data (long form)
# -----------------------------
# Covariates (N x K). If you want an intercept, put it as column 1:
X <- matrix(rnorm(N_SITES * K), nrow = N_SITES, ncol = K)
colnames(X) <- colnames(B_true)

# species assignment per row
species_labels <- tr$tip.label
sp <- sample(species_labels, size = N_SITES, replace = TRUE)
species_index <- match(sp, species_labels) # 1..J

# site RE
f_true <- rnorm(N_SITES, 0, SIGMA_F_TRUE)

offset <- rep(OFFSET_CONST, N_SITES)

eta <- rowSums(X * B_true[species_index, , drop = FALSE]) + f_true + offset
y <- rpois(N_SITES, lambda = exp(eta))

# -----------------------------
# 4) Induce missing species + data-deficient species
# -----------------------------
missing_species <- pick_species(species_labels, N_MISSING_SPECIES, MISSING_PICK, seed = SEED + 10)
remaining_species <- setdiff(species_labels, missing_species)

def_species <- pick_species(remaining_species, N_DEFICIENT_SPECIES, DEF_PICK, seed = SEED + 20)

# Drop all rows for missing species (fully missing)
keep_idx <- !(sp %in% missing_species)

# For deficient species, keep only DEF_KEEP_N_PER_SPECIES rows (others dropped)
if (length(def_species) > 0 && DEF_KEEP_N_PER_SPECIES > 0) {
    for (ds in def_species) {
        ds_idx_all <- which(sp == ds & keep_idx)
        if (length(ds_idx_all) > DEF_KEEP_N_PER_SPECIES) {
            drop_extra <- sample(ds_idx_all, size = length(ds_idx_all) - DEF_KEEP_N_PER_SPECIES, replace = FALSE)
            keep_idx[drop_extra] <- FALSE
        }
    }
}

# Apply subsetting
X_obs <- X[keep_idx, , drop = FALSE]
y_obs <- y[keep_idx]
offset_obs <- offset[keep_idx]
sp_obs <- sp[keep_idx]
species_obs_index <- match(sp_obs, species_labels)

N_obs <- length(y_obs)

cat("\n--- SIM SUMMARY ---\n")
cat("J tips:", J, "\n")
cat("N total rows:", N_SITES, "\n")
cat("N observed rows after missingness/deficiency:", N_obs, "\n")
cat("Missing species (fully absent):", length(missing_species), "\n")
cat("Deficient species (thinned):", length(def_species), "\n\n")

# -----------------------------
# 5) Fit model via TMB / RTMB-style Laplace
# -----------------------------
stopifnot(file.exists(CPP_FILE))

model_base <- sub("\\.cpp$", "", basename(CPP_FILE))

# compile / load
so_path <- file.path(dirname(CPP_FILE), paste0(model_base, .Platform$dynlib.ext))
if (!file.exists(so_path) || file.info(CPP_FILE)$mtime > file.info(so_path)$mtime) {
    message("Compiling TMB model: ", CPP_FILE)
    old <- getwd()
    setwd(dirname(CPP_FILE))
    TMB::compile(basename(CPP_FILE))
    setwd(old)
}
dyn.load(TMB::dynlib(file.path(dirname(CPP_FILE), model_base)))

# data for fitted model (matches your LGCP background structure)
dat <- list(
    N = N_obs,
    J = J,
    K = K,
    species = as.integer(species_obs_index), # 1-based
    X = X_obs,
    y = as.numeric(y_obs),
    D_phylo = D_phylo,
    offset = as.numeric(offset_obs)
)

# init (OPT2-style)
par <- list(
    B = matrix(0, J, K),
    alpha = 1.0,
    rho = median(D_phylo),
    sigma_f = 0.5,
    f = rep(0, N_obs)
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
    control = list(iter.max = TMB_ITER, eval.max = TMB_ITER * 2)
)
t2 <- Sys.time()

par_hat <- obj$env$parList(opt$par)
B_hat <- par_hat$B
alpha_hat <- par_hat$alpha
rho_hat <- par_hat$rho
sigma_f_hat <- par_hat$sigma_f

cat("\n--- FIT SUMMARY ---\n")
cat("Convergence code:", opt$convergence, "\n")
cat("Runtime (mins):", round(as.numeric(difftime(t2, t1, units = "mins")), 2), "\n")
cat("alpha_hat:", round(alpha_hat, 3), " rho_hat:", round(rho_hat, 3), " sigma_f_hat:", round(sigma_f_hat, 3), "\n")

# -----------------------------
# 6) Conditional prediction of B for missing species (and optional deficient species)
#     Under fitted GP prior: vec over species for each covariate k:
#       B[,k] ~ MVN(0, Kphy(alpha_hat, rho_hat))
#     Predict B[missing, k] | B[observed, k]
# -----------------------------
make_Kphy <- function(D, alpha, rho, jitter = 1e-6) {
    K <- (alpha^2) * exp(-(D^2) / (2 * rho^2))
    diag(K) <- diag(K) + jitter
    K
}

cond_pred_B <- function(B_hat, D_phylo, alpha_hat, rho_hat,
                        target_species, conditioned_species,
                        jitter = 1e-8) {
    # returns list with mean (|target| x K) and cov list per k
    Kphy <- make_Kphy(D_phylo, alpha_hat, rho_hat)

    idx_t <- match(target_species, rownames(B_hat))
    idx_c <- match(conditioned_species, rownames(B_hat))
    stopifnot(all(!is.na(idx_t)), all(!is.na(idx_c)))

    K_cc <- Kphy[idx_c, idx_c, drop = FALSE]
    K_tc <- Kphy[idx_t, idx_c, drop = FALSE]
    K_tt <- Kphy[idx_t, idx_t, drop = FALSE]

    # stable solve via chol
    K_cc <- K_cc + diag(jitter, nrow(K_cc))
    L <- chol(K_cc)
    invKcc_times <- function(M) {
        # solve(K_cc, M) using chol
        backsolve(L, forwardsolve(t(L), M))
    }

    Kcc_inv_Bc <- invKcc_times(B_hat[idx_c, , drop = FALSE])

    mu_t <- K_tc %*% Kcc_inv_Bc # |target| x K
    cov_t <- K_tt - K_tc %*% invKcc_times(t(K_tc)) # |target| x |target|
    list(mean = mu_t, cov = cov_t)
}

# Decide what you want to predict:
# - "missing_species" are fully absent from data (hardest)
# - you may also want to include def_species as "targets" to see shrinkage
observed_species <- setdiff(species_labels, missing_species)

pred_missing <- NULL
if (length(missing_species) > 0) {
    pred_missing <- cond_pred_B(
        B_hat = B_hat,
        D_phylo = D_phylo,
        alpha_hat = alpha_hat,
        rho_hat = rho_hat,
        target_species = missing_species,
        conditioned_species = observed_species
    )
}

# -----------------------------
# 7) Diagnostics for betas (ballpark evaluation)
# -----------------------------
beta_summary_over_species <- function(Bmat) {
    data.frame(
        k = seq_len(ncol(Bmat)),
        mean_over_species = colMeans(Bmat),
        sd_over_species = apply(Bmat, 2, sd)
    )
}

true_sum <- beta_summary_over_species(B_true)
hat_sum <- beta_summary_over_species(B_hat)

out <- merge(true_sum, hat_sum, by = "k", suffixes = c("_true", "_hat"))

cat("\n--- B SUMMARY (across ALL species) ---\n")
print(out)

# Evaluate missing-species prediction quality (only possible in simulation)
if (!is.null(pred_missing)) {
    B_true_missing <- B_true[missing_species, , drop = FALSE]
    B_pred_missing <- pred_missing$mean
    colnames(B_pred_missing) <- colnames(B_true_missing)
    rownames(B_pred_missing) <- missing_species

    rmse <- sqrt(mean((B_true_missing - B_pred_missing)^2))
    corv <- cor(as.numeric(B_true_missing), as.numeric(B_pred_missing))

    cat("\n--- MISSING SPECIES B PREDICTION (posterior mean under fitted GP) ---\n")
    cat("RMSE:", round(rmse, 3), "\n")
    cat("Correlation:", round(corv, 3), "\n")

    # optional: per-covariate RMSE
    rmse_k <- sapply(seq_len(K), function(k) sqrt(mean((B_true_missing[, k] - B_pred_missing[, k])^2)))
    names(rmse_k) <- colnames(B_true_missing)
    cat("RMSE by covariate:\n")
    print(round(rmse_k, 3))
}

# -----------------------------
# 8) Save everything (optional)
# -----------------------------
res <- list(
    settings = list(
        seed = SEED, J = J, K = K, N_sites = N_SITES, N_obs = N_obs,
        missing_species = missing_species,
        deficient_species = def_species,
        DEF_KEEP_N_PER_SPECIES = DEF_KEEP_N_PER_SPECIES
    ),
    tree = tr,
    traits = Z,
    D_phylo = D_phylo,
    truth = list(B_true = B_true, sigma_f_true = SIGMA_F_TRUE),
    observed_data = list(
        X = X_obs, y = y_obs, species = sp_obs, offset = offset_obs
    ),
    fit = list(
        opt = opt,
        B_hat = B_hat,
        alpha_hat = alpha_hat,
        rho_hat = rho_hat,
        sigma_f_hat = sigma_f_hat
    ),
    pred_missing = pred_missing
)

save(res, file = "sim_evolutionary_rtmb_condpred.RData")
cat("\nSaved: sim_evolutionary_rtmb_condpred.RData\n")
# ============================================================
