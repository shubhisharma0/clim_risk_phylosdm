# ====================== CLIM_RISK_PHYLOSDM | COMPUTE SENSITIVITY ======================
# Computes climate sensitivity based on model coefficients and species distributions
# Sensitivity = partial derivative of habitat suitability with respect to climate
# For quadratic terms: d/dx(β1*x + β2*x²) = β1 + 2*β2*x
# For linear terms: d/dx(β*x) = β
# ========================================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
    library(terra)
    library(rstan)
})

# ---- Set-up paths ----
HPC <- Sys.getenv("HPC")
if (HPC != "FALSE") {
    root <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm"
    epath <- "/vast/palmer/pi/jetz/ss4224/env"
    message("Running on HPC")

    # Job-specific temp directory
    slurm_job_id <- Sys.getenv("SLURM_JOB_ID")
    slurm_array_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
    if (slurm_job_id != "") {
        job_id <- if (slurm_array_id != "") paste0(slurm_job_id, "_", slurm_array_id) else slurm_job_id
    } else {
        job_id <- paste0("local_", Sys.getpid())
    }
    temp_dir <- file.path("/vast/palmer/scratch/jetz/ss4224/clim_risk/temp", job_id)
    dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

    Sys.setenv(TMPDIR = temp_dir)
    terra::terraOptions(
        tempdir = temp_dir,
        memfrac = 0.1,
        todisk = TRUE
    )
} else {
    root <- "~/clim_risk_phylosdm"
    epath <- "~/env"
    temp_dir <- tempdir()
    message("Running locally")
}

# ---- Command line arguments or manual set ----
# ---- Command line arguments or manual set ----
if (interactive()) {
    EXP_ROOT <- "v0"
    EXP_ID <- "sub1000"
    DATASET <- "amphibians"
    CLUSTER <- "Kass1"
    FSP <- "ALL"
    REPNO <- 1
    NREP <- 1
    MODEL_TYPE <- "STAN"
    MODEL_NAME <- "LGCP_background"
} else {
    args <- commandArgs(trailingOnly = TRUE)
    EXP_ROOT <- args[1]
    EXP_ID <- args[2]
    DATASET <- args[3]
    CLUSTER <- args[4]
    FSP <- args[5]
    REPNO <- as.integer(args[6])
    NREP <- as.integer(args[7])
    MODEL_TYPE <- args[8]
    MODEL_NAME <- args[9]
}
message(sprintf("Computing climate sensitivity for cluster: %s, rep: %d", CLUSTER, REPNO))

# ---- Paths ----
analysis_dir <- file.path(root, "analysis")
raw_dir <- file.path(root, "raw_data")
res_dir <- file.path(root, "res")
dpath <- file.path(raw_dir, DATASET)
data_dir <- file.path(dpath, CLUSTER)

# Input/output directories
if (HPC != "FALSE") {
    binary_dir <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm/analysis/v0/binary_maps"
    sensitivity_dir <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm/analysis/v0/sensitivity"
} else {
    binary_dir <- file.path(analysis_dir, EXP_ROOT, "binary_maps")
    sensitivity_dir <- file.path(analysis_dir, EXP_ROOT, "sensitivity")
}

cond_pred_dir <- file.path(analysis_dir, EXP_ROOT, "cond_pred")
eval_dir <- file.path(analysis_dir, EXP_ROOT, "eval")

dir.create(sensitivity_dir, recursive = TRUE, showWarnings = FALSE)

# ---- File naming ----
base_prefix <- paste0(EXP_ROOT, "_", EXP_ID, "_", CLUSTER, "_", FSP, "_rep_", REPNO)

# ---- Helper functions ----
cleanup_terra_temps <- function(temp_dir = NULL) {
    if (is.null(temp_dir)) {
        temp_dir <- terra::terraOptions()$tempdir
    }
    if (!is.null(temp_dir) && dir.exists(temp_dir)) {
        temp_files <- list.files(temp_dir, pattern = "^spat_", full.names = TRUE)
        if (length(temp_files) > 0) {
            suppressWarnings(unlink(temp_files))
        }
    }
    gc()
}

# ---- Load TRAINED species list from model_data ----
model_data_file <- file.path(data_dir, paste0(base_prefix, "_model_data.Rdata"))
if (!file.exists(model_data_file)) {
    stop(sprintf("Model data file not found: %s", model_data_file))
}
load(model_data_file) # loads 'model_data'

species_names <- model_data$species_names
J <- length(species_names)
message(sprintf("Found %d trained species in cluster %s", J, CLUSTER))

# ---- Load scaling stats ----
scales_file <- file.path(data_dir, paste0(CLUSTER, "_env_scales.csv"))
if (!file.exists(scales_file)) {
    stop(sprintf("Scales file not found: %s", scales_file))
}
scales_df <- read.csv(scales_file, stringsAsFactors = FALSE)

scales_df$variable <- c(
    "meanTemp", "tempSeason", "precipWet", "precipSeason",
    "cloudCover", "EVI", "TRI", "elevation"
)
message("Loaded scaling statistics")

# ---- Load model fit ----
result_file <- file.path(res_dir, paste0(base_prefix, "_", MODEL_TYPE, "_fit.Rdata"))
if (!file.exists(result_file)) {
    stop(sprintf("Result file not found: %s", result_file))
}
load(result_file) # loads 'result'
stan_fit <- result$fit
message("Loaded model fit")

# ---- Extract posterior samples ----
posterior_B <- rstan::extract(stan_fit, pars = "B")$B # (S, J, K)
S <- dim(posterior_B)[1]
J_model <- dim(posterior_B)[2]
K <- dim(posterior_B)[3]
message(sprintf("Posterior: S=%d samples, J=%d species, K=%d predictors", S, J_model, K))

if (J != J_model) {
    stop(sprintf("Species count mismatch: model_data has %d, posterior has %d", J, J_model))
}

# ---- Load conditional predictions ----
cond_pred_file <- file.path(cond_pred_dir, paste0(base_prefix, "_cond_pred.Rdata"))
has_cond_pred <- file.exists(cond_pred_file)
if (has_cond_pred) {
    load(cond_pred_file) # loads 'B_cond_list'
    message(sprintf("Loaded conditional predictions: %d species", length(B_cond_list)))
} else {
    message("Conditional predictions not found, will use model coefficients only")
    B_cond_list <- NULL
}

# ---- Load best model selection ----
best_model_file <- file.path(eval_dir, paste0(base_prefix, "_best_model.csv"))
if (file.exists(best_model_file)) {
    best_model_df <- read.csv(best_model_file, stringsAsFactors = FALSE)
    message(sprintf("Loaded best model selection: %d species", nrow(best_model_df)))
} else {
    message("Best model file not found, defaulting to 'model' for all species")
    best_model_df <- data.frame(
        species_name = species_names,
        best_model = rep("model", length(species_names)),
        stringsAsFactors = FALSE
    )
}

# ---- Define coefficient indices ----
COEF_IDX <- list(
    meanTemp = 2,
    meanTemp2 = 3,
    tempSeason = 4,
    precipWet = 5,
    precipWet2 = 6,
    precipSeason = 7
)

# ---- Climate variables for sensitivity ----
clim_vars <- list(
    bio1_temp = list(
        linear_coef = "meanTemp",
        quad_coef = "meanTemp2",
        env_file = "CHELSA_bio_1.tif",
        scale_var = "meanTemp"
    ),
    bio4_tempSeas = list(
        linear_coef = "tempSeason",
        quad_coef = NULL,
        env_file = "CHELSA_bio_4.tif",
        scale_var = "tempSeason"
    ),
    bio13_precipWet = list(
        linear_coef = "precipWet",
        quad_coef = "precipWet2",
        env_file = "CHELSA_bio_13.tif",
        scale_var = "precipWet"
    ),
    bio15_precipSeas = list(
        linear_coef = "precipSeason",
        quad_coef = NULL,
        env_file = "CHELSA_bio_15.tif",
        scale_var = "precipSeason"
    )
)

# ---- Load climate rasters ----
message("\nLoading climate rasters...")
clim_rasters <- list()
for (var_name in names(clim_vars)) {
    f <- file.path(epath, clim_vars[[var_name]]$env_file)
    if (!file.exists(f)) {
        stop(sprintf("Climate file not found: %s", f))
    }
    clim_rasters[[var_name]] <- terra::rast(f)
    message(sprintf("  Loaded %s", clim_vars[[var_name]]$env_file))
}

# ---- Load world vector for ocean masking ----
message("\nLoading world vector for ocean masking...")
world_file <- file.path(root, "raw_data", "world.Rdata")
if (file.exists(world_file)) {
    load(world_file) # loads 'world' (sf object)
    world_vect <- terra::vect(world)
    message("  Loaded world vector")
} else {
    message("  Warning: world.Rdata not found, skipping ocean masking")
    world_vect <- NULL
}

# ---- Initialize results storage ----
sensitivity_results <- data.frame(
    species_name = character(),
    range_cells = integer(),
    bio1_temp_beta_linear = numeric(),
    bio1_temp_beta_quad = numeric(),
    bio1_temp_sens_mean = numeric(),
    bio1_temp_sens_sd = numeric(),
    bio4_tempSeas_beta = numeric(),
    bio4_tempSeas_sens_mean = numeric(),
    bio4_tempSeas_sens_sd = numeric(),
    bio13_precipWet_beta_linear = numeric(),
    bio13_precipWet_beta_quad = numeric(),
    bio13_precipWet_sens_mean = numeric(),
    bio13_precipWet_sens_sd = numeric(),
    bio15_precipSeas_beta = numeric(),
    bio15_precipSeas_sens_mean = numeric(),
    bio15_precipSeas_sens_sd = numeric(),
    stringsAsFactors = FALSE
)

# ---- Track processing status ----
success_species <- character()
missing_binary <- character()
error_species <- character()

# ---- Main loop over species ----
message(sprintf("\n========== Computing sensitivity for %d species ==========", J))

for (j in seq_len(J)) {
    sp_name <- species_names[j]
    message(sprintf("\n[%d/%d] Processing species: %s", j, J, sp_name))

    tryCatch(
        {
            # ---- Load binary map ----
            binary_file <- file.path(binary_dir, paste0(base_prefix, "_", sp_name, "_binary.tif"))
            if (!file.exists(binary_file)) {
                message(sprintf("  Binary map not found: %s", basename(binary_file)))
                missing_binary <- c(missing_binary, sp_name)
                next
            }

            binary_rast <- terra::rast(binary_file)
            sp_extent <- terra::ext(binary_rast)

            # Count range cells
            range_cells <- terra::global(binary_rast, "sum", na.rm = TRUE)[[1]]
            message(sprintf("  Range cells: %d", range_cells))

            if (range_cells == 0 || is.na(range_cells)) {
                message("  No presence cells in binary map, skipping")
                next
            }

            # ---- Crop world vector to species extent for faster masking ----
            if (!is.null(world_vect)) {
                world_crop <- terra::crop(world_vect, sp_extent)
            } else {
                world_crop <- NULL
            }

            # ---- Get coefficients for this species ----
            best_row <- which(best_model_df$species_name == sp_name)
            if (length(best_row) == 0) {
                best_model <- "model"
            } else {
                best_model <- best_model_df$best_model[best_row]
            }

            if (best_model == "cond_pred" && has_cond_pred && sp_name %in% names(B_cond_list)) {
                B_cond_sp <- as.matrix(B_cond_list[[sp_name]])
                beta_bar <- as.vector(colMeans(B_cond_sp))
                message("  Using conditional prediction coefficients")
            } else {
                beta_bar <- as.vector(colMeans(posterior_B[, j, , drop = FALSE]))
                message("  Using model coefficients")
            }

            # ---- Initialize results for this species ----
            sp_results <- list(
                species_name = sp_name,
                range_cells = as.integer(range_cells)
            )

            # ---- Compute sensitivity for each climate variable ----
            for (var_name in names(clim_vars)) {
                var_info <- clim_vars[[var_name]]

                # Get coefficients
                beta_linear <- beta_bar[COEF_IDX[[var_info$linear_coef]]]

                if (!is.null(var_info$quad_coef)) {
                    beta_quad <- beta_bar[COEF_IDX[[var_info$quad_coef]]]
                } else {
                    beta_quad <- NULL
                }

                # Crop and resample climate raster to match binary map
                clim_crop <- terra::crop(clim_rasters[[var_name]], sp_extent, snap = "out")
                clim_resamp <- terra::resample(clim_crop, binary_rast, method = "bilinear")

                # Scale the climate raster
                scale_row <- which(scales_df$variable == var_info$scale_var)
                mu <- scales_df$mean[scale_row]
                sd_val <- scales_df$sd[scale_row]
                if (sd_val == 0) sd_val <- 1
                clim_scaled <- (clim_resamp - mu) / sd_val

                # Compute sensitivity
                if (!is.null(beta_quad)) {
                    sensitivity_rast <- beta_linear + 2 * beta_quad * clim_scaled
                    sp_results[[paste0(var_name, "_beta_linear")]] <- beta_linear
                    sp_results[[paste0(var_name, "_beta_quad")]] <- beta_quad
                } else {
                    sensitivity_rast <- terra::rast(clim_scaled)
                    terra::values(sensitivity_rast) <- beta_linear
                    sp_results[[paste0(var_name, "_beta")]] <- beta_linear
                }

                # Mask to species range
                sensitivity_masked <- terra::mask(sensitivity_rast, binary_rast, maskvalues = 0)

                # Mask out ocean areas
                if (!is.null(world_crop)) {
                    sensitivity_masked <- terra::mask(sensitivity_masked, world_crop)
                }

                # Compute statistics within range (after ocean masking)
                sens_stats <- terra::global(sensitivity_masked, c("mean", "sd"), na.rm = TRUE)
                sp_results[[paste0(var_name, "_sens_mean")]] <- sens_stats$mean
                sp_results[[paste0(var_name, "_sens_sd")]] <- sens_stats$sd

                # Write sensitivity raster
                out_file <- file.path(
                    sensitivity_dir,
                    paste0(base_prefix, "_", sp_name, "_sensitivity_", var_name, ".tif")
                )
                terra::writeRaster(sensitivity_masked, out_file, overwrite = TRUE)

                # Clean up
                rm(clim_crop, clim_resamp, clim_scaled, sensitivity_rast, sensitivity_masked)
            }

            message("  Wrote 4 sensitivity rasters")

            # Print summary
            message(sprintf(
                "  Bio1 sensitivity (mean): %.4f ± %.4f",
                sp_results$bio1_temp_sens_mean, sp_results$bio1_temp_sens_sd
            ))
            message(sprintf(
                "  Bio13 sensitivity (mean): %.4f ± %.4f",
                sp_results$bio13_precipWet_sens_mean, sp_results$bio13_precipWet_sens_sd
            ))

            # ---- Add to results ----
            sensitivity_results <- rbind(sensitivity_results, as.data.frame(sp_results))

            success_species <- c(success_species, sp_name)

            # Clean up
            rm(binary_rast, world_crop)
            cleanup_terra_temps(temp_dir)
        },
        error = function(e) {
            message(sprintf("  ERROR: %s", e$message))
            error_species <- c(error_species, sp_name)
            cleanup_terra_temps(temp_dir)
        }
    )
}

# ---- Save results ----
message("\n========== Saving results ==========")

results_csv <- file.path(sensitivity_dir, paste0(base_prefix, "_sensitivity_summary.csv"))
write.csv(sensitivity_results, results_csv, row.names = FALSE)
message(sprintf("Saved sensitivity summary: %s", results_csv))

results_rdata <- file.path(sensitivity_dir, paste0(base_prefix, "_sensitivity_summary.Rdata"))
save(sensitivity_results, file = results_rdata)
message(sprintf("Saved sensitivity Rdata: %s", results_rdata))

# ---- Print summary ----
message("\n========== Sensitivity Summary ==========")
message(sprintf("Species processed: %d", length(success_species)))
message(sprintf("Missing binary maps: %d", length(missing_binary)))
message(sprintf("Errors: %d", length(error_species)))

if (nrow(sensitivity_results) > 0) {
    message(sprintf(
        "\nMean bio1 (temp) sensitivity: %.4f ± %.4f",
        mean(sensitivity_results$bio1_temp_sens_mean, na.rm = TRUE),
        sd(sensitivity_results$bio1_temp_sens_mean, na.rm = TRUE)
    ))
    message(sprintf(
        "Mean bio13 (precip) sensitivity: %.4f ± %.4f",
        mean(sensitivity_results$bio13_precipWet_sens_mean, na.rm = TRUE),
        sd(sensitivity_results$bio13_precipWet_sens_mean, na.rm = TRUE)
    ))
}

if (length(missing_binary) > 0 && length(missing_binary) <= 10) {
    message("\nMissing binary maps:")
    message(paste("  ", missing_binary, collapse = "\n"))
}

if (length(error_species) > 0) {
    message("\nSpecies with errors:")
    message(paste("  ", error_species, collapse = "\n"))
}

# ---- Final cleanup ----
message("\nFinal cleanup...")
if (HPC != "FALSE" && exists("temp_dir") && dir.exists(temp_dir)) {
    unlink(temp_dir, recursive = TRUE)
    message(sprintf("Removed job temp directory: %s", temp_dir))
}

message(sprintf("\n========== Done! Sensitivity results saved to: %s ==========", sensitivity_dir))
