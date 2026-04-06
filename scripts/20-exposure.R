# ====================== CLIM_RISK_PHYLOSDM | COMPUTE EXPOSURE ======================
# Computes climate exposure based on binary predictions of the species distribution maps
# Exposure = climate change within species range (masked by binary SDM)
# ========================================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
    library(terra)
})

# ---- Set-up paths ----
HPC <- Sys.getenv("HPC")
if (HPC != "FALSE") {
    root <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm"
    epath <- "/vast/palmer/pi/jetz/ss4224/env"
    cc_dir <- "/vast/palmer/pi/jetz/ss4224/env_ssp370_2041-2060"
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
    cc_dir <- "~/env_ssp370_2041-2060"
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
message(sprintf("Computing climate exposure for cluster: %s, rep: %d", CLUSTER, REPNO))

# ---- Paths ----
analysis_dir <- file.path(root, "analysis")
raw_dir <- file.path(root, "raw_data")

# Input/output directories
if (HPC != "FALSE") {
    binary_dir <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm/analysis/v0/binary_maps"
    exposure_dir <- "/vast/palmer/pi/jetz/ss4224/clim_risk_phylosdm/analysis/v0/clim_exposure"
} else {
    binary_dir <- file.path(analysis_dir, EXP_ROOT, "binary_maps")
    exposure_dir <- file.path(analysis_dir, EXP_ROOT, "clim_exposure")
}

dir.create(exposure_dir, recursive = TRUE, showWarnings = FALSE)

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

# ---- Load species list ----
splist_file <- file.path(raw_dir, DATASET, "spList.Rdata")
if (!file.exists(splist_file)) {
    stop(sprintf("Species list not found: %s", splist_file))
}
load(splist_file) # loads 'spList'

species_names <- spList[[CLUSTER]]
if (is.null(species_names)) {
    stop(sprintf("Cluster '%s' not found in spList", CLUSTER))
}
J <- length(species_names)
message(sprintf("Found %d species in cluster %s", J, CLUSTER))

# ---- Load climate layers ----
message("\nLoading climate layers...")

clim_layers <- c(
    "CHELSA_bio_1.tif", # Mean annual temperature
    "CHELSA_bio_4.tif", # Temperature seasonality
    "CHELSA_bio_13.tif", # Precipitation of wettest month
    "CHELSA_bio_15.tif" # Precipitation seasonality
)

clim_names <- c("bio1_temp", "bio4_tempSeas", "bio13_precipWet", "bio15_precipSeas")

# Load baseline climate
r_clim_baseline <- list()
for (i in seq_along(clim_layers)) {
    f <- file.path(epath, clim_layers[i])
    if (!file.exists(f)) {
        stop(sprintf("Baseline climate file not found: %s", f))
    }
    r_clim_baseline[[i]] <- terra::rast(f)
}
names(r_clim_baseline) <- clim_names
message(sprintf("  Loaded %d baseline climate layers", length(r_clim_baseline)))

# Load future climate (SSP370 2041-2060)
clim_cc_layers <- c(
    "CHELSA_bio_1_cc.tif",
    "CHELSA_bio_4_cc.tif",
    "CHELSA_bio_13_cc.tif",
    "CHELSA_bio_15_cc.tif"
)

r_clim_future <- list()
for (i in seq_along(clim_cc_layers)) {
    f <- file.path(cc_dir, clim_cc_layers[i])
    if (!file.exists(f)) {
        stop(sprintf("Future climate file not found: %s", f))
    }
    r_clim_future[[i]] <- terra::rast(f)
}
names(r_clim_future) <- clim_names
message(sprintf("  Loaded %d future climate layers", length(r_clim_future)))

# ---- Initialize results storage ----
exposure_results <- data.frame(
    species_name = character(),
    range_cells = integer(),
    bio1_temp_mean_baseline = numeric(),
    bio1_temp_mean_future = numeric(),
    bio1_temp_mean_change = numeric(),
    bio1_temp_sd_change = numeric(),
    bio4_tempSeas_mean_baseline = numeric(),
    bio4_tempSeas_mean_future = numeric(),
    bio4_tempSeas_mean_change = numeric(),
    bio4_tempSeas_sd_change = numeric(),
    bio13_precipWet_mean_baseline = numeric(),
    bio13_precipWet_mean_future = numeric(),
    bio13_precipWet_mean_change = numeric(),
    bio13_precipWet_sd_change = numeric(),
    bio15_precipSeas_mean_baseline = numeric(),
    bio15_precipSeas_mean_future = numeric(),
    bio15_precipSeas_mean_change = numeric(),
    bio15_precipSeas_sd_change = numeric(),
    stringsAsFactors = FALSE
)

# ---- Track processing status ----
success_species <- character()
missing_binary <- character()
error_species <- character()

# ---- Main loop over species ----
message(sprintf("\n========== Computing exposure for %d species ==========", J))

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

            # Get extent of binary map
            sp_extent <- terra::ext(binary_rast)

            # Count range cells
            range_cells <- terra::global(binary_rast, "sum", na.rm = TRUE)[[1]]
            message(sprintf("  Range cells: %d", range_cells))

            if (range_cells == 0 || is.na(range_cells)) {
                message("  No presence cells in binary map, skipping")
                next
            }

            # ---- Initialize results for this species ----
            sp_results <- list(
                species_name = sp_name,
                range_cells = as.integer(range_cells)
            )

            # ---- Process each climate variable ----
            for (i in seq_along(clim_names)) {
                var_name <- clim_names[i]

                # Crop climate layers to species extent
                baseline_crop <- terra::crop(r_clim_baseline[[i]], sp_extent, snap = "out")
                future_crop <- terra::crop(r_clim_future[[i]], sp_extent, snap = "out")

                # Resample to match binary raster
                baseline_resamp <- terra::resample(baseline_crop, binary_rast, method = "bilinear")
                future_resamp <- terra::resample(future_crop, binary_rast, method = "bilinear")

                # Compute climate change (difference)
                clim_change <- future_resamp - baseline_resamp

                # Mask to species range (where binary == 1)
                baseline_masked <- terra::mask(baseline_resamp, binary_rast, maskvalues = 0)
                future_masked <- terra::mask(future_resamp, binary_rast, maskvalues = 0)
                clim_change_masked <- terra::mask(clim_change, binary_rast, maskvalues = 0)

                # Compute statistics within range
                baseline_stats <- terra::global(baseline_masked, c("mean", "sd"), na.rm = TRUE)
                future_stats <- terra::global(future_masked, c("mean", "sd"), na.rm = TRUE)
                change_stats <- terra::global(clim_change_masked, c("mean", "sd"), na.rm = TRUE)

                # Store results
                sp_results[[paste0(var_name, "_mean_baseline")]] <- baseline_stats$mean
                sp_results[[paste0(var_name, "_mean_future")]] <- future_stats$mean
                sp_results[[paste0(var_name, "_mean_change")]] <- change_stats$mean
                sp_results[[paste0(var_name, "_sd_change")]] <- change_stats$sd

                # Clean up
                rm(
                    baseline_crop, future_crop, baseline_resamp, future_resamp,
                    clim_change, baseline_masked, future_masked, clim_change_masked
                )
            }

            # ---- Write exposure rasters for this species ----
            message("  Writing exposure rasters...")

            for (i in seq_along(clim_names)) {
                var_name <- clim_names[i]

                # Crop and resample again for output
                baseline_crop <- terra::crop(r_clim_baseline[[i]], sp_extent, snap = "out")
                future_crop <- terra::crop(r_clim_future[[i]], sp_extent, snap = "out")
                baseline_resamp <- terra::resample(baseline_crop, binary_rast, method = "bilinear")
                future_resamp <- terra::resample(future_crop, binary_rast, method = "bilinear")

                # Compute change and mask
                clim_change <- future_resamp - baseline_resamp
                clim_change_masked <- terra::mask(clim_change, binary_rast, maskvalues = 0)

                # Write output
                out_file <- file.path(
                    exposure_dir,
                    paste0(base_prefix, "_", sp_name, "_exposure_", var_name, ".tif")
                )
                terra::writeRaster(clim_change_masked, out_file, overwrite = TRUE)

                rm(
                    baseline_crop, future_crop, baseline_resamp, future_resamp,
                    clim_change, clim_change_masked
                )
            }
            message("  Wrote 4 exposure rasters")

            # ---- Add to results ----
            exposure_results <- rbind(exposure_results, as.data.frame(sp_results))

            # Print summary
            message(sprintf(
                "  Bio1 (temp) change: %.2f ± %.2f",
                sp_results$bio1_temp_mean_change,
                sp_results$bio1_temp_sd_change
            ))
            message(sprintf(
                "  Bio13 (precip) change: %.2f ± %.2f",
                sp_results$bio13_precipWet_mean_change,
                sp_results$bio13_precipWet_sd_change
            ))

            success_species <- c(success_species, sp_name)

            # Clean up
            rm(binary_rast)
            cleanup_terra_temps(temp_dir)
        },
        error = function(e) {
            message(sprintf("  ERROR: %s", e$message))
            error_species <- c(error_species, sp_name)
            cleanup_terra_temps(temp_dir)
        }
    )
}

# ---- Compute summary metrics ----
message("\n========== Computing summary exposure metrics ==========")

if (nrow(exposure_results) > 0) {
    # Add composite exposure score (normalized sum of absolute changes)
    exposure_results$exposure_score <- NA

    if (nrow(exposure_results) > 1) {
        # Z-score normalize each change variable
        bio1_z <- scale(abs(exposure_results$bio1_temp_mean_change))
        bio4_z <- scale(abs(exposure_results$bio4_tempSeas_mean_change))
        bio13_z <- scale(abs(exposure_results$bio13_precipWet_mean_change))
        bio15_z <- scale(abs(exposure_results$bio15_precipSeas_mean_change))

        # Handle NAs
        bio1_z[is.na(bio1_z)] <- 0
        bio4_z[is.na(bio4_z)] <- 0
        bio13_z[is.na(bio13_z)] <- 0
        bio15_z[is.na(bio15_z)] <- 0

        # Composite score (mean of z-scores)
        exposure_results$exposure_score <- rowMeans(cbind(bio1_z, bio4_z, bio13_z, bio15_z))
    }
}

# ---- Save results ----
message("\n========== Saving results ==========")

# Save CSV
results_csv <- file.path(exposure_dir, paste0(base_prefix, "_exposure_summary.csv"))
write.csv(exposure_results, results_csv, row.names = FALSE)
message(sprintf("Saved exposure summary: %s", results_csv))

# Save Rdata
results_rdata <- file.path(exposure_dir, paste0(base_prefix, "_exposure_summary.Rdata"))
save(exposure_results, file = results_rdata)
message(sprintf("Saved exposure Rdata: %s", results_rdata))

# ---- Print summary ----
message("\n========== Exposure Summary ==========")
message(sprintf("Species processed: %d", length(success_species)))
message(sprintf("Missing binary maps: %d", length(missing_binary)))
message(sprintf("Errors: %d", length(error_species)))

if (nrow(exposure_results) > 0) {
    message(sprintf(
        "\nMean temperature change (bio1): %.2f ± %.2f",
        mean(exposure_results$bio1_temp_mean_change, na.rm = TRUE),
        sd(exposure_results$bio1_temp_mean_change, na.rm = TRUE)
    ))
    message(sprintf(
        "Mean precip change (bio13): %.2f ± %.2f",
        mean(exposure_results$bio13_precipWet_mean_change, na.rm = TRUE),
        sd(exposure_results$bio13_precipWet_mean_change, na.rm = TRUE)
    ))

    # Top 5 most exposed species
    if (nrow(exposure_results) >= 5 && !all(is.na(exposure_results$exposure_score))) {
        top_exposed <- exposure_results[order(-exposure_results$exposure_score), ][1:5, ]
        message("\nTop 5 most exposed species:")
        for (i in 1:5) {
            message(sprintf(
                "  %d. %s (score: %.2f, temp change: %.2f)",
                i, top_exposed$species_name[i],
                top_exposed$exposure_score[i],
                top_exposed$bio1_temp_mean_change[i]
            ))
        }
    }
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

message(sprintf("\n========== Done! Exposure results saved to: %s ==========", exposure_dir))
