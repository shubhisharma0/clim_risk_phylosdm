#!/usr/bin/env Rscript

# ============================================================
# dsq_retry_builder.R
#
# Build a DSQ retry job-file by subsetting an existing job-file.
#
# Two modes:
#   (A) From scheduler output (recommended):
#       Rscript scripts/handlers/dsq_retry_builder.R \
#         --jobid 64658193 \
#         --states FAILED,OUT_OF_MEMORY \
#         --in  jobs/job_array_make_datasets_for_03-gen_data.txt \
#         --out jobs/job_array_make_datasets_for_03-gen_data_retry.txt
#
#   (B) From an explicit 0-based index file:
#       Rscript scripts/handlers/dsq_retry_builder.R \
#         --indices checks/dsq_reruns/retry_indices.txt \
#         --in  jobs/job_array_make_datasets_for_03-gen_data.txt \
#         --out jobs/job_array_make_datasets_for_03-gen_data_retry.txt
#
# Notes:
# - DSQ indices are 0-based (as reported by dsqa).
# - DSQ index i corresponds to line (i + 1) in the job file.
# ============================================================

suppressWarnings(suppressPackageStartupMessages({
  # No external packages required
}))

# -----------------------------
# Utilities
# -----------------------------
die <- function(...) {
  msg <- paste0(..., collapse = "")
  stop(msg, call. = FALSE)
}

trim <- function(x) trimws(x)

split_csv <- function(x) {
  if (is.null(x) || !nzchar(trim(x))) {
    return(character())
  }
  parts <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
  parts <- trim(parts)
  parts[nzchar(parts)]
}

read_nonempty_lines <- function(path) {
  if (is.null(path)) {
    return(character())
  }
  if (!file.exists(path)) die("ERROR: File not found: ", path, "\n")
  x <- readLines(path, warn = FALSE)
  x <- gsub("\r", "", x) # strip Windows CR if present
  x <- x[nzchar(trim(x))] # drop blank lines
  x
}

parse_indices_vector <- function(idx_raw, context_label = "indices") {
  if (length(idx_raw) == 0) {
    return(integer())
  }
  idx <- suppressWarnings(as.integer(idx_raw))
  if (any(is.na(idx))) {
    bad <- idx_raw[is.na(idx)]
    die(
      "ERROR: Some lines in ", context_label, " are not integers:\n",
      paste(bad, collapse = "\n"),
      "\n"
    )
  }
  sort(unique(idx))
}

# Extract 0-based indices for specified states from dsqa output text
extract_indices_from_dsqa <- function(dsqa_text, states_keep) {
  # Expect lines like:
  # OUT_OF_MEMORY      21    8,21-23,29,...
  # FAILED              2     53,215
  lines <- unlist(strsplit(dsqa_text, "\n", fixed = TRUE), use.names = FALSE)
  lines <- trim(lines)
  lines <- lines[nzchar(lines)]

  # Identify the "State Summary" section if present; otherwise parse all lines.
  # We parse any line that starts with a state token and has an indices field.
  parsed <- list()

  for (ln in lines) {
    # Skip header-ish lines
    if (grepl("^State\\s+Num_Jobs\\s+Indices", ln)) next
    if (grepl("^-{2,}$", ln)) next
    if (grepl("^State Summary", ln)) next

    # Split on whitespace; expect >= 3 tokens: STATE, NUM, INDICES
    tok <- strsplit(ln, "\\s+")[[1]]
    if (length(tok) < 3) next

    state <- tok[1]
    if (!(state %in% states_keep)) next

    # Indices may contain commas and hyphens; may be token 3 or later depending on spacing.
    # dsqa prints: STATE <spaces> NUM <spaces> INDICES
    # Combine everything from token 3 onward (safer).
    idx_field <- paste(tok[3:length(tok)], collapse = " ")
    idx_field <- gsub("\\s+", "", idx_field) # remove spaces inside indices field

    if (!nzchar(idx_field)) next

    parsed[[state]] <- idx_field
  }

  if (length(parsed) == 0) {
    return(integer())
  }

  # Expand an indices field like: "0-7,9-20,24,30"
  expand_idx_field <- function(s) {
    if (!nzchar(s)) {
      return(integer())
    }
    parts <- unlist(strsplit(s, ",", fixed = TRUE), use.names = FALSE)
    out <- integer()

    for (p in parts) {
      if (!nzchar(p)) next
      if (grepl("-", p, fixed = TRUE)) {
        ab <- unlist(strsplit(p, "-", fixed = TRUE), use.names = FALSE)
        if (length(ab) != 2) next
        a <- suppressWarnings(as.integer(ab[1]))
        b <- suppressWarnings(as.integer(ab[2]))
        if (is.na(a) || is.na(b)) next
        if (b < a) {
          # Rare, but handle gracefully
          out <- c(out, seq.int(b, a))
        } else {
          out <- c(out, seq.int(a, b))
        }
      } else {
        v <- suppressWarnings(as.integer(p))
        if (!is.na(v)) out <- c(out, v)
      }
    }
    sort(unique(out))
  }

  all_fields <- unname(unlist(parsed, use.names = FALSE))
  idx <- unlist(lapply(all_fields, expand_idx_field), use.names = FALSE)
  sort(unique(as.integer(idx)))
}

run_dsqa <- function(jobid) {
  if (is.null(jobid) || !nzchar(trim(jobid))) {
    die("ERROR: --jobid was not provided.\n")
  }
  cmd <- sprintf("dsqa -j %s 2>&1", shQuote(jobid))
  out <- tryCatch(
    system(cmd, intern = TRUE),
    error = function(e) character()
  )
  if (length(out) == 0) {
    return("")
  }
  result <- paste(out, collapse = "\n")
  result
}

usage <- function() {
  cat(
    "\n",
    "dsq_retry_builder.R\n",
    "\n",
    "Build a DSQ retry job-file by subsetting an existing job-file.\n",
    "\n",
    "Required:\n",
    "  --in   <path>   Input job file (one command per line)\n",
    "  --out  <path>   Output retry job file\n",
    "\n",
    "Choose ONE source of indices:\n",
    "  (A) From dsqa:\n",
    "      --jobid  <jobid>\n",
    "      --states <comma-separated states, default: FAILED,OUT_OF_MEMORY>\n",
    "\n",
    "  (B) From explicit index file:\n",
    "      --indices <path to file with 0-based indices, one per line>\n",
    "\n",
    "Optional:\n",
    "  --dry_run   Print summary but do not write output\n",
    "  --verbose   Print dsqa output and parsing details\n",
    "\n",
    "Examples:\n",
    "  Rscript scripts/handlers/dsq_retry_builder.R \\\n",
    "    --jobid 64658193 --states FAILED,OUT_OF_MEMORY \\\n",
    "    --in jobs/job_array_make_datasets.txt \\\n",
    "    --out jobs/job_array_make_datasets_retry.txt\n",
    "\n",
    "  Rscript scripts/handlers/dsq_retry_builder.R \\\n",
    "    --indices checks/dsq_reruns/retry_indices.txt \\\n",
    "    --in jobs/job_array_make_datasets.txt \\\n",
    "    --out jobs/job_array_make_datasets_retry.txt\n",
    "\n",
    sep = ""
  )
}

# -----------------------------
# Parse args (simple flag parser)
# -----------------------------
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) == 0 || any(argv %in% c("-h", "--help"))) {
  usage()
  quit(status = 0)
}

get_flag <- function(flag) flag %in% argv

get_value <- function(flag) {
  w <- which(argv == flag)
  if (length(w) == 0) {
    return(NULL)
  }
  if (w[1] == length(argv)) die("ERROR: Missing value after ", flag, "\n")
  argv[w[1] + 1]
}

in_path <- get_value("--in")
out_path <- get_value("--out")
jobid <- get_value("--jobid")
states_str <- get_value("--states")
idx_path <- get_value("--indices")
dry_run <- get_flag("--dry_run")
verbose <- get_flag("--verbose")

if (is.null(in_path) || is.null(out_path)) {
  usage()
  die("ERROR: --in and --out are required.\n")
}

if (!file.exists(in_path)) {
  die("ERROR: Input job file not found: ", in_path, "\n")
}

# Enforce exactly one index source
has_jobid <- !is.null(jobid)
has_idxfile <- !is.null(idx_path)

if (has_jobid && has_idxfile) {
  die("ERROR: Provide only one of --jobid or --indices (not both).\n")
}
if (!has_jobid && !has_idxfile) {
  die("ERROR: You must provide either --jobid (dsqa mode) or --indices (index-file mode).\n")
}

# Default states
if (is.null(states_str) || !nzchar(trim(states_str))) {
  states_keep <- c("FAILED", "OUT_OF_MEMORY")
} else {
  states_keep <- split_csv(states_str)
  if (length(states_keep) == 0) die("ERROR: --states was provided but parsed to an empty list.\n")
}

# -----------------------------
# Read job file
# -----------------------------
job_lines <- readLines(in_path, warn = FALSE)
job_lines <- gsub("\r", "", job_lines)
n_jobs <- length(job_lines)

if (n_jobs == 0) {
  die("ERROR: Input job file has 0 lines: ", in_path, "\n")
}

# -----------------------------
# Determine indices to rerun
# -----------------------------
idx <- integer()
dsqa_text <- NULL

if (has_idxfile) {
  idx_raw <- read_nonempty_lines(idx_path)
  idx <- parse_indices_vector(idx_raw, context_label = "--indices file")
} else {
  dsqa_text <- run_dsqa(jobid)
  if (is.null(dsqa_text) || !nzchar(trim(dsqa_text))) {
    die("ERROR: dsqa returned no output for jobid=", jobid, "\n")
  }
  idx <- extract_indices_from_dsqa(dsqa_text, states_keep)

  if (verbose) {
    cat("---- dsqa raw output ----\n")
    cat(dsqa_text, "\n")
    cat("-------------------------\n")
  }
}

if (length(idx) == 0) {
  cat("---- DSQ retry rebuild summary ----\n")
  if (has_jobid) {
    cat("JobID:", jobid, "\n")
    cat("States:", paste(states_keep, collapse = ","), "\n")
  } else {
    cat("Indices file:", idx_path, "\n")
  }
  cat("Job file lines:", n_jobs, "\n")
  cat("Selected indices: 0\n")
  cat("Nothing to do: no matching indices found.\n")
  quit(status = 0)
}

# -----------------------------
# Validate bounds
# -----------------------------
oob <- idx[idx < 0 | idx >= n_jobs]
if (length(oob) > 0) {
  die(
    "ERROR: Some indices are out of bounds for job file of length ", n_jobs, ":\n",
    paste(oob, collapse = ", "),
    "\n"
  )
}

# -----------------------------
# Subset job lines (DSQ i -> line i+1)
# -----------------------------
retry_lines <- job_lines[idx + 1]

# -----------------------------
# Write output
# -----------------------------
if (!dry_run) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(retry_lines, out_path)
}

# -----------------------------
# Verification (strict)
# -----------------------------
# This is strictly correct because we directly subset by line-number.
# Still, warn about duplicates because they can confuse humans.
dup_ct <- sum(duplicated(job_lines))
dup_retry_ct <- sum(duplicated(retry_lines))

# -----------------------------
# Summary
# -----------------------------
cat("---- DSQ retry rebuild summary ----\n")
if (has_jobid) {
  cat("JobID:", jobid, "\n")
  cat("States:", paste(states_keep, collapse = ","), "\n")
} else {
  cat("Indices file:", idx_path, "\n")
}
cat("Job file lines:", n_jobs, "\n")
cat("Selected indices:", length(idx), "\n")
cat("Wrote retry job lines:", length(retry_lines), "\n")
cat("Output file:", out_path, "\n")
cat("Verification: OK — output is exactly job_lines[idx + 1].\n")

if (dup_ct > 0) {
  cat("NOTE: Input job file contains duplicated lines (count duplicates:", dup_ct, ").\n")
}
if (dup_retry_ct > 0) {
  cat("NOTE: Retry job file contains duplicated lines (count duplicates:", dup_retry_ct, ").\n")
}

if (verbose) {
  cat("Indices:\n")
  cat(paste(idx, collapse = ", "), "\n")
}

if (dry_run) {
  cat("Dry run: did not write output.\n")
}

# ------------------------------------------------------------
# Example DSQ submission:
# dsq --job-file jobs/job_array_make_datasets_retry.txt \
#     --batch-file jobs/dsq_make_datasets_retry.sbatch \
#     --mem=64G -t 06:00:00 -p day
# ------------------------------------------------------------
