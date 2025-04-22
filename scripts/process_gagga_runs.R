#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readxl)
  library(writexl)
})

# ----------------- Utility Functions -----------------
get_option <- function(arg_name, default) {
  args <- commandArgs(trailingOnly = TRUE)
  m    <- grep(paste0("^--", arg_name, "="), args)
  if (length(m)) {
    val <- sub(paste0("^--", arg_name, "="), "", args[m])
    if (arg_name %in% c("top_percent", "n_cliques")) return(as.numeric(val))
    return(val)
  }
  default
}

show_help <- function() {
  cat(
"Usage: Rscript process_gagga_runs.R [options]\n\nOptions:\n",
"  --file_path=PATH             Aggregated output (default: 'all_gagga_data.xlsx')\n",
"  --top_percent=NUM            Top % cutoff for coarse filter (required)\n",
"  --n_cliques=NUM              Number to pass as k when thinning (required)\n",
"  --debug                      Limit to first 200 files (flag)\n",
"  --o=PATH                     Final hingesets output (default: 'hingesets.xlsx')\n",
"  --data_directory=DIR         Directory of gagga_*.xlsx (default: 'out/gagga')\n",
"  --external_overhangs=STRING  External overhangs (default: 'AATG, GCTT')\n",
"  --no-thin                    Skip thinning in Python (flag)\n",
"  --help                       Show help\n"
  )
}

# ---- parse_args ----
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) { show_help(); q(status = 0) }

file_path          <- get_option("file_path",         "all_gagga_data.xlsx")
top_percent        <- get_option("top_percent",       NA)
n_cliques          <- get_option("n_cliques",         NA)
debug_mode         <- "--debug"         %in% args
output_file        <- get_option("o",                 "hingesets.xlsx")
data_directory     <- get_option("data_directory",    "out/gagga")
external_overhangs <- get_option("external_overhangs","AATG, GCTT")
do_thin_flag       <- !("--no-thin"     %in% args)

if (is.na(top_percent) || is.na(n_cliques)) {
  stop("Error: --top_percent and --n_cliques are required. Use --help.")
}

# ----------------- AGGREGATION STEP -----------------
setwd(data_directory)
all_files <- list.files(pattern="^gagga_.*\\.xlsx$")
if (debug_mode) {
  all_files <- head(all_files, 200)
  message("[DEBUG] Using ", length(all_files), " files")
}

read_file <- function(fn) {
  sheets <- readxl::excel_sheets(fn)
  sh     <- if ("Results" %in% sheets) "Results" else sheets[[1]]
  readxl::read_excel(fn, sheet = sh, col_names = TRUE)
}

all_data_list <- lapply(all_files, function(fn) {
  tryCatch({
    dat <- read_file(fn)
    # unify raw_fidelity
    if ("fidelity" %in% names(dat)) {
      dat <- dat %>%
        rename(raw_fidelity = fidelity) %>%
        mutate(raw_fidelity = as.numeric(str_remove_all(raw_fidelity, "[^0-9\\.]")))
    } else if ("Score" %in% names(dat)) {
      dat <- dat %>%
        mutate(raw_fidelity = as.numeric(Score) * 100)
    } else {
      stop("No fidelity column in ", fn)
    }
    # unify setA
    if ("hf_oh_set" %in% names(dat)) {
      dat <- rename(dat, setA = hf_oh_set)
    } else if ("Set" %in% names(dat)) {
      dat <- rename(dat, setA = Set)
    }
    # clean setA, compute size
    dat %>%
      mutate(
        setA = str_remove_all(setA, "\\[|\\]|'"),
        setA = str_squish(str_replace_all(setA, ",", ", "))
      ) %>%
      mutate(
        set_size = lengths(str_split(setA, ",\\s*")),
        file     = fn
      )
  }, error = function(e) {
    message("[ERROR] ", fn, ": ", e$message)
    NULL
  })
})

all_data <- bind_rows(all_data_list) %>%
  filter(!is.na(raw_fidelity))

message("Total aggregated rows: ", nrow(all_data))
message("Writing aggregated data to: ", file_path)
writexl::write_xlsx(all_data, file_path)

# ----------------- COARSE FILTER -----------------
coarse_p <- 2 * top_percent / 100
message(sprintf("Coarse filtering: top %.1f%% per set_size...", 2 * top_percent))
working <- all_data %>%
  group_by(set_size) %>%
  filter(
    raw_fidelity >= quantile(
      raw_fidelity,
      probs = 1 - coarse_p,
      na.rm = TRUE
    )
  ) %>%
  ungroup()
message("Rows after coarse filter: ", nrow(working))

# ----------------- PALINDROME FILTER -----------------
pal_list <- c('ATAT','TATA','CGCG','GCGC','GATC','CTAG','GTAC','CATG',
              'AATT','TTAA','CCGG','GGCC','AGCT','TCGA','ACGT','TGCA')
message("Filtering palindromes...")
working <- working %>%
  rowwise() %>%
  filter(!any(str_split(setA, ",\\s*")[[1]] %in% pal_list)) %>%
  ungroup()
message("Rows after palindrome filter: ", nrow(working))

# ----------------- WRITE & CALL PYTHON RESCORER -----------------
thin_file <- paste0(tools::file_path_sans_ext(file_path), "_top",
                    top_percent, "pct.xlsx")
message("Writing pre‑thin data to: ", thin_file)
writexl::write_xlsx(working, thin_file)

py_cmd <- paste(
  "python ../../scripts/rescore.py",
  paste0("--i=", thin_file),
  paste0("--k=", n_cliques),
  paste0("--external='", external_overhangs, "'"),
  if (!do_thin_flag) "--no-thin" else "",
  sep = " "
)
message("Running Python rescoring/thinning:\n  ", py_cmd)
system(py_cmd)

message("Done. Outputs:\n  • top_rescored.xlsx\n  • hingests.xlsx")
