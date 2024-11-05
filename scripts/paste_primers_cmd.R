# The main function of this script is to take the chiseled/hinged oligo 
# sequences and append indexing primers to them so that the fragments can 
# be amplified in pools prior to Golden Gate assembly. The script takes as 
# input a file containing the fragmented sequences and a file of indexing 
# primers. Depending on the provided command-line options, it can handle 
# various configurations, including multi-step assemblies and different 
# types of sequences.
#
# The script:
# - Assigns indexing primers to each sequence fragment.
# - Splits larger genes into fragment sets, each amplified by different 
#   primers.
# - Supports two-step assembly processes for more complex sequences.
# - Handles potential sequence duplicates by renaming them and ensuring 
#   unique identifiers/barcodes.
# - Outputs the final oligo sequences in FASTA format ready for ordering.
#
# Command-Line Arguments:
# -i, --input_pops     : Path to the input pops file (Excel), 
#                        default="pops.xlsx".
# -p, --input_primers  : Path to the input primers file (CSV), 
#                        default="data/10K_primers_renamed.csv".
# -o, --output_file    : Prefix for output file names, 
#                        default="popped_output".
# -n, --primer_index   : Starting index for assigning primers, 
#                        default=1.
# -l, --pcr_5p_cut     : Sequence to add to the 5' end for cutting, 
#                        default="CGTCTCA".
# -r, --pcr_3p_cut     : Sequence to add to the 3' end for cutting, 
#                        default="AGAGACG".
# -m, --max_fragments  : Max fragments per PCR pool, 
#                        default=7.
# -t, --run_type       : Type of sequences ("cds" for FASTA or "gb" for 
#                        GenBank), default="cds".
# -e, --two_step       : Enables two-step assembly mode ("on" or "off"), 
#                        default="off".

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(optparse)
  library(openxlsx)
})

option_list <- list(
  make_option(c("-i", "--input_pops"), default = "pops.xlsx", 
              type = "character", 
              help = "Path to the input pops file [default = %default]."),
  make_option(c("-p", "--input_primers"), 
              default = "data/10K_primers_renamed.csv", type = "character", 
              help = "Path to the input primers file [default = %default]."),
  make_option(c("-o", "--output_file"), default = "popped_output", 
              type = "character", 
              help = "Name of the output file prefix [default = %default]."),
  make_option(c("-n", "--primer_index"), default = 1, type = "numeric", 
              help = "Where to start adding primers from [default = %default]."),
  make_option(c("-l", "--pcr_5p_cut"), default = "CGTCTCA", type = "character", 
              help = "5' sequence for cutting [default = %default]."),
  make_option(c("-r", "--pcr_3p_cut"), default = "AGAGACG", type = "character", 
              help = "3' sequence for cutting [default = %default]."),
  make_option(c("-m", "--max_fragments"), default = 7, type = "numeric", 
              help = "Max number of fragments to include per PCR pool [default = %default]."),
  make_option(c("-t", "--run_type"), default = "cds", type = "character", 
              help = "Run type - gb or cds (i.e., fasta)."),
  make_option(c("-e", "--two_step"), default = "off", type = "character", 
              help = "Two-step setting.")
)



# Parse the command-line options
opt <- parse_args(OptionParser(option_list = option_list))

# Load the pre-indexed fragment data
input_pops <- tryCatch({ 
    readxl::read_xlsx(opt$input_pops)
  },
  error = function(e) {
    cat("An error with the xls file read occurred:", 
        conditionMessage(e), "\n")
    NULL  # Return NULL or any other appropriate value if needed
  }
)

# Load the two-step fragment data if present
if (opt$two_step == "on") {
  two_step_data <- tryCatch({ 
      readxl::read_xlsx(opt$input_pops, sheet = "step_1_fragments")
    },
    error = function(e) {
      cat("An error with the two-step data read occurred:", 
          conditionMessage(e), "\n")
      NULL  # Return NULL or any other appropriate value if needed
    }
  )
  two_step_data <- two_step_data %>%
    select(c(base_id, second_step_assembly_seq = Chiseled_sequence))
}

# Split larger genes into sets of fragments amplified by different primers
# Set max_fragments to a high number if you want genes to always be amplified 
# by the same index primer set
input_pops <- input_pops %>%
  mutate(original_row = row_number()) %>%
  group_by(Seq_ID) %>%
  mutate(total_fragments = n()) %>%
  ungroup() %>%
  mutate(Subgroups = ceiling(total_fragments / opt$max_fragments)) %>%
  group_by(Seq_ID) %>%
  mutate(Subgroup_ID = (row_number() - 1) %/% (total_fragments / Subgroups) + 1) %>%
  ungroup() %>%
  mutate(New_Seq_ID = if_else(Subgroups > 1, 
                              paste(Seq_ID, Subgroup_ID, sep = "-"), 
                              as.character(Seq_ID))) %>%
  select(-Subgroups, -Subgroup_ID) %>%
  # Replace the original Seq_ID with New_Seq_ID
  mutate(base_id = Seq_ID) %>%
  mutate(Seq_ID = New_Seq_ID) %>%
  select(-New_Seq_ID) %>%
  select(-original_row)

primer_list  <- read.csv(opt$input_primers)
pcr_5p_cut  <- opt$pcr_5p_cut
pcr_3p_cut <- opt$pcr_3p_cut

# Initialize the dataframe to store the modified fragments
modified_fragments_df <- tibble()

# Group by Seq_ID and Fragment_n, add a row number within each group
input_pops <- input_pops %>%
  group_by(Seq_ID, Fragment_n) %>%
  mutate(row_num = row_number()) %>%
  ungroup()  # Ungroup here

# Identify duplicates
duplicates <- input_pops %>%
  group_by(Seq_ID, Fragment_n) %>%
  filter(n() > 1) %>%
  ungroup()  # Ungroup here

# If duplicates are found, print a warning
if (nrow(duplicates) > 0) {
  # Rename Seq_ID for duplicates
  input_pops <- input_pops %>%
    mutate(new_name = ifelse(Seq_ID %in% unique(duplicates$Seq_ID), 
                             paste(Seq_ID, row_num, sep = "_"), 
                             Seq_ID))
  
  # Extract the old and new names for the warning message
  name_changes <- input_pops %>%
    filter(Seq_ID %in% unique(duplicates$Seq_ID)) %>%
    select(Seq_ID, new_name)
  
  # Create the warning message
  cat(str_c("The following Seq_IDs are duplicated and were renamed.",
            paste0(name_changes$Seq_ID, " => ", name_changes$new_name, 
                   collapse = "\n"), sep = "\n"))
  
  # Update the Seq_ID column with the new names
  input_pops$Seq_ID <- input_pops$new_name
  
  # Drop the row_num and new_name columns
  input_pops <- select(input_pops, -row_num, -new_name)
}

# Get the unique Seq_IDs
unique_seq_ids <- unique(input_pops$Seq_ID)

# Iterate over unique Seq_IDs and corresponding primers
for (i in 1:length(unique_seq_ids)) {
  # Get the fragments for the specific Seq_ID
  fragments <- filter(input_pops, Seq_ID == unique_seq_ids[i])
  
  # Adjust the index for selecting the primer, considering the starting 
  # primer_index
  primer_selector <- (i + opt$primer_index - 2) %% nrow(primer_list) + 1
  
  if (primer_selector == 1 && i != 1) {
    warning("Primer index has exceeded the number of available primers. 
             Recycling from the start.")
  }
  
  # Use the adjusted primer_selector to select the primer
  F_seq <- primer_list$F_seq[primer_selector]
  R_seq <- primer_list$R_seq[primer_selector]
  R_seq_rc <- primer_list$R_seq_rc[primer_selector]
  primer_id <- primer_list$name[primer_selector]
  
  # Modify the fragments based on their position within each Seq_ID group
  fragments_with_additions <- fragments %>%
    group_by(Seq_ID) %>%
    mutate(
      oligo = case_when(
        # Modify the first two lines in this mutate if you need to have any 
        # special sequences on the 5' and 3' ends of the final assembled 
        # fragment -- currently they are identical primer + BsmBI site 
        # (default) added to each end
        row_number() == 1 ~ paste0(F_seq, pcr_5p_cut, Fragment, pcr_3p_cut, 
                                   R_seq_rc), # First fragment
        row_number() == n() ~ paste0(F_seq, pcr_5p_cut, Fragment, pcr_3p_cut, 
                                     R_seq_rc), # Last fragment
        TRUE ~ paste0(F_seq, pcr_5p_cut, Fragment, pcr_3p_cut, R_seq_rc) # 
        # Middle fragments
      ),
      f_primer = F_seq, 
      r_primer = R_seq,
      primer_index = primer_id,
      n_fragments = length(oligo)
    ) %>%
    ungroup()
  
  # Append the modified fragments to the dataframe
  modified_fragments_df <- bind_rows(modified_fragments_df, 
                                     fragments_with_additions)
}

modified_fragments_df <- modified_fragments_df %>% 
  mutate(id = str_c(Seq_ID, "_", Fragment_n), oligo, primer_index, 
         f_primer, r_primer, frag_length = nchar(oligo))

if (opt$run_type == "gb") {
  out <- modified_fragments_df %>% 
    select(accession, base_id, total_fragments, frag_id = id, oligo, 
           primer_index, f_primer, r_primer, n_fragments_in_set = n_fragments, 
           frag_length, left_oh = Left_Overhang, right_oh = Right_Overhang, 
           internal_ohs = Internal_overhangs, fidelity = Fidelity, 
           input_seq = Original_sequence, chiseled_seq = Chiseled_sequence, 
           Seq_ID)
}

if (opt$run_type == "cds") {
  out <- modified_fragments_df %>% 
    select(accession, base_id, total_fragments, frag_id = id, oligo, 
           primer_index, f_primer, r_primer, n_fragments_in_set = n_fragments, 
           frag_length, left_oh = Left_Overhang, right_oh = Right_Overhang, 
           internal_ohs = Internal_overhangs, fidelity = Fidelity, 
           input_seq = Original_sequence, chiseled_seq = Chiseled_sequence, 
           Seq_ID, input_cai = Original_CAI, chiseled_cai = Chiseled_CAI)
}

oligos_out <- out %>%
  mutate(frag_id = str_replace(frag_id, "(.*\\.\\d+_)(\\d+)", "\\1FRAG_\\2")) %>%
  transmute(oligos = str_c(">", frag_id, "_PRI_", primer_index, "\n", oligo))

primers_out <- out %>% 
  transmute(primers = str_c(">", primer_index, "_F", "\n", f_primer, "\n", ">", 
                            primer_index, "_R", "\n", r_primer)) %>% 
  distinct()

chisels_out <- out %>% 
  group_by(base_id) %>% 
  slice_head(n = 1) %>%
  transmute(chisels = str_c(">", base_id, " ", primer_index, "\n", 
                            chiseled_seq, "\n")) %>% 
  ungroup() %>%
  select(chisels) 

input_seqs_out <- out %>% 
  group_by(base_id) %>% 
  slice_head(n = 1) %>%
  transmute(input_seqs = str_c(">", base_id, "\n", input_seq, "\n")) %>% 
  ungroup() %>%
  select(input_seqs)

if (opt$two_step == "on") {
  out <- out %>% mutate(full_base_id = gsub("_L0_\\d+$", "", base_id))
  two_step_data <- two_step_data %>% 
    mutate(full_base_id = gsub("_L0_\\d+$", "", base_id))
} else {
  out <- out %>% mutate(full_base_id = base_id)
}

if (opt$run_type == "gb") {

  data_out <- out %>% 
    select(c(accession, pop_id = base_id, frag_id, primer_index, oligo, 
             left_oh, right_oh, frag_length))

  nano_out <- out %>% 
    group_by(base_id) %>%
    select(c(accession, pop_id = base_id, primer_index, total_fragments, 
             chiseled_seq)) %>%
    slice(1) %>%
    arrange(primer_index)

  seqs_out <- out %>% 
    group_by(base_id) %>%
    select(c(accession, pop_id = base_id, primer_index, f_primer, r_primer, 
             total_fragments, n_fragments_in_set, internal_ohs, fidelity, input_seq, 
             chiseled_seq)) %>%
    slice(1) %>%
    arrange(primer_index)  

  if (opt$two_step == "on") {
    seqs_out <- out %>%
      left_join(two_step_data, by = "full_base_id", 
                relationship = "many-to-many") %>%
      mutate(base_id = base_id.x) %>%
      group_by(base_id) %>%
      mutate(total_fragments = n()) %>%    
      select(c(accession, pop_id = base_id, primer_index, total_fragments, 
               n_fragments_in_set, internal_ohs, fidelity, f_primer, r_primer, input_seq, 
               chiseled_seq, second_step_assembly_seq)) %>%
      slice(1) %>%
      arrange(primer_index)
  }
}

if (opt$run_type == "cds") {

  data_out <- out %>% 
    select(c(accession, pop_id = base_id, frag_id, primer_index, oligo, 
             left_oh, right_oh, frag_length))

  nano_out <- out %>% 
    group_by(base_id) %>%
    select(c(accession, pop_id = base_id, primer_index, total_fragments, 
             chiseled_seq)) %>%
    slice(1) %>%
    arrange(primer_index)

  seqs_out <- out %>% 
    group_by(base_id) %>%
    select(c(accession, pop_id = base_id, primer_index, f_primer, r_primer, 
             total_fragments, n_fragments_in_set, internal_ohs, fidelity, f_primer, r_primer, 
             input_cai, chiseled_cai, input_seq, chiseled_seq)) %>%
    slice(1) %>%
    arrange(primer_index)

  if (opt$two_step == "on") {
    seqs_out <- out %>%
      left_join(two_step_data, by = "full_base_id", 
                relationship = "many-to-many") %>%
      mutate(base_id = base_id.x) %>%
      group_by(base_id) %>%
      mutate(total_fragments = n()) %>%
      select(c(accession, pop_id = base_id, primer_index, total_fragments, 
               n_fragments_in_set, internal_ohs, fidelity, f_primer, r_primer, input_cai, 
               chiseled_cai, input_seq, chiseled_seq, 
               second_step_assembly_seq)) %>%
      slice(1) %>%
      arrange(primer_index)
  }
}

# Construct the file paths
oligos_file_path <- file.path(
  paste0(opt$output_file, "_oligo_pool.fasta"))
primers_file_path <- file.path(
  paste0(opt$output_file, "_index_primers_required.fasta"))
seqs_file_path <- file.path(
  paste0(opt$output_file, "_designed_seqs.fasta"))
nano_file_path <- file.path(
  paste0(opt$output_file, "_nanoseq.csv"))

write_csv2(oligos_out, oligos_file_path, col_names = FALSE, quote = "none")
write_csv2(primers_out, primers_file_path, col_names = FALSE, quote = "none")
write_csv2(chisels_out, seqs_file_path, col_names = FALSE, quote = "none")
write_csv2(nano_out, nano_file_path, col_names = TRUE, quote = "none")

excel_file_path <- file.path(paste0(opt$output_file, "_all_data.xlsx"))
wb <- loadWorkbook(excel_file_path)

# Uncomment this if you don't want the input data from pop included in the 
# final excel file that has all of the data
removeWorksheet(wb, sheet = "pre_indexed_data")
if (opt$two_step == "on") {
  removeWorksheet(wb, sheet = "step_1_fragments")
}

addWorksheet(wb, "fragment_data")
writeData(wb, sheet = "fragment_data", data_out)
addWorksheet(wb, "seq_data")
writeData(wb, sheet = "seq_data", seqs_out)

saveWorkbook(wb, excel_file_path, overwrite = TRUE)

cat(str_c("The barcodes used start at index ", opt$primer_index, 
          " and continue through to ", primer_selector, 
          "\nUse: --primer_index ", primer_selector + 1, 
          " to add new barcodes in your next run.\n\n"))

cat(str_c(n_distinct(modified_fragments_df$base_id), 
          " sequence(s) poppified.\n\nResults saved to:\n\n ", 
          paste0(opt$output_file, "_all_data.xlsx\n "), 
          paste0(opt$output_file, "_oligo_pool.fasta\n "), 
          paste0(opt$output_file, "_index_primers_required.fasta\n "))) 
