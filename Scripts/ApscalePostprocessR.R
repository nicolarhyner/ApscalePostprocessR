## =============================================================================
## Create and filter ESV tables from apscale metabarcoding pipelines
## using assigned taxonomies from apscale_blast or BOLDigger3
##
## Tested with: 
## - apscale = 4.3.0
## - apscale_blast = 1.3.2
## - BOLDigger3 = 2.2.0
##
## Author: Nicola Rhyner
## Year: 2026
## nicola.rhyner@zhaw.ch
## https://github.com/nicolarhyner
##
## CHECK BEFORE RUNNING
## - Set paths to input files
## - Dashes will be changed to underscores, make sure this is ok for your data
## - Define parameters and variables
## - Specify negative sample columns for blank correction
## - Inspect controls in the exported *.raw table
## =============================================================================


# load libraries ---------------------------------------------------------------
library(dplyr)
library(readxl)
library(writexl)
library(stringr)

# Create required folder structure --------------------------------
project_dir <- getwd()

subdirs <- c("Intermediate", "Output", "Scripts")

for (dir in subdirs) {
  if (!dir.exists(file.path(project_dir, dir))) {
    dir.create(file.path(project_dir, dir))
    print(paste("Required subdirectory", dir, "did not exist and was created"))
  } else {
    print(paste("Required subdirectory", dir, "already exists"))
  }
}


# read input files --------------------------------
abundance_table <- read_excel("")
#taxonomy_table_apscale_blast <- read_excel("") 
taxonomy_table_boldigger <- read_excel("")


# Define columns in the taxonomy tables that are not samples. Will be used to select sample columns later on
tax_cols_apscale_blast <- c("unique ID", "Kingdom", "Phylum", "Class", 
                      "Order", "Family", "Genus", "Species", "Similarity", "evalue", "Flag", 
                      "Ambiguous taxa", "Status")

tax_cols_boldigger <- c("id", "phylum", "class", "order",
                        "family", "genus", "species", "pct_identity", "status", "records", "records_ratio", "selected_level",
                        "BIN","flags")

# --- Configuration variables -----------------------------------------------------
# Project & labeling
project_name      <- "test" # give your output tables a project name
tax_assignment    <- "boldigger"    # select tax assignment method "apscale_blast" or "boldigger"
outputfile_suffix <- "BOLD"     # Can be i.e. the used database (Midori/BOLD)

# [Optional] Sample name prefix or suffix to remove 
samplename_prefix <- "no"
samplename_suffix <- "no"

# Filter parameters
tag_switch    <- 0.000   # proportion (0–1. i.e. 1% = 0.01) of reads to remove from each ESV
min_abundance <- 0.000  # proporiton (0–1) of reads to remove from each sample

# Tax assignment columns that you do not need in the final tables
unwanted_cols_apscale_blast <- c("evalue", "Status")
unwanted_cols_boldigger     <- c("status", "records", "records_ratio")

# Create a vector with your negative controls for the blank correction  --------------------------------
# Use sample names after any prefix/suffix removal
# Uncomment and adapt when used:
negatives <- c("NIC27_27", "NIC28_28", "NIC29_29")


# --- 0: Load and harmonize taxonomy table ------------------------------------------
if (tax_assignment == "apscale_blast") {
  
  taxonomy_table <- taxonomy_table_apscale_blast %>%
    # first select all relevant columns (before removing anything)
    select(all_of(tax_cols_apscale_blast)) %>%
    # remove unwanted columns
    select(-all_of(unwanted_cols_apscale_blast)) %>%
    # rename columns after selection
    rename(
      hash = `unique ID`,
      Ambiguous_taxa = `Ambiguous taxa`
    )
  
} else if (tax_assignment == "boldigger") {
  
  taxonomy_table <- taxonomy_table_boldigger %>%
    # first select all relevant columns
    select(all_of(tax_cols_boldigger)) %>%
    # remove unwanted columns
    select(-all_of(unwanted_cols_boldigger)) %>%
    # rename id column
    rename(hash = id)
  
} else {
  stop("tax_assignment must be 'apscale_blast' or 'boldigger'")
}

# overwrite metadata_cols to reflect renamed columns
metadata_cols <- colnames(taxonomy_table)

# Merge the files by ESV has to a taxonomy abundance table 
db <- merge(taxonomy_table, abundance_table, by = "hash", all = TRUE)

# remove unwanted prefix
if (samplename_prefix != "no") {
  colnames(db) <- gsub(samplename_prefix, "", colnames(db))
}

# remove unwanted suffix
if (samplename_suffix != "no") {
  colnames(db) <- gsub(samplename_suffix, "", colnames(db))
}

# Change all sequences to uppercase 
db$sequence <- toupper(db$sequence)

# reorder columns - I like the sequence at the end --------------------------------
db <- db %>%
   relocate(c(sequence), .after = last_col())


# --- 0: Identify sample columns automatically --------------------------------------
# substitute every dash to underscore to prevent excel from reading colnames as dates or numbers
colnames(db) <- gsub("-", "_", colnames(db))

# Sample columns = all columns that are NOT metadata
sample_cols <- setdiff(colnames(db), metadata_cols)

# Make sure sequence stays metadata, not a sample
sample_cols <- sample_cols[sample_cols != "sequence"]

# Create to to dataframes for samples and metadata for more convenience
sample_df <- db %>% select(all_of(sample_cols))
metadata_df <- db %>%select(all_of(metadata_cols))

# Create a column for the total number of reads per ESV
db$total_reads <- rowSums(sample_df, na.rm = TRUE)

# Sort hits in a decreasing way by first best identity and then total reads, order total_reads columns
if (tax_assignment == "apscale_blast") {
  db <- arrange(db, desc(Similarity), desc(total_reads))
} else if (tax_assignment == "boldigger") {
  db <- arrange(db, desc(pct_identity), desc(total_reads))
}

# Relocate the total reads column before the sequence
db <- db %>%
  relocate("total_reads", .before = "sequence")


# Export first raw ESV table --------------------------------
# Export raw taxonomy abundance table to the intermediate folder
subfolder <- "Intermediate"
raw <- paste0(project_name, "_table_raw", "_", outputfile_suffix, ".csv")
filepath <- file.path(subfolder,raw)
write.csv2(db, filepath, row.names=F) 


#  Once we have our raw OTU/ESV table we can start to curate/filter ---------------------
#  You can skip step 2 & 4 by adjusting setting min_abundance and tag_switch in the beginning to 0

# 2: Correct tag switching/index hopping in your samples
# 3: Perform a blank correction
# 4: Minimal abundance filtering of each sample 
# 5: [Optional] collapse identical species ESVs 


# --- 2: Tag-switch / Index-hopping correction -----------------------------------

# Applies tag-switch correction by subtracting OTU/ESV wise 'tag_switch' × total_reads 
# From the abundance. Negative values become zero; all values are rounded.

db.corr <- db   # copy raw db for these manipulations

if (tag_switch > 0) {
  
  message("Applying tag-switch correction (fraction = ", tag_switch, ")")
  
  # For each sample column: subtract proportionally to total reads
  for (sc in sample_cols) {
    db.corr[[sc]] <- db.corr[[sc]] - (tag_switch * db.corr$total_reads)
    db.corr[[sc]] <- pmax(db.corr[[sc]], 0)   # negative values to 0
    db.corr[[sc]] <- round(db.corr[[sc]])     # round values
  }
  
  # Calculate corrected total read sum per ESV
  db.corr$corrected_reads <- rowSums(db.corr[, sample_cols])
  
  # Adjust column order
  db.corr <- db.corr %>%
    relocate(corrected_reads, total_reads, .before = sequence)
  
  # Remove ESVs with no reads
  removed <- sum(db.corr$corrected_reads == 0)
  message(removed, " ESVs removed after tag-switch correction")
  
  db.corr <- db.corr %>% filter(corrected_reads > 0)
  
} else {
  
  message("Tag-switch correction skipped (tag_switch = 0)")
  
  # Just calculate read sums, no modification of abundances
  db.corr <- db %>%
    mutate(corrected_reads = rowSums(across(all_of(sample_cols))))
}

# Export
subfolder <- "Intermediate"
tagswitch <- paste0(project_name, "_table_tagswitch_", tag_switch, "_", outputfile_suffix, ".csv")
write.csv2(db.corr, file.path(subfolder, tagswitch), row.names = FALSE)


# --- 3: Minimal abundance filtering per sample ----------------------------------

# This step filters based on sample-level abundance thresholds using a 
# relative minimum abundance (proportion of total sample reads) 
# Be aware of any minimum abundance filter in apscale!
  
db.corr2 <- db.corr   # copy db.corr for these manipulations

if (min_abundance > 0) {
  
  message("Applying minimal relative abundance filter = ", min_abundance)
  
  for (sc in sample_cols) {
    db.corr2[[sc]] <- db.corr2[[sc]] - (sum(db.corr2[[sc]]) * min_abundance)
    db.corr2[[sc]] <- pmax(db.corr2[[sc]], 0)
    db.corr2[[sc]] <- round(db.corr2[[sc]])
  }
  
  db.corr2$corrected_reads <- rowSums(db.corr2[, sample_cols])
  
  message(sum(db.corr2$corrected_reads == 0),
          " ESVs will be removed (all reads lost after min-abundance filter)")
  
  db.corr2 <- db.corr2 %>% filter(corrected_reads > 0)
  
} else {
  message("Min-abundance filter skipped (min_abundance = 0)")
}

# Export
subfolder <- "Intermediate"
minabund <- paste0(project_name, "_table_tagswitch_", tag_switch, "_min_abundance_",
                   min_abundance, "_", outputfile_suffix, ".csv")
write.csv2(db.corr2, file.path(subfolder, minabund), row.names = FALSE)


# --- 4: Blank correction (strict) ----------------------------------------------

# For every ESV, the maximum read count detected in the specified negative controls is
# determined and subtracted from all sample abundances. Values below zero are set to zero.

db.corr3 <- db.corr2   # copy db.corr2 for these manipulations

if (exists("negatives") && length(negatives) > 0) {
  
  neg_cols <- intersect(negatives, sample_cols)
  
  if (length(neg_cols) == 0) {
    warning("No matching negative control columns found. Skipping blank correction.")
  } else {
    
    message("Applying blank correction using negatives: ",
            paste(neg_cols, collapse=", "))
    
    # max per ESV across negatives
    row_max <- apply(db.corr3[, neg_cols, drop = FALSE], 1,
                     max, na.rm = TRUE)
    
    # subtract from each sample column
    for (sc in sample_cols) {
      db.corr3[[sc]] <- pmax(db.corr3[[sc]] - row_max, 0)
    }
    
    # recompute reads
    db.corr3$corrected_reads <- rowSums(db.corr3[, sample_cols, drop = FALSE])
    
    message(sum(db.corr3$corrected_reads == 0),
            " ESVs removed after blank correction.")
    
    # remove ESVs with zero reads
    db.corr3 <- db.corr3 %>% filter(corrected_reads > 0)
  }
  
} else {
  message("Blank correction skipped (no negatives specified).")
}

# remove sample columns with zero reads
sample_sums <- colSums(db.corr3[, sample_cols, drop = FALSE])
zero_samples <- names(sample_sums[sample_sums == 0])

if (length(zero_samples) > 0) {
  
  message(length(zero_samples),
          " sample(s) with 0 reads removed after blank correction:")
  message(paste(zero_samples, collapse = ", "))
  
  db.corr3 <- db.corr3 %>% select(-all_of(zero_samples))
  
} else {
  message("No samples with 0 reads removed after blank correction.")
}


# Recalculate sample_cols after removing zero-sample columns
sample_cols <- intersect(sample_cols, names(db.corr3))
message("Updated sample_cols after blank correction:")
print(sample_cols)

# Export
subfolder <- "Output"
blankcorr <- paste0(project_name, "_filtered_minabund_", min_abundance,
                    "_blanks_", outputfile_suffix, ".csv")
write.csv2(db.corr3, file.path(subfolder, blankcorr), row.names = FALSE)


# [OPTIONAL]Collapse identical species ESVs ----------------------------------------
# Collapse species-level ESVs belonging to the same species
# Collapse / summarise their flags & reads. 
# Should work for apscale_blast/midori and boldigger v3

# ⚠️
# This code aims to provide a pragmatic solution and is certainly not the only
# possible approach. ESVs sharing the same taxonomic assignment above species 
# level  are kept separate to avoid biasing richness estimates too much.


message("\n Extended sample column check before collapsing ")

# Missing sample columns
missing_in_db <- setdiff(sample_cols, names(db.corr3))
if (length(missing_in_db) > 0) {
  message(" These sample_cols do NOT exist in db.corr3:")
  print(missing_in_db)
} else {
  message("All sample_cols exist in db.corr3.")
}

# Check numeric
existing <- intersect(sample_cols, names(db.corr3))
non_numeric <- existing[!sapply(db.corr3[existing], is.numeric)]

if (length(non_numeric) > 0) {
  message(" These sample_cols exist but are NOT numeric:")
  print(non_numeric)
} else {
  message("All sample_cols are numeric.")
}

# Final usable sample columns
valid_sample_cols <- setdiff(existing, non_numeric)

message("\n final sample_cols used for collapsing:")
print(valid_sample_cols)

# Auto detect and check relevant column names for the collapsing
species_col <- if ("species" %in% names(db.corr3)) "species" else 
  if ("Species" %in% names(db.corr3)) "Species" else 
    stop("No species column found!")

identity_col <- if ("pct_identity" %in% names(db.corr3)) "pct_identity" else 
  if ("Similarity" %in% names(db.corr3)) "Similarity" else NA

flags_col <- if ("flags" %in% names(db.corr3)) "flags" else 
  if ("Flag" %in% names(db.corr3)) "Flag" else NA

message("Detected species column: ", species_col)
message("Detected identity column: ", identity_col)
message("Detected flags column: ", flags_col, "\n")


# Select the collapsable rows 
if ("selected_level" %in% names(db.corr3)) {
  mergeable_rows <- db.corr3 %>%
    filter(selected_level == "species", !is.na(.data[[species_col]]))
} else {
  mergeable_rows <- db.corr3 %>%
    filter(!is.na(.data[[species_col]]))
}

message("Collapsalpe rows detected: ", nrow(mergeable_rows))


# collapse the flags as well 
collapse_flag_values <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  paste(sort(unique(unlist(strsplit(x, "\\|")))), collapse = "|")
}

# Define metatada columns
read_cols <- intersect(c("total_reads", "corrected_reads"), names(db.corr3))

metadata_cols <- setdiff(
  names(mergeable_rows),
  c(sample_cols, read_cols, species_col)
)

# Remove old flags column from metadata
metadata_cols <- metadata_cols[metadata_cols != flags_col]


# Perform the collapsing 
message("\n perform species collapsing ")

summarized_db <- mergeable_rows %>%
  group_by(.data[[species_col]]) %>%
  summarise(
    across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)),
    across(all_of(read_cols), ~ sum(.x, na.rm = TRUE)),
    flags = collapse_flag_values(.data[[flags_col]]),
    across(all_of(metadata_cols), ~ dplyr::first(.x)),
    row_count = n(),
    .groups = "drop"
  )

colnames(summarized_db)[colnames(summarized_db) == species_col] <- species_col

# select non-collapsable rows
remaining_db <- db.corr3 %>%
  filter(!hash %in% mergeable_rows$hash) %>%
  mutate(row_count = NA_integer_)

# combine the dataframes
collapsed_db <- bind_rows(remaining_db, summarized_db)

# reorder original + row_count after flags
original_columns <- colnames(db.corr3)
insert_after <- match(flags_col, original_columns)

if (!is.na(insert_after)) {
  new_order <- append(original_columns, "row_count", after = insert_after)
} else {
  new_order <- c(original_columns, "row_count")
}

collapsed_db <- collapsed_db %>% select(any_of(new_order))

# Sort by identity and then descending by total reads
if (!is.na(identity_col)) {
  collapsed_db <- collapsed_db %>%
    arrange(desc(.data[[identity_col]]), desc(total_reads))
}

# Print the summary of these steps
message("\n Post collapsing summary ")
message("Eligible rows:   ", nrow(mergeable_rows))
message("Species merged:  ", n_distinct(summarized_db[[species_col]]))
message("Final rows:      ", nrow(collapsed_db))

# export the table
subfolder <- "Output"
collapsed_file <- paste0(
  project_name, "_ESV_table_",
  outputfile_suffix, "_filtered_",
  min_abundance, "_collapsed.csv"
)

if (!dir.exists(subfolder)) dir.create(subfolder, recursive = TRUE)

filepath <- file.path(subfolder, collapsed_file)
write.csv2(collapsed_db, filepath, row.names = FALSE)

message("✔ Collapsed table written to: ", filepath)

