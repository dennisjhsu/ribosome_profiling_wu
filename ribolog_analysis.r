## ribolog_analysis.R
## Ribolog-based analysis of ribosome profiling data
## CELP bias correction computed separately per condition
##
## Output: per-gene text files containing codon-level observed counts,
## CELP bias coefficients, and corrected counts for each sample.

library(ggplot2)
library(Ribolog)
library(dplyr)

# ---- CONFIGURATION ----
# Update all paths below to match your directory structure

annotation_file   = "path/to/annotation.txt"
longest_cds_file  = "path/to/cDNA_longest_CDS.txt"
rpf_bam_folder    = "path/to/rpf_deduplicated_bams/"
output_base_dir   = "path/to/output/"

# Sample indices for each condition (adjust to match your BAM file order)
condition_1_idx = 1:3   # e.g., fed
condition_2_idx = 4:6   # e.g., starved

# Read length range for P-site analysis
read_lengths = c(20:32)

# ---- MODULE 1: P-SITE OFFSET AND CODON COUNTS ----

annotation = Ribolog::read_annotation(annotation_file)

# Load RPF BAMs
# refseq_sep = "." strips Ensembl version suffixes (e.g., ENST00000000233.10 -> ENST00000000233)
# This is handled automatically by bamtolist_rW but NOT by bam2count
reads_list = Ribolog::bamtolist_rW(
  bamfolder  = rpf_bam_folder,
  annotation = annotation,
  refseq_sep = "."
)

# Compute P-site offsets and assign P-site positions
psite_offset = Ribolog::psite_rW(reads_list)
reads_psite  = Ribolog::psite_info_rW(reads_list, psite_offset)

# QC plots
Ribolog::print_read_ldist(reads_list, "rpf_read_length_distribution.pdf")
Ribolog::print_period_region(reads_psite, "periodicity_by_region.pdf")
Ribolog::print_period_region_length(reads_psite, "periodicity_by_length.pdf")

# ---- MODULE 1: CELP BIAS CORRECTION (per condition) ----

# Generate codon-level read counts per condition
codon_counts_cond1 = Ribolog::psite_to_codon_count(
  reads_psite[condition_1_idx], read_lengths, annotation, longest_cds_file
)
codon_counts_cond2 = Ribolog::psite_to_codon_count(
  reads_psite[condition_2_idx], read_lengths, annotation, longest_cds_file
)

# Run CELP bias correction per condition
celp_cond1 = Ribolog::CELP_bias(codon_counts_cond1)
celp_cond2 = Ribolog::CELP_bias(codon_counts_cond2)

# ---- SAVE PER-GENE CODON-LEVEL OUTPUT ----
# Each file contains: codon_number, codon_type, aa_type,
#                      observed_count, bias_coefficient, corrected_count
#
# Output structure:
#   output_base_dir/
#     sample_name_1/
#       ENST00000000233.txt
#       ENST00000000412.txt
#       ...
#     sample_name_2/
#       ...

save_celp_per_gene = function(celp_result, out_dir) {
  corrected = celp_result$tr_codon_read_count_loess_corrected
  for (sample_name in names(corrected)) {
    sample_dir = file.path(out_dir, sample_name)
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    transcripts = names(corrected[[sample_name]])
    for (tr in transcripts) {
      write.table(
        corrected[[sample_name]][tr],
        file = file.path(sample_dir, paste0(tr, ".txt")),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
    }
    cat("Saved", length(transcripts), "gene files for", sample_name, "\n")
  }
}

save_celp_per_gene(celp_cond1, output_base_dir)
save_celp_per_gene(celp_cond2, output_base_dir)
