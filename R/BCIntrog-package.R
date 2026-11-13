#' BCIntrog: Identify donor-derived introgression segments
#'
#' BCIntrog identifies donor-parent-derived genomic segments in plant backcross
#' populations from whole-genome resequencing VCF data. The package filters
#' parental informative markers, runs a three-state hidden Markov model
#' (D = donor homozygous, R = recurrent homozygous, H = heterozygous), detects
#' introgressed segments, summarizes donor and recurrent genome proportions, and
#' generates genome-wide visualizations.
#'
#' @keywords internal
#' @importFrom rlang .data
"_PACKAGE"

utils::globalVariables(c(
  ".data", "CHROM", "POS", "REF", "ALT", "Genotype", "State", "Viterbi",
  "Start", "End", "Length_bp", "SNP_number", "Observed_SNP_number",
  "segment_id", "Prob_D", "Prob_R", "Prob_H", "Confidence",
  "Mean_Prob_D", "Mean_Prob_R", "Mean_Prob_H", "Mean_Confidence"
))
