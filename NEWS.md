# BCIntrog 0.1.0

## Initial release

- Implemented `read_vcf()` for reading VCF/VCF.GZ and extracting CHROM, POS, REF, ALT, GT and GQ.
- Implemented `select_site()` for selecting biallelic markers with opposite homozygous parents, filtering by GQ, and calling offspring D/R/H genotypes.
- Implemented `emission_probability()` and `transition_probability()`.
- Implemented `run_hmm()` with Viterbi and Forward-Backward; missing observations are treated as uninformative.
- Implemented `segment_detection()` with 1-based closed intervals, small-segment processing, and at most two iterations. Its `segments` table contains `CHROM`, `segment_id`, `Start`, `End`, `Length_bp`, `SNP_number`, `Observed_SNP_number`, `Viterbi` and `State`.
- Implemented `summary.donor_segments()` for donor genome equivalents, recurrent parent recovery, and per-state/per-chromosome summaries.
- Implemented `visualization()` for genome-wide marker and segment plots.
- Implemented `run_analysis()` as a one-step pipeline with optional result files.
- Added five vignettes, unit tests, and a full Demo integration test.
- `R CMD check`: 0 errors, 0 warnings, 0 notes.
