# Per-target coverage histograms
Create histograms of per-target read counts using per-target coverage stats generated with Picard. A list of *per_target_coverage.metrics files for one or more samples needs to be provided as input. 

The script produces:
- histograms of read counts per target for each sample
- violin plots of read counts for targets in each sample
- violin plots of normalized coverage values are produced
- csv with low covered regions based on a specified mean_coverage threshold shared across samples
- excel with normalized coverage and mean_coverage of the shared lowly covered regions across samples (as provided in the per target coverage metrics file of Picard).

# Requirements
- ggplot2
- dplyr
- openxlsx
