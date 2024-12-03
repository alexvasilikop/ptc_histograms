# Per-target coverage histograms
Create histograms of per-target read counts using per-target coverage stats generated with Picard. A list of *per_target_coverage.metrics files for one or more samples needs to be provided as input. The script produces histograms of read counts per target for all targets and violin plots of read counts for targets in each sample. Also violin plots of normalized coverage values are produced.

# Requirements
- ggplot2
- dplyr
- openxlsx
