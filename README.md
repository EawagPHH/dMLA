#README

##Multiplex-testing on *E. coli* Samples Analysis

This project analyzes the performance of a multiplex assay on two sets of *E. coli* isolates. The script ```ecoli_samples.R``` processes sequencing data, calculates thresholds for detecting true positives/negatives, and visualizes the results.

###Prerequisites
The analysis requires R and the following packages:

+ ggplot2
+ dplyr
+ gridExtra
+ ggpubr
+ tidyr
+ seqinr
+ readr

###Data Files
Ensure the following CSV files are in the data directory:

+ ```ecoli_set1_output.csv```: Contains the read counts for set 1.
+ ```ecoli_set2_output.csv```: Contains the read counts for set 2.
+ ```wgs_results_ecoli_set1.csv```: Contains whole genome sequencing results for set 1.
+ ```wgs_results_ecoli_set2.csv```: Contains whole genome sequencing results for set 2.
+ ```ecoli_merged.csv```: Merged results of set 1 and set 2 for final visualization.

###Analysis Workflow
1. **Data Import and Formatting:** The script reads and formats read count data for E. coli sets 1 and 2, mapping sample IDs to their actual names.
2. **False Positive Distribution:** Filters negative controls and generates histograms showing the distribution of false positives across different probe pairs.
3. **Threshold Calculation:** Calculates the 99.9% quantile of the gamma distribution of false positives for each target to set detection thresholds. These thresholds are used to adjust the read counts, ensuring that counts below the threshold are considered zero.
4. **WGS Comparison:** Compares adjusted read counts with whole genome sequencing (WGS) data to categorize results as true positives, true negatives, false positives, or false negatives.
5. **Visualization:** 
	+ Generates heatmaps showing the detection status of each sample-probe pair.
	+ Combines results from both sets to provide a comprehensive view of the assay performance.
6. **Duplicate Analysis:** Assesses consistency in probe detection between duplicate samples to evaluate assay reliability.
7. **Results Export:** Exports cleaned and annotated results to CSV files for further review.

### Output
**Heatmaps:** Visual representations of the assay performance, highlighting true/false positives and negatives.
**Summary Tables:** Counts of true/false positives/negatives for quality control.

### How to Run
1. Set your working directory to the location of your data files.
2. Run the script in R to perform the analysis.
3. Review the plots and exported CSV files for insights into assay performance.

### Saving Results
The script includes commented-out lines to save processed data (```write.csv``` functions). Uncomment these lines to save the results as CSV files.