#README

## Overview

This project involves analyzing multiplex assay performance on *E. coli* isolates and synthetic DNA templates using two R scripts: ```positive_samples.R``` for synthetic DNA templates and ```ecoli_samples.R ```for *E. coli* samples. The analysis includes data processing, threshold calculations for detecting true positives/negatives, and result visualization.

## Prerequisites

The analysis requires R and the following packages:

+ ggplot2
+ dplyr
+ gridExtra
+ ggpubr
+ tidyr
+ seqinr
+ readr

## Scripts and Data Files

### 1. ```positive_samples.R``` - Synthetic DNA Templates
This script analyzes the read counts from two replicates of synthetic DNA templates to distinguish true and false positives/negatives.

#### Data Files:

+ ```positive_output_tube1.csv```: Data for Replicate 1.
+ ```positive_output_tube2.csv```: Data for Replicate 2.

#### Workflow:

1. **Data Import and Formatting:** Loads data and maps sample IDs to corresponding codes.
2. **Initial Visualization:** Creates heatmaps showing the distribution of read counts before filtering.
3. **False Positive Filtering:** Identifies false positives, calculates a 99.9% quantile threshold using a gamma distribution, and adjusts read counts.
4. **Error Categorization:** Counts true/false positives/negatives based on the adjusted read counts.
5. **Filtered Visualization:** Generates heatmaps of adjusted read counts.
6. **Aggregate Analysis:** Merges replicates, calculates mean counts, and visualizes in a bar plot.

#### Output:

+ Adjusted read counts (e.g., ```positive_adjusted_tube1.csv```).
+ Plots comparing read counts before and after filtering, with aggregate analysis results.

### 2. ```ecoli_samples.R``` - *E. coli* Samples
This script analyzes multiplex assay performance on two sets of *E. coli* isolates.

#### Data Files:

+ ```ecoli_set1_output.csv```: Read counts for set 1.
+ ```ecoli_set2_output.csv```: Read counts for set 2.
+ ```wgs_results_ecoli_set1.csv```: Whole genome sequencing results for set 1.
+ ```wgs_results_ecoli_set2.csv```: Whole genome sequencing results for set 2.
+ ```ecoli_merged.csv```: Merged results for final visualization.

#### Workflow:

1. **Data Import and Formatting:** Reads and formats data, mapping sample IDs to actual names.
2. **False Positive Analysis:** Filters negative controls and visualizes the distribution of false positives.
3. **Threshold Calculation:** Sets detection thresholds using a 99.9% quantile of the gamma distribution and adjusts read counts accordingly.
4. **WGS Comparison:** Compares adjusted counts with WGS data to classify results as true positives/negatives or false positives/negatives.
5. **Visualization:** Generates heatmaps and summary tables to assess assay performance and reliability.
6. **Duplicate Analysis:** Evaluates consistency between duplicates.

#### Output:

+ Heatmaps and summary tables of assay performance.
+ Exported cleaned and annotated results.

#### How to Run

1. Set your working directory to the location containing the data files.
2. Run each script in R with the required packages installed.
3. Review generated plots and CSV files for insights into assay performance.

#### Saving Results

Uncomment the ```write.csv``` lines in the scripts to save processed data.

This project provides a comprehensive analysis of multiplex assay results by filtering false positives and comparing to WGS data, aiding in the evaluation of assay accuracy and reliability.