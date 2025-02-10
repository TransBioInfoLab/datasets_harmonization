# Harmonization of DNA methylation datasets

A repository for harmonizing 5 datasets: AIBL, ADNI, FHS, MESA, and HRS

Process:
Step 1 - Process each Dataset
1. Read in the Raw data
2. QC samples
3. QC probes
4. Identify low confidence probes using detection P values, and impute them
5. Normalize data

Step 2 - Harmonize with each Dataset
We batch-correct data using Harman method
