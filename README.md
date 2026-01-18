# Harmonization of DNA methylation datasets

This repository contains scripts for harmonizing DNA methylation (DNAm) datasets from four studies: **ADNI**, **FHS**, **MESA**, and **HRS**. The harmonized DNAm data are intended for evaluating methylation-based predictors and currently include **baseline** DNAm from **cognitively unimpaired** participants only.

## Overview of processing steps

1. **Sample-level QC**
   - Removed samples with discordant predicted vs. recorded sex
   - Removed samples with bisulfite conversion rate < 85%
   - Removed samples with detection *P*-values > 0.01 for > 5% of autosomal probes

2. **Probe-level QC**
   - Removed probes that did not start with `"cg"`
   - Removed probes on the mitochondrial chromosome (`chrM`)
   - Removed probes missing from the annotation file

3. **Imputation**
   - A moderate number of probes had missing values.
   - Missing probe values (and values with detection *P* > 0.01) were imputed using the `methyLImp2` R package.

4. **Platform harmonization**
   - All datasets (ADNI, FHS, MESA, HRS) used Illumina EPIC arrays; therefore, platform harmonization was not required.

5. **Normalization**
   - Applied BMIQ normalization within each dataset.

6. **Batch effect correction**
   - Performed using the `harman` R package.

## Computational environment

This pipeline was executed on a workstation running **Ubuntu 24.04.3 LTS** with **16 CPU cores** and **1.5 TB RAM**. Peak observed memory usage was approximately **102 GB** (for the FHS dataset).
