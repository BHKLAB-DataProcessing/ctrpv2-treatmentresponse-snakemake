# Data Sources

## Cancer Therapeutics Response Portal (CTRP) Dataset

### Overview

The Cancer Therapeutics Response Portal (CTRP) is a public resource providing information about the response of cancer cell lines to small-molecule probes and drugs. This project focuses specifically on CTRPv2.0, which includes a comprehensive dataset of drug sensitivity measurements across numerous cancer cell lines.

### External Data Sources

#### Cell Line Information

- **Name**: CTRP Cell Line Information
- **Version/Date**: CTRPv2.0
- **URL**: https://ctd2-data.nci.nih.gov/Public/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/
- **Access Method**: Direct download
- **Data Format**: Tab-separated values (TSV)
- **Citation**: Seashore-Ludlow, B., Rees, M.G., Cheah, J.H., et al. (2015). Harnessing Connectivity in a Large-Scale Small-Molecule Sensitivity Dataset. Cancer Discovery, 5(11), 1210-1223. https://doi.org/10.1158/2159-8290.CD-15-0235
- **License**: [CTD² Data Portal Terms of Use](https://ctd2-data.nci.nih.gov/Public/termsofuse/)

#### Treatment Information

- **Name**: CTRP Compound Information
- **Version/Date**: CTRPv2.0
- **URL**: https://ctd2-data.nci.nih.gov/Public/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/
- **Access Method**: Direct download
- **Data Format**: Tab-separated values (TSV)
- **Citation**: Same as cell line information
- **License**: [CTD² Data Portal Terms of Use](https://ctd2-data.nci.nih.gov/Public/termsofuse/)

#### Drug Response Data

- **Name**: CTRP Drug Sensitivity Data
- **Version/Date**: CTRPv2.0
- **URL**: https://ctd2-data.nci.nih.gov/Public/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/
- **Access Method**: Direct download
- **Data Format**: Tab-separated values (TSV)
- **Citation**: Same as cell line information
- **License**: [CTD² Data Portal Terms of Use](https://ctd2-data.nci.nih.gov/Public/termsofuse/)

### Processed Data

#### Sample Metadata (CTRPv2.0_sampleMetadata.tsv)

- **Name**: Preprocessed Sample Metadata
- **Input Data**: CTRP Cell Line Information
- **Processing Scripts**: `workflow/scripts/R/preprocessMetadata.R`
- **Output Location**: `metadata/CTRPv2.0_sampleMetadata.tsv`

| Column Name         | Description                         | Example Values                                    |
| ------------------- | ----------------------------------- | ------------------------------------------------- |
| sampleid            | Unique identifier for the cell line | 697, 5637, 2313287                                |
| tissueid            | Tissue of origin                    | haematopoietic_and_lymphoid_tissue, urinary_tract |
| master_ccl_id       | Master cell line identifier         | 1, 3, 4                                           |
| sample_availability | Where the sample is available       | ccle;public, ccle                                 |
| ccle_primary_hist   | Primary histology                   | lymphoid_neoplasm, carcinoma, glioma              |

#### Treatment Metadata (CTRPv2.0_treatmentMetadata.tsv)

- **Name**: Preprocessed Treatment Metadata
- **Input Data**: CTRP Compound Information
- **Processing Scripts**: `workflow/scripts/R/preprocessMetadata.R`
- **Output Location**: `metadata/CTRPv2.0_treatmentMetadata.tsv`

| Column Name                    | Description                          | Example Values                                  |
| ------------------------------ | ------------------------------------ | ----------------------------------------------- |
| master_cpd_id                  | Master compound identifier           | 1788, 3588, 12877                               |
| cpd_name                       | Compound name                        | CIL55, BRD4132, BRD6340                         |
| broad_cpd_id                   | Broad Institute compound identifier  | BRD-K46556387, BRD-K86574132, BRD-K35716340     |
| top_test_conc_umol             | Top test concentration in μM         | 10, 160, 33                                     |
| cpd_status                     | Compound status (probe, FDA)         | probe, FDA                                      |
| inclusion_rationale            | Rationale for inclusion in the study | pilot-set, chromatin;pilot-set                  |
| gene_symbol_of_protein_target  | Gene symbol for protein target       | S1PR3, BAX, RARA;RARB;RARG                      |
| target_or_activity_of_compound | Target or activity description       | agonist of sphingosine 1-phosphate receptor 3   |
| source_name                    | Source of the compound               | Columbia University, ChemDiv Inc., Enamine Ltd. |
| source_catalog_id              | Catalog identifier from source       | 4998-1380, 1988-0090                            |
| cpd_smiles                     | SMILES chemical notation             | Chemical structure in SMILES format             |
| treatmentid                    | Standardized treatment ID            | CIL55, BRD4132, BRD6340                         |

#### Treatment Response: Raw

- **Name**: Preprocessed Raw Treatment Response
- **Input Data**: CTRP Drug Sensitivity Data
- **Processing Scripts**: `workflow/scripts/R/preprocessTreatmentResponse.R`
- **Output Location**: `data/procdata/CTRPv2.0_preprocessed_TreatmentResponse_raw.csv`

| Column Name   | Description                  | Units/Notes                            |
| ------------- | ---------------------------- | -------------------------------------- |
| sampleid      | Cell line identifier         | Same as `sampleid` in metadata         |
| culture_media | Culture media used           | DMEM001, RPMI001, etc.                 |
| treatmentid   | Drug identifier              | Same as `treatmentid` in metadata      |
| experiment_id | Experiment identifier        | Numerical identifier for experiment    |
| dose          | Drug concentration           | μM (micromolar)                        |
| viability     | Cell viability at given dose | Percentage, values relative to control |
| tech_rep      | Technical replicate number   | 1, 2, 3, etc.                          |
