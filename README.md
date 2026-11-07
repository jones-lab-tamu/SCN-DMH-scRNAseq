
# SCN-DMH-scRNAseq

Analysis code for the single-cell RNA-sequencing characterization of dorsomedial hypothalamus (DMH)-projecting suprachiasmatic nucleus (SCN) neurons, conducted in the Jones Lab at Texas A&M University.

## Contents

- `scn_pipeline_Whole_Droplet.R`: R analysis pipeline used to generate the manuscript's SCN transcriptomic figures and supplementary analyses.
- `sessionInfo.txt`: R session information documenting the software environment and package versions used for analysis.

The analysis script was developed by Mustafa Ozturgut in the McCreedy Lab at Texas A&M University.

## Sequencing data

Raw sequencing reads and processed Cell Ranger count matrices have been submitted to NCBI Gene Expression Omnibus (GEO).

**GEO accession:** Pending assignment.

The dataset comprises pooled SCN tissue from ten mice following retrograde labeling of DMH-projecting neurons, analyzed using 10x Genomics Flex Gene Expression sequencing with custom EGFP probes.

## Running the analysis

The script was developed using R and the package versions documented in `sessionInfo.txt`.

It expects the following directories within its working folder:

```text
Jones Lab/
├── Jones Lab data/
│   ├── matrix.mtx.gz
│   ├── features.tsv.gz
│   └── barcodes.tsv.gz
├── reference_data/
│   ├── GSE167927_CT7.5_SCN_Cell_Group_Genes.xlsx
│   ├── GSE167927_CT7.5_Cell_Groups_All_Sequenced_Cells.xlsx
│   └── Wen2020_NatNeurosci_Supplementary_Tables_1-6.xlsx
└── scn_pipeline_Whole_Droplet.R
```

The count matrix can be obtained from the GEO dataset. The script uses the Matrix Market representation of the Cell Ranger raw feature-barcode matrix. The GEO deposit provides an HDF5 representation of these counts, which can be converted to Matrix Market format or obtained in that format from the corresponding Cell Ranger output.

The three reference workbooks contain published SCN gene lists from Morris et al. (2021; GEO accession GSE167927) and Wen et al. (2020). Obtain these files from the original publications and associated supplementary materials, retaining the filenames shown above.

Before running the script, update the `setwd()` statement near the beginning to point to the local analysis directory.

The script produces figures and analysis outputs in the `whole droplet paper figures` directory.

## Reproducibility

The pipeline includes cell-level quality control, neuronal and putative SCN identification, whole-droplet EGFP transcript recovery, clustering, differential-expression analysis, EGFP detection-threshold sensitivity analyses, and comparisons with published SCN transcriptional signatures.

## Code attribution

Copyright © Mustafa Ozturgut. Analysis conducted for the Jones Lab, Texas A&M University.
