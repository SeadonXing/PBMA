# Code adapted from https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-gex-count

#!/bin/bash

export PATH=$PATH:/Software/cellranger-6.1.2

cellranger count --id=BrM02_count \

                 --transcriptome=/GenomeFiles/Human/Ensembl/GRCh38/CellRanger/refdata-cellranger-GRCh38-1.2.0 \

                 --fastqs=/Project/BrM/CellRanger_fastq/BrM02 \

                 --expect-cells=10000 \

                 --localcores=16 \

                 --localmem=32
